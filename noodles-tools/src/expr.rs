//! samtools filter expressions (`view -e`): HTSlib's expression language (`hts_expr.c`) over the
//! record fields and tags of `bam_sym_lookup` (`sam.c`), as of HTSlib 1.24.
//!
//! HTSlib parses and evaluates the expression together, once per record. Here it is parsed once
//! into a tree, and evaluating the tree repeats HTSlib's value handling step by step, including
//! its quirks: no short-circuiting, "undefined" values (absent tags, `NaN`), a tag's presence
//! counting as true even when its value is zero, and values that keep fields from earlier steps.
//! Syntax errors are reported when the expression is parsed, before any record is read; errors
//! that depend on the data (arithmetic on strings, array tags) are reported for the first
//! record that causes one, as in HTSlib.

use std::{cell::OnceCell, cmp::Ordering, fmt::Write as _};

use anyhow::{anyhow, bail};
use noodles_sam::{self as sam, header::record::value::map::read_group::tag::LIBRARY};
use regex::bytes::{Regex, RegexBuilder};

use crate::bam_raw;

const FLAG_UNMAPPED: u16 = 0x4;
const CIGAR_OPS: &[u8; 16] = b"MIDNSHP=XB??????";
const CIGAR_SOFT_CLIP: u32 = 4;
const CIGAR_HARD_CLIP: u32 = 5;

/// A parsed filter expression.
#[derive(Debug)]
pub struct Expression {
    root: Node,
}

impl Expression {
    /// Parses an expression.
    ///
    /// # Errors
    ///
    /// Returns an error if the expression is malformed or names an unknown field or function.
    pub fn parse(text: &str) -> anyhow::Result<Self> {
        let parser = Parser {
            src: text.as_bytes(),
        };

        let (root, end) = parser
            .expression(0)
            .map_err(|e| anyhow!("invalid filter expression {text:?}: {e}"))?;

        if parser.at(parser.ws(end)) != 0 {
            bail!("invalid filter expression {text:?}: unable to parse expression at {text}");
        }

        Ok(Self { root })
    }

    /// Returns whether an encoded BAM record passes (`sam_passes_filter`).
    ///
    /// # Errors
    ///
    /// Returns an error if the expression cannot be applied to the record, e.g., arithmetic on a
    /// string or an array tag.
    pub fn matches(
        &self,
        record: bam_raw::Record<'_>,
        header: &sam::Header,
    ) -> anyhow::Result<bool> {
        let cx = Context {
            record,
            header,
            cigar: OnceCell::new(),
        };

        let mut res = Value::default();
        eval(&self.root, &cx, &mut res).map_err(|e| {
            anyhow!(
                "couldn't process filter expression for record {}: {e}",
                String::from_utf8_lossy(record.name())
            )
        })?;

        // `hts_filter_eval_`: strings are true if present (even empty); numbers if nonzero.
        if res.is_str {
            res.is_true |= res.s.is_some();
        } else if res.exists() {
            res.is_true |= res.d != 0.0;
        }

        Ok(res.is_true)
    }
}

/// `hts_expr_val_t`.
#[derive(Clone, Debug, Default)]
struct Value {
    /// Whether `s` (rather than `d`) is the value.
    is_str: bool,
    /// Forces the value to be true (set by tags that are present).
    is_true: bool,
    /// The string; `None` is HTSlib's NULL string, which is undefined.
    s: Option<Vec<u8>>,
    /// The number; `NaN` is undefined.
    d: f64,
}

impl Value {
    /// `hts_expr_val_exists`
    fn exists(&self) -> bool {
        !(self.is_str && self.s.is_none()) && !(!self.is_str && self.d.is_nan())
    }

    /// `hts_expr_val_existsT`: defined, or undefined but forced true.
    fn exists_t(&self) -> bool {
        self.is_true || self.exists()
    }

    /// `hts_expr_val_undef`. As in HTSlib, an allocated string stays allocated (but empty).
    fn set_undefined(&mut self) {
        self.clear_string();
        self.is_true = false;
        self.is_str = false;
        self.d = f64::NAN;
    }

    /// `ks_clear`: empties the string without freeing it.
    fn clear_string(&mut self) {
        if let Some(s) = &mut self.s {
            s.clear();
        }
    }

    fn set_string(&mut self, value: &[u8]) {
        self.is_str = true;
        let s = self.s.get_or_insert_with(Vec::new);
        s.clear();
        s.extend_from_slice(value);
    }

    fn set_number(&mut self, d: f64) {
        self.is_str = false;
        self.d = d;
    }

    fn set_bool(&mut self, b: bool) {
        self.is_true = b;
        self.d = f64::from(u8::from(b));
    }

    /// The truth test of `&&` and `||`: forced true, a defined string, or a nonzero (or `NaN`)
    /// number.
    fn truth(&self) -> bool {
        self.is_true || (self.is_str && self.s.is_some()) || self.d != 0.0
    }

    /// The string as C sees it, up to the first NUL.
    fn c_str(&self) -> &[u8] {
        c_str(self.s.as_deref().unwrap_or_default())
    }
}

fn c_str(s: &[u8]) -> &[u8] {
    memchr::memchr(0, s).map_or(s, |i| &s[..i])
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Symbol {
    Cigar,
    EndPos,
    Flag(Option<u16>),
    HardClipLength,
    Library,
    MappingQuality,
    MatePosition,
    MateReferenceName,
    MateReferenceId,
    CigarOpCount,
    Position,
    QueryLength,
    Name,
    Quality,
    ReferenceLength,
    ReferenceName,
    ReferenceId,
    Sequence,
    SoftClipLength,
    TemplateLength,
    Tag([u8; 2]),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Function {
    Average,
    Exists,
    Exp,
    Length,
    Log,
    Max,
    Min,
    Sqrt,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum UnaryOp {
    Plus,
    Minus,
    Not,
    BitNot,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ArithOp {
    Mul,
    Div,
    Rem,
    Add,
    Sub,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum BitOp {
    And,
    Xor,
    Or,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum CmpOp {
    Ge,
    Gt,
    Le,
    Lt,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum EqOp {
    Eq,
    Ne,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum LogicOp {
    And,
    Or,
}

#[derive(Debug)]
enum Node {
    Number(f64),
    String(Vec<u8>),
    Symbol(Symbol),
    Function(Function, Box<Node>),
    Default(Box<Node>, Box<Node>),
    Pow(Box<Node>, Box<Node>),
    Unary(UnaryOp, Box<Node>),
    /// `*`, `/`, `%` or `+`, `-` chains; HTSlib reuses one value for every right operand.
    Arith(Box<Node>, Vec<(ArithOp, Node)>),
    /// One bitwise precedence level.
    Bit(BitOp, Box<Node>, Vec<Node>),
    Compare(CmpOp, Box<Node>, Box<Node>),
    Equal(EqOp, Box<Node>, Box<Node>),
    Match {
        negated: bool,
        lhs: Box<Node>,
        rhs: Box<Node>,
        /// The compiled pattern when `rhs` is a string literal.
        regex: Option<Regex>,
    },
    Logic(Box<Node>, Vec<(LogicOp, Node)>),
}

type ParseResult = Result<(Node, usize), String>;

/// A recursive-descent parser with HTSlib's grammar, precedence, and whitespace rules: only
/// spaces and tabs separate tokens, and a function name must be followed directly by `(`.
struct Parser<'a> {
    src: &'a [u8],
}

impl Parser<'_> {
    /// The byte at `i`, or NUL past the end (as in a C string).
    fn at(&self, i: usize) -> u8 {
        self.src.get(i).copied().unwrap_or(0)
    }

    fn starts_with(&self, i: usize, prefix: &[u8]) -> bool {
        self.src.get(i..).is_some_and(|s| s.starts_with(prefix))
    }

    /// `ws`: skips spaces and tabs.
    fn ws(&self, mut i: usize) -> usize {
        while matches!(self.at(i), b' ' | b'\t') {
            i += 1;
        }
        i
    }

    fn expression(&self, i: usize) -> ParseResult {
        self.and_expr(i)
    }

    /// `&&` and `||`, at one precedence level, left to right.
    fn and_expr(&self, i: usize) -> ParseResult {
        let (first, mut end) = self.eq_expr(i)?;
        let mut rest = Vec::new();

        loop {
            let j = self.ws(end);

            let op = if self.starts_with(j, b"&&") {
                LogicOp::And
            } else if self.starts_with(j, b"||") {
                LogicOp::Or
            } else {
                break;
            };

            let (node, e) = self.eq_expr(j + 2)?;
            rest.push((op, node));
            end = e;
        }

        if rest.is_empty() {
            Ok((first, end))
        } else {
            Ok((Node::Logic(Box::new(first), rest), end))
        }
    }

    /// `==`, `!=`, `=~`, `!~`, right-associative.
    fn eq_expr(&self, i: usize) -> ParseResult {
        let (lhs, end) = self.cmp_expr(i)?;
        let j = self.ws(end);

        let (a, b) = (self.at(j), self.at(j + 1));

        let op = match (a, b) {
            (b'=', b'=') => Some(EqOp::Eq),
            (b'!', b'=') => Some(EqOp::Ne),
            (b'=' | b'!', b'~') => None,
            _ => return Ok((lhs, end)),
        };

        let (rhs, end) = self.eq_expr(j + 2)?;

        let node = match op {
            Some(op) => Node::Equal(op, Box::new(lhs), Box::new(rhs)),
            None => {
                let regex = match &rhs {
                    Node::String(pattern) => Some(compile_regex(pattern)?),
                    _ => None,
                };

                Node::Match {
                    negated: a == b'!',
                    lhs: Box::new(lhs),
                    rhs: Box::new(rhs),
                    regex,
                }
            }
        };

        Ok((node, end))
    }

    /// `>=`, `>`, `<=`, `<`, right-associative.
    fn cmp_expr(&self, i: usize) -> ParseResult {
        let (lhs, end) = self.bitor_expr(i)?;
        let j = self.ws(end);

        let (op, len) = match (self.at(j), self.at(j + 1)) {
            (b'>', b'=') => (CmpOp::Ge, 2),
            (b'>', _) => (CmpOp::Gt, 1),
            (b'<', b'=') => (CmpOp::Le, 2),
            (b'<', _) => (CmpOp::Lt, 1),
            _ => return Ok((lhs, end)),
        };

        let (rhs, end) = self.cmp_expr(j + len)?;
        Ok((Node::Compare(op, Box::new(lhs), Box::new(rhs)), end))
    }

    fn bitor_expr(&self, i: usize) -> ParseResult {
        self.bit_level(i, BitOp::Or)
    }

    /// One bitwise level (`|` over `^` over `&`), left to right.
    fn bit_level(&self, i: usize, op: BitOp) -> ParseResult {
        let operand = |p| match op {
            BitOp::Or => self.bit_level(p, BitOp::Xor),
            BitOp::Xor => self.bit_level(p, BitOp::And),
            BitOp::And => self.add_expr(p),
        };

        let (first, mut end) = operand(i)?;
        let mut rest = Vec::new();

        loop {
            let j = self.ws(end);

            let matched = match op {
                BitOp::Or => self.at(j) == b'|' && self.at(j + 1) != b'|',
                BitOp::Xor => self.at(j) == b'^',
                BitOp::And => self.at(j) == b'&' && self.at(j + 1) != b'&',
            };

            if !matched {
                break;
            }

            let (node, e) = operand(j + 1)?;
            rest.push(node);
            end = e;
        }

        if rest.is_empty() {
            Ok((first, end))
        } else {
            Ok((Node::Bit(op, Box::new(first), rest), end))
        }
    }

    fn add_expr(&self, i: usize) -> ParseResult {
        self.arith_level(
            i,
            &[(b'+', ArithOp::Add), (b'-', ArithOp::Sub)],
            Self::mul_expr,
        )
    }

    fn mul_expr(&self, i: usize) -> ParseResult {
        self.arith_level(
            i,
            &[
                (b'*', ArithOp::Mul),
                (b'/', ArithOp::Div),
                (b'%', ArithOp::Rem),
            ],
            Self::unary_expr,
        )
    }

    fn arith_level(
        &self,
        i: usize,
        ops: &[(u8, ArithOp)],
        operand: fn(&Self, usize) -> ParseResult,
    ) -> ParseResult {
        let (first, mut end) = operand(self, i)?;
        let mut rest = Vec::new();

        loop {
            let j = self.ws(end);
            let c = self.at(j);

            let Some(&(_, op)) = ops.iter().find(|(b, _)| *b == c) else {
                break;
            };

            let (node, e) = operand(self, j + 1)?;
            rest.push((op, node));
            end = e;
        }

        if rest.is_empty() {
            Ok((first, end))
        } else {
            Ok((Node::Arith(Box::new(first), rest), end))
        }
    }

    /// Unary `+` and `-` take a simple expression; `!` and `~` take a unary one.
    fn unary_expr(&self, i: usize) -> ParseResult {
        let j = self.ws(i);

        let (op, (node, end)) = match self.at(j) {
            b'+' => (UnaryOp::Plus, self.simple_expr(j + 1)?),
            b'-' => (UnaryOp::Minus, self.simple_expr(j + 1)?),
            b'!' => (UnaryOp::Not, self.unary_expr(j + 1)?),
            b'~' => (UnaryOp::BitNot, self.unary_expr(j + 1)?),
            _ => return self.simple_expr(j),
        };

        Ok((Node::Unary(op, Box::new(node)), end))
    }

    /// A parenthesized expression, number, string, field, or function call.
    fn simple_expr(&self, i: usize) -> ParseResult {
        let i = self.ws(i);

        if self.at(i) == b'(' {
            let (node, end) = self.expression(i + 1)?;
            let j = self.ws(end);

            if self.at(j) != b')' {
                return Err("missing ')'".into());
            }

            return Ok((node, j + 1));
        }

        if let Some((d, end)) = strtod(self.src, i) {
            return Ok((Node::Number(d), end));
        }

        if self.at(i) == b'"' {
            return self.string(i);
        }

        if let Some((symbol, end)) = self.symbol(i)? {
            return Ok((Node::Symbol(symbol), end));
        }

        self.function(i)
    }

    /// A string literal, with `\"`, `\\`, `\t`, `\n`, and `\r` escapes (a backslash before any
    /// other character is kept).
    fn string(&self, i: usize) -> ParseResult {
        let mut e = i + 1;

        while !matches!(self.at(e), 0 | b'"') {
            if self.at(e) == b'\\' {
                e += 1 + usize::from(self.at(e + 1) != 0);
            } else {
                e += 1;
            }
        }

        if self.at(e) != b'"' {
            return Err("unterminated string".into());
        }

        let raw = &self.src[i + 1..e];
        let mut s = Vec::with_capacity(raw.len());
        let mut chars = raw.iter().copied();

        while let Some(c) = chars.next() {
            if c != b'\\' {
                s.push(c);
                continue;
            }

            match chars.next() {
                Some(b'"') => s.push(b'"'),
                Some(b'\\') => s.push(b'\\'),
                Some(b't') => s.push(b'\t'),
                Some(b'n') => s.push(b'\n'),
                Some(b'r') => s.push(b'\r'),
                Some(other) => s.extend([b'\\', other]),
                None => s.push(b'\\'),
            }
        }

        Ok((Node::String(s), e + 1))
    }

    /// `bam_sym_lookup`'s names, matched as prefixes as HTSlib does (so `posx` is `pos`
    /// followed by a syntax error).
    fn symbol(&self, i: usize) -> Result<Option<(Symbol, usize)>, String> {
        const NAMES: &[(&[u8], Symbol)] = &[
            (b"cigar", Symbol::Cigar),
            (b"endpos", Symbol::EndPos),
            (b"hclen", Symbol::HardClipLength),
            (b"library", Symbol::Library),
            (b"mapq", Symbol::MappingQuality),
            (b"mpos", Symbol::MatePosition),
            (b"mrname", Symbol::MateReferenceName),
            (b"mrefid", Symbol::MateReferenceId),
            (b"ncigar", Symbol::CigarOpCount),
            (b"pos", Symbol::Position),
            (b"pnext", Symbol::MatePosition),
            (b"qlen", Symbol::QueryLength),
            (b"qname", Symbol::Name),
            (b"qual", Symbol::Quality),
            (b"rlen", Symbol::ReferenceLength),
            (b"rname", Symbol::ReferenceName),
            (b"rnext", Symbol::MateReferenceName),
            (b"refid", Symbol::ReferenceId),
            (b"seq", Symbol::Sequence),
            (b"sclen", Symbol::SoftClipLength),
            (b"tlen", Symbol::TemplateLength),
        ];

        const FLAGS: &[(&[u8], u16)] = &[
            (b"paired", 0x1),
            (b"proper_pair", 0x2),
            (b"unmap", 0x4),
            (b"munmap", 0x8),
            (b"reverse", 0x10),
            (b"mreverse", 0x20),
            (b"read1", 0x40),
            (b"read2", 0x80),
            (b"secondary", 0x100),
            (b"qcfail", 0x200),
            (b"dup", 0x400),
            (b"supplementary", 0x800),
        ];

        if self.starts_with(i, b"flag") {
            let j = i + 4;

            if self.at(j) != b'.' {
                return Ok(Some((Symbol::Flag(None), j)));
            }

            return FLAGS
                .iter()
                .find(|(name, _)| self.starts_with(j + 1, name))
                .map(|&(name, bit)| Some((Symbol::Flag(Some(bit)), j + 1 + name.len())))
                .ok_or_else(|| "unrecognised flag string".into());
        }

        if self.at(i) == b'['
            && self.at(i + 1) != 0
            && self.at(i + 2) != 0
            && self.at(i + 3) == b']'
        {
            let tag = [self.at(i + 1), self.at(i + 2)];
            return Ok(Some((Symbol::Tag(tag), i + 4)));
        }

        Ok(NAMES
            .iter()
            .find(|(name, _)| self.starts_with(i, name))
            .map(|&(name, symbol)| (symbol, i + name.len())))
    }

    /// `func_expr`: `NAME(expr)`, or `default(a, b)` and `pow(a, b)`, whose comma must follow
    /// the first argument directly.
    fn function(&self, i: usize) -> ParseResult {
        const UNARY: &[(&[u8], Function)] = &[
            (b"avg(", Function::Average),
            (b"exists(", Function::Exists),
            (b"exp(", Function::Exp),
            (b"length(", Function::Length),
            (b"log(", Function::Log),
            (b"min(", Function::Min),
            (b"max(", Function::Max),
            (b"sqrt(", Function::Sqrt),
        ];

        let (node, end) = if let Some(&(name, function)) =
            UNARY.iter().find(|(name, _)| self.starts_with(i, name))
        {
            let (arg, end) = self.expression(i + name.len())?;
            (Node::Function(function, Box::new(arg)), end)
        } else if let Some((name, is_default)) = [(b"default(".as_slice(), true), (b"pow(", false)]
            .into_iter()
            .find(|(name, _)| self.starts_with(i, name))
        {
            let (a, end) = self.expression(i + name.len())?;

            if self.at(end) != b',' {
                return Err("expected ',' after the first argument".into());
            }

            let (b, end) = self.expression(self.ws(end + 1))?;
            let (a, b) = (Box::new(a), Box::new(b));
            let node = if is_default {
                Node::Default(a, b)
            } else {
                Node::Pow(a, b)
            };

            (node, end)
        } else {
            let rest = String::from_utf8_lossy(&self.src[i.min(self.src.len())..]);
            return Err(format!("unknown field or function at {rest:?}"));
        };

        let j = self.ws(end);

        if self.at(j) != b')' {
            return Err("missing ')'".into());
        }

        Ok((node, j + 1))
    }
}

/// C's `strtod` (as used by HTSlib's `hts_str2dbl`): leading whitespace, a sign, then a decimal
/// or hexadecimal number, `inf`, `infinity`, or `nan`. Returns the value and the end of the
/// number, or `None` if there is no number.
fn strtod(src: &[u8], start: usize) -> Option<(f64, usize)> {
    let at = |i: usize| src.get(i).copied().unwrap_or(0);
    let ci = |i: usize, word: &[u8]| {
        src.get(i..i + word.len())
            .is_some_and(|s| s.eq_ignore_ascii_case(word))
    };

    let mut i = start;
    while matches!(at(i), b' ' | b'\t' | b'\n' | b'\x0b' | b'\x0c' | b'\r') {
        i += 1;
    }

    let negative = at(i) == b'-';
    if matches!(at(i), b'+' | b'-') {
        i += 1;
    }

    let sign = if negative { -1.0 } else { 1.0 };

    if ci(i, b"inf") {
        let end = if ci(i, b"infinity") { i + 8 } else { i + 3 };
        return Some((sign * f64::INFINITY, end));
    }

    if ci(i, b"nan") {
        let mut end = i + 3;

        if at(end) == b'(' {
            let mut j = end + 1;
            while at(j).is_ascii_alphanumeric() || at(j) == b'_' {
                j += 1;
            }
            if at(j) == b')' {
                end = j + 1;
            }
        }

        return Some((f64::NAN, end));
    }

    if at(i) == b'0'
        && matches!(at(i + 1), b'x' | b'X')
        && (at(i + 2).is_ascii_hexdigit() || (at(i + 2) == b'.' && at(i + 3).is_ascii_hexdigit()))
    {
        let (d, end) = hex_float(src, i + 2);
        return Some((sign * d, end));
    }

    let digits_start = i;
    let mut digits = 0;

    while at(i).is_ascii_digit() {
        i += 1;
        digits += 1;
    }

    if at(i) == b'.' {
        i += 1;
        while at(i).is_ascii_digit() {
            i += 1;
            digits += 1;
        }
    }

    if digits == 0 {
        return None;
    }

    if matches!(at(i), b'e' | b'E') {
        let mut j = i + 1;
        if matches!(at(j), b'+' | b'-') {
            j += 1;
        }
        if at(j).is_ascii_digit() {
            while at(j).is_ascii_digit() {
                j += 1;
            }
            i = j;
        }
    }

    let text = std::str::from_utf8(&src[digits_start..i]).ok()?;
    let d: f64 = text.parse().ok()?;
    Some((sign * d, i))
}

/// The value of a hexadecimal floating-point number's digits (after `0x`) and optional binary
/// exponent, and its end.
fn hex_float(src: &[u8], start: usize) -> (f64, usize) {
    let at = |i: usize| src.get(i).copied().unwrap_or(0);
    let digit = |c: u8| f64::from((c as char).to_digit(16).unwrap_or(0));

    let mut i = start;
    let mut value = 0.0;
    let mut scale = 0i32;

    while at(i).is_ascii_hexdigit() {
        value = value * 16.0 + digit(at(i));
        i += 1;
    }

    if at(i) == b'.' {
        i += 1;
        while at(i).is_ascii_hexdigit() {
            value = value * 16.0 + digit(at(i));
            scale -= 4;
            i += 1;
        }
    }

    if matches!(at(i), b'p' | b'P') {
        let mut j = i + 1;
        let negative = at(j) == b'-';
        if matches!(at(j), b'+' | b'-') {
            j += 1;
        }
        if at(j).is_ascii_digit() {
            let mut exponent = 0i32;
            while at(j).is_ascii_digit() {
                exponent = exponent
                    .saturating_mul(10)
                    .saturating_add(i32::from(at(j) - b'0'));
                j += 1;
            }
            scale = scale.saturating_add(if negative { -exponent } else { exponent });
            i = j;
        }
    }

    (value * 2f64.powi(scale), i)
}

/// Compiles a POSIX extended regular expression (`regcomp` with `REG_EXTENDED | REG_NOSUB`) as
/// a byte regex.
fn compile_regex(pattern: &[u8]) -> Result<Regex, String> {
    let translated = translate_ere(c_str(pattern))?;

    RegexBuilder::new(&translated)
        .unicode(false)
        .dot_matches_new_line(true)
        .build()
        .map_err(|e| format!("failed regex: {e}"))
}

/// Rewrites a POSIX ERE in the Rust regex syntax: a backslash is literal inside a bracket
/// expression, characters Rust treats as class operators are escaped there, GNU's `` \` `` and
/// `\'` become `\A` and `\z`, a backslash before another letter or digit is the letter itself
/// (GNU's `\w`, `\W`, `\s`, `\S`, `\b`, `\B`, `\<`, and `\>` keep their meaning), and bytes
/// outside ASCII are written as escapes.
fn translate_ere(pattern: &[u8]) -> Result<String, String> {
    let mut out = String::with_capacity(pattern.len() * 2);
    let mut i = 0;

    let push_literal = |out: &mut String, c: u8| {
        if c.is_ascii() && !c.is_ascii_control() {
            out.push_str(&regex::escape(&char::from(c).to_string()));
        } else {
            let _ = write!(out, "\\x{c:02X}");
        }
    };

    while i < pattern.len() {
        let c = pattern[i];

        match c {
            b'\\' => {
                let Some(&next) = pattern.get(i + 1) else {
                    return Err("trailing backslash".into());
                };

                match next {
                    b'w' | b'W' | b's' | b'S' | b'b' | b'B' | b'<' | b'>' => {
                        out.push('\\');
                        out.push(char::from(next));
                    }
                    b'`' => out.push_str("\\A"),
                    b'\'' => out.push_str("\\z"),
                    b'1'..=b'9' => return Err("back-references are not supported".into()),
                    _ => push_literal(&mut out, next),
                }

                i += 2;
            }
            b'[' => {
                i = translate_bracket(pattern, i, &mut out)?;
            }
            _ if !c.is_ascii() || c.is_ascii_control() => {
                push_literal(&mut out, c);
                i += 1;
            }
            _ => {
                out.push(char::from(c));
                i += 1;
            }
        }
    }

    Ok(out)
}

/// Translates the bracket expression starting at `pattern[start]` (`[`), returning the index
/// after it.
fn translate_bracket(pattern: &[u8], start: usize, out: &mut String) -> Result<usize, String> {
    let mut i = start + 1;
    out.push('[');

    if pattern.get(i) == Some(&b'^') {
        out.push('^');
        i += 1;
    }

    // A `]` first in the list is literal.
    if pattern.get(i) == Some(&b']') {
        out.push_str("\\]");
        i += 1;
    }

    loop {
        let Some(&c) = pattern.get(i) else {
            return Err("unmatched [".into());
        };

        match c {
            b']' => {
                out.push(']');
                return Ok(i + 1);
            }
            b'[' if matches!(pattern.get(i + 1), Some(b':' | b'=' | b'.')) => {
                let kind = pattern[i + 1];
                let close = [kind, b']'];
                let body_start = i + 2;
                let len = pattern[body_start..]
                    .windows(2)
                    .position(|w| w == close)
                    .ok_or("unterminated character class")?;
                let body = &pattern[body_start..body_start + len];

                if kind == b':' {
                    out.push_str("[:");
                    out.push_str(std::str::from_utf8(body).map_err(|e| e.to_string())?);
                    out.push_str(":]");
                } else {
                    // Collating elements and equivalence classes of one character.
                    for &b in body {
                        push_class_literal(out, b);
                    }
                }

                i = body_start + len + 2;
            }
            b'-' => {
                out.push('-');
                i += 1;
            }
            _ => {
                push_class_literal(out, c);
                i += 1;
            }
        }
    }
}

fn push_class_literal(out: &mut String, c: u8) {
    if c.is_ascii_alphanumeric() {
        out.push(char::from(c));
    } else if c.is_ascii() && !c.is_ascii_control() && c != b' ' {
        out.push('\\');
        out.push(char::from(c));
    } else {
        let _ = write!(out, "\\x{c:02X}");
    }
}

/// The record a filter is evaluated on.
struct Context<'r> {
    record: bam_raw::Record<'r>,
    header: &'r sam::Header,
    /// The CIGAR, read on first use (with a `CG` tag's operations, as HTSlib reads them).
    cigar: OnceCell<Result<Vec<u32>, String>>,
}

impl Context<'_> {
    fn cigar(&self) -> Result<&[u32], String> {
        self.cigar
            .get_or_init(|| {
                let mut ops = Vec::new();
                self.record
                    .read_cigar(&mut ops)
                    .map(|()| ops)
                    .map_err(|e| e.to_string())
            })
            .as_deref()
            .map_err(Clone::clone)
    }

    fn reference_name(&self, id: i32) -> &[u8] {
        usize::try_from(id)
            .ok()
            .and_then(|i| self.header.reference_sequences().get_index(i))
            .map_or(b"*".as_slice(), |(name, _)| &name[..])
    }

    /// `bam_get_library`: the `LB` of the record's read group, or the empty string.
    fn library(&self) -> Result<&[u8], String> {
        let Some(range) = self.record.find_field(*b"RG").map_err(|e| e.to_string())? else {
            return Ok(&[]);
        };

        let field = &self.record.as_bytes()[range];

        if field[2] != b'Z' {
            return Ok(&[]);
        }

        let id = c_str(&field[3..]);

        Ok(self
            .header
            .read_groups()
            .get::<[u8]>(id)
            .and_then(|read_group| read_group.other_fields().get(&LIBRARY))
            .map_or(&[][..], |lb| &lb[..]))
    }
}

/// `bam_cigar_type`: bit 1 consumes the query, bit 2 the reference.
fn cigar_type(op: u32) -> u32 {
    (0x3c1a7 >> ((op & 0xf) << 1)) & 3
}

fn eval(node: &Node, cx: &Context<'_>, res: &mut Value) -> Result<(), String> {
    match node {
        Node::Number(d) => res.set_number(*d),
        Node::String(s) => res.set_string(s),
        Node::Symbol(symbol) => lookup(*symbol, cx, res)?,
        Node::Function(function, arg) => eval_function(*function, arg, cx, res)?,
        Node::Default(a, b) => {
            eval(a, cx, res)?;
            let mut val = Value::default();
            eval(b, cx, &mut val)?;

            if !res.exists_t() {
                *res = val;
            }
        }
        Node::Pow(a, b) => {
            eval(a, cx, res)?;
            let mut val = Value::default();
            eval(b, cx, &mut val)?;

            if !res.exists() || !val.exists() {
                res.set_undefined();
            } else if res.is_str || val.is_str {
                return Err("arithmetic on strings".into());
            } else {
                res.set_number(res.d.powf(val.d));
            }

            if res.d.is_nan() {
                res.set_undefined();
            }
        }
        Node::Unary(op, operand) => eval_unary(*op, operand, cx, res)?,
        Node::Arith(first, rest) => eval_arith(first, rest, cx, res)?,
        Node::Bit(op, first, rest) => eval_bit(*op, first, rest, cx, res)?,
        Node::Compare(op, lhs, rhs) => {
            eval(lhs, cx, res)?;
            let mut val = Value::default();
            eval(rhs, cx, &mut val)?;

            if !res.exists() || !val.exists() {
                res.set_undefined();
            } else {
                let b = if res.is_str && res.s.is_some() && val.is_str && val.s.is_some() {
                    let ord = res.c_str().cmp(val.c_str());
                    match op {
                        CmpOp::Ge => ord != Ordering::Less,
                        CmpOp::Gt => ord == Ordering::Greater,
                        CmpOp::Le => ord != Ordering::Greater,
                        CmpOp::Lt => ord == Ordering::Less,
                    }
                } else {
                    !res.is_str
                        && !val.is_str
                        && match op {
                            CmpOp::Ge => res.d >= val.d,
                            CmpOp::Gt => res.d > val.d,
                            CmpOp::Le => res.d <= val.d,
                            CmpOp::Lt => res.d < val.d,
                        }
                };

                res.set_bool(b);
                res.is_str = false;
            }
        }
        Node::Equal(op, lhs, rhs) => {
            eval(lhs, cx, res)?;
            let mut val = Value::default();
            eval(rhs, cx, &mut val)?;

            if !res.exists() || !val.exists() {
                res.set_undefined();
            } else {
                let both_strings = res.s.is_some() && val.s.is_some();
                let b = match op {
                    EqOp::Eq if res.is_str => both_strings && res.c_str() == val.c_str(),
                    EqOp::Eq => !val.is_str && res.d == val.d,
                    EqOp::Ne if res.is_str => !both_strings || res.c_str() != val.c_str(),
                    EqOp::Ne => val.is_str || res.d != val.d,
                };
                res.set_bool(b);
            }

            res.is_str = false;

            if !val.exists() || !res.exists() {
                res.set_undefined();
            }
        }
        Node::Match {
            negated,
            lhs,
            rhs,
            regex,
        } => {
            eval(lhs, cx, res)?;
            let mut val = Value::default();
            eval(rhs, cx, &mut val)?;

            if !val.is_str || !res.is_str {
                return Err("regular expression matching needs strings".into());
            }

            if let (Some(pattern), Some(_)) = (&val.s, &res.s) {
                let compiled;
                let regex = match regex {
                    Some(regex) => regex,
                    None => {
                        compiled = compile_regex(pattern)?;
                        &compiled
                    }
                };

                let matched = regex.is_match(res.c_str());
                res.set_bool(matched != *negated);
            } else {
                // A missing pattern or string is false.
                res.is_true = false;
            }

            res.is_str = false;

            if !val.exists() || !res.exists() {
                res.set_undefined();
            }
        }
        Node::Logic(first, rest) => {
            eval(first, cx, res)?;

            for (op, node) in rest {
                let mut val = Value::default();
                eval(node, cx, &mut val)?;

                let (res_t, val_t) = (res.exists_t(), val.exists_t());

                let undefined = match op {
                    LogicOp::And => !res_t || !val_t,
                    LogicOp::Or => {
                        (!res_t && !val_t) || (!res_t && !val.truth()) || (!val_t && !res.truth())
                    }
                };

                if undefined {
                    res.set_undefined();
                    res.d = 0.0;
                } else {
                    let b = match op {
                        LogicOp::And => res.truth() && val.truth(),
                        LogicOp::Or => res.truth() || val.truth(),
                    };
                    res.set_bool(b);
                    res.is_str = false;
                }
            }
        }
    }

    Ok(())
}

fn eval_function(
    function: Function,
    arg: &Node,
    cx: &Context<'_>,
    res: &mut Value,
) -> Result<(), String> {
    eval(arg, cx, res)?;

    match function {
        Function::Length | Function::Min | Function::Max | Function::Average => {
            if !res.is_str {
                return Err("function needs a string".into());
            }

            let s: &[u8] = res.s.as_deref().unwrap_or_default();

            let d = match function {
                Function::Length => s.len() as f64,
                Function::Min => s.iter().min().map_or(f64::NAN, |&b| f64::from(b)),
                Function::Max => s.iter().max().map_or(f64::NAN, |&b| f64::from(b)),
                _ if s.is_empty() => 0.0,
                _ => s.iter().map(|&b| f64::from(b)).sum::<f64>() / s.len() as f64,
            };

            res.set_number(d);
        }
        Function::Exists => {
            let b = res.exists_t();
            res.set_bool(b);
            res.is_str = false;
        }
        Function::Exp | Function::Log | Function::Sqrt => {
            let d = match function {
                Function::Exp => res.d.exp(),
                Function::Log => res.d.ln(),
                _ => res.d.sqrt(),
            };

            res.set_number(d);

            if d.is_nan() {
                res.set_undefined();
            }
        }
    }

    Ok(())
}

fn eval_unary(
    op: UnaryOp,
    operand: &Node,
    cx: &Context<'_>,
    res: &mut Value,
) -> Result<(), String> {
    eval(operand, cx, res)?;

    match op {
        UnaryOp::Plus | UnaryOp::Minus => {
            if !res.exists() {
                res.set_undefined();
            } else if res.is_str {
                return Err("arithmetic on strings".into());
            } else {
                if op == UnaryOp::Minus {
                    res.d = -res.d;
                }
                res.is_true = res.d != 0.0;
            }
        }
        UnaryOp::Not => {
            if res.is_true {
                res.set_bool(false);
            } else if !res.exists() {
                // Undefined values negate to true.
                res.set_bool(true);
            } else if res.is_str {
                // A defined string, even empty, is true, so its negation is false.
                res.set_bool(false);
            } else {
                let b = res.d as i64 == 0;
                res.set_bool(b);
            }

            res.is_str = false;
        }
        UnaryOp::BitNot => {
            if !res.exists() {
                res.set_undefined();
            } else if res.is_str {
                return Err("bitwise operation on a string".into());
            } else {
                let d = !(res.d as i64) as f64;
                res.d = d;
                res.is_true = d != 0.0;
            }
        }
    }

    Ok(())
}

fn eval_arith(
    first: &Node,
    rest: &[(ArithOp, Node)],
    cx: &Context<'_>,
    res: &mut Value,
) -> Result<(), String> {
    eval(first, cx, res)?;

    // HTSlib reuses one value for the right operands, so fields it does not set carry over.
    let mut val = Value::default();

    for (op, node) in rest {
        eval(node, cx, &mut val)?;

        let undefined = !val.exists() || !res.exists();

        if !undefined && (val.is_str || res.is_str) {
            return Err("arithmetic on strings".into());
        }

        match op {
            ArithOp::Mul | ArithOp::Div => {
                if undefined {
                    res.set_undefined();
                }

                if *op == ArithOp::Mul {
                    res.d *= val.d;
                } else {
                    res.d /= val.d;
                }

                res.is_true = res.exists() && res.d != 0.0;
            }
            ArithOp::Rem => {
                // HTSlib computes an undefined left operand as a number here (giving 0 rather
                // than undefined) and crashes on divisors in (-1, 1); both stay undefined here.
                let divisor = val.d as i64;

                if undefined || divisor == 0 {
                    res.set_undefined();
                } else {
                    res.d = (res.d as i64).wrapping_rem(divisor) as f64;
                    res.is_true = res.d != 0.0;
                }
            }
            ArithOp::Add | ArithOp::Sub => {
                if *op == ArithOp::Add {
                    res.d += val.d;
                } else {
                    res.d -= val.d;
                }

                if undefined {
                    res.set_undefined();
                } else {
                    res.is_true = res.d != 0.0;
                }
            }
        }
    }

    Ok(())
}

fn eval_bit(
    op: BitOp,
    first: &Node,
    rest: &[Node],
    cx: &Context<'_>,
    res: &mut Value,
) -> Result<(), String> {
    eval(first, cx, res)?;

    let mut val = Value::default();
    let mut undefined = false;

    for node in rest {
        eval(node, cx, &mut val)?;

        if !val.exists() || !res.exists() {
            undefined = true;
        } else if res.is_str || val.is_str {
            return Err("bitwise operation on a string".into());
        } else {
            let (a, b) = (res.d as i64, val.d as i64);
            let d = match op {
                BitOp::And => a & b,
                BitOp::Xor => a ^ b,
                BitOp::Or => a | b,
            } as f64;
            res.d = d;
            res.is_true = d != 0.0;
        }
    }

    if undefined {
        res.set_undefined();
    }

    Ok(())
}

/// `bam_sym_lookup`: numeric fields set the number, string fields the string, and neither
/// touches the other fields of the value.
fn lookup(symbol: Symbol, cx: &Context<'_>, res: &mut Value) -> Result<(), String> {
    let record = cx.record;
    res.is_str = false;

    let number = |n: i64| n as f64;

    match symbol {
        Symbol::Cigar => {
            let ops = cx.cigar()?;
            let mut s = String::new();

            if ops.is_empty() {
                s.push('*');
            }

            for &op in ops {
                let _ = write!(
                    s,
                    "{}{}",
                    op >> 4,
                    char::from(CIGAR_OPS[(op & 0xf) as usize])
                );
            }

            res.set_string(s.as_bytes());
        }
        Symbol::EndPos => {
            let rlen = if record.flags() & FLAG_UNMAPPED != 0 {
                0
            } else {
                reference_length(cx.cigar()?)
            };

            res.set_number(number(i64::from(record.position()) + rlen.max(1)));
        }
        Symbol::Flag(mask) => {
            let flags = record.flags();
            res.set_number(f64::from(mask.map_or(flags, |bit| flags & bit)));
        }
        Symbol::HardClipLength => {
            let ops = cx.cigar()?;
            let mut n = 0;

            if let Some(&first) = ops.first()
                && first & 0xf == CIGAR_HARD_CLIP
            {
                n += first >> 4;
            }

            if ops.len() > 1
                && let Some(&last) = ops.last()
                && last & 0xf == CIGAR_HARD_CLIP
            {
                n += last >> 4;
            }

            res.set_number(f64::from(n));
        }
        Symbol::SoftClipLength => {
            res.set_number(f64::from(soft_clip_length(cx.cigar()?)));
        }
        Symbol::Library => {
            let library = cx.library()?;
            res.set_string(library);
        }
        Symbol::MappingQuality => res.set_number(f64::from(record.mapping_quality())),
        Symbol::MatePosition => {
            res.set_number(number(i64::from(record.mate_position()) + 1));
        }
        Symbol::MateReferenceName => {
            let name = cx.reference_name(record.mate_reference_sequence_id());
            res.set_string(name);
        }
        Symbol::MateReferenceId => {
            res.set_number(f64::from(record.mate_reference_sequence_id()));
        }
        Symbol::CigarOpCount => res.set_number(number(cx.cigar()?.len() as i64)),
        Symbol::Position => res.set_number(number(i64::from(record.position()) + 1)),
        Symbol::QueryLength => {
            let n: i64 = cx
                .cigar()?
                .iter()
                .filter(|&&op| cigar_type(op) & 1 != 0)
                .map(|&op| i64::from(op >> 4))
                .sum();
            res.set_number(number(n));
        }
        Symbol::Name => res.set_string(record.name()),
        Symbol::Quality => res.set_string(record.quality_scores()),
        Symbol::ReferenceLength => res.set_number(number(reference_length(cx.cigar()?))),
        Symbol::ReferenceName => {
            let name = cx.reference_name(record.reference_sequence_id());
            res.set_string(name);
        }
        Symbol::ReferenceId => res.set_number(f64::from(record.reference_sequence_id())),
        Symbol::Sequence => {
            let mut s = res.s.take().unwrap_or_default();
            s.clear();
            record.decode_sequence(&mut s);
            res.s = Some(s);
            res.is_str = true;
        }
        Symbol::TemplateLength => res.set_number(f64::from(record.template_length())),
        Symbol::Tag(tag) => lookup_tag(tag, cx, res)?,
    }

    Ok(())
}

fn reference_length(ops: &[u32]) -> i64 {
    ops.iter()
        .filter(|&&op| cigar_type(op) & 2 != 0)
        .map(|&op| i64::from(op >> 4))
        .sum()
}

/// `sclen`: soft clips at either end, inside any hard clips.
fn soft_clip_length(ops: &[u32]) -> u32 {
    let op = |i: usize| ops[i] & 0xf;
    let len = |i: usize| ops[i] >> 4;
    let n = ops.len();

    let mut sclen = 0;
    let mut left = 0;

    if n > 0 && op(0) == CIGAR_SOFT_CLIP {
        sclen += len(0);
    } else if n > 1 && op(0) == CIGAR_HARD_CLIP && op(1) == CIGAR_SOFT_CLIP {
        left = 1;
        sclen += len(1);
    }

    if n > left + 1 && op(n - 1) == CIGAR_SOFT_CLIP {
        sclen += len(n - 1);
    } else if n > left + 2 && op(n - 1) == CIGAR_HARD_CLIP && op(n - 2) == CIGAR_SOFT_CLIP {
        sclen += len(n - 2);
    }

    sclen
}

/// `[XX]`: a present tag is true even if zero; an absent one is an undefined string.
fn lookup_tag(tag: [u8; 2], cx: &Context<'_>, res: &mut Value) -> Result<(), String> {
    let record = cx.record;

    let mut range = record.find_field(tag).map_err(|e| e.to_string())?;

    // HTSlib moves a long CIGAR out of its `CG` tag when it reads the record.
    if tag == *b"CG"
        && record
            .long_cigar_field()
            .map_err(|e| e.to_string())?
            .is_some()
    {
        range = None;
    }

    let Some(range) = range else {
        res.is_str = true;
        res.clear_string();
        res.d = 0.0;
        res.is_true = false;
        return Ok(());
    };

    let field = &record.as_bytes()[range];
    let (ty, value) = (field[2], &field[3..]);
    res.is_true = true;

    let int = |n: i64| n as f64;

    match ty {
        b'Z' | b'H' => res.set_string(c_str(value)),
        b'A' => res.set_string(&value[..1]),
        b'c' => res.set_number(f64::from(value[0].cast_signed())),
        b'C' => res.set_number(f64::from(value[0])),
        b's' => res.set_number(f64::from(i16::from_le_bytes([value[0], value[1]]))),
        b'S' => res.set_number(f64::from(u16::from_le_bytes([value[0], value[1]]))),
        b'i' => res.set_number(f64::from(i32::from_le_bytes(word(value)))),
        b'I' => res.set_number(int(i64::from(u32::from_le_bytes(word(value))))),
        b'f' => res.set_number(f64::from(f32::from_le_bytes(word(value)))),
        b'd' => {
            let mut b = [0; 8];
            b.copy_from_slice(&value[..8]);
            res.set_number(f64::from_le_bytes(b));
        }
        other => {
            return Err(format!(
                "aux type '{}' not yet supported by filters",
                char::from(other)
            ));
        }
    }

    Ok(())
}

fn word(value: &[u8]) -> [u8; 4] {
    [value[0], value[1], value[2], value[3]]
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::alignment::io::Write as _;

    fn header() -> sam::Header {
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\n@RG\tID:rg1\tLB:lib1\n@RG\tID:rg2\n"
            .parse()
            .expect("valid header")
    }

    fn encode(header: &sam::Header, line: &str) -> Vec<u8> {
        let record: sam::Record = sam::Record::try_from(line.as_bytes()).expect("valid record");
        let mut writer = noodles_bam::io::Writer::from(Vec::new());
        writer
            .write_alignment_record(header, &record)
            .expect("encodable record");
        writer.get_ref()[4..].to_vec()
    }

    const READ: &str = "r1\t99\tchr1\t100\t30\t2H3S10M2D5M4S\t=\t300\t250\tACGTACGTACGTACGTACGTAC\tIIIIIIIIIIIIIIIIIIIII#\tRG:Z:rg1\tNM:i:3\tXA:A:q\tXF:f:0.5\tZZ:i:0\tBA:B:c,1,2";
    const UNMAPPED: &str = "u1\t4\t*\t0\t0\t*\t*\t0\t0\tAC\t*";

    fn check(expr: &str, line: &str) -> anyhow::Result<bool> {
        let header = header();
        let raw = encode(&header, line);
        let record = bam_raw::Record::new(&raw)?;
        Expression::parse(expr)?.matches(record, &header)
    }

    fn holds(expr: &str) -> bool {
        check(expr, READ).unwrap_or_else(|e| panic!("{expr}: {e}"))
    }

    #[test]
    fn test_fields() {
        assert!(holds(
            "mapq == 30 && pos == 100 && mpos == 300 && pnext == 300"
        ));
        assert!(holds(
            "flag == 99 && flag.paired && flag.read1 && !flag.read2"
        ));
        assert!(holds("flag.proper_pair == 2 && flag.mreverse == 32"));
        assert!(holds(
            "rname == \"chr1\" && rnext == \"chr1\" && mrname == \"chr1\""
        ));
        assert!(holds("refid == 0 && mrefid == 0 && tlen == 250"));
        assert!(holds("qname == \"r1\" && library == \"lib1\""));
        assert!(holds("cigar == \"2H3S10M2D5M4S\" && ncigar == 6"));
        assert!(holds("qlen == 22 && rlen == 17 && endpos == 116"));
        assert!(holds("sclen == 7 && hclen == 2"));
        assert!(holds("seq =~ \"^ACGT\" && length(seq) == 22"));
        assert!(holds("min(qual) == 2 && max(qual) == 40"));
    }

    #[test]
    fn test_tags() {
        assert!(holds("[NM] == 3 && [XA] == \"q\" && [XF] == 0.5"));
        // A present tag is true even when zero; an absent one is false.
        assert!(holds("[ZZ]"));
        assert!(!holds("[YY]"));
        assert!(holds("![YY]"));
        assert!(holds("exists([ZZ]) && !exists([YY])"));
        assert!(holds("default([YY], 7) == 7 && default([NM], 7) == 3"));
        assert!(check("[BA] == 1", READ).is_err());
        assert!(check("[NM] + qname", READ).is_err());
    }

    #[test]
    fn test_unmapped() -> anyhow::Result<()> {
        assert!(check(
            "rname == \"*\" && pos == 0 && endpos == 0 && cigar == \"*\" && library == \"\"",
            UNMAPPED
        )?);
        // Missing qualities are 0xff bytes.
        assert!(check("min(qual) == 255", UNMAPPED)?);
        Ok(())
    }

    #[test]
    fn test_operators() {
        assert!(holds("1 + 2 * 3 == 7"));
        assert!(holds("(1 + 2) * 3 == 9"));
        assert!(holds("7 % 3 == 1 && -7 % 3 == -1"));
        assert!(holds(
            "(5 & 3) == 1 && (5 | 3) == 7 && (5 ^ 3) == 6 && ~0 == -1"
        ));
        // Bitwise operators bind tighter than comparisons, unlike C.
        assert!(holds("flag & 64 == 64"));
        assert!(holds("2 > 1 && 1 >= 1 && 1 < 2 && 2 <= 2"));
        assert!(holds("\"abc\" < \"abd\" && \"b\" > \"abc\""));
        assert!(holds("1 != \"1\" && !(1 == \"1\")"));
        assert!(holds("0 || 1"));
        assert!(!holds("1 && 0"));
        // && and || share a precedence level, left to right.
        assert!(!holds("1 || 1 && 0"));
        assert!(holds(
            "sqrt(16) == 4 && pow(2, 10) == 1024 && exp(0) == 1 && log(1) == 0"
        ));
        assert!(holds("avg(\"AC\") == 66"));
        assert!(holds("0x10 == 16 && 1e2 == 100 && .5 == 0.5"));
        assert!(holds("qname !~ \"^x\" && qname =~ \"r[0-9]\""));
        assert!(holds("\"a\\\"b\" == \"a\\\"b\" && length(\"a\\tb\") == 3"));
    }

    #[test]
    fn test_undefined() {
        // Undefined values propagate through arithmetic and comparisons.
        assert!(!holds("[YY] + 1 == 1"));
        assert!(!holds("[YY] % 2 == 0"));
        assert!(!holds("[YY] > 0"));
        assert!(!holds("[YY] || 0"));
        assert!(holds("[YY] || 1"));
        assert!(!holds("[YY] =~ \"x\""));
        assert!(!holds("[YY] !~ \"x\""));
        // A missing function value (min of an empty string) is undefined.
        assert!(!holds("min(\"\") == 0"));
        assert!(holds("length([YY]) == 0"));
    }

    #[test]
    fn test_parse_errors() {
        for bad in [
            "mapq >",
            "(mapq",
            "foo",
            "flag.bogus",
            "posx",
            "\"abc",
            "default(1 ,2)",
            "length (seq)",
            "mapq = 1",
            "qname =~ \"(\"",
        ] {
            assert!(Expression::parse(bad).is_err(), "{bad}");
        }

        assert!(Expression::parse("default(1, 2)").is_ok());
    }

    #[test]
    fn test_strtod() {
        assert_eq!(strtod(b"12.5x", 0), Some((12.5, 4)));
        assert_eq!(strtod(b"0x1p4", 0), Some((16.0, 5)));
        assert_eq!(strtod(b"0x", 0), Some((0.0, 1)));
        assert_eq!(strtod(b"1e", 0), Some((1.0, 1)));
        assert_eq!(strtod(b"-inf", 0), Some((f64::NEG_INFINITY, 4)));
        assert_eq!(strtod(b"abc", 0), None);
        assert_eq!(strtod(b".", 0), None);
    }

    #[test]
    fn test_translate_ere() -> Result<(), String> {
        assert_eq!(translate_ere(b"a.c")?, "a.c");
        assert_eq!(translate_ere(br"[\.]")?, r"[\\\.]");
        assert_eq!(translate_ere(b"[]a]")?, r"[\]a]");
        assert_eq!(translate_ere(b"[[:digit:]]+")?, "[[:digit:]]+");
        assert_eq!(translate_ere(br"\d")?, "d");
        assert!(translate_ere(br"(a)\1").is_err());
        Ok(())
    }
}
