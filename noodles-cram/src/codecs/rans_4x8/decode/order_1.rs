use std::io;

use super::{
    SLOT_COUNT, order_0, read_states, renormalize, state_cumulative_frequency, state_renormalize,
    state_step,
};
use crate::{
    codecs::rans_4x8::{ALPHABET_SIZE, STATE_COUNT},
    io::reader::num::read_u8,
};

type Frequencies = [[u16; ALPHABET_SIZE]; ALPHABET_SIZE]; // F
type CumulativeFrequencies = Frequencies; // C
type CumulativeFrequenciesSymbolsTable = [[u8; 4096]; ALPHABET_SIZE];

pub fn decode(src: &mut &[u8], dst: &mut [u8]) -> io::Result<()> {
    let frequencies = read_frequencies(src)?;
    let cumulative_frequencies = build_cumulative_frequencies(&frequencies);

    match SymbolTables::new(&frequencies, &cumulative_frequencies) {
        Some(tables) => decode_with_context_tables(src, dst, &tables),
        None => decode_with_symbols_table(src, dst, &frequencies, &cumulative_frequencies),
    }
}

/// Decodes with the cumulative frequencies symbols tables of all 256 contexts (the general
/// method).
fn decode_with_symbols_table(
    src: &mut &[u8],
    dst: &mut [u8],
    frequencies: &Frequencies,
    cumulative_frequencies: &CumulativeFrequencies,
) -> io::Result<()> {
    let cumulative_frequencies_symbols_table =
        build_cumulative_frequencies_symbols_table(cumulative_frequencies);

    let mut states = read_states(src)?;
    let mut prev_syms = [0; STATE_COUNT];

    let [chunk_0, chunk_1, chunk_2, chunk_3, chunk_4] = split_chunks(dst);
    let chunks = chunk_0.iter_mut().zip(chunk_1).zip(chunk_2).zip(chunk_3);

    for (((d0, d1), d2), d3) in chunks {
        let dsts = [d0, d1, d2, d3];

        for (state, (prev_sym, d)) in states.iter_mut().zip(prev_syms.iter_mut().zip(dsts)) {
            let i = usize::from(*prev_sym);
            let f = state_cumulative_frequency(*state);
            let sym = cumulative_frequencies_symbols_table[i][usize::from(f)];

            *d = sym;

            let j = usize::from(sym);
            *state = state_step(*state, frequencies[i][j], cumulative_frequencies[i][j]);
            *state = state_renormalize(*state, src)?;

            *prev_sym = sym;
        }
    }

    let mut state = states[3];
    let mut prev_sym = prev_syms[3];

    for d in chunk_4 {
        let i = usize::from(prev_sym);
        let f = state_cumulative_frequency(state);
        let sym = cumulative_frequencies_symbols_table[i][usize::from(f)];

        *d = sym;

        let j = usize::from(sym);
        state = state_step(state, frequencies[i][j], cumulative_frequencies[i][j]);
        state = state_renormalize(state, src)?;

        prev_sym = sym;
    }

    Ok(())
}

/// Symbol tables of the contexts that occur: 4 KiB per context, so that the tables of all
/// contexts in use stay in cache.
struct SymbolTables<'f> {
    /// The start of each context's table in `symbols`, or `None` if the context does not occur.
    offsets: [Option<usize>; ALPHABET_SIZE],
    symbols: Vec<u8>,
    frequencies: &'f Frequencies,
    cumulative_frequencies: &'f CumulativeFrequencies,
}

impl<'f> SymbolTables<'f> {
    fn new(
        frequencies: &'f Frequencies,
        cumulative_frequencies: &'f CumulativeFrequencies,
    ) -> Option<Self> {
        let n = frequencies
            .iter()
            .filter(|f| f.iter().any(|&g| g > 0))
            .count();
        let mut symbols = vec![0; n * SLOT_COUNT];
        let mut offsets = [None; ALPHABET_SIZE];
        let mut next = 0;

        for ((f, g), offset) in frequencies
            .iter()
            .zip(cumulative_frequencies)
            .zip(&mut offsets)
        {
            if f.iter().all(|&x| x == 0) {
                continue;
            }

            let total: usize = f.iter().map(|&x| usize::from(x)).sum();

            if total > SLOT_COUNT {
                return None;
            }

            let table = &mut symbols[next..next + SLOT_COUNT];

            for (sym, (&x, &y)) in f.iter().zip(g).enumerate() {
                let start = usize::from(y);
                table[start..start + usize::from(x)].fill(sym as u8);
            }

            table[total..].fill(u8::MAX);
            *offset = Some(next);
            next += SLOT_COUNT;
        }

        Some(Self {
            offsets,
            symbols,
            frequencies,
            cumulative_frequencies,
        })
    }

    /// Decodes one symbol in the context of the previous symbol, as the general method does;
    /// returns the symbol and the state before renormalization.
    #[inline(always)]
    fn decode(&self, s: u32, prev_sym: u8) -> (u8, u32) {
        let i = usize::from(prev_sym);
        let sym = match self.offsets[i] {
            Some(start) => self.symbols[start + usize::from(state_cumulative_frequency(s))],
            // A context without frequencies: the general method maps every slot to 255.
            None => u8::MAX,
        };
        let j = usize::from(sym);
        (
            sym,
            state_step(s, self.frequencies[i][j], self.cumulative_frequencies[i][j]),
        )
    }
}

/// Decodes with the symbol tables of the contexts that occur: the same result as
/// [`decode_with_symbols_table`], without building tables for unused contexts.
fn decode_with_context_tables(
    src: &mut &[u8],
    dst: &mut [u8],
    tables: &SymbolTables<'_>,
) -> io::Result<()> {
    let mut states = read_states(src)?;
    let mut prev_syms = [0; STATE_COUNT];

    let [chunk_0, chunk_1, chunk_2, chunk_3, chunk_4] = split_chunks(dst);
    let chunks = chunk_0.iter_mut().zip(chunk_1).zip(chunk_2).zip(chunk_3);

    for (((d0, d1), d2), d3) in chunks {
        let dsts = [d0, d1, d2, d3];

        for (state, (prev_sym, d)) in states.iter_mut().zip(prev_syms.iter_mut().zip(dsts)) {
            let (sym, s) = tables.decode(*state, *prev_sym);
            *d = sym;
            *state = renormalize(s, src)?;
            *prev_sym = sym;
        }
    }

    let mut state = states[3];
    let mut prev_sym = prev_syms[3];

    for d in chunk_4 {
        let (sym, s) = tables.decode(state, prev_sym);
        *d = sym;
        state = renormalize(s, src)?;
        prev_sym = sym;
    }

    Ok(())
}

fn read_frequencies(src: &mut &[u8]) -> io::Result<Frequencies> {
    let mut frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    let mut sym = read_u8(src)?;
    let mut prev_sym = sym;

    loop {
        let f = order_0::read_frequencies(src)?;
        frequencies[usize::from(sym)] = f;

        sym = read_u8(src)?;

        if sym == 0 {
            break;
        }

        if sym - 1 == prev_sym {
            let len = read_u8(src)?;

            for _ in 0..len {
                let f = order_0::read_frequencies(src)?;
                frequencies[usize::from(sym)] = f;
                sym += 1;
            }
        }

        prev_sym = sym;
    }

    Ok(frequencies)
}

fn build_cumulative_frequencies(frequencies: &Frequencies) -> CumulativeFrequencies {
    let mut cumulative_frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    for (f, g) in frequencies.iter().zip(&mut cumulative_frequencies) {
        *g = order_0::build_cumulative_frequencies(f);
    }

    cumulative_frequencies
}

pub fn build_cumulative_frequencies_symbols_table(
    cumulative_freqs: &CumulativeFrequencies,
) -> Box<CumulativeFrequenciesSymbolsTable> {
    let mut tables = Box::new([[0; 4096]; 256]);

    for (table, cumulative_freqs) in tables.iter_mut().zip(cumulative_freqs) {
        *table = order_0::build_cumulative_frequencies_symbols_table(cumulative_freqs);
    }

    tables
}

fn split_chunks(dst: &mut [u8]) -> [&mut [u8]; 5] {
    let chunk_size = dst.len() / STATE_COUNT;

    let (left_chunk, right_chunk) = dst.split_at_mut(2 * chunk_size);
    let (chunk_0, chunk_1) = left_chunk.split_at_mut(chunk_size);
    let (chunk_2, chunk_3_4) = right_chunk.split_at_mut(chunk_size);
    let (chunk_3, chunk_4) = chunk_3_4.split_at_mut(chunk_size);

    [chunk_0, chunk_1, chunk_2, chunk_3, chunk_4]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::codecs::rans_4x8::decode::test_support::{Rng, random_frequencies};

    #[test]
    fn test_decode_with_context_tables_matches_symbols_table() {
        let mut rng = Rng(0x2545_f491_4f6c_dd1d);

        for _ in 0..300 {
            // A few contexts with frequencies (some summing to less than 4095), others unused.
            let mut frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];

            for _ in 0..1 + rng.below(6) {
                let ctx = if rng.below(3) == 0 { 0 } else { rng.below(256) } as usize;
                frequencies[ctx] = random_frequencies(&mut rng, 4096);
            }

            let cumulative_frequencies = build_cumulative_frequencies(&frequencies);
            let tables =
                SymbolTables::new(&frequencies, &cumulative_frequencies).expect("valid tables");

            let src: Vec<u8> = (0..16 + rng.below(96)).map(|_| rng.next() as u8).collect();
            let len = rng.below(64) as usize;
            let (mut expected, mut actual) = (vec![0; len], vec![0; len]);

            let expected_result = decode_with_symbols_table(
                &mut &src[..],
                &mut expected,
                &frequencies,
                &cumulative_frequencies,
            );
            let actual_result = decode_with_context_tables(&mut &src[..], &mut actual, &tables);

            match (expected_result, actual_result) {
                (Ok(()), Ok(())) => assert_eq!(actual, expected),
                (Err(e), Err(a)) => assert_eq!(a.kind(), e.kind()),
                (e, a) => panic!("results differ: {e:?} != {a:?}"),
            }
        }
    }

    #[test]
    fn test_read_frequencies() -> io::Result<()> {
        const NUL: u8 = 0x00;

        // § 2.1.2 "Frequency table: Order-1 encoding" (563e8ab 2025-04-07):
        // "abracadabraabracadabraabracadabraabracadabrad".
        let src = [
            0x00, // symbols[0] = '\0' {
            b'a', //   symbols[1] = 'a'
            0x04, //     frequencies['\0']['a'] = 4
            0x00, // }
            b'a', // symbols[0] = 'a' {
            b'a', //   symbols[1] = 'a'
            0x03, //     frequencies['a']['a'] = 3
            b'b', //   symbols[1] = 'b'
            0x02, //     run length = 2
            0x08, //     frequencies['a']['b'] = 8
            0x04, //     frequencies['a']['c'] = 4
            0x05, //     frequencies['a']['d'] = 5
            0x00, // }
            b'b', // symbols[0] = 'b' {
            0x02, //   run length = 2
            b'r', //   symbols[1] = 'r'
            0x08, //     frequencies['b']['r'] = 8
            0x00, // }
            //    // symbols[0] = 'c' {
            b'a', //   symbols[1] = 'a'
            0x04, //     frequencies['c']['a'] = 4
            0x00, // }
            //    // symbols[0] = 'd' {
            b'a', //   symbols[1] = 'a'
            0x04, //     frequencies['d']['a'] = 4
            0x00, // }
            b'r', // symbols[0] = 'r' {
            b'a', //   symbols[1] = 'a'
            0x08, //     frequencies['r']['a'] = 8
            0x00, // }
            0x00, // EOF
        ];

        let mut expected = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];
        expected[usize::from(NUL)][usize::from(b'a')] = 4;
        expected[usize::from(b'a')][usize::from(b'a')] = 3;
        expected[usize::from(b'a')][usize::from(b'b')] = 8;
        expected[usize::from(b'a')][usize::from(b'c')] = 4;
        expected[usize::from(b'a')][usize::from(b'd')] = 5;
        expected[usize::from(b'b')][usize::from(b'r')] = 8;
        expected[usize::from(b'c')][usize::from(b'a')] = 4;
        expected[usize::from(b'd')][usize::from(b'a')] = 4;
        expected[usize::from(b'r')][usize::from(b'a')] = 8;

        assert_eq!(read_frequencies(&mut &src[..])?, expected);

        Ok(())
    }
}
