mod header;
mod order_0;
mod order_1;

use std::io;

use self::header::read_header;
use super::{ALPHABET_SIZE, LOWER_BOUND, Order, STATE_COUNT};
use crate::io::reader::num::{read_u8, read_u32_le};

pub fn decode(mut src: &[u8]) -> io::Result<Vec<u8>> {
    let (order, _, uncompressed_size) = read_header(&mut src)?;

    let mut dst = vec![0; uncompressed_size];

    match order {
        Order::Zero => order_0::decode(&mut src, &mut dst)?,
        Order::One => order_1::decode(&mut src, &mut dst)?,
    }

    Ok(dst)
}

fn read_states(src: &mut &[u8]) -> io::Result<[u32; STATE_COUNT]> {
    let mut states = [0; STATE_COUNT];

    for state in &mut states {
        *state = read_u32_le(src)?;
    }

    Ok(states)
}

pub fn state_cumulative_frequency(s: u32) -> u16 {
    (s & 0x0fff) as u16
}

pub fn state_step(s: u32, f: u16, g: u16) -> u32 {
    u32::from(f) * (s >> 12) + (s & 0x0fff) - u32::from(g)
}

pub fn state_renormalize(mut s: u32, src: &mut &[u8]) -> io::Result<u32> {
    while s < LOWER_BOUND {
        let b = read_u8(src).map(u32::from)?;
        s = (s << 8) | b;
    }

    Ok(s)
}

/// The number of cumulative frequency slots (2^12).
const SLOT_COUNT: usize = 4096;

/// Builds the slot lookup of a frequency table, mapping each slot to what the general method
/// (the cumulative frequencies symbols table and [`state_step`]) decodes there.
///
/// Each slot maps to `sym | f << 8 | (slot - g) << 20`, where `sym` is the slot's symbol, `f` its
/// frequency and `g` its cumulative frequency. The next state is then `f * (s >> 12) + (slot -
/// g)`, which equals [`state_step`]. Slots past the total (frequencies are normalized to 4095)
/// map to symbol 255, as in the general method.
///
/// Returns `false` (and leaves `dst` unspecified) if a frequency does not fit in 12 bits or the
/// frequencies sum to more than [`SLOT_COUNT`]; such tables are decoded with the general method.
fn build_slot_lookup(
    frequencies: &[u16; ALPHABET_SIZE],
    cumulative_frequencies: &[u16; ALPHABET_SIZE],
    dst: &mut [u32],
) -> bool {
    const MAX_FREQUENCY: u16 = 0x0fff;

    let total: usize = frequencies.iter().map(|&f| usize::from(f)).sum();

    if total > SLOT_COUNT || frequencies.iter().any(|&f| f > MAX_FREQUENCY) {
        return false;
    }

    for (sym, (&f, &g)) in frequencies.iter().zip(cumulative_frequencies).enumerate() {
        let (f, g) = (u32::from(f), usize::from(g));

        for (offset, slot) in (0..f).zip(&mut dst[g..]) {
            *slot = sym as u32 | (f << 8) | (offset << 20);
        }
    }

    let last = ALPHABET_SIZE - 1;
    let (f, g) = (
        u32::from(frequencies[last]),
        usize::from(cumulative_frequencies[last]),
    );

    for (i, slot) in dst.iter_mut().enumerate().take(SLOT_COUNT).skip(total) {
        *slot = last as u32 | (f << 8) | (((i - g) as u32) << 20);
    }

    true
}

/// Decodes one symbol with a slot lookup (see [`build_slot_lookup`]); returns the symbol and
/// the state before renormalization.
#[inline(always)]
fn decode_slot(s: u32, lookup: &[u32]) -> (u8, u32) {
    let entry = lookup[(s & 0x0fff) as usize];
    let f = (entry >> 8) & 0x0fff;
    (entry as u8, f * (s >> 12) + (entry >> 20))
}

/// [`state_renormalize`] for the slot lookup path. A state needs at most two bytes unless its
/// symbol's frequency is 0, so with two bytes left this reads without branching on the data.
#[inline(always)]
fn renormalize(s: u32, src: &mut &[u8]) -> io::Result<u32> {
    if let [b0, b1, ..] = **src {
        let n0 = u32::from(s < LOWER_BOUND);
        let s = (s << (8 * n0)) | (u32::from(b0) & n0.wrapping_neg());
        let n1 = u32::from(s < LOWER_BOUND);
        let b = if n0 == 1 { b1 } else { b0 };
        let s = (s << (8 * n1)) | (u32::from(b) & n1.wrapping_neg());
        *src = &src[(n0 + n1) as usize..];

        if s >= LOWER_BOUND {
            return Ok(s);
        }

        return renormalize_bytewise(s, src);
    }

    renormalize_bytewise(s, src)
}

fn renormalize_bytewise(mut s: u32, src: &mut &[u8]) -> io::Result<u32> {
    while s < LOWER_BOUND {
        let Some((&b, rest)) = src.split_first() else {
            return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
        };

        *src = rest;
        s = (s << 8) | u32::from(b);
    }

    Ok(s)
}

/// Helpers for randomized decoder tests.
#[cfg(test)]
pub(super) mod test_support {
    use super::ALPHABET_SIZE;

    /// A deterministic xorshift generator for randomized tests.
    pub(super) struct Rng(pub(super) u64);

    impl Rng {
        pub(super) fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }

        pub(super) fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    /// Random frequencies over a few symbols, summing to at most `max_total`.
    pub(super) fn random_frequencies(rng: &mut Rng, max_total: u64) -> [u16; ALPHABET_SIZE] {
        let mut frequencies = [0; ALPHABET_SIZE];
        let mut left = 1 + rng.below(max_total);

        for _ in 0..1 + rng.below(8) {
            let sym = if rng.below(4) == 0 {
                255
            } else {
                rng.below(256)
            } as usize;
            let f = 1 + rng.below(left.min(4095));
            frequencies[sym] = f as u16;
            left -= f;

            if left == 0 {
                break;
            }
        }

        frequencies
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_states() -> io::Result<()> {
        let src = [
            0x05, 0x00, 0x00, 0x00, // states[0] = 5
            0x08, 0x00, 0x00, 0x00, // states[1] = 8
            0x0d, 0x00, 0x00, 0x00, // states[2] = 13
            0x15, 0x00, 0x00, 0x00, // states[3] = 21
        ];

        assert_eq!(read_states(&mut &src[..])?, [5, 8, 13, 21]);

        Ok(())
    }

    #[test]
    fn test_decode_with_order_0() -> io::Result<()> {
        let expected = b"noodles";

        let src = [
            0x00, 0x25, 0x00, 0x00, 0x00, 0x07, 0x00, 0x00, 0x00, 0x64, 0x82, 0x49, 0x65, 0x00,
            0x82, 0x49, 0x6c, 0x82, 0x49, 0x6e, 0x82, 0x49, 0x6f, 0x00, 0x84, 0x92, 0x73, 0x82,
            0x49, 0x00, 0xe2, 0x06, 0x83, 0x18, 0x74, 0x7b, 0x41, 0x0c, 0x2b, 0xa9, 0x41, 0x0c,
            0x25, 0x31, 0x80, 0x03,
        ];

        let actual = decode(&src)?;
        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_decode_with_order_1() -> io::Result<()> {
        let expected = b"noodles";

        let src = [
            0x01, 0x3b, 0x00, 0x00, 0x00, 0x07, 0x00, 0x00, 0x00, 0x00, 0x64, 0x84, 0x00, 0x6e,
            0x84, 0x00, 0x6f, 0x00, 0x87, 0xff, 0x00, 0x64, 0x6c, 0x8f, 0xff, 0x00, 0x65, 0x00,
            0x73, 0x8f, 0xff, 0x00, 0x6c, 0x65, 0x8f, 0xff, 0x00, 0x6e, 0x6f, 0x8f, 0xff, 0x00,
            0x6f, 0x00, 0x64, 0x87, 0xff, 0x6f, 0x88, 0x00, 0x00, 0x00, 0x00, 0x04, 0x00, 0x02,
            0x02, 0x28, 0x00, 0x01, 0x02, 0x28, 0x00, 0x01, 0x02, 0x60, 0x00, 0x02,
        ];

        let actual = decode(&src)?;
        assert_eq!(actual, expected);

        Ok(())
    }
}
