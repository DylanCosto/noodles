use std::io;

use super::{
    SLOT_COUNT, build_slot_lookup, decode_slot, read_states, renormalize,
    state_cumulative_frequency, state_renormalize, state_step,
};
use crate::{
    codecs::rans_4x8::{ALPHABET_SIZE, STATE_COUNT},
    io::reader::num::{read_itf8_as, read_u8},
};

type Frequencies = [u16; ALPHABET_SIZE]; // F
type CumulativeFrequencies = Frequencies; // C
type CumulativeFrequenciesSymbolsTable = [u8; 4096];

pub fn decode(src: &mut &[u8], dst: &mut [u8]) -> io::Result<()> {
    let frequencies = read_frequencies(src)?;
    let cumulative_frequencies = build_cumulative_frequencies(&frequencies);

    let mut lookup = [0; SLOT_COUNT];

    if build_slot_lookup(&frequencies, &cumulative_frequencies, &mut lookup) {
        decode_with_slot_lookup(src, dst, &lookup)
    } else {
        decode_with_symbols_table(src, dst, &frequencies, &cumulative_frequencies)
    }
}

/// Decodes with a slot lookup (see [`build_slot_lookup`]): the same result as
/// [`decode_with_symbols_table`], with one table read per symbol.
fn decode_with_slot_lookup(src: &mut &[u8], dst: &mut [u8], lookup: &[u32]) -> io::Result<()> {
    let mut states = read_states(src)?;
    let (chunks, remainder) = dst.as_chunks_mut::<STATE_COUNT>();

    for chunk in chunks {
        for (d, state) in chunk.iter_mut().zip(&mut states) {
            let (sym, s) = decode_slot(*state, lookup);
            *d = sym;
            *state = renormalize(s, src)?;
        }
    }

    for (d, state) in remainder.iter_mut().zip(&mut states) {
        let (sym, s) = decode_slot(*state, lookup);
        *d = sym;
        *state = renormalize(s, src)?;
    }

    Ok(())
}

/// Decodes with the cumulative frequencies symbols table (the general method).
fn decode_with_symbols_table(
    src: &mut &[u8],
    dst: &mut [u8],
    frequencies: &Frequencies,
    cumulative_frequencies: &CumulativeFrequencies,
) -> io::Result<()> {
    let cumulative_frequencies_symbols_table =
        build_cumulative_frequencies_symbols_table(cumulative_frequencies);

    let mut states = read_states(src)?;
    let (chunks, remainder) = dst.as_chunks_mut::<STATE_COUNT>();

    for chunk in chunks {
        for (d, state) in chunk.iter_mut().zip(&mut states) {
            let f = state_cumulative_frequency(*state);
            let sym = cumulative_frequencies_symbols_table[usize::from(f)];

            *d = sym;

            let i = usize::from(sym);
            *state = state_step(*state, frequencies[i], cumulative_frequencies[i]);
            *state = state_renormalize(*state, src)?;
        }
    }

    for (d, state) in remainder.iter_mut().zip(&mut states) {
        let f = state_cumulative_frequency(*state);
        let sym = cumulative_frequencies_symbols_table[usize::from(f)];

        *d = sym;

        let i = usize::from(sym);
        *state = state_step(*state, frequencies[i], cumulative_frequencies[i]);
        *state = state_renormalize(*state, src)?;
    }

    Ok(())
}

pub(super) fn read_frequencies(src: &mut &[u8]) -> io::Result<Frequencies> {
    const NUL: u8 = 0x00;

    let mut frequencies = [0; ALPHABET_SIZE];

    let mut sym = read_u8(src)?;
    let mut prev_sym = sym;

    loop {
        let f = read_itf8_as(src)?;
        frequencies[usize::from(sym)] = f;

        sym = read_u8(src)?;

        if sym == NUL {
            break;
        }

        if sym - 1 == prev_sym {
            let len = read_u8(src)?;

            for _ in 0..len {
                let f = read_itf8_as(src)?;
                frequencies[usize::from(sym)] = f;
                sym += 1;
            }
        }

        prev_sym = sym;
    }

    Ok(frequencies)
}

pub(super) fn build_cumulative_frequencies(frequencies: &Frequencies) -> CumulativeFrequencies {
    let mut cumulative_frequencies = [0; ALPHABET_SIZE];

    let mut f = cumulative_frequencies[0];

    for (next_f, g) in cumulative_frequencies[1..].iter_mut().zip(frequencies) {
        *next_f = f + g;
        f = *next_f;
    }

    cumulative_frequencies
}

pub(super) fn build_cumulative_frequencies_symbols_table(
    cumulative_freqs: &CumulativeFrequencies,
) -> CumulativeFrequenciesSymbolsTable {
    // Each symbol covers [C[sym], C[sym + 1]); the last symbol and any slots past the total map
    // to 255.
    let mut table = [u8::MAX; 4096];

    for (sym, w) in cumulative_freqs.windows(2).enumerate() {
        let (start, end) = (usize::from(w[0]).min(4096), usize::from(w[1]).min(4096));

        if start < end {
            table[start..end].fill(sym as u8);
        }
    }

    table
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::codecs::rans_4x8::decode::test_support::{Rng, random_frequencies};

    #[test]
    fn test_decode_with_slot_lookup_matches_symbols_table() {
        let mut rng = Rng(0x9e37_79b9_7f4a_7c15);

        for _ in 0..2000 {
            let frequencies = random_frequencies(&mut rng, 4096);
            let cumulative_frequencies = build_cumulative_frequencies(&frequencies);
            let mut lookup = [0; SLOT_COUNT];
            assert!(build_slot_lookup(
                &frequencies,
                &cumulative_frequencies,
                &mut lookup
            ));

            // Random states and stream bytes, so slots past the total are reached too.
            let src: Vec<u8> = (0..16 + rng.below(64)).map(|_| rng.next() as u8).collect();
            let len = rng.below(48) as usize;
            let (mut expected, mut actual) = (vec![0; len], vec![0; len]);

            let expected_result = decode_with_symbols_table(
                &mut &src[..],
                &mut expected,
                &frequencies,
                &cumulative_frequencies,
            );
            let actual_result = decode_with_slot_lookup(&mut &src[..], &mut actual, &lookup);

            match (expected_result, actual_result) {
                (Ok(()), Ok(())) => assert_eq!(actual, expected),
                (Err(e), Err(a)) => assert_eq!(a.kind(), e.kind()),
                (e, a) => panic!("results differ: {e:?} != {a:?}"),
            }
        }
    }

    #[test]
    fn test_build_cumulative_frequencies_symbols_table_matches_slot_scan() {
        // The table as originally built: scan the slots, advancing past each symbol's range.
        fn scan(cumulative_freqs: &CumulativeFrequencies) -> CumulativeFrequenciesSymbolsTable {
            let mut table = [0; 4096];
            let mut sym = 0;

            for (f, g) in (0u16..).zip(&mut table) {
                while sym < u8::MAX && f >= cumulative_freqs[usize::from(sym + 1)] {
                    sym += 1;
                }

                *g = sym;
            }

            table
        }

        let mut rng = Rng(0x1234_5678_9abc_def1);

        for _ in 0..2000 {
            let frequencies = random_frequencies(&mut rng, 4096);
            let cumulative_frequencies = build_cumulative_frequencies(&frequencies);
            assert_eq!(
                build_cumulative_frequencies_symbols_table(&cumulative_frequencies),
                scan(&cumulative_frequencies)
            );
        }

        let empty = build_cumulative_frequencies(&[0; ALPHABET_SIZE]);
        assert_eq!(
            build_cumulative_frequencies_symbols_table(&empty),
            scan(&empty)
        );
    }

    #[test]
    fn test_build_slot_lookup_rejects_wide_frequencies() {
        let mut frequencies = [0; ALPHABET_SIZE];
        frequencies[usize::from(b'a')] = 4096;
        let cumulative_frequencies = build_cumulative_frequencies(&frequencies);
        let mut lookup = [0; SLOT_COUNT];
        assert!(!build_slot_lookup(
            &frequencies,
            &cumulative_frequencies,
            &mut lookup
        ));
    }

    #[test]
    fn test_read_frequencies() -> io::Result<()> {
        let src = [
            b'a', // symbol = 'a'
            0x05, // frequencies['a'] = 5
            b'b', // symbol = 'b'
            0x02, // run length = 2
            0x02, // frequencies['b'] = 2
            0x01, // frequencies['c'] = 1
            0x01, // frequencies['d'] = 1
            b'r', // symbol = 'r'
            0x02, // frequencies['r'] = 2
            0x00, // EOF
        ];

        let mut expected = [0; ALPHABET_SIZE];
        expected[usize::from(b'a')] = 5;
        expected[usize::from(b'b')] = 2;
        expected[usize::from(b'c')] = 1;
        expected[usize::from(b'd')] = 1;
        expected[usize::from(b'r')] = 2;

        assert_eq!(read_frequencies(&mut &src[..])?, expected);

        Ok(())
    }

    #[test]
    fn test_build_cumulative_frequencies() {
        let mut frequencies = [0; ALPHABET_SIZE];
        frequencies[usize::from(b'a')] = 5;
        frequencies[usize::from(b'b')] = 2;
        frequencies[usize::from(b'c')] = 1;
        frequencies[usize::from(b'd')] = 1;
        frequencies[usize::from(b'r')] = 2;

        let mut expected = [0; ALPHABET_SIZE];
        expected[..usize::from(b'b')].fill(0);
        expected[usize::from(b'b')..usize::from(b'c')].fill(5);
        expected[usize::from(b'c')..usize::from(b'd')].fill(7);
        expected[usize::from(b'd')..usize::from(b'e')].fill(8);
        expected[usize::from(b'e')..usize::from(b's')].fill(9);
        expected[usize::from(b's')..].fill(11);

        assert_eq!(build_cumulative_frequencies(&frequencies), expected);
    }
}
