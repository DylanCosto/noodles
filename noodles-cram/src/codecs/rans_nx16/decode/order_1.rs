#[cfg(target_arch = "x86_64")]
mod avx2;

use std::io;

use super::{
    order_0, read_alphabet, read_states, renormalize_states, split_off, state_renormalize,
};
use crate::{
    codecs::rans_nx16::ALPHABET_SIZE,
    io::reader::num::{read_u8, read_uint7, read_uint7_as},
};

/// The largest frequency table shift.
///
/// htscodecs writes 10 or 12. Every context with a nonzero row gets a `1 << bits` table, so a
/// larger shift would only let a few bytes of malformed input allocate and fill tens of MiB.
const MAX_BITS: u32 = 12;

/// Decoding tables for every context, stored contiguously so that a symbol costs one context
/// offset lookup rather than a pointer chase. Context `c` has the `1 << bits` entries starting at
/// `offsets[c]`. Contexts absent from the alphabet, and contexts whose frequencies are all zero,
/// share an empty table at offset 0, so `offsets[c] + (1 << bits) <= entries.len()` for every
/// context.
///
/// Each entry packs, for cumulative frequency slot `f`, the symbol whose range contains `f` (bits
/// 0..8), its frequency minus 1 (bits 8..20), and `f` minus its cumulative frequency (bits
/// 20..32). One load then gives everything a decoding step needs. A frequency of 0, which only
/// malformed tables have, wraps to 4096; the output of a malformed stream is meaningless either
/// way.
struct Tables {
    entries: Vec<u32>,
    offsets: [u32; ALPHABET_SIZE],
}

impl Tables {
    #[inline(always)]
    fn decode_symbol(&self, prev_sym: u8, state: &mut u32, bits: u32) -> u8 {
        let i =
            self.offsets[usize::from(prev_sym)] as usize + (*state & ((1 << bits) - 1)) as usize;
        let entry = self.entries[i];
        *state = apply_entry(*state, entry, bits);
        entry as u8
    }
}

/// Packs a symbol, its frequency, and a slot's offset into its range as a [`Tables`] entry.
fn pack_entry(sym: u8, frequency: u32, offset: u32) -> u32 {
    u32::from(sym) | ((frequency.wrapping_sub(1) & 0xfff) << 8) | (offset << 20)
}

/// Fills `dst` (`1 << bits` entries) with the packed entries of [`super::build_symbol_table`]'s table,
/// without building its two arrays first.
fn build_entries(frequencies: &[u32; ALPHABET_SIZE], bits: u32, dst: &mut [u32]) -> io::Result<()> {
    let total_frequency = 1 << bits;
    let mut cumulative_frequency: u32 = 0;

    for (sym, &f) in (0..=u8::MAX).zip(frequencies) {
        let next_cumulative_frequency = cumulative_frequency
            .checked_add(f)
            .filter(|&n| n <= total_frequency)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid rANS Nx16 frequency table",
                )
            })?;

        let range = cumulative_frequency as usize..next_cumulative_frequency as usize;

        for (offset, entry) in (0..).zip(&mut dst[range]) {
            *entry = pack_entry(sym, f, offset);
        }

        cumulative_frequency = next_cumulative_frequency;
    }

    // As in `build_symbol_table`, slots at or past the total frequency resolve to the last
    // symbol.
    let last_frequency = frequencies[usize::from(u8::MAX)];
    let last_cumulative_frequency = cumulative_frequency - last_frequency;

    for (slot, entry) in
        (cumulative_frequency..).zip(dst.iter_mut().skip(cumulative_frequency as usize))
    {
        *entry = pack_entry(u8::MAX, last_frequency, slot - last_cumulative_frequency);
    }

    Ok(())
}

/// Advances a state with a [`Tables`] entry, before renormalization.
#[inline(always)]
fn apply_entry(state: u32, entry: u32, bits: u32) -> u32 {
    let frequency = ((entry >> 8) & 0xfff) + 1;
    let offset = entry >> 20;
    // Valid tables cannot overflow here; wrapping keeps malformed initial states panic-free.
    frequency.wrapping_mul(state >> bits).wrapping_add(offset)
}

pub fn decode(src: &mut &[u8], dst: &mut [u8], state_count: usize) -> io::Result<()> {
    let (tables, bits) = read_frequencies(src)?;

    match state_count {
        4 => decode_with_states::<4, false>(src, dst, &tables, bits),
        32 => decode_with_states::<32, true>(src, dst, &tables, bits),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid rANS Nx16 state count",
        )),
    }
}

/// Decodes with `N` states, using SIMD where the CPU supports it if `SIMD` is set.
fn decode_with_states<const N: usize, const SIMD: bool>(
    src: &mut &[u8],
    dst: &mut [u8],
    tables: &Tables,
    bits: u32,
) -> io::Result<()> {
    let mut states = read_states::<N>(src)?;
    let mut prev_syms = [0u8; N];

    let chunk_size = dst.len() / N;
    let (body, last_chunk) = dst.split_at_mut(chunk_size * N);

    let simd_decoded = SIMD
        && chunk_size > 0
        && try_decode_simd(
            tables,
            bits,
            &mut states,
            &mut prev_syms,
            src,
            body,
            chunk_size,
        )?;

    if chunk_size > 0 && !simd_decoded {
        let mut chunks = body.chunks_exact_mut(chunk_size);
        // `body` has exactly `N` chunks of `chunk_size`.
        let mut dsts: [&mut [u8]; N] = std::array::from_fn(|_| chunks.next().unwrap());

        for i in 0..chunk_size {
            for ((d, state), prev_sym) in dsts.iter_mut().zip(&mut states).zip(&mut prev_syms) {
                let sym = tables.decode_symbol(*prev_sym, state, bits);
                d[i] = sym;
                *prev_sym = sym;
            }

            // The state steps are independent, so they run ahead of the ordered renormalization
            // reads.
            renormalize_states(&mut states, src)?;
        }
    }

    let mut state = states[N - 1];
    let mut prev_sym = prev_syms[N - 1];

    for d in last_chunk {
        let sym = tables.decode_symbol(prev_sym, &mut state, bits);
        *d = sym;
        state = state_renormalize(state, src)?;
        prev_sym = sym;
    }

    Ok(())
}

/// Decodes every round of the body with SIMD, if the CPU supports it for `N` states, and returns
/// whether it did.
#[cfg(target_arch = "x86_64")]
fn try_decode_simd<const N: usize>(
    tables: &Tables,
    bits: u32,
    states: &mut [u32; N],
    prev_syms: &mut [u8; N],
    src: &mut &[u8],
    body: &mut [u8],
    chunk_size: usize,
) -> io::Result<bool> {
    let (Ok(states), Ok(prev_syms)) = (
        <&mut [u32; 32]>::try_from(states.as_mut_slice()),
        <&mut [u8; 32]>::try_from(prev_syms.as_mut_slice()),
    ) else {
        return Ok(false);
    };

    if !std::arch::is_x86_feature_detected!("avx2") {
        return Ok(false);
    }

    // SAFETY: AVX2 is available.
    unsafe { avx2::decode(tables, bits, states, prev_syms, src, body, chunk_size)? };

    Ok(true)
}

#[cfg(not(target_arch = "x86_64"))]
fn try_decode_simd<const N: usize>(
    _: &Tables,
    _: u32,
    _: &mut [u32; N],
    _: &mut [u8; N],
    _: &mut &[u8],
    _: &mut [u8],
    _: usize,
) -> io::Result<bool> {
    Ok(false)
}

fn read_frequencies(src: &mut &[u8]) -> io::Result<(Tables, u32)> {
    const STATE_COUNT: usize = 4;

    let n = read_u8(src)?;
    let bits = u32::from(n >> 4);
    let is_compressed = (n & 0x01) != 0;

    if bits > MAX_BITS {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid rANS Nx16 frequency table",
        ));
    }

    let tables = if is_compressed {
        let uncompressed_size = read_uint7_as(src)?;
        let compressed_size = read_uint7_as(src)?;
        let mut compressed_data = split_off(src, compressed_size)?;
        let mut dst = vec![0; uncompressed_size];
        order_0::decode(&mut compressed_data, &mut dst, STATE_COUNT)?;
        read_frequencies_inner(&mut &dst[..], bits)?
    } else {
        read_frequencies_inner(src, bits)?
    };

    Ok((tables, bits))
}

fn read_frequencies_inner(src: &mut &[u8], bits: u32) -> io::Result<Tables> {
    let alphabet = read_alphabet(src)?;

    // Every row is read before any table is allocated, so all-zero rows take no space.
    let mut rows = Vec::new();

    for (context, _) in alphabet.iter().enumerate().filter(|(_, a)| **a) {
        let mut fs = [0; ALPHABET_SIZE];
        let mut iter = alphabet.iter().zip(fs.iter_mut()).filter(|(b, _)| **b);

        while let Some((_, f)) = iter.next() {
            *f = read_uint7(src)?;

            if *f == 0 {
                let n = read_u8(src).map(usize::from)?;
                for _ in iter.by_ref().take(n) {}
            }
        }

        if fs.iter().any(|&f| f > 0) {
            order_0::normalize_frequencies(&mut fs, bits);
            rows.push((context, fs));
        }
    }

    let size = 1 << bits;

    let mut entries = vec![0; (rows.len() + 1) * size];
    let mut offsets = [0; ALPHABET_SIZE];

    let empty_row = [0; ALPHABET_SIZE];
    let tables = [(None, &empty_row)]
        .into_iter()
        .chain(rows.iter().map(|(context, fs)| (Some(*context), fs)));

    for (i, ((context, fs), dst)) in tables.zip(entries.chunks_exact_mut(size)).enumerate() {
        build_entries(fs, bits, dst)?;

        if let Some(context) = context {
            // SAFETY: `i * size < (ALPHABET_SIZE + 1) << MAX_BITS`, which fits in a `u32`.
            offsets[context] = (i * size) as u32;
        }
    }

    Ok(Tables { entries, offsets })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::codecs::rans_nx16::{Flags, encode::encode};

    #[test]
    fn test_decode_round_trips_32_states() -> io::Result<()> {
        // Quality-like runs, and uniform bytes that renormalize nearly every state every round.
        let mut state = 1u32;
        let mut next = move || {
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            state >> 8
        };

        let qualities: Vec<u8> = (0..200_003)
            .scan(30u8, |q, _| {
                let r = next();
                if r % 7 == 0 {
                    *q = 2 + (r % 40) as u8;
                }
                Some(*q + 33)
            })
            .collect();
        let uniform: Vec<u8> = (0..100_001).map(|_| next() as u8).collect();

        for src in [&qualities[..], &uniform[..]] {
            for len in [32, 33, 63, 64, 65, 1_000, 4_097, src.len()] {
                let src = &src[..len];
                let encoded = encode(Flags::ORDER | Flags::N32, src)?;
                assert_eq!(super::super::decode(&encoded, len)?, src, "len = {len}");
                assert_eq!(decode_scalar(&encoded, len)?, src, "len = {len}");

                // The decoder reads the stream to its end, so a truncated one runs out of input.
                let truncated = &encoded[..encoded.len() - 2];
                assert!(
                    super::super::decode(truncated, len)
                        .is_err_and(|e| e.kind() == io::ErrorKind::UnexpectedEof),
                    "len = {len}"
                );
            }
        }

        Ok(())
    }

    /// Decodes a 32-state order-1 stream without SIMD.
    fn decode_scalar(mut src: &[u8], len: usize) -> io::Result<Vec<u8>> {
        let flags = read_u8(&mut src)?;
        assert_eq!(Flags::from(flags), Flags::ORDER | Flags::N32);
        let size: usize = read_uint7_as(&mut src)?;
        assert_eq!(size, len);

        let (tables, bits) = read_frequencies(&mut src)?;
        let mut dst = vec![0; len];
        decode_with_states::<32, false>(&mut src, &mut dst, &tables, bits)?;

        Ok(dst)
    }

    #[test]
    fn test_read_frequencies_shares_the_empty_table_and_rejects_shifts_above_12() -> io::Result<()>
    {
        // Shift 12, uncompressed; alphabet {A, C}; row A: A=1, C=1; row C: all zero (0, then a
        // run of 1 more zero).
        let mut src = [0xc0, b'A', b'C', 0x00, 1, 1, 0, 1];
        let (tables, bits) = read_frequencies(&mut &src[..])?;

        assert_eq!(bits, 12);
        assert_eq!(tables.entries.len(), 2 << 12);
        assert_eq!(tables.offsets[usize::from(b'A')], 1 << 12);
        assert_eq!(tables.offsets[usize::from(b'C')], 0);

        src[0] = 0xd0; // shift 13
        let result = read_frequencies(&mut &src[..]);
        assert!(result.is_err_and(|e| e.kind() == io::ErrorKind::InvalidData));

        Ok(())
    }
}
