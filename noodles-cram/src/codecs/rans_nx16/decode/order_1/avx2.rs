//! AVX2 decoding of 32-state order-1 rANS Nx16 streams.
//!
//! Each round decodes one symbol per state, 8 states per vector: the context offsets and table
//! entries are gathered, the states advance together, and renormalization expands the next input
//! words into the lanes that need them, in state order. This is the same computation as the
//! scalar loop, so the two can hand off between rounds.

use std::{
    arch::x86_64::{
        __m256i, _mm_cvtsi32_si128, _mm_loadu_si128, _mm256_add_epi32, _mm256_and_si256,
        _mm256_blendv_epi8, _mm256_castsi256_ps, _mm256_cmpeq_epi32, _mm256_cvtepu16_epi32,
        _mm256_i32gather_epi32, _mm256_loadu_si256, _mm256_min_epu32, _mm256_movemask_ps,
        _mm256_mullo_epi32, _mm256_or_si256, _mm256_permutevar8x32_epi32, _mm256_set1_epi32,
        _mm256_setzero_si256, _mm256_slli_epi32, _mm256_srl_epi32, _mm256_srli_epi32,
        _mm256_storeu_si256,
    },
    io,
};

use super::Tables;

const STATE_COUNT: usize = 32;
const LANES: usize = 8;
const VECTORS: usize = STATE_COUNT / LANES;

/// The input a round may read: two bytes per state.
const ROUND_WINDOW: usize = 2 * STATE_COUNT;

/// For each 8-bit mask of lanes that renormalize, the input word each lane takes: lane `k` takes
/// word `popcount(mask & ((1 << k) - 1))`. Lanes that do not renormalize take word 0, unused.
static EXPAND: [[u32; LANES]; 256] = build_expand();

const fn build_expand() -> [[u32; LANES]; 256] {
    let mut table = [[0; LANES]; 256];
    let mut mask = 0;

    while mask < 256 {
        let mut word = 0;
        let mut lane = 0;

        while lane < LANES {
            if mask & (1 << lane) != 0 {
                table[mask][lane] = word;
                word += 1;
            }

            lane += 1;
        }

        mask += 1;
    }

    table
}

/// Decodes rounds `0..chunk_size` of the 32 interleaved streams into `dst` (state `j` writes
/// `dst[j * chunk_size + i]` in round `i`), updating `states`, `prev_syms`, and `src` as the scalar
/// loop does.
///
/// # Errors
///
/// Returns an [`io::ErrorKind::UnexpectedEof`] error if `src` ends before the last round, as the
/// scalar loop does.
///
/// # Safety
///
/// The CPU must support AVX2.
#[target_feature(enable = "avx2")]
pub(super) unsafe fn decode(
    tables: &Tables,
    bits: u32,
    states: &mut [u32; STATE_COUNT],
    prev_syms: &mut [u8; STATE_COUNT],
    src: &mut &[u8],
    dst: &mut [u8],
    chunk_size: usize,
) -> io::Result<()> {
    assert_eq!(dst.len(), STATE_COUNT * chunk_size);

    let slot_mask = _mm256_set1_epi32((1 << bits) - 1);
    let shift = _mm_cvtsi32_si128(bits as i32);
    let sym_mask = _mm256_set1_epi32(0xff);
    let frequency_mask = _mm256_set1_epi32(0xfff);
    let one = _mm256_set1_epi32(1);
    let renormalize_max = _mm256_set1_epi32((1 << 15) - 1);

    let offsets = tables.offsets.as_ptr().cast::<i32>();
    let entries = tables.entries.as_ptr().cast::<i32>();

    let mut s: [__m256i; VECTORS] = std::array::from_fn(|v| {
        // SAFETY: `states` has `VECTORS * LANES` values.
        unsafe { _mm256_loadu_si256(states.as_ptr().add(v * LANES).cast()) }
    });

    let mut prev: [__m256i; VECTORS] = std::array::from_fn(|v| {
        let syms: [u32; LANES] = std::array::from_fn(|k| u32::from(prev_syms[v * LANES + k]));
        // SAFETY: `syms` has `LANES` values.
        unsafe { _mm256_loadu_si256(syms.as_ptr().cast()) }
    });

    let mut syms = [0u32; STATE_COUNT];
    let mut tail = [0; ROUND_WINDOW];

    for round in 0..chunk_size {
        for v in 0..VECTORS {
            // SAFETY: `prev` lanes are symbols (< 256 = `offsets.len()`), and every offset plus a
            // slot (< `1 << bits`) is within `entries` (see `Tables`). Both tables are far smaller
            // than `i32::MAX` elements.
            let entry = unsafe {
                let offset = _mm256_i32gather_epi32::<4>(offsets, prev[v]);
                let i = _mm256_add_epi32(offset, _mm256_and_si256(s[v], slot_mask));
                _mm256_i32gather_epi32::<4>(entries, i)
            };

            let frequency = _mm256_add_epi32(
                _mm256_and_si256(_mm256_srli_epi32::<8>(entry), frequency_mask),
                one,
            );
            let offset = _mm256_srli_epi32::<20>(entry);
            let x = _mm256_srl_epi32(s[v], shift);
            s[v] = _mm256_add_epi32(_mm256_mullo_epi32(frequency, x), offset);

            prev[v] = _mm256_and_si256(entry, sym_mask);

            // SAFETY: `syms` has `LANES` values from `v * LANES`.
            unsafe { _mm256_storeu_si256(syms.as_mut_ptr().add(v * LANES).cast(), prev[v]) };
        }

        for (j, &sym) in syms.iter().enumerate() {
            dst[j * chunk_size + round] = sym as u8;
        }

        // States renormalize in order, each taking the next input word if it is below 2^15.
        let mut renormalize = [_mm256_setzero_si256(); VECTORS];
        let mut masks = [0; VECTORS];
        let mut needed = 0;

        for ((state, renormalize), mask) in s.iter().zip(&mut renormalize).zip(&mut masks) {
            // `min(s, 2^15 - 1) == s` is an unsigned `s < 2^15`.
            *renormalize = _mm256_cmpeq_epi32(_mm256_min_epu32(*state, renormalize_max), *state);
            *mask = _mm256_movemask_ps(_mm256_castsi256_ps(*renormalize)) as usize;
            needed += 2 * mask.count_ones() as usize;
        }

        // Near the end of the input, the loads read a padded copy. The padding only reaches
        // lanes that do not renormalize.
        let window: &[u8] = if src.len() >= ROUND_WINDOW {
            &src[..ROUND_WINDOW]
        } else if needed <= src.len() {
            tail[..src.len()].copy_from_slice(src);
            &tail
        } else {
            return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
        };

        let mut n = 0;

        for ((state, renormalize), mask) in s.iter_mut().zip(renormalize).zip(masks) {
            // SAFETY: `n <= 2 * LANES * v` before vector `v`, so the 16-byte load ends within
            // `window` (`2 * STATE_COUNT` bytes), and `EXPAND[mask]` has `LANES` values.
            let words = unsafe {
                let words = _mm256_cvtepu16_epi32(_mm_loadu_si128(window.as_ptr().add(n).cast()));
                let expand = _mm256_loadu_si256(EXPAND[mask].as_ptr().cast());
                _mm256_permutevar8x32_epi32(words, expand)
            };

            let renormalized = _mm256_or_si256(_mm256_slli_epi32::<16>(*state), words);
            *state = _mm256_blendv_epi8(*state, renormalized, renormalize);

            n += 2 * mask.count_ones() as usize;
        }

        *src = &src[n..];
    }

    for v in 0..VECTORS {
        // SAFETY: `states` and `syms` have `LANES` values from `v * LANES`.
        unsafe {
            _mm256_storeu_si256(states.as_mut_ptr().add(v * LANES).cast(), s[v]);
            _mm256_storeu_si256(syms.as_mut_ptr().add(v * LANES).cast(), prev[v]);
        }
    }

    for (prev_sym, &sym) in prev_syms.iter_mut().zip(&syms) {
        *prev_sym = sym as u8;
    }

    Ok(())
}
