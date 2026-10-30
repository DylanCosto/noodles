//! Nucleotide sequence helpers.

/// [`complement`] for every byte.
static COMPLEMENTS: [u8; 256] = {
    let mut table = [0; 256];
    let mut i = 0;
    while i < 256 {
        table[i] = complement(i as u8);
        i += 1;
    }
    table
};

/// Reverse complements IUPAC nucleotides in place, preserving case, as samtools does.
pub fn reverse_complement(sequence: &mut [u8]) {
    sequence.reverse();

    for base in sequence {
        *base = COMPLEMENTS[usize::from(*base)];
    }
}

const fn complement(base: u8) -> u8 {
    let c = match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        other => other,
    };

    if base.is_ascii_lowercase() {
        c.to_ascii_lowercase()
    } else {
        c
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        let mut s = b"ACGTNacgtnRY".to_vec();
        reverse_complement(&mut s);
        assert_eq!(s, b"RYnacgtNACGT");
    }
}
