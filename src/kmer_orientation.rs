use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt;

pub struct ReferenceKmers {
    kmers_forward: HashSet<String>,
    kmers_reverse: HashSet<String>,
    kmer_size: usize,
}

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
enum ReadOrientation {
    FWD,
    REV,
    UNKNOWN
}
impl fmt::Display for ReadOrientation {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}
impl ReferenceKmers {
    /// generate a set of unique forward and reverse kmers for a given reference
    pub fn generate_kmers(sequence: &str, kmer_size: &usize) -> ReferenceKmers {
        let forward_kmers = ReferenceKmers::sequence_to_kmers(sequence,kmer_size);
        let reverse_kmers = ReferenceKmers::sequence_to_kmers(&ReferenceKmers::reverse_complement_sequence(sequence),kmer_size);
        ReferenceKmers{kmers_forward: forward_kmers, kmers_reverse: reverse_kmers, kmer_size: *kmer_size}
    }

    pub fn sequence_to_kmers(sequence: &str, kmer_size: &usize) -> HashSet<String> {
        let mut kmer_set = HashSet::new();
        for seq in sequence.chars().collect::<Vec<char>>().windows(*kmer_size) {
            kmer_set.insert(seq.into_iter().collect());
        }
        kmer_set
    }

    pub fn reverse_complement_sequence(sequence: &str) -> String {
        sequence.chars().map(|letter| ReferenceKmers::complement(&letter)).rev().collect::<String>()
    }

    pub fn complement(base: &char) -> char {
        match base {
            'A' | 'a' => 'T',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            'T' | 't' => 'A',
            _ => 'N',
        }
    }

    pub fn vote_orientation(&self, sequence: &str, minRatio: &f32, min_count: &usize, kmer_size: &usize) -> ReadOrientation {
        let mut counts : HashMap<ReadOrientation, usize> = HashMap::new();
        let test_kmers = ReferenceKmers::sequence_to_kmers(sequence, kmer_size);

        counts.insert(ReadOrientation::FWD,self.kmers_forward.intersection(&test_kmers).collect::<Vec<&String>>().len());
        counts.insert(ReadOrientation::REV,self.kmers_reverse.intersection(&test_kmers).collect::<Vec<&String>>().len());
        match ReferenceKmers::max_key_by_value(&counts).cloned() {
            Some(p) => {
                let max_key_count = counts[&p];
                let total: usize = counts.iter().map(|(k, v)| v).sum();
                if max_key_count >= *min_count && (max_key_count as f32)/(total  as f32) >= *minRatio {
                    p
                } else {
                    ReadOrientation::UNKNOWN
                }
            }
            None => ReadOrientation::UNKNOWN
        }

    }

    // https://stackoverflow.com/questions/62525693/how-do-i-get-the-key-associated-with-the-maximum-value-of-a-rust-hashmap
    fn max_key_by_value<K, V>(a_hash_map: &HashMap<K, V>) -> Option<&K>
        where
            V: Ord,
    {
        a_hash_map
            .iter()
            .max_by(|a, b| a.1.cmp(&b.1))
            .map(|(k, _v)| k)
    }

}


#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_complement() {
        assert_eq!(ReferenceKmers::complement(&'A'), 'T');
        assert_eq!(ReferenceKmers::complement(&'c'), 'G');
        assert_eq!(ReferenceKmers::complement(&'Y'), 'N');
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(ReferenceKmers::reverse_complement_sequence(&"ACGGT"), "ACCGT");
        assert_eq!(ReferenceKmers::reverse_complement_sequence(&"TTTTA"), "TAAAA");
        assert_eq!(ReferenceKmers::reverse_complement_sequence(&"CCGAA"), "TTCGG");
    }

    #[test]
    fn test_sequence_to_kmers() {
        let kmers = ReferenceKmers::sequence_to_kmers(&"ACGGT", &3);
        assert!(kmers.contains("ACG"));
        assert!(kmers.contains("CGG"));
        assert!(kmers.contains("GGT"));
    }

    #[test]
    fn test_kmer_orientation_basic() {
        let kmers = ReferenceKmers::generate_kmers(&"ACGGTAATTGGCC", &5);
        let orientation = kmers.vote_orientation("ACGGT", &0.5, &1,&5);
        assert_eq!(orientation,ReadOrientation::FWD)
    }

    #[test]
    fn test_kmer_orientation_multi() {
        let kmers = ReferenceKmers::generate_kmers(&"ACGGTCCGGTTTAATTAGAGATTTTT", &5);
        let orientation = kmers.vote_orientation("ACGGTCCGGTTTAATTAGAGATTTTT", &0.5, &20,&5);
        assert_eq!(orientation,ReadOrientation::FWD)
    }

    #[test]
    fn test_kmer_orientation_below_threshold() {
        let kmers = ReferenceKmers::generate_kmers(&"ACGGTCCGGTTTAATTAGAGATTTTT", &5);
        let orientation = kmers.vote_orientation("ACGGTCCGGTTTAATTAGAGATTTTT", &0.5, &23,&5);
        assert_eq!(orientation,ReadOrientation::UNKNOWN)
    }
}