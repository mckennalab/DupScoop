use std::collections::HashMap;


pub struct PairedScores {
    counts: HashMap<String, f64>,
    name: &'static str,
}

impl PairedScores {
    pub fn default_scores(match_score: &f64, mismatch_score: &f64) -> PairedScores {
        let mut my_counts: HashMap<String, f64> = HashMap::new();

        let bases = vec!['A', 'C', 'G', 'T', 'U', 'R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N'];
        for base1_index in 0..bases.len() {
            for base2_index in 0..bases.len() {
                let key = bases.get(base1_index).unwrap().to_string() + "-" + &bases.get(base2_index).unwrap().to_string();
                if bases.get(base1_index).unwrap() == bases.get(base2_index).unwrap() {
                    my_counts.insert(key,*match_score);
                } else {
                    my_counts.insert(key,*mismatch_score);
                }
            }
        }
        PairedScores{counts: my_counts,name: "Default"}
    }
}