use std::collections::HashMap;

#[allow(dead_code)]
pub struct PairedScores {
    counts: HashMap<String, f64>,
    name: &'static str,
}

#[allow(dead_code)]
impl PairedScores {
    pub fn default_scores(match_score: &f64, mismatch_score: &f64) -> PairedScores {
        let mut my_counts: HashMap<String, f64> = HashMap::new();

        let bases = vec!['A', 'C', 'G', 'T', 'U', 'R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N'];
        for base1_index in 0..bases.len() {
            for base2_index in 0..bases.len() {
                let key = bases.get(base1_index).unwrap().to_string() + "-" + &bases.get(base2_index).unwrap().to_string();
                my_counts.insert(key, match (bases.get(base1_index).unwrap(),bases.get(base2_index).unwrap()) {
                    (_x, _y) if _x == _y => {*match_score} // they match regardless of what they are
                    ('R','A') | ('R','G') => {*match_score} // they match regardless of what they are
                    ('Y','C') | ('Y','T') => {*match_score} // they match regardless of what they are
                    ('M','A') | ('M','C') => {*match_score} // they match regardless of what they are
                    ('S','C') | ('S','G') => {*match_score} // they match regardless of what they are
                    ('W','A') | ('W','T') => {*match_score} // they match regardless of what they are
                    ('B','C') | ('B','G') | ('B','T') => {*match_score} // they match regardless of what they are
                    ('D','A') | ('D','G') | ('D','T') => {*match_score} // they match regardless of what they are
                    ('H','A') | ('H','C') | ('H','T') => {*match_score} // they match regardless of what they are
                    ('V','A') | ('V','C') | ('V','G') => {*match_score} // they match regardless of what they are
                    ('N',_) => {*match_score} // they match regardless of what they are
                    (_,_) => {*mismatch_score}
                });
                let key2 = bases.get(base1_index).unwrap().to_string() + "-" + &bases.get(base2_index).unwrap().to_string();
                my_counts.insert(key2, match (bases.get(base1_index).unwrap(),bases.get(base2_index).unwrap()) {
                    (_x, _y) if _x == _y => {*match_score} // they match regardless of what they are
                    ('R','A') | ('R','G') => {*match_score} // they match regardless of what they are
                    ('Y','C') | ('Y','T') => {*match_score} // they match regardless of what they are
                    ('M','A') | ('M','C') => {*match_score} // they match regardless of what they are
                    ('S','C') | ('S','G') => {*match_score} // they match regardless of what they are
                    ('W','A') | ('W','T') => {*match_score} // they match regardless of what they are
                    ('B','C') | ('B','G') | ('B','T') => {*match_score} // they match regardless of what they are
                    ('D','A') | ('D','G') | ('D','T') => {*match_score} // they match regardless of what they are
                    ('H','A') | ('H','C') | ('H','T') => {*match_score} // they match regardless of what they are
                    ('V','A') | ('V','C') | ('V','G') => {*match_score} // they match regardless of what they are
                    ('N',_) => {*match_score} // they match regardless of what they are
                    (_,_) => {*mismatch_score}
                });
            }
        }
        PairedScores{counts: my_counts,name: "Default"}
    }
}