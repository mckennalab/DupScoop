mod convex;
mod kmer_orientation;
mod score_matrix;
pub mod mymatrix;
pub mod needleman;
pub mod smith_waterman_no_diag;

extern crate bio;
extern crate clap;
extern crate csv;
extern crate matrix;
extern crate string_builder;
extern crate indicatif;

use std::fs::File;
use std::io::prelude::*;
use std::io::{Write};
use needleman::{Scores, Alignment};
use clap::{Arg, App};
use std::cmp::min;
use std::iter::FromIterator;


fn main() -> std::io::Result<()> {
    let matches = App::new("DupScoop")
        .version("1.0")
        .author("Aaron M. <aaronatwpi@gmail.com>")
        .about("Deduplicate plasmid assemblies")
        .arg(Arg::with_name("reference")
            .short("r")
            .long("ref")
            .value_name("FILE")
            .help("The reference we will align to")
            .takes_value(true))
        .arg(Arg::with_name("output")
             .short("o")
             .long("output")
             .value_name("FILE")
             .help("the alignment output file")
             .takes_value(true))
        .arg(Arg::with_name("minLength")
            .short("m")
            .long("min")
            .value_name("INT")
            .help("the minimum length for a segment to be considered a duplication")
            .takes_value(true))
        .arg(Arg::with_name("minScoreProportion")
            .short("s")
            .long("score")
            .value_name("FLOAT")
            .help("the score ")
            .takes_value(true))
        .arg(Arg::with_name("minDiagDistance")
            .short("d")
            .long("diagonaldist")
            .value_name("INT")
            .help("the score ")
            .takes_value(true))
        .get_matches();

    let output_file = matches.value_of("output").unwrap_or("output.fa");
    let reference_file = matches.value_of("reference").unwrap_or("reference.fa");

    let min_score_prop: f64 = matches.value_of("minScoreProportion").unwrap_or("0.9").parse::<f64>().unwrap();
    let min_length: u64 = matches.value_of("minLength").unwrap_or("1000").parse::<u64>().unwrap();
    let diag_dist: i32 = matches.value_of("minDiagDistance").unwrap_or("10").parse::<i32>().unwrap();

    let reference_as_chars = reference_to_sequence(reference_file).unwrap();
    let mut reference_as_chars_duplicated = reference_to_sequence(&reference_file).unwrap();
    let mut refClone = reference_as_chars_duplicated.clone();
    reference_as_chars_duplicated.append(&mut refClone);

    let scores = Scores::plasmid_aligmment_scores();

    let mut output = File::create(output_file).unwrap();

    let check_dups = check_for_duplicate_region(&reference_as_chars, &reference_as_chars_duplicated, min_score_prop, min_length, &scores, diag_dist);
    if check_dups.0 {
        let rotated_reference = rotate_reference(&reference_as_chars,  check_dups.1);
        let resulting_reference = String::from_iter(align_and_remove_dup(&rotated_reference, min_score_prop, min_length, &scores, diag_dist));
        println!("Dup!");
        writeln!(output,">reference\n{}\n",resulting_reference);
    } else {
        println!("No dups found!");
        writeln!(output,">reference\n{}\n",String::from_iter(reference_as_chars));
    }

    Ok(())
}

fn align_and_remove_dup(reference: &Vec<char>, min_score_prop: f64, min_length: u64, scores: &Scores, diag_dist: i32) -> Vec<char> {
    let alignment = smith_waterman_no_diag::smith_waterman_no_diag(&reference, &reference, &scores, diag_dist);
    let seq_one_aligned= String::from_iter(alignment.seq_one_aligned.clone().into_iter().filter(|&x| x != '-'));
    let seq_two_aligned= String::from_iter(alignment.seq_two_aligned.clone().into_iter().filter(|&x| x != '-'));
    let min_size = min(seq_one_aligned.len(), seq_two_aligned.len());

    let mut start_del = alignment.start_x;
    let mut end_del = alignment.end_x;
    println!("Alignment starts and stops {},{} with score {}, and {},{} with lengths {} and {}, sequences {} and {} from {} and {}",alignment.start_x,
             alignment.end_x,
             alignment.score,
             alignment.start_y,
             alignment.end_y,
             seq_one_aligned.len(),
             seq_two_aligned.len(),
             seq_one_aligned,
             seq_two_aligned,
             String::from_iter(alignment.seq_one_aligned),
             String::from_iter(alignment.seq_two_aligned));

    if seq_two_aligned.len() == min_size {
        start_del = alignment.start_y;
        end_del = alignment.end_y;
    }
    let split_at_start = reference.split_at(start_del);
    let mut first_half: Vec<char> = split_at_start.0.to_vec();
    let mut second_half: Vec<char> = split_at_start.1[(end_del-start_del)..].to_vec();
    first_half.append(&mut second_half);
    first_half

}

fn rotate_reference(reference: &Vec<char>,offset: usize) -> Vec<char> {
    let mut new_ref = reference.clone();
    new_ref.rotate_right(offset);
    new_ref
}

fn aligned_distance(alignment: &Alignment) -> u32 {
    let it = alignment.seq_one_aligned.iter().zip(alignment.seq_two_aligned.iter());
    let mut differences = 0;
    for (i, (x, y)) in it.enumerate() {
        if x.to_uppercase().to_string() != y.to_uppercase().to_string() {
            differences += 1
        }
    }
    differences
}

fn check_for_duplicate_region(reference: &Vec<char>, reference_dup: &Vec<char>, min_score_prop: f64, min_length: u64, scores: &Scores, diag_dist: i32) -> (bool, usize) {
    let alignment = smith_waterman_no_diag::smith_waterman_no_diag(&reference, &reference_dup, &scores, diag_dist);
    let length_one = alignment.end_x - alignment.start_x;
    let length_two = alignment.end_y - alignment.start_y;
    let min_size = min(length_one, length_two);
    let seq_one_aligned= String::from_iter(alignment.seq_one_aligned.clone().into_iter().filter(|&x| x != '-'));
    let seq_two_aligned= String::from_iter(alignment.seq_two_aligned.clone().into_iter().filter(|&x| x != '-'));

    let seq1_aligned_len = alignment.seq_one_aligned.clone().len() as f64;
    let seq1_aligned = alignment.seq_one_aligned.clone();
    let seq2_aligned = alignment.seq_two_aligned.clone();
    let start_y = alignment.start_y;
    let differences = aligned_distance(&alignment);
    let matching_prop = 1.0 - (differences as f64/ seq1_aligned_len);


    println!("PRE starts and stops {},{} with score {}, and {},{} with lengths {} and {}, sequences {} and {} from {} and {}",alignment.start_x,
             alignment.end_x,
             alignment.score,
             alignment.start_y,
             alignment.end_y,
             seq_one_aligned.len(),
             seq_two_aligned.len(),
             seq_one_aligned,
             seq_two_aligned,
             String::from_iter(seq1_aligned),
             String::from_iter(seq2_aligned));

    (min_size > min_length as usize && min_score_prop < matching_prop,start_y)
}

fn reference_to_sequence(reference_file: &str) -> Result< Vec<char>, std::io::Error> {
    let mut file = File::open(reference_file)?;
    let mut reference = String::new();
    file.read_to_string(&mut reference)?;

    // slice off the first line (name), and combine the rest into a single string of characters
    let first_line_marker = reference.find("\n");
    println!("first endline: {}", first_line_marker.unwrap());

    match first_line_marker {
        // The division was valid
        Some(x) => {
            reference = reference.split_off(x).replace("\n", "");
        }
        // The division was invalid
        None => panic!("We couldn't find the reference name in your input file {}", reference_file),
    }

    let reference_as_chars: Vec<char> = reference.to_string().to_uppercase().chars().collect();
    Ok(reference_as_chars)
}

