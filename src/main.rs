mod convex;
mod kmer_orientation;
mod score_matrix;
pub mod mymatrix;
pub mod needleman_no_diag;
pub mod smith_waterman;

extern crate bio;
extern crate clap;
extern crate csv;
extern crate matrix;
extern crate string_builder;

use std::fs::File;
use std::io::prelude::*;
use std::io::{Write};
use needleman_no_diag::{Scores};
use clap::{Arg, App};


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
        .arg(Arg::with_name("aligner")
             .short("a")
             .long("aligner")
             .value_name("STRING")
             .help("the alignment method")
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
        .get_matches();

    // pull out the command line arguments and setup input and output
    // ----------------------------------------------------------------
    let input_file = matches.value_of("input").unwrap_or("input.fq");
    let output_file = matches.value_of("output").unwrap_or("output.fa");
    let reference_file = matches.value_of("reference").unwrap_or("reference.fa");
    let read_file = matches.value_of("input").unwrap_or("reference.fa");
    let aligner = matches.value_of("aligner").unwrap_or("swa");
    println!("input : {}, output : {} aligner : {}",input_file, output_file, aligner);

    // read in the reference, and covert to a string
    // ----------------------------------------------------------------
    let mut file = File::open(reference_file)?;
    let mut reference = String::new();
    file.read_to_string(&mut reference)?;

    // slice off the first line (name), and combine the rest into a single string of characters
    let first_line_marker = reference.find("\n");
    println!("first endline: {}" ,first_line_marker.unwrap());

    match first_line_marker {
        // The division was valid
        Some(x) => {
            reference = reference.split_off(x).replace("\n", "");
        }
        // The division was invalid
        None    => panic!("We couldn't find the reference name in your input file {}", input_file),
    }

    let reference_as_chars: Vec<char> = reference.to_string().chars().collect();

    //let kmer_orientation = ReferenceKmers::generate_kmers(&reference,&20);


    // read in sequences and align them
    // ----------------------------------------------------------------
    let scores = Scores::default_scores();
    
    let mut file = File::open(read_file)?;
    let mut read = String::new();
    file.read_to_string(&mut read)?;

    // slice off the first line (name), and combine the rest into a single string of characters
    let first_line_marker = read.find("\n");
    println!("first endline: {}" ,first_line_marker.unwrap());

    match first_line_marker {
        // The division was valid
        Some(x) => {
            read = read.split_off(x).replace("\n", "");
        }
        // The division was invalid
        None    => panic!("We couldn't find the reference name in your input file {}", input_file),
    }

    let read_as_chars: Vec<char> = read.to_string().chars().collect();

    let mut output = File::create(output_file).unwrap();

    let mut count = 0;

    let alignment = smith_waterman::smith_waterman_no_diag(&reference_as_chars, &read_as_chars, &scores, 10);
                
    let str1align: String = alignment.seq_one_aligned.into_iter().collect();
    let str2align: String = alignment.seq_two_aligned.into_iter().collect();
    println!("From {},{} to {},{}",alignment.start_x,alignment.start_y,alignment.end_x,alignment.end_y);
    writeln!(output, ">{}", "ref")?;
    writeln!(output, "{}", str1align)?;
    writeln!(output, ">{}", "read")?;
    writeln!(output, "{}", str2align)?;

    count += 1;
    if count % 1000 == 0 {
        println!("Processed {} reads", count);
    }

    
    //output.close();

    // we're ok!
    Ok(())
}

