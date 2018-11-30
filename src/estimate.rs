use std::io;
use std::error::Error;
use std::fs::File;

use serde_json;
use clap;
use csv;
use bio::stats::{Prob, LogProb};

use libprosic;
use libprosic::model::AlleleFreq;
use libprosic::estimation;
use libprosic::model;

pub fn effective_mutation_rate(matches: &clap::ArgMatches) -> Result<(), Box<Error>> {
    let min_af = value_t!(matches, "min-af", f64).unwrap_or(0.12);
    let max_af = value_t!(matches, "max-af", f64).unwrap_or(0.25);
    let mut reader = csv::Reader::from_reader(io::stdin());
    let freqs = try!(reader.decode().collect::<Result<Vec<f64>, _>>());
    let estimate = estimation::effective_mutation_rate::estimate(freqs.into_iter().filter(|&f| {
        f >= min_af && f <= max_af
    }).map(|f| AlleleFreq(f)));

    // print estimated mutation rate to stdout
    println!("{}", estimate.effective_mutation_rate());

    // if --fit is given, print data visualizing model fit
    if let Some(path) = matches.value_of("fit") {
        let mut f = try!(File::create(path));
        serde_json::to_writer(&mut f, &estimate)?;
    }
    Ok(())
}
