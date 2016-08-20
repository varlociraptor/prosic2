use std::io;
use std::error::Error;
use std::io::prelude::*;
use std::fs::File;

use rustc_serialize::json;
use clap;
use csv;
use itertools::Itertools;

use libprosic;
use libprosic::model::AlleleFreq;
use libprosic::estimation;

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
        let json = json::encode(&estimate).unwrap();
        let mut f = try!(File::create(path));
        try!(f.write_all(json.as_bytes()));
    }
    Ok(())
}


struct DummyEvent {
    pub name: String
}


impl libprosic::Event for DummyEvent {
    fn name(&self) -> &str {
        &self.name
    }
}


pub fn fdr(matches: &clap::ArgMatches) -> Result<(), Box<Error>> {
    let inbcf = matches.value_of("calls").unwrap();
    let outbcf = matches.value_of("output").unwrap_or("-");
    let events = matches.values_of("events").unwrap().map(|e| {
        DummyEvent { name: e.to_owned() }
    }).collect_vec();

    try!(estimation::fdr::annotate(&inbcf, &outbcf, &events));

    Ok(())
}
