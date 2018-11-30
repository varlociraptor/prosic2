use std::error::Error;

use clap;
use bio::stats::{Prob, LogProb};
use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use itertools::Itertools;

use libprosic;
use libprosic::filtration;
use libprosic::model;

struct DummyEvent {
    pub name: String
}


impl libprosic::Event for DummyEvent {
    fn name(&self) -> &str {
        &self.name
    }
}


/// Parse `VariantType` from command line arguments.
pub fn parse_vartype(vartype: &str, min_len: Option<u32>, max_len: Option<u32>) -> Result<model::VariantType, Box<Error>> {
    Ok(match (vartype, min_len, max_len) {
        ("SNV", _, _) => model::VariantType::SNV,
        ("INS", Some(min_len), Some(max_len)) => model::VariantType::Insertion(Some(min_len..max_len)),
        ("DEL", Some(min_len), Some(max_len)) => model::VariantType::Deletion(Some(min_len..max_len)),
        ("INS", _, _) => model::VariantType::Insertion(None),
        ("DEL", _, _) => model::VariantType::Deletion(None),
        _ => {
            return Err(Box::new(clap::Error {
                message: "unsupported variant type (supported: SNV, INS, DEL)".to_owned(),
                kind: clap::ErrorKind::InvalidValue,
                info: None
            }));
        }
    })
}


pub fn filter(matches: &clap::ArgMatches) -> Result<(), Box<Error>> {
    let events = matches.values_of("events").unwrap();
    let events = events.into_iter().map(|event| DummyEvent { name: event.to_owned() }).collect_vec();

    if let Some(matches) = matches.subcommand_matches("control-fdr") {
        let call_bcf = matches.value_of("calls").unwrap();

        let vartype = matches.value_of("vartype").unwrap();
        let min_len = value_t!(matches, "min-len", u32).ok();
        let max_len = value_t!(matches, "max-len", u32).ok();
        let vartype = parse_vartype(vartype, min_len, max_len)?;
        let alpha = value_t!(matches, "fdr", f64)?;
        let alpha = LogProb::from(Prob::checked(alpha)?);
        filtration::fdr::control_fdr::<_, _, &str>(
            call_bcf, None, &events, &vartype, alpha
        )?;
    } else if let Some(matches) = matches.subcommand_matches("posterior-odds") {
        let score = matches.value_of("odds").unwrap();
        let min_evidence = match score {
            "none" => KassRaftery::None,
            "barely" => KassRaftery::Barely,
            "positive" => KassRaftery::Positive,
            "strong" => KassRaftery::Strong,
            "very-strong" => KassRaftery::VeryStrong,
            _ => panic!("bug: unexpected KassRaftery score"),
        };
        filtration::posterior_odds::filter_by_odds::<_, &str, &str>(
            None, None, &events, min_evidence
        )?;
    }

    Ok(())
}
