use std::error::Error;

use clap;
use libprosic;
use libprosic::model::AlleleFreq;
use rust_htslib::bam;

pub fn tumor_normal(matches: &clap::ArgMatches) -> Result<(), Box<Error>> {
    // read command line parameters
    let tumor_mean_insert_size = value_t!(matches, "insert-size-mean", f64).unwrap();
    let tumor_sd_insert_size = value_t!(matches, "insert-size-sd", f64).unwrap();
    let normal_mean_insert_size = value_t!(matches, "normal-insert-size-mean", f64).unwrap_or(tumor_mean_insert_size);
    let normal_sd_insert_size = value_t!(matches, "normal-insert-size-sd", f64).unwrap_or(tumor_sd_insert_size);
    let normal_heterozygosity = value_t!(matches, "heterozygosity", f64).unwrap_or(1.25E-4);
    let ploidy = value_t!(matches, "ploidy", u32).unwrap_or(2);
    let tumor_effective_mutation_rate = value_t!(matches, "effective-mutation-rate", f64).unwrap();
    let deletion_factor = value_t!(matches, "deletion-factor", f64).unwrap_or(0.03);
    let insertion_factor = value_t!(matches, "insertion-factor", f64).unwrap_or(0.01);
    let tumor_purity = value_t!(matches, "purity", f64).unwrap_or(1.0);
    let min_somatic_af = value_t!(matches, "min-somatic-af", f64).map(|af| AlleleFreq(af)).unwrap_or(AlleleFreq(0.05));
    let pileup_window = value_t!(matches, "pileup-window", u32).unwrap_or(2500);
    let no_fragment_evidence = matches.is_present("omit-fragment-evidence");
    let no_secondary = matches.is_present("omit-secondary-alignments");
    let normal = matches.value_of("normal").unwrap();
    let tumor = matches.value_of("tumor").unwrap();
    let candidates = matches.value_of("candidates").unwrap_or("-");
    let output = matches.value_of("output").unwrap_or("-");

    let tumor_bam = try!(bam::IndexedReader::new(&tumor));
    let normal_bam = try!(bam::IndexedReader::new(&normal));
    let genome_size = (0..tumor_bam.header.target_count()).fold(0, |s, tid| {
        s + tumor_bam.header.target_len(tid).unwrap() as u64
    });

    // init tumor sample
    let tumor_sample = libprosic::Sample::new(
        tumor_bam,
        pileup_window,
        !no_fragment_evidence,
        !no_secondary,
        libprosic::InsertSize {
            mean: tumor_mean_insert_size,
            sd: tumor_sd_insert_size
        },
        libprosic::likelihood::LatentVariableModel::new(tumor_purity)
    );

    // init normal sample
    let normal_sample = libprosic::Sample::new(
        normal_bam,
        pileup_window,
        !no_fragment_evidence,
        !no_secondary,
        libprosic::InsertSize {
            mean: normal_mean_insert_size,
            sd: normal_sd_insert_size
        },
        libprosic::likelihood::LatentVariableModel::new(1.0)
    );


    let prior_model = libprosic::priors::TumorNormalModel::new(
        ploidy,
        tumor_effective_mutation_rate,
        deletion_factor,
        insertion_factor,
        genome_size,
        tumor_purity,
        normal_heterozygosity
    );

    // init joint model
    let mut joint_model = libprosic::model::ContinuousVsDiscreteModel::new(
        tumor_sample,
        normal_sample,
        prior_model
    );

    // setup events
    let events = [
        libprosic::case_control::CaseControlEvent {
            name: "germline".to_owned(),
            af_case: AlleleFreq(0.0)..AlleleFreq(1.0),
            af_control: vec![AlleleFreq(0.5), AlleleFreq(1.0)]
        },
        libprosic::case_control::CaseControlEvent {
            name: "somatic".to_owned(),
            af_case: min_somatic_af..AlleleFreq(1.0),
            af_control: vec![AlleleFreq(0.0)]
        }
    ];
    // call absent variants as the complement of the other events
    let absent_event = libprosic::ComplementEvent { name: "absent".to_owned() };

    // perform calling
    try!(libprosic::case_control::call(
        &candidates,
        &output,
        &events,
        Some(&absent_event),
        &mut joint_model
    ));

    Ok(())
}
