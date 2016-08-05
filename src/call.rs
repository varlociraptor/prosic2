use std::error::Error;

use clap;
use libprosic;
use rust_htslib::bam;

pub fn tumor_normal(matches: &clap::ArgMatches) -> Result<(), Box<Error+Send+Sync>> {
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
    let min_somatic_af = value_t!(matches, "min-somatic-af", f64).unwrap_or(0.05);
    let pileup_window = value_t!(matches, "pileup-window", u32).unwrap_or(2500);
    let normal = matches.value_of("normal").unwrap();
    let tumor = matches.value_of("tumor").unwrap();

    let tumor_bam = bam::IndexedReader::new(&tumor).expect(&format!("Error reading BAM file {}.", tumor));
    let normal_bam = bam::IndexedReader::new(&normal).expect(&format!("Error reading BAM file {}.", normal));
    let genome_size = (0..tumor_bam.header.target_count()).fold(0, |s, tid| {
        s + tumor_bam.header.target_len(tid).unwrap() as u64
    });

    // init tumor sample
    let tumor_sample = libprosic::Sample::new(
        tumor_bam,
        pileup_window,
        libprosic::InsertSize {
            mean: tumor_mean_insert_size,
            sd: tumor_sd_insert_size
        },
        libprosic::priors::TumorModel::new(
            ploidy,
            tumor_effective_mutation_rate,
            deletion_factor,
            insertion_factor,
            genome_size,
            tumor_purity,
            normal_heterozygosity
        ),
        libprosic::likelihood::LatentVariableModel::new(tumor_purity)
    );

    // init normal sample
    let normal_sample = libprosic::Sample::new(
        normal_bam,
        pileup_window,
        libprosic::InsertSize {
            mean: normal_mean_insert_size,
            sd: normal_sd_insert_size
        },
        libprosic::priors::InfiniteSitesNeutralVariationModel::new(
            ploidy,
            normal_heterozygosity
        ),
        libprosic::likelihood::LatentVariableModel::new(1.0)
    );

    // init joint model
    let mut joint_model = libprosic::model::ContinuousVsDiscreteModel::new(
        tumor_sample,
        normal_sample
    );

    // setup events
    let events = [
        libprosic::case_control::Event{
            name: "GERMLINE".to_owned(),
            af_case: 0.0..1.0,
            af_control: vec![0.5, 1.0]
        },
        libprosic::case_control::Event{
            name: "SOMATIC".to_owned(),
            af_case: min_somatic_af..1.0,
            af_control: vec![0.0]
        }
    ];

    // perform calling
    try!(libprosic::case_control::call(&"-", &"-", &events, &mut joint_model));

    Ok(())
}
