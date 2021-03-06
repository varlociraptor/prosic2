use std::error::Error;

use clap;
use libprosic;
use libprosic::model::{ContinuousAlleleFreqs, DiscreteAlleleFreqs};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use bio::stats::Prob;


fn path_or_pipe(arg: Option<&str>) -> Option<&str> {
    arg.map_or(None, |f| if f == "-" { None } else { Some(f) })
}


fn alignment_properties(path: &str) -> Result<libprosic::AlignmentProperties, Box<Error>> {
    let mut bam = bam::Reader::from_path(path)?;
    Ok(libprosic::AlignmentProperties::estimate(&mut bam)?)
}


pub fn tumor_normal(matches: &clap::ArgMatches) -> Result<(), Box<Error>> {
    let normal_heterozygosity = Prob::checked(value_t!(matches, "heterozygosity", f64).unwrap())?;
    let ploidy = value_t!(matches, "ploidy", u32).unwrap();
    let tumor_effective_mutation_rate = value_t!(matches, "effective-mutation-rate", f64).unwrap();
    let deletion_factor = value_t!(matches, "deletion-factor", f64).unwrap();
    let insertion_factor = value_t!(matches, "insertion-factor", f64).unwrap();
    let tumor_purity = value_t!(matches, "purity", f64).unwrap();
    let pileup_window = value_t!(matches, "pileup-window", u32).unwrap();
    let no_fragment_evidence = matches.is_present("omit-fragment-evidence");
    let omit_snvs = matches.is_present("omit-snvs");
    let omit_indels = matches.is_present("omit-indels");
    let normal = matches.value_of("normal").unwrap();
    let tumor = matches.value_of("tumor").unwrap();
    let candidates = path_or_pipe(matches.value_of("candidates"));
    let output = path_or_pipe(matches.value_of("output"));
    let reference = matches.value_of("reference").unwrap();
    let observations = matches.value_of("observations");
    let flat_priors = matches.is_present("flat-priors");
    let exclusive_end = matches.is_present("exclusive-end");
    let indel_haplotype_window = value_t!(matches, "indel-window", u32).unwrap();

    let prob_spurious_ins = Prob::checked(value_t_or_exit!(matches, "prob-spurious-ins", f64))?;
    let prob_spurious_del = Prob::checked(value_t_or_exit!(matches, "prob-spurious-del", f64))?;
    let prob_ins_extend = Prob::checked(value_t_or_exit!(matches, "prob-ins-extend", f64))?;
    let prob_del_extend = Prob::checked(value_t_or_exit!(matches, "prob-del-extend", f64))?;

    let max_indel_len = value_t!(matches, "max-indel-len", u32).unwrap();


    let tumor_alignment_properties = alignment_properties(&tumor)?;
    let normal_alignment_properties = alignment_properties(&normal)?;

    info!("estimated tumor properties: {:?}", tumor_alignment_properties);
    info!("estimated normal properties: {:?}", normal_alignment_properties);


    let tumor_bam = bam::IndexedReader::from_path(&tumor)?;
    let normal_bam = bam::IndexedReader::from_path(&normal)?;
    let genome_size = (0..tumor_bam.header().target_count()).fold(0, |s, tid| {
        s + tumor_bam.header().target_len(tid).unwrap() as u64
    });

    // init tumor sample
    let tumor_sample = libprosic::Sample::new(
        tumor_bam,
        pileup_window,
        !no_fragment_evidence,
        false,
        false,
        tumor_alignment_properties,
        libprosic::likelihood::LatentVariableModel::new(tumor_purity),
        prob_spurious_ins,
        prob_spurious_del,
        prob_ins_extend,
        prob_del_extend,
        indel_haplotype_window
    );

    // init normal sample
    let normal_sample = libprosic::Sample::new(
        normal_bam,
        pileup_window,
        !no_fragment_evidence,
        false,
        false,
        normal_alignment_properties,
        libprosic::likelihood::LatentVariableModel::new(1.0),
        prob_spurious_ins,
        prob_spurious_del,
        prob_ins_extend,
        prob_del_extend,
        indel_haplotype_window
    );

    // setup events
    let events = [
        libprosic::call::pairwise::PairEvent {
            name: "germline".to_owned(),
            af_case: ContinuousAlleleFreqs::inclusive(0.0..1.0),
            af_control: DiscreteAlleleFreqs::feasible(ploidy).not_absent()
        },
        libprosic::call::pairwise::PairEvent {
            name: "somatic".to_owned(),
            af_case: ContinuousAlleleFreqs::left_exclusive(0.0..1.0),
            af_control: DiscreteAlleleFreqs::absent()
        },
        libprosic::call::pairwise::PairEvent {
            name: "absent".to_owned(),
            af_case: ContinuousAlleleFreqs::inclusive(0.0..0.0),
            af_control: DiscreteAlleleFreqs::absent()
        }
    ];

    if !flat_priors {
        let prior_model = libprosic::priors::TumorNormalModel::new(
            ploidy,
            tumor_effective_mutation_rate,
            deletion_factor,
            insertion_factor,
            genome_size,
            normal_heterozygosity
        );

        // init joint model
        let mut joint_model = libprosic::model::PairCaller::new(
            tumor_sample,
            normal_sample,
            prior_model
        );

        // perform calling
        libprosic::call::pairwise::call::<
            _, _, _,
            libprosic::model::PairCaller<
                libprosic::model::ContinuousAlleleFreqs,
                libprosic::model::DiscreteAlleleFreqs,
                libprosic::model::priors::TumorNormalModel
            >, _, _, _, _>
        (
            candidates,
            output,
            &reference,
            &events,
            &mut joint_model,
            omit_snvs,
            omit_indels,
            Some(max_indel_len),
            observations.as_ref(),
            exclusive_end
        )
    } else {
        let prior_model = libprosic::priors::FlatTumorNormalModel::new(ploidy);

        // init joint model
        let mut joint_model = libprosic::model::PairCaller::new(
            tumor_sample,
            normal_sample,
            prior_model
        );

        // perform calling
        libprosic::call::pairwise::call::<
            _, _, _,
            libprosic::model::PairCaller<
                libprosic::model::ContinuousAlleleFreqs,
                libprosic::model::DiscreteAlleleFreqs,
                libprosic::model::priors::FlatTumorNormalModel
            >, _, _, _, _>
        (
            candidates,
            output,
            &reference,
            &events,
            &mut joint_model,
            omit_snvs,
            omit_indels,
            Some(max_indel_len),
            observations.as_ref(),
            exclusive_end
        )
    }
}
