use std::error::Error;

use clap;
use libprosic;
use libprosic::model::AlleleFreq;
use rust_htslib::bam;
use bio::stats::Prob;

pub fn tumor_normal(matches: &clap::ArgMatches) -> Result<(), Box<Error>> {
    // read command line parameters
    let tumor_mean_insert_size = value_t!(matches, "insert-size-mean", f64).unwrap();
    let tumor_sd_insert_size = value_t!(matches, "insert-size-sd", f64).unwrap();
    let normal_mean_insert_size = value_t!(matches, "normal-insert-size-mean", f64).unwrap_or(tumor_mean_insert_size);
    let normal_sd_insert_size = value_t!(matches, "normal-insert-size-sd", f64).unwrap_or(tumor_sd_insert_size);
    let normal_heterozygosity = try!(Prob::checked(value_t!(matches, "heterozygosity", f64).unwrap_or(1.25E-4)));
    let ploidy = value_t!(matches, "ploidy", u32).unwrap_or(2);
    let tumor_effective_mutation_rate = value_t!(matches, "effective-mutation-rate", f64).unwrap_or(2000.0);
    let deletion_factor = value_t!(matches, "deletion-factor", f64).unwrap_or(0.03);
    let insertion_factor = value_t!(matches, "insertion-factor", f64).unwrap_or(0.01);
    let tumor_purity = value_t!(matches, "purity", f64).unwrap_or(1.0);
    let min_somatic_af = value_t!(matches, "min-somatic-af", f64).map(|af| AlleleFreq(af)).unwrap_or(AlleleFreq(0.05));
    let pileup_window = value_t!(matches, "pileup-window", u32).unwrap_or(2500);
    let no_fragment_evidence = matches.is_present("omit-fragment-evidence");
    let no_secondary = matches.is_present("omit-secondary-alignments");
    let no_mapq = matches.is_present("omit-mapq");
    let adjust_mapq = matches.is_present("adjust-mapq");
    let omit_snvs = matches.is_present("omit-snvs");
    let omit_indels = matches.is_present("omit-indels");
    let normal = matches.value_of("normal").unwrap();
    let tumor = matches.value_of("tumor").unwrap();
    let candidates = matches.value_of("candidates").unwrap_or("-");
    let output = matches.value_of("output").unwrap_or("-");
    let reference = matches.value_of("reference").unwrap();
    let observations = matches.value_of("observations");
    let flat_priors = matches.is_present("flat-priors");
    let exclusive_end = matches.is_present("exclusive-end");

    let prob_spurious_isize = try!(Prob::checked(value_t!(matches, "prob-spurious-isize", f64).unwrap_or(0.0)));

    let max_indel_dist = value_t!(matches, "max-indel-dist", u32).unwrap_or(50);
    let max_indel_len_diff = value_t!(matches, "max-indel-len-diff", u32).unwrap_or(20);
    let max_indel_len = value_t!(matches, "max-indel-len", u32).unwrap_or(1000);

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
        !no_mapq,
        adjust_mapq,
        libprosic::InsertSize {
            mean: tumor_mean_insert_size,
            sd: tumor_sd_insert_size
        },
        libprosic::likelihood::LatentVariableModel::new(tumor_purity),
        prob_spurious_isize
    ).max_indel_dist(max_indel_dist).max_indel_len_diff(max_indel_len_diff);

    // init normal sample
    let normal_sample = libprosic::Sample::new(
        normal_bam,
        pileup_window,
        !no_fragment_evidence,
        !no_secondary,
        !no_mapq,
        adjust_mapq,
        libprosic::InsertSize {
            mean: normal_mean_insert_size,
            sd: normal_sd_insert_size
        },
        libprosic::likelihood::LatentVariableModel::new(1.0),
        prob_spurious_isize
    ).max_indel_dist(max_indel_dist).max_indel_len_diff(max_indel_len_diff);

    // setup events
    let events = [
        libprosic::call::pairwise::PairEvent {
            name: "germline".to_owned(),
            af_case: AlleleFreq(0.0)..AlleleFreq(1.0),
            af_control: vec![AlleleFreq(0.5), AlleleFreq(1.0)]
        },
        libprosic::call::pairwise::PairEvent {
            name: "somatic".to_owned(),
            af_case: min_somatic_af..AlleleFreq(1.0),
            af_control: vec![AlleleFreq(0.0)]
        }
    ];
    // call absent variants as the complement of the other events
    let absent_event = libprosic::ComplementEvent { name: "absent".to_owned() };

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
            &candidates,
            &output,
            &reference,
            &events,
            Some(&absent_event),
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
            &candidates,
            &output,
            &reference,
            &events,
            Some(&absent_event),
            &mut joint_model,
            omit_snvs,
            omit_indels,
            Some(max_indel_len),
            observations.as_ref(),
            exclusive_end
        )
    }
}


pub fn normal_normal(matches: &clap::ArgMatches) -> Result<(), Box<Error>> {
    // read command line parameters
    let mean_insert_size = value_t!(matches, "insert-size-mean", f64).unwrap();
    let sd_insert_size = value_t!(matches, "insert-size-sd", f64).unwrap();
    let heterozygosity = try!(Prob::checked(value_t!(matches, "heterozygosity", f64).unwrap_or(1.25E-4)));
    let ploidy = value_t!(matches, "ploidy", u32).unwrap_or(2);
    let pileup_window = value_t!(matches, "pileup-window", u32).unwrap_or(2500);
    let no_fragment_evidence = matches.is_present("omit-fragment-evidence");
    let no_secondary = matches.is_present("omit-secondary-alignments");
    let no_mapq = matches.is_present("omit-mapq");
    let adjust_mapq = matches.is_present("adjust-mapq");
    let omit_snvs = matches.is_present("omit-snvs");
    let omit_indels = matches.is_present("omit-indels");
    let first = matches.value_of("first").unwrap();
    let second = matches.value_of("second").unwrap();
    let candidates = matches.value_of("candidates").unwrap_or("-");
    let output = matches.value_of("output").unwrap_or("-");
    let reference = matches.value_of("reference").unwrap();
    let observations = matches.value_of("observations");
    let flat_priors = matches.is_present("flat-priors");
    let exclusive_end = matches.is_present("exclusive-end");

    let prob_spurious_isize = try!(Prob::checked(value_t!(matches, "prob-spurious-isize", f64).unwrap_or(0.0)));

    let max_indel_dist = value_t!(matches, "max-indel-dist", u32).unwrap_or(50);
    let max_indel_len_diff = value_t!(matches, "max-indel-len-diff", u32).unwrap_or(20);
    let max_indel_len = value_t!(matches, "max-indel-len", u32).unwrap_or(1000);

    let first_bam = try!(bam::IndexedReader::new(&first));
    let second_bam = try!(bam::IndexedReader::new(&second));
    let genome_size = (0..first_bam.header.target_count()).fold(0, |s, tid| {
        s + first_bam.header.target_len(tid).unwrap() as u64
    });
    let insert_size = libprosic::InsertSize {
        mean: mean_insert_size,
        sd: sd_insert_size
    };

    // init normal sample
    let init_sample = |bam| {
        libprosic::Sample::new(
            bam,
            pileup_window,
            !no_fragment_evidence,
            !no_secondary,
            !no_mapq,
            adjust_mapq,
            insert_size,
            libprosic::likelihood::LatentVariableModel::new(1.0),
            prob_spurious_isize,
        ).max_indel_dist(max_indel_dist).max_indel_len_diff(max_indel_len_diff)
    };

    let first_sample = init_sample(first_bam);
    let second_sample = init_sample(second_bam);

    // setup events
    let events = [
        libprosic::call::pairwise::PairEvent {
            name: "both".to_owned(),
            af_case: vec![AlleleFreq(0.5), AlleleFreq(1.0)],
            af_control: vec![AlleleFreq(0.5), AlleleFreq(1.0)]
        },
        libprosic::call::pairwise::PairEvent {
            name: "first".to_owned(),
            af_case: vec![AlleleFreq(0.5), AlleleFreq(1.0)],
            af_control: vec![AlleleFreq(0.0)]
        },
        libprosic::call::pairwise::PairEvent {
            name: "second".to_owned(),
            af_case: vec![AlleleFreq(0.0)],
            af_control: vec![AlleleFreq(0.5), AlleleFreq(1.0)]
        }
    ];
    // call absent variants as the complement of the other events
    let absent_event = libprosic::ComplementEvent { name: "absent".to_owned() };

    if !flat_priors {
        panic!("only flat priors are supported for now (use --flat-priors)");
    } else {
        let prior_model = libprosic::priors::FlatNormalNormalModel::new(ploidy);

        // init joint model
        let mut joint_model = libprosic::model::PairCaller::new(
            first_sample,
            second_sample,
            prior_model
        );

        // perform calling
        libprosic::call::pairwise::call::<
            _, _, _,
            libprosic::model::PairCaller<
                libprosic::model::DiscreteAlleleFreqs,
                libprosic::model::DiscreteAlleleFreqs,
                libprosic::model::priors::FlatNormalNormalModel
            >, _, _, _, _>
        (
            &candidates,
            &output,
            &reference,
            &events,
            Some(&absent_event),
            &mut joint_model,
            omit_snvs,
            omit_indels,
            Some(max_indel_len),
            observations.as_ref(),
            exclusive_end
        )
    }
}
