extern crate libprosic;
extern crate rust_htslib;
#[macro_use]
extern crate log;
extern crate fern;
#[macro_use]
extern crate clap;

use clap::{App,AppSettings};

use std::process;

use rust_htslib::bam;

fn main() {
    // parse command line
    let yaml = load_yaml!("cli.yaml");
    let matches = App::from_yaml(yaml)
                      .version(env!("CARGO_PKG_VERSION"))
                      .global_settings(&[AppSettings::SubcommandRequired,
                                         AppSettings::ColoredHelp])
                      .get_matches();

    // setup logger
    let logger_config = fern::DispatchConfig {
        format: Box::new(|msg: &str, level: &log::LogLevel, _: &log::LogLocation| {
          match level {
              &log::LogLevel::Debug => format!("DEBUG: {}", msg),
              _ => msg.to_owned()
          }
        }),
        output: vec![fern::OutputConfig::stderr()],
        level: log::LogLevelFilter::Trace,
    };
    if let Err(e) = fern::init_global_logger(
        logger_config,
        if matches.is_present("verbose") { log::LogLevelFilter::Debug } else { log::LogLevelFilter::Info }
    ) {
        panic!("Failed to initialize logger: {}", e);
    }

    if let Some(matches) = matches.subcommand_matches("tumor-normal") {
        // read command line parameters
        let normal_mean_insert_size = value_t!(matches, "normal_mean_insert_size", f64).unwrap();
        let normal_sd_insert_size = value_t!(matches, "normal_sd_insert_size", f64).unwrap();
        let normal_heterozygosity = value_t!(matches, "normal_heterozygosity", f64).unwrap_or(0.001);
        let ploidy = value_t!(matches, "ploidy", u32).unwrap_or(2);
        let tumor_mean_insert_size = value_t!(matches, "tumor_mean_insert_size", f64).unwrap();
        let tumor_sd_insert_size = value_t!(matches, "tumor_sd_insert_size", f64).unwrap();
        let tumor_effective_mutation_rate = value_t!(matches, "tumor_effective_mutation_rate", f64).unwrap();
        let tumor_purity = value_t!(matches, "tumor_purity", f64).unwrap_or(1.0);
        let min_somatic_af = value_t!(matches, "min_somatic_af", f64).unwrap_or(0.05);
        let pileup_window = value_t!(matches, "pileup_window", u32).unwrap_or(2500);
        let normal = matches.value_of("normal").unwrap();
        let tumor = matches.value_of("tumor").unwrap();

        if let Ok(tumor_bam) = bam::IndexedReader::new(&tumor) {
            if let Ok(normal_bam) = bam::IndexedReader::new(&normal) {
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
                        genome_size,
                        tumor_purity,
                        normal_heterozygosity
                    )
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
                    )
                );

                // init joint model
                let mut joint_model = libprosic::JointModel::new(
                    libprosic::LatentVariableModel::new(tumor_purity),
                    libprosic::LatentVariableModel::new(1.0),
                    tumor_sample,
                    normal_sample
                );

                // setup events
                let events = [
                    libprosic::Event{ name: "GERMLINE".to_owned(), af_case: 0.0..1.0, af_control: vec![0.5, 1.0] },
                    libprosic::Event{ name: "SOMATIC".to_owned(), af_case: min_somatic_af..1.0, af_control: vec![0.0] }
                ];

                // perform calling
                if let Err(msg) = libprosic::call(&"-", &"-", &events, &mut joint_model) {
                    error!("{}", msg);
                    process::exit(1);
                }
            } else {
                error!("Error reading normal BAM.");
                process::exit(1);
            }
        } else {
            error!("Error reading tumor BAM.");
            process::exit(1);
        }
    }
}
