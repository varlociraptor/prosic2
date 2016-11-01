# Postprocessing for Somatic Mutation Calling (PROSIC)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/prosic/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/prosic/badges/downloads.svg)](http://bioconda.github.io/recipes/prosic/README.html)

PROSIC is a caller for somatic variants in tumor-normal sample pairs, sequenced with any next-generation sequencing technology.
It provides a novel latent variable model that integrates various levels of uncertainty, and thereby allows to properly asses the probability of having a somatic variant while controlling the false discovery rate.

## Installation

PROSIC is available via [Bioconda](https://bioconda.github.io), a distribution
of bioinformatics software for the conda package manager.
Bioconda can be set up in any Linux environment, even without admin rights.
With Bioconda set up, PROSIC can be installed via

	$ conda install prosic

## Usage

The purpose of PROSIC is to call somatic insertions and deletions (indels) on tumor/normal sample pairs.
For this, PROSIC requires a VCF file with preliminary indel calls, e.g. obtained with [Delly](https://github.com/tobiasrausch/delly) or [Lancet](https://github.com/nygenome/lancet).
Then, calling with PROSIC consists of two steps.

### Step 1: Calling

Variants are called by applying PROSIC to the preliminary calls, i.e.

    $ prosic tumor-normal tumor.bam normal.bam < pre-calls.vcf > prosic-calls.bcf

PROSIC then annotates the initial calls with probabilities for the events somatic, germline and absent (`PROB_SOMATIC`, `PROB_GERMLINE`, `PROB_ABSENT`).
Issue `prosic tumor-normal --help` for information about additional parameters.

### Step 2: Controlling FDR

Second, the FDR (here for somatic deletions) can be controlled by

	$ prosic control-fdr --event SOMATIC --var DEL < prosic-calls.bcf > thresholds.txt

The resulting tab-separated table `thresholds.txt` contains thresholds that can be applied to the corresponding `PROB_SOMATIC` field in `prosic-calls.bcf`, in order to control the FDR at different levels.

# Authors

* Original model: [Louis Dijkstra](https://github.com/louisdijkstra)
* Extended model and implementation: [Johannes Köster](https://johanneskoester.bitbucket.org)