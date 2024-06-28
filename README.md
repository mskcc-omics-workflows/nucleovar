[![GitHub Actions CI Status](https://github.com/msk/nucleovar/workflows/nf-core%20CI/badge.svg)](https://github.com/msk/nucleovar/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/msk/nucleovar/workflows/nf-core%20linting/badge.svg)](https://github.com/msk/nucleovar/actions?query=workflow%3A%22nf-core+linting%22)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)

[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)


## Introduction

**msk/nucleovar** is a bioinformatics pipeline that ...

<!-- TODO nf-core:
   Processes a variety of sample BAM files through three variant callers (Mutect (v.1.1.5), VarDict, and GATK Mutect2). Output VCF Files are normalized, sorted, and concatenated, proceeding to be annotated and converted into a MAF format file. The following MAF file is tagged with the presence/absence of specific variant criteria, resulting in a final output MAF file containing variants filtered by criteria set forth by the ACCESS pipeline. 
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Read in core samplesheet (containing case and control samples) and auxillary bams samplesheet
2. Run case and control samples through variant callers (Mutect v.1.1.5, VarDict, and GATK Mutect2)
3. Take output VCF Files from each variant caller and normalize, sort, concatenate and annotate using BCFtools suite.
4. Convert output VCF File into MAF file format and annotate using Genome Nexus (option provided to invoke PERL VCF2MAF script as well)
5. Tag output MAF file using MSK ACCESS pipeline criteria (presence of hotspots and removal of specific variant annotations)
6. Run tagged MAF file through a traceback subworkflow which tags the file with presence of genotypes and performs specific tagging and concatentation.
7. Tag output file is run through filtering based on criteria set forth by ACCESS_filters script. 


## Usage

:::note
If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
with `-profile test` before running the workflow on actual data.
:::

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`core_samplesheet.csv`:

```csv
patient_id,sample_id,type,maf,duplex_bam,duplex_bai,simplex_bam,simplex_bai
PATIENT1,SAMPLE1,case,null,path/to/duplex.bam,path/to/duplex.bai,path/to/simplex.bam,path/to/simplex.bai
```
Each row represents an individual case and control sample.

`aux_bams_samplesheet.csv`:

```csv
sample_id,normal_path,duplex_path,simplex_path,type
SAMPLE1,/path/to/normal.bam,path/to/duplex.bam,path/to/simplex.bam,curated
```

Each row represents an individual sample which may contain a standard bam (if an unmatched or matched normal sample), or an individual sample which contains a simplex and duplex bam (if a curated or plasma sample)

-->

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run msk/nucleovar/main.nf \
   --input core_samplesheet.csv \
   --aux_bams aux_bams.csv \
   --rules_json rules.json \
   -profile <docker/singularity/.../institute> \
   --outdir <OUTDIR>
```

:::warning
Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).
:::

## Credits

msk/nucleovar was originally written by @rnaidu and @buehlere.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use msk/nucleovar for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
