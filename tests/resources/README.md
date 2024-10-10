[![GitHub Actions CI Status](https://github.com/msk/nucleovar/workflows/nf-core%20CI/badge.svg)](https://github.com/msk/nucleovar/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/msk/nucleovar/workflows/nf-core%20linting/badge.svg)](https://github.com/msk/nucleovar/actions?query=workflow%3A%22nf-core+linting%22)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

## Resources

**stub** houses files relevant to stub/dummy testing of nucleovar.
This includes testing framework for basic functionality of nucleovar.

**v1.0** houses files relevant to end to end testing of nucleovar.
This includes testing framework for basic functionality of nucleovar plus archival patient data obtained
from previous voyager runs. Due to data privacy, this data is located on terra (filepaths provided in the test.config) and only non-patient data files are stored in this folder.

## Files in **stub**

1.`stub_input_samplesheet.csv`: 1.`stub_aux_bams.csv`:

## Files in **v1.0**

1.`input_samplesheet.csv`: 1.`stub_aux_bams.csv`:

## Usage

invoke the **stub** profile data using `-profile stub -stub`
invoke the **v1.0** profile data using `-profile test`

## RUN COMMAND FOR STUB

```bash
nextflow run msk/nucleovar/main.nf \
   --input stub_input_samplesheet.csv \
   --aux_bams stub_aux_bams.csv \
   --rules_json stub_rules.json \
   --fasta stub_ref.fasta \
   --fai stub_ref.fasta.fai \
   --dict stub_ref.dict \
   --canonical_bed stub_canonical.bed \
   --target_bed stub_target.bed \
   --blocklist stub_blocklist.txt \
   --canonical_tx_ref stub_canonical_tx_ref.tsv \
   --hotspots stub_hotspots.maf \
   --annotator genomenexus \
   -profile <docker/singularity/.../institute> \
   --outdir <OUTDIR>
```

## RUN COMMAND FOR V1.0

```bash
nextflow run msk/nucleovar/main.nf \
   --input input_samplesheet.csv \
   --aux_bams aux_bams.csv \
   --rules_json rules.json \
   --fasta ref.fasta \
   --fai ref.fasta.fai \
   --dict ref.dict \
   --canonical_bed canonical.bed \
   --target_bed target.bed \
   --blocklist blocklist.txt \
   --canonical_tx_ref canonical_tx_ref.tsv \
   --hotspots hotspots.maf \
   --annotator genomenexus \
   -profile <docker/singularity/.../institute> \
   --outdir <OUTDIR>
```
