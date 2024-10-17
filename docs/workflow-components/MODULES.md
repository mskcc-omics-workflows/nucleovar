## Bcftools concat (v.1.15.1)

### description:

BCFtools is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. All commands work transparently with both VCFs and BCFs, both uncompressed and BGZF-compressed. Most commands accept VCF, bgzipped VCF and BCF with filetype detected automatically even when streaming from a pipe. Indexed VCF and BCF will work in all situations. Un-indexed VCF and BCF and streams will work in most, but not all situations. BCFtools Concat is a module that will concatenate multiple VCF file together in one format, even if individual VCFs belong to different variant callers.

_BCFtools concat is used in nucleovar to..._ concatenate the output VCF files from each variant caller into a final file that is converted and annotated into a MAF. BCFtools concat is also invoked to concatenate the standard and complexvariants VCF files together (produced by VarDict variant caller).

### run command:

```bash
bcftools concat \\
    --output ${prefix}_vardict_concat.vcf.gz \\
    --threads $task.cpus \\
    ${vcf1} ${vcf2}
```

### parameters:

**vcf1**: VCF file that is being concatenated
**vcf2**: VCF file that is being concatenated
**threads**: Use multithreading with INT worker threads. The option is currently used only for the compression of the output stream, only when --output-type is b or z. Default: 0.
**output**: output file name

## Bcftools annotate (v.1.15.1)

### description:

BCFtools is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. All commands work transparently with both VCFs and BCFs, both uncompressed and BGZF-compressed. Most commands accept VCF, bgzipped VCF and BCF with filetype detected automatically even when streaming from a pipe. Indexed VCF and BCF will work in all situations. Un-indexed VCF and BCF and streams will work in most, but not all situations. BCFtools Annotate is used to add and/or remove annotations from a given VCF file using a source "annotations" VCF file.

_BCFtools annotate is used in nucleovar to..._ annotate the output VCF files from each variant caller against the source variant caller VCF.

### run command:

```bash
bcftools annotate \\
--header-lines ${header_file} \\
--annotations ${singular_vcf} \\
-c CHROM,POS,REF,ALT,INFO/TYPE --mark-sites +MUTECT --output-type v \\
--output ${prefix}.vcf ${combined_vcf}

```

### parameters:

**header-lines**: Text file that includes information about the variant caller used to create the source VCF file to annotate with.
**annotations**: VCF file that is being used to annotate against the concatenated or core VCF
**c**: columns to include in the annotated final VCF
**--mark-sites**: Variant Caller to mark sites with
**--output-type**: Output file format, default is VCF
**--output**: Name of output file
**combined_vcf**: Source VCF file that is being annotated

## Bcftools norm (v.1.15.1)

### description:

BCFtools is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. All commands work transparently with both VCFs and BCFs, both uncompressed and BGZF-compressed. Most commands accept VCF, bgzipped VCF and BCF with filetype detected automatically even when streaming from a pipe. Indexed VCF and BCF will work in all situations. Un-indexed VCF and BCF and streams will work in most, but not all situations. BCFtools norm is used to normalize VCF files.

_BCFtools norm is used in nucleovar to..._ normalize output VCF files from the VarDict (standard variants), MuTect, and MuTect2 variant callers.

### run command:

```bash
bcftools norm \\
   --fasta-ref ${fasta} \\
   --output ${prefix}_norm.vcf.gz \\
   --threads $task.cpus \\
   ${vcf}
```

### parameters:

**--fasta-ref**: Reference Fasta file
**--output**: Name for output normalized VCF file
**--threads**: Use multithreading with INT worker threads. The option is currently used only for the compression of the output stream, only when --output-type is b or z. Default: 0.
**vcf**: Input VCF file that user wants normalized.

## Bcftools sort (v.1.15.1)

### description:

BCFtools is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. All commands work transparently with both VCFs and BCFs, both uncompressed and BGZF-compressed. Most commands accept VCF, bgzipped VCF and BCF with filetype detected automatically even when streaming from a pipe. Indexed VCF and BCF will work in all situations. Un-indexed VCF and BCF and streams will work in most, but not all situations. BCFTools sort is used to sort VCF files after normalization.

_BCFtools sort is used in nucleovar to..._ sort output VCF files from the VarDict (standard variants), MuTect, and MuTect2 variant callers.

### run command:

```bash
bcftools sort \\
   --output ${prefix}.vcf.gz \\
   --temp-dir . \\
   $vcf
```

### parameters:

**vcf**: Input VCF file that user wants sorted.
**--temp-dir**: Temporary directory to hold cached files.
**--output**: Output filename

## Bedtools genomecov (v.2.31.1)

### description:

Bedtools is a set of tools for genomic analysis tasks, specifically enabling genome arithmetic (merge, count, complement) on various file types. Bedtools genomecov computes histograms (default), per-base reports (-d) and BEDGRAPH (-bg) summaries of feature coverage (e.g., aligned sequences) for a given genome.

_Bedtools genomecov is used in nucleovar to..._ generate an intermediary bedgraph file for an input Tumor BAM, so that a BED file may be created for use in variant calling, if user doesn't specify a target BED file upon running nucleovar.

### run command:

```bash
bedtools genomecov -ibam ${sorted_tumor_bam} -bg > ${prefix}.bedgraph
```

### parameters:

**-ibam**: input BAM file
**-bg**: output bedgraph file

## Bedtools merge (v.2.31.1)

### description:

Bedtools is a set of tools for genomic analysis tasks, specifically enabling genome arithmetic (merge, count, complement) on various file types. Bedtools merge combines overlapping or “book-ended” features in an interval file into a single feature which spans all of the combined features.

_Bedtools merge is used in nucleovar to..._ convert an intermediary bedgraph file to BED file so that a BED file may be created for use in variant calling, if user doesn't specify a target BED file upon running nucleovar.

### run command:

```bash
bedtools merge -i ${prefix}.bedgraph > ${prefix}_target.bed
```

### parameters:

**-i**: input bedgraph file

## VardictJava (v1.5.1)

### description:

Java port of the VarDict variant discovery program Variant Caller.

_vardictjava is used in nucleovar to..._ generate a standard and complex variants VCF file for a pair of tumor and normal BAM files.

### run command:

```bash
export JAVA_OPTS='"-Xms${task.memory.toMega()/4}m" "-Xmx${task.memory.toGiga()}g" "-Dsamjdk.reference_fasta=${fasta}"'
   vardict-java \\
      -c 1 -S 2 -E 3 -g 4 \\
      ${input} \\
      -th ${task.cpus} \\
      -G ${fasta} \\
      ${bed} \\
      ${filter} \\
      ${convert_to_vcf} \\
      -f .0002 -r 1 -N \\
      "${meta.case_id}|${meta.control_id}" \\
      > ${prefix}.vcf
```

### parameters:

**input**: input set of tumor and normal BAM files with their respective index.
**-th**: thread count
**-G**: input reference FASTA file
**bed**: input BED file
**filter**: flag to determine if VarDict Java is to be run in single or paired-sample mode. default is paired-sample.
**convert_to_vcf**: flag to determine if VarDict Java is to be run in single or paired-sample mode. default is paired-sample.

## MuTect (v.1.5.1)

### description:

This module wraps MuTect v.1.1.5 software, which generates VCF files and standard output from case-control sample BAM files, BED file coordinates, and a reference fasta.

_mutect1 is used in nucleovar to..._ generate a VCF file for a pair of tumor and normal BAM files.

### run command:

```bash
java -Xmx28g -Xms256m -XX:-UseGCOverheadLimit -jar /opt/mutect/muTect-1.1.5.jar -T MuTect \
   --downsample_to_coverage 50000 --downsampling_type NONE \
   --enable_extended_output --fraction_contamination 0.0005 \
   --minimum_mutation_cell_fraction 0.0005  \
   --read_filter BadCigar \
   --input_file:tumor ${case_bam} \
   --input_file:normal ${control_bam} \
   --intervals ${bed_file} \
   --tumor_sample_name ${case_sample_name} \
   --reference_sequence ${fasta_file} \
   --normal_sample_name ${control_sample_name} \
   ${args2} \
   --out ${case_sample_name}.${control_sample_name}.mutect.txt \
   --vcf ${case_sample_name}.${control_sample_name}.mutect.vcf
```

### parameters:

**input_file** specifies input BAM file, either deriving from the tumor or normal sample.
**intervals** specifies input BED file
**reference_sequence**: specifies input reference FASTA file
**tumor_sample_name**: tumor sample name
**normal_sample_name**: normal sample name
**out**: filename for output MuTect text file
**vcf**: filename for output MuTect VCF file

## MuTect2 (called by GATK v.4.5.0)

### description:

This module wraps MuTect2 software called by GATK v.4.5.0, which generates VCF files and standard output from case-control sample BAM files, BED file coordinates, and a reference fasta.

_mutect2 is used in nucleovar to..._ generate a VCF file for a pair of tumor and normal BAM files.

### run command:

```bash
gatk "Mutect2" \\
   -R ${fasta_file} \\
   -I ${case_bam} \\
   -I ${control_bam} \\
   -tumor ${case_sample_name} \\
   -normal ${control_sample_name} \\
   --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \\
   --minimum-allele-fraction .0002 \\
   --output ${meta.id}.mutect2.vcf
```

### parameters:

**-R**: input reference FASTA
**-I**: Tumor sample BAM file
**-I**: Normal sample BAM file
**-tumor**: Tumor sample name to match BAM file with
**-normal**: Normal sample name to match BAM file with
**--output**: Output filename

## VarDict Filter (called by pv v.0.2.7)

### description:

Post-processing command to filter a case-controlled VarDict version VCF file."

_vardict filter is used in nucleovar to..._ filter a VarDict VCF file after generation.

### run command:

```bash
 pv vardict case-control filter \
   --inputVcf ${vardict_vcf_file} \
   --tsampleName ${meta.case_id} \
   --alleledepth 1 --minQual 0 \
   --totalDepth 20 --tnRatio 1 \
   --variantFraction 5e-05
```

### parameters:

**inputVcf**: input VarDict VCF file
**tsampleName**: Tumor Sample Name

## MuTect Filter (called by pv v.0.2.7)

### description:

Post-processing command to filter a case-controlled MuTect1 VCF file.

_mutect1 filter is used in nucleovar to..._ filter a MuTect1 VCF file after generation.

### run command:

```bash
pv mutect1 case-control filter \
--inputVcf ${mutect_vcf_file} \
--inputTxt ${mutect_txt_file} \
--refFasta ${reference_fasta} \
--tsampleName ${meta.case_id} \
--outDir .
```

### parameters:

**inputVcf**: input MuTect1 VCF file
**inputTxt**: input MuTect1 Text file
**refFasta**: input reference FASTA file
**tsampleName**: tumor sample name
**outDir**: filepath to output directory

## MuTect2 Filter (called by pv v.0.2.7)

### description:

Post-processing command to filter a case-controlled MuTect2 VCF file.

_mutect2 filter is used in nucleovar to..._ filter a MuTect2 VCF file after generation.

### run command:

```bash
pv mutect2 case-control filter \
--inputVcf ${mutect_vcf_file} \
--tsampleName ${meta.case_id} \
--outDir .
```

### parameters:

**inputVcf**: input MuTect2 VCF file
**tsampleName**: tumor sample name
**outDir**: filepath to output directory

## BCFtools index (v.1.15.1)

### description:

BCFtools is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. All commands work transparently with both VCFs and BCFs, both uncompressed and BGZF-compressed. Most commands accept VCF, bgzipped VCF and BCF with filetype detected automatically even when streaming from a pipe. Indexed VCF and BCF will work in all situations. Un-indexed VCF and BCF and streams will work in most, but not all situations. BCFtools index is a module that will create an index file for any input VCF.

_bcftools index is used in nucleovar to..._ Create index files for an input VCF.

### run command:

```bash
bcftools index \\
   --threads $task.cpus \\
   ${vcf}.gz
```

### parameters:

**threads**: Use multithreading with INT worker threads. The option is currently used only for the compression of the output stream, only when --output-type is b or z. Default: 0.
**vcf**: input VCF file. Must be in compressed file format.

## GenomeNexus VCF2MAF-lite (v.0.0.1)

### description:

DSL2 module to perform a native conversion of VCF file type to MAF file type using the genomenexus suite of conversion and annotation tools.

_genomenexus vcf2maf is used in nucleovar to..._ Convert an input VCF file to MAF file format.

### run command:

```python3
/vcf2maf-lite/vcf2maf_lite.py -i ${vcf} \
--retain-info AC,AF,AC_nfe_seu,AC_afr,AF_afr,MUTECT,set,TYPE,FAILURE_REASON
```

### parameters:

**-i**: input VCF file

## GenomeNexus AnnotationPipeline (v.0.0.1)

### description:

DSL2 module to perform annotation of a MAF file type using the genomenexus suite of conversion and annotation tools.

_genomenexus annotationpipeline is used in nucleovar to..._ Annotate an input MAF file.

### run command:

```bash
java -jar /genome-nexus-annotation-pipeline/annotationPipeline/target/annotationPipeline.jar \\
--filename ${input_maf} \\
--output-filename ${meta.id}_annotated.maf
```

### parameters:

**filename**: input pre-annotated MAF file
**output-filename**: output post-annotated MAF file

## pv maf concat (v.0.2.7)

### description:

a flexible command for concatenating maf files

_pv maf concat is used in nucleovar to..._ Concatenate MAF files together in the traceback subworkflow.

### run command:

```bash
pv maf concat \\
   $flagFiles \\
   --output $output \\
```

### parameters:

**$flagFiles**: variable corresponding to a list of MAF files which need to be concatenated.
**--output**: output concatenated MAF filename.

## genotype_variants_all (v.0.3.9)

### description:

module supports genotyping and merging small variants (SNV and INDELS).

_genotypevariantsall is used in nucleovar to..._ Genotype individual MAF files in the traceback subworkflow.

### run command:

```bash
genotype_variants small_variants all \\
   -i ${maf} \\
   -r ${fasta} \\
   -g /usr/local/bin/GetBaseCountsMultiSample \\
   $patient \\
   $bams_standard \\
   $bam_liquid \\
   $sample \\
   --tumor_name_override
```

### parameters:

**-i**: input maf file
**-r**: input reference fasta file
**-g**: filepath to GBCMS script
**patient**: patient ID
**bams_standard**: standard BAMs (tumor and normal sample)
**bam_liquid**: auxiliary BAMS (curated,plasma,matched_normal,unmatched_normal,standard)
**sample**: sample ID
**tumor_name_override**: flag to make sure correct sample IDs are being written in output MAF for normal and tumor sample.

## pv maf tagtraceback (v.0.2.7)

### description:

a flexible command for tagging maf files

_pv maf tagtraceback is used in nucleovar to..._ Assign a sample type to each output MAF from traceback subworkflow, so they can be identified in downstream modules.

### run command:

```bash
 pv maf tag traceback \\
   -m $maf \\
   $sampleFiles \\
   --output $output
```

### parameters:

**-m**: input MAF file to be tagged
**sampleFiles**: samplesheets to pull sample type metadata from
**output**: output tagged MAF file name

## maf_processing (v.0.2.7)

### description:

Module that is run downstream of traceback subworkflow which tags an input MAF file based on a separate "rules" JSON input file, as well as a separate hotspots input MAF file.

### run command:

```bash
pv maf tag access \\
--maf ${genotyped_maf} \\
--rules ${rules_file} \\
--hotspots ${hotspots} \\
--output ${meta.id}_tagged.maf
```

### parameters:

**maf**: input MAF file
**rules**: input rules JSON file to tag MAF by
**hotspots**: input hotspots MAF file to tag MAF by
**output**: output filename

## access_filters (v.0.2.7)

### description:

Module run to tag input MAF file with a series of filtering criteria set forth by ACCESS.

### run command:

```bash
pv maf filter access_filters \
-f ${traceback_maf} \
-a ${anno_maf} \
-ts ${meta.case_id}  \
-ns ${meta.control_id} \
-bl ${blocklist} \
--output ${meta.id}
```

### parameters:

**-f**: Fillout MAF file to subset (direct output from traceback subworkflow)
**-a**: Annotated MAF file to subset (direct input file from beginning of traceback subworkflow)
**-ts**: tumor sample name
**-ns**: normal sample name
**-bl**: Optional input blocklist file for access filtering criteria.
**output**: output filename

## tag_by_variant_annotation (v.0.2.7)

### description:

Module to tag output MAF file by criteria set forth in mPATH. Also splits MAF file up into output text files for each category of variants (annotated,silent,dropped,nonpanel_annotated,nonpanel_silent)

### run command:

```bash
pv maf tag by_variant_classification \
--maf ${access_filtered_maf} \
--canonical_tx_ref ${canonical_tx_ref} \
--output_dir .
```

### parameters:

**maf**:input MAF file (direct output from access filters)
**canonical_tx_ref**: Text file containing canonical transcript references
**output_dir**: final output directory
