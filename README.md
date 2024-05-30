# mutect2ConsensusTumorOnly

The Mutect2Consensus workflow will process umiConsensus outputs for the tumour data through mutect2 in tumour only mode to call variants then use information from the matched normal to identify likely germline variants.

## Overview

## Dependencies

* [gatk 3.6-0](https://gatk.broadinstitute.org)
* [python 3.9](https://www.python.org/downloads/)
* [vep 105.0](https://useast.ensembl.org/info/docs/tools/vep/)
* [gatk 4.1.6.0](https://gatk.broadinstitute.org/)
* [tabix 0.2.6](https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2/download)
* [vcf2maf 1.6](https://github.com/mskcc/vcf2maf)
* [pandas 1.4.2](https://pandas.pydata.org/)


## Usage

### Cromwell
```
java -jar cromwell.jar run mutect2ConsensusTumorOnly.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`tumorInputGroup`|InputGroup|partitioned bam files from umiConsensus outputs for tumor sample
`outputFileNamePrefix`|String|Prefix to use for output file
`intervalFile`|String|interval file to subset variant calls
`inputIntervalsToParalellizeBy`|String|intervals for parallelization
`tumorName`|String|Name of the tumor sample
`reference`|String|reference version
`combineVariants.workflows`|Array[String]|array of ids of producer workflows


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`mutect2.filter_timeout`|Int|12|Hours before task timeout
`mutect2.filter_memory`|Int|16|Memory allocated for job
`mutect2.filter_filterExtraArgs`|String?|None|Extra arguments
`mutect2.mergeStats_timeout`|Int|5|Hours before task timeout
`mutect2.mergeStats_memory`|Int|4|Memory allocated for job
`mutect2.mergeStats_modules`|String|"gatk/4.1.6.0"|Names and versions of modules to load
`mutect2.mergeVCFs_timeout`|Int|12|Hours before task timeout
`mutect2.mergeVCFs_memory`|Int|4|Memory allocated for job
`mutect2.runMutect2_timeout`|Int|24|Hours before task timeout
`mutect2.runMutect2_memory`|Int|32|Memory allocated for job
`mutect2.runMutect2_threads`|Int|4|Number of threads to request
`mutect2.runMutect2_mutect2ExtraArgs`|String?|None|Extra arguments
`mutect2.runMutect2_mutectTag`|String|"mutect2"|Tag
`mutect2.splitStringToArray_modules`|String|""|Names and versions of modules to load
`mutect2.splitStringToArray_timeout`|Int|1|Hours before task timeout
`mutect2.splitStringToArray_memory`|Int|1|Memory allocated for job
`mutect2.splitStringToArray_lineSeparator`|String|","|line separator
`mutect2.normalBam`|File?|None|Input normal file (bam or sam)
`mutect2.normalBai`|File?|None|Index file for normal bam
`mutect2.pon`|File?|None|pon
`mutect2.ponIdx`|File?|None|pon ID
`mutect2.gnomad`|File?|None|gnomad
`mutect2.gnomadIdx`|File?|None|gnomad ID
`getFileName.jobMemory`|Int|4|memory allocated to preprocessing, in GB
`getFileName.timeout`|Int|1|timeout in hours
`getFileName.threads`|Int|1|number of cpu threads to be used
`combineVariants.jobMemory`|Int|24|memory allocated to preprocessing, in GB
`combineVariants.timeout`|Int|20|timeout in hours
`combineVariants.threads`|Int|8|number of cpu threads to be used
`annotation.modules`|String|"samtools/1.9 bcftools/1.9 htslib/1.9 tabix/1.9"|module for running preprocessing
`annotation.jobMemory`|Int|24|memory allocated to preprocessing, in GB
`annotation.timeout`|Int|20|timeout in hours
`annotation.threads`|Int|8|number of cpu threads to be used
`variantEffectPredictor.mergeVcfs_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`variantEffectPredictor.mergeVcfs_threads`|Int|4|Requested CPU threads.
`variantEffectPredictor.mergeVcfs_overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`variantEffectPredictor.mergeVcfs_jobMemory`|Int|24|Memory allocated to job (in GB).
`variantEffectPredictor.mergeVcfs_extraArgs`|String?|None|Additional arguments to be passed directly to the command.
`variantEffectPredictor.mergeVcfs_modules`|String|"gatk/4.1.7.0"|Required environment modules.
`variantEffectPredictor.mergeMafs_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`variantEffectPredictor.mergeMafs_threads`|Int|4|Requested CPU threads.
`variantEffectPredictor.mergeMafs_jobMemory`|Int|24|Memory allocated to job (in GB).
`variantEffectPredictor.mergeMafs_modules`|String|"tabix/0.2.6"|Required environment modules
`variantEffectPredictor.vcf2maf_timeout`|Int|48|Hours before task timeout
`variantEffectPredictor.vcf2maf_threads`|Int|4|Requested CPU threads
`variantEffectPredictor.vcf2maf_jobMemory`|Int|32|Memory allocated for this job (GB)
`variantEffectPredictor.vcf2maf_bufferSize`|Int|200|The buffer size
`variantEffectPredictor.vcf2maf_minHomVaf`|Float|0.7|The minimum vaf for homozygous calls
`variantEffectPredictor.vcf2maf_vepStats`|Boolean|true|If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'
`variantEffectPredictor.vcf2maf_species`|String|"homo_sapiens"|Species name
`variantEffectPredictor.vcf2maf_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.tumorOnlyAlign_timeout`|Int|6|Hours before task timeout
`variantEffectPredictor.tumorOnlyAlign_threads`|Int|4|Requested CPU threads
`variantEffectPredictor.tumorOnlyAlign_jobMemory`|Int|32|Memory allocated for this job (GB)
`variantEffectPredictor.tumorOnlyAlign_modules`|String|"bcftools/1.9 tabix/0.2.6"|Required environment modules
`variantEffectPredictor.tumorOnlyAlign_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.vep_timeout`|Int|16|Hours before task timeout
`variantEffectPredictor.vep_threads`|Int|4|Requested CPU threads
`variantEffectPredictor.vep_jobMemory`|Int|32|Memory allocated for this job (GB)
`variantEffectPredictor.vep_vepStats`|Boolean|true|If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'
`variantEffectPredictor.vep_species`|String|"homo_sapiens"|Species name
`variantEffectPredictor.vep_addParam`|String?|None|Additional vep parameters
`variantEffectPredictor.vep_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.subsetVcf_timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`variantEffectPredictor.subsetVcf_threads`|Int|4|Requested CPU threads.
`variantEffectPredictor.subsetVcf_jobMemory`|Int|32|Memory allocated to job (in GB).
`variantEffectPredictor.subsetVcf_modules`|String|"bcftools/1.9"|Required environment modules
`variantEffectPredictor.subsetVcf_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.chromosomeArray_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`variantEffectPredictor.chromosomeArray_threads`|Int|4|Requested CPU threads.
`variantEffectPredictor.chromosomeArray_jobMemory`|Int|1|Memory allocated to job (in GB).
`variantEffectPredictor.getSampleNames_timeout`|Int|1|Hours before task timeout
`variantEffectPredictor.getSampleNames_threads`|Int|4|Requested CPU threads
`variantEffectPredictor.getSampleNames_jobMemory`|Int|1|Memory allocated for this job (GB)
`variantEffectPredictor.targetBedTask_timeout`|Int|6|Hours before task timeout
`variantEffectPredictor.targetBedTask_threads`|Int|4|Requested CPU threads
`variantEffectPredictor.targetBedTask_jobMemory`|Int|32|Memory allocated for this job (GB)
`variantEffectPredictor.targetBedTask_modules`|String|"bedtools/2.27 tabix/0.2.6"|Required environment modules
`variantEffectPredictor.targetBedTask_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.normalName`|String?|None|Name of the normal sample
`filterMaf.mafNormalFile`|File?|None|input file for normal sample
`filterMaf.freqList`|String|"$MAF_FILTERING_ROOT/TGL.frequency.20210609.annot.txt"|frequency list used in maf annotation
`filterMaf.genesToKeep`|String|"$MAF_FILTERING_ROOT/genes_to_keep.txt"|gene list in maf filtering
`filterMaf.modules`|String|"python/3.9 pandas/1.4.2 maf-filtering/2023-10-06"|module for running preprocessing
`filterMaf.jobMemory`|Int|8|memory allocated to preprocessing, in GB
`filterMaf.timeout`|Int|1|timeout in hours
`filterMaf.threads`|Int|1|number of cpu threads to be used


### Outputs

Output | Type | Description
---|---|---
`tumorDcsScVcf`|File|DCS vcf for tumor sample
`tumorDcsScVcfIndex`|File|DCS vcf index for tumor sample
`tumorSscsScVcf`|File|SSCS vcf for tumor sample
`tumorSscsScVcfIndex`|File|SSCS vcf index for tumor sample
`tumorAllUniqueVcf`|File|vcf of DCS + singletons for tumor sample
`tumorAllUniqueVcfIndex`|File|vcf index for DCS + singletons for tumor sample
`tumorVepVcf`|File|vep vcf for tumor sample
`tumorVepVcfIndex`|File|vep vcf index for tumor sample
`tumorMafOutput`|File?|maf output for tumor sample
`filterredMaf`|File?|maf file after filtering

## Commands
This section lists command(s) run by mutect2ConsensusTumorOnly workflow

* Running mutect2ConsensusTumorOnly

```
    basename ~{fileName} | cut -d. -f1 
```
```
  python3<<CODE
  import subprocess
  import sys
  inputStrings = []
  v = "~{sep=' ' inputVcfs}"
  vcfFiles = v.split()
  w = "~{sep=' ' workflows}"
  workflowIds = w.split()
  priority = "~{priority}"
  
  if len(vcfFiles) != len(workflowIds):
      print("The arrays with input files and their respective workflow names are not of equal size!")
  else:
      for f in range(0, len(vcfFiles)):
          inputStrings.append("--variant:" + workflowIds[f] + " " + vcfFiles[f])

  javaMemory = ~{jobMemory} - 6 
  gatkCommand  = "$JAVA_ROOT/bin/java -Xmx" + str(javaMemory) + "G -jar $GATK_ROOT/GenomeAnalysisTK.jar "
  gatkCommand += "-T CombineVariants "
  gatkCommand += " ".join(inputStrings)
  gatkCommand += " -R ~{referenceFasta} "
  gatkCommand += "-o ~{outputPrefix}_combined.vcf.gz "
  gatkCommand += "-genotypeMergeOptions PRIORITIZE "
  gatkCommand += "-priority " + priority
  gatkCommand += " 2>&1"

  result_output = subprocess.run(gatkCommand, shell=True)
  sys.exit(result_output.returncode)
  
```
```
  bcftools annotate -a ~{uniqueVcf} \
 -c FMT/AD,FMT/DP ~{mergedVcf} -Oz \
 -o "~{outputPrefix}.merged.vcf.gz"

 tabix -p vcf "~{outputPrefix}.merged.vcf.gz"
```
```
    python3<<CODE
    ## Adapted from https://github.com/oicr-gsi/djerba/blob/GCGI-806_v1.0.0-dev/src/lib/djerba/plugins/tar/snv_indel/plugin.py
    ## this code will filter a maf file, generated from tumor-only mutect2 calls to identify likely germline calls generated from a mutect2 calls from the matched normal
    import pandas as pd
    maf_file_path = "~{mafFile}"
    maf_normal_path = "~{mafNormalFile}"
    freq_list_path = "~{freqList}"
    output_path_prefix = "~{outputPrefix}"
    genes_to_keep_path = "~{genesToKeep}"

    if maf_normal_path:
      df_bc = pd.read_csv(maf_normal_path,
                      sep = "\t",
                      on_bad_lines="error",
                      compression='gzip',
                      skiprows=[0])

    df_pl = pd.read_csv(maf_file_path,
                    sep = "\t",
                    on_bad_lines="error",
                    compression='gzip',
                    skiprows=[0])
    df_freq = pd.read_csv(freq_list_path,
                  sep = "\t")
    with open(genes_to_keep_path) as f:
      GENES_TO_KEEP = f.read()


    for row in df_pl.iterrows():
      hugo_symbol = row[1]['Hugo_Symbol']
      chromosome = row[1]['Chromosome']
      start_position = row[1]['Start_Position']
      reference_allele = row[1]['Reference_Allele']
      allele = row[1]['Allele']

      # If there is normal input, annotate rows with information from the matched normal and from the frequency table
      if maf_normal_path:
        # Lookup the entry in the BC and annotate the tumour maf with
        #   n_depth, n_ref_count, n_alt_count

        row_lookup = df_bc[(df_bc['Hugo_Symbol'] == hugo_symbol) & 
                    (df_bc['Chromosome'] == chromosome) & 
                    (df_bc['Start_Position'] == start_position) &
                    (df_bc['Reference_Allele'] == reference_allele) &
                    (df_bc['Allele'] == allele)]


        # If there's only one entry, take its normal values
        if len(row_lookup) == 1:
            df_pl.at[row[0], "n_depth"] = row_lookup['n_depth'].item()
            df_pl.at[row[0], "n_ref_count"] = row_lookup['n_ref_count'].item()
            df_pl.at[row[0], "n_alt_count"] = row_lookup['n_alt_count'].item()
      
        # If the entry isn't in the table, 
        # or if there is more than one value and so you can't choose which normal values to take, 
        # set them as 0
        else:
            df_pl.at[row[0], "n_depth"] = 0
            df_pl.at[row[0], "n_ref_count"] = 0
            df_pl.at[row[0], "n_alt_count"] = 0
            
      # Lookup the entry in the frequency table and annotate the tumour maf with Freq
    
      row_lookup = df_freq[(df_freq['Start_Position'] == row[1]['Start_Position']) &
                          (df_freq['Reference_Allele'] == row[1]['Reference_Allele']) &
                          ((df_freq['Tumor_Seq_Allele'] == row[1]['Tumor_Seq_Allele1']) |
                          (df_freq['Tumor_Seq_Allele'] == row[1]['Tumor_Seq_Allele2']))]

      if len(row_lookup) > 0:
          df_pl.at[row[0], 'Freq'] = row_lookup['Freq'].item()
      else:
          df_pl.at[row[0], 'Freq'] = 0

    # Filter the maf to remove rows based on various criteria, but always maintaining genes in the GENES_TO_KEEP list  
    for row in df_pl.iterrows():
        hugo_symbol = row[1]['Hugo_Symbol']
        frequency = row[1]['Freq']
        gnomAD_AF = row[1]['gnomAD_AF']
        n_alt_count = row[1]['n_alt_count']
        if hugo_symbol not in GENES_TO_KEEP or frequency > 0.1 or n_alt_count > 4 or gnomAD_AF > 0.001:
            df_pl = df_pl.drop(row[0])   

    df_pl.to_csv(output_path_prefix + '_filtered_maf_for_tar.maf.gz', sep = "\t", compression='gzip', index=False)
    CODE
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
