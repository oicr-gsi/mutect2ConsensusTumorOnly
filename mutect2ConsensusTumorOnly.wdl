version 1.0

import "imports/pull_mutect2.wdl" as mutect2
import "imports/pull_variantEffectPredictor.wdl" as vep


struct BamAndBamIndex {
  File bam
  File bamIndex
}

struct InputGroup {
  BamAndBamIndex dcsScBamAndIndex
  BamAndBamIndex sscsScBamAndIndex
  BamAndBamIndex allUniqueBamAndIndex
}

struct GenomeResources {
  String inputRefFasta
  String combineVariants_modules
}

workflow mutect2ConsensusTumorOnly {
  input {
    InputGroup tumorInputGroup
    String reference
    String intervalFile
    String inputIntervalsToParalellizeBy
    String tumorName
    String outputFileNamePrefix
    String gatk
    Boolean filterMafFile
  }

  Map[String,GenomeResources] run_resources = {
    "hg19": {
      "inputRefFasta": "$HG19_ROOT/hg19_random.fa",
      "combineVariants_modules": "gatk/3.6-0 tabix/0.2.6 hg19/p13",
      },
    "hg38": {
      "inputRefFasta": "$HG38_ROOT/hg38_random.fa",
      "combineVariants_modules": "gatk/3.6-0 tabix/0.2.6 hg38/p12",
      }
  }
  

  parameter_meta {
    tumorInputGroup: "partitioned bam files from umiConsensus outputs for tumor sample"
    intervalFile: "interval file to subset variant calls"
    inputIntervalsToParalellizeBy: "intervals for parallelization"
    tumorName: "Name of the tumor sample"
    reference: "reference version"
    outputFileNamePrefix: "Prefix to use for output file"
    gatk: "gatk version to be used"
    filterMafFile: "whether filter the maf file"
  }

  Array[BamAndBamIndex]partitionedBams = [tumorInputGroup.dcsScBamAndIndex, tumorInputGroup.sscsScBamAndIndex, tumorInputGroup.allUniqueBamAndIndex]
  scatter ( bamAndIndex in partitionedBams ) {
    call mutect2.mutect2 {
      input:
        tumorBam = bamAndIndex.bam,
        tumorBai = bamAndIndex.bamIndex,
        intervalFile = intervalFile,
        intervalsToParallelizeBy = inputIntervalsToParalellizeBy,
        reference = reference,
        gatk = gatk,
        outputFileNamePrefix = outputFileNamePrefix
      }
    }

  Array[File] mutect2FilteredVcfFiles = mutect2.filteredVcfFile
  Array[File] mutect2FilteredVcfIndexes = mutect2.filteredVcfIndex

  call combineVariants {
    input: 
      inputVcfs = [mutect2FilteredVcfFiles[0],mutect2FilteredVcfFiles[1]],
      inputIndexes = [mutect2FilteredVcfIndexes[0],mutect2FilteredVcfIndexes[1]],
      priority = "mutect2-dcsSc,mutect2-sscsSc",
      outputPrefix = outputFileNamePrefix,
      referenceFasta = run_resources[reference].inputRefFasta,
      modules = run_resources[reference].combineVariants_modules
  }

  call annotation {
    input: 
      uniqueVcf = mutect2FilteredVcfFiles[2],
      uniqueVcfIndex = mutect2FilteredVcfIndexes[2],
      mergedVcf = combineVariants.combinedVcf,
      mergedVcfIndex = combineVariants.combinedIndex,
      outputPrefix = outputFileNamePrefix
  }

  call vep.variantEffectPredictor {
    input: 
      vcfFile = annotation.annotatedCombinedVcf,
      vcfIndex = annotation.annotatedCombinedIndex,
      toMAF = true,
      onlyTumor = true,
      tumorOnlyAlign_updateTagValue = true,
      vcf2maf_retainInfoProvided = true,
      tumorName = outputFileNamePrefix,
      reference = reference
  }

  File? tumor_Maf = variantEffectPredictor.outputMaf

  if (filterMafFile && defined(tumor_Maf)) {
    call filterMaf {
      input:
      mafFile = tumor_Maf,
      outputPrefix = outputFileNamePrefix
    }
  }

  meta {
    author: "Gavin Peng"
    email: "gpeng@oicr.on.ca"
    description: "The Mutect2Consensus workflow will process umiConsensus outputs for the tumour data through mutect2 in tumour only mode to call variants and annotation."
    dependencies: [
      {
        name: "gatk/3.6-0",
        url: "https://gatk.broadinstitute.org"
      },
      {
        name: "python/3.9",
        url: "https://www.python.org/downloads/"
      },
      {
        name: "vep/105.0",
        url: "https://useast.ensembl.org/info/docs/tools/vep/"
      },
      {
        name: "gatk/4.1.6.0",
        url: "https://gatk.broadinstitute.org/"
      },
      {
        name: "tabix/0.2.6",
        url: "https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2/download"
      },
      {
        name: "vcf2maf/1.6",
        url: "https://github.com/mskcc/vcf2maf"
      },
      {
        name: "pandas/1.4.2",
        url: "https://pandas.pydata.org/"
      }
    ]

    output_meta: {
      tumorVcf: {
          description: "vep vcf for tumor sample",
          vidarr_label: "tumorVcf"
      },
      tumorVcfIndex: {
          description: "vep vcf index for tumor sample",
          vidarr_label: "tumorVcfIndex"
      },
      tumorMaf: {
          description: "maf output for tumor sample",
          vidarr_label: "tumorMaf"
      },
      filteredMaf: {
          description: "maf file after filtering",
          vidarr_label: "filterredMaf"
      }
    }
  }

  output {
    File tumorVcf = variantEffectPredictor.outputVcf
    File tumorVcfIndex = variantEffectPredictor.outputTbi
    File? tumorMaf = tumor_Maf
    File? filteredMaf = filterMaf.filteredMaf
  }
}

task combineVariants {
input {
 Array[File] inputVcfs
 Array[File] inputIndexes
 Array[String] workflows
 String referenceFasta
 String outputPrefix 
 String modules
 String priority
 Int jobMemory = 24
 Int timeout = 20
 Int threads = 8
}

parameter_meta {
 inputVcfs: "array of input vcf files"
 inputIndexes: "array of tabix indexes for vcf files"
 workflows: "array of ids of producer workflows"
 referenceFasta: "path to the reference FASTA file"
 outputPrefix: "prefix for output file"
 modules: "modules for running preprocessing"
 priority: "Comma-separated list defining priority of workflows when combining variants"
 jobMemory: "memory allocated to preprocessing, in GB"
 timeout: "timeout in hours"
 threads: "number of cpu threads to be used"
}

command <<<
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
  CODE
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  cpu:     "~{threads}"
  timeout: "~{timeout}"
}

output {
  File combinedVcf = "~{outputPrefix}_combined.vcf.gz"
  File combinedIndex = "~{outputPrefix}_combined.vcf.gz.tbi"
}
}

task annotation {
input {
 File uniqueVcf 
 File uniqueVcfIndex
 File mergedVcf
 File mergedVcfIndex
 String outputPrefix
 String modules = "samtools/1.9 bcftools/1.9 htslib/1.9 tabix/1.9"
 Int jobMemory = 24
 Int timeout = 20
 Int threads = 8
}

parameter_meta {
 uniqueVcf: "input unique vcf files"
 uniqueVcfIndex: "input unique tabix indexes for vcf files"
 mergedVcf: "input merged vcf"
 mergedVcfIndex: "input merged vcf index"
 outputPrefix: "prefix for output file"
 modules: "module for running preprocessing"
 jobMemory: "memory allocated to preprocessing, in GB"
 timeout: "timeout in hours"
 threads: "number of cpu threads to be used"
}

command <<<
  bcftools annotate -a ~{uniqueVcf} \
 -c FMT/AD,FMT/DP ~{mergedVcf} -Oz \
 -o "~{outputPrefix}.merged.vcf.gz"

 tabix -p vcf "~{outputPrefix}.merged.vcf.gz"
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  cpu:     "~{threads}"
  timeout: "~{timeout}"
}

output {
  File annotatedCombinedVcf = "~{outputPrefix}.merged.vcf.gz"
  File annotatedCombinedIndex = "~{outputPrefix}.merged.vcf.gz.tbi"
}
}

task filterMaf {
  input {
    File? mafFile
    File? mafNormalFile
    String freqList ="$MAF_FILTERING_ROOT/TGL.frequency.20210609.annot.txt"
    String outputPrefix 
    String modules = "python/3.9 pandas/1.4.2 maf-filtering/2024-07-10"
    Int jobMemory = 8
    Int timeout = 1
    Int threads = 1
  }

  parameter_meta {
    mafFile: "input maf file for tumor sample"
    mafNormalFile: "input file for normal sample"
    freqList: "frequency list used in maf annotation"
    outputPrefix: "prefix for output file"
    modules: "module for running preprocessing"
    jobMemory: "memory allocated to preprocessing, in GB"
    timeout: "timeout in hours"
    threads: "number of cpu threads to be used"
  }


  command <<<
    python3<<CODE
    ## Adapted from https://github.com/oicr-gsi/djerba/blob/GCGI-806_v1.0.0-dev/src/lib/djerba/plugins/tar/snv_indel/plugin.py
    ## this code will filter a maf file, generated from tumor-only mutect2 calls to identify likely germline calls generated from a mutect2 calls from the matched normal
    import pandas as pd
    maf_file_path = "~{mafFile}"
    maf_normal_path = "~{mafNormalFile}"
    freq_list_path = "~{freqList}"
    output_path_prefix = "~{outputPrefix}"
    clean_columns = ["t_depth", "t_ref_count", "t_alt_count", "n_depth", "n_ref_count", "n_alt_count", "gnomAD_AF"]

    if maf_normal_path:
      df_bc = pd.read_csv(maf_normal_path,
                      sep = "\t",
                      on_bad_lines="error",
                      compression='gzip',
                      skiprows=[0])

      # Clean up df_bc if normal maf available
      for column in clean_columns:
        # Convert to numeric, setting errors='coerce' to turn non-numeric values into NaN
        df_bc[column] = pd.to_numeric(df_bc[column], errors='coerce')
        # Replace NaN with 0
        df_bc[column] = df_bc[column].fillna(0)

    df_pl = pd.read_csv(maf_file_path,
                    sep = "\t",
                    on_bad_lines="error",
                    compression='gzip',
                    skiprows=[0])
    
    # Clean up df_pl
    for column in clean_columns:
      # Convert to numeric, setting errors='coerce' to turn non-numeric values into NaN
      df_pl[column] = pd.to_numeric(df_pl[column], errors='coerce')
      # Replace NaN with 0
      df_pl[column] = df_pl[column].fillna(0)

    df_freq = pd.read_csv(freq_list_path,
                  sep = "\t")


    for row in df_pl.iterrows():
      chromosome = row[1]['Chromosome']
      start_position = row[1]['Start_Position']
      reference_allele = row[1]['Reference_Allele']
      allele = row[1]['Allele']

      # If there is normal input, annotate rows with information from the matched normal and from the frequency table
      if maf_normal_path:
        # Lookup the entry in the BC and annotate the tumour maf with
        #   n_depth, n_ref_count, n_alt_count

        row_lookup = df_bc[
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

    # Filter the maf to remove rows based on various criteria
    for row in df_pl.iterrows():
        frequency = row[1]['Freq']
        gnomAD_AF = row[1]['gnomAD_AF']
        n_alt_count = row[1]['n_alt_count']
        if  frequency > 0.1 or n_alt_count > 4 or gnomAD_AF > 0.001:
            df_pl = df_pl.drop(row[0])   

    df_pl.to_csv(output_path_prefix + '_filtered_maf.gz', sep = "\t", compression='gzip', index=False)
    CODE
  >>>

  runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  cpu:     "~{threads}"
  timeout: "~{timeout}"
  }

  output {
    File filteredMaf = "~{outputPrefix}_filtered_maf.gz"
  }
}