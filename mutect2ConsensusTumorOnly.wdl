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
  String inputRefDict
  String inputRefFasta
  String inputRefFai
  String inputMutectModules
  String combineVariants_modules
  String variantEffectPredictor_vcf2maf_modules
  String variantEffectPredictor_vcf2maf_ncbiBuild
  String variantEffectPredictor_vcf2maf_vepCacheDir
  String variantEffectPredictor_vcf2maf_vepPath
  String variantEffectPredictor_vep_modules
  String variantEffectPredictor_vep_ncbiBuild
  String variantEffectPredictor_vep_vepCacheDir
}

workflow mutect2ConsensusTumorOnly {
  input {
    InputGroup tumorInputGroup
    String outputFileNamePrefix
    String intervalFile
    String inputIntervalsToParalellizeBy
    String tumorName
    String reference
  }

  Map[String,GenomeResources] resources = {
    "hg19": {
      "inputRefDict": "$HG19_ROOT/hg19_random.dict",
      "inputRefFai": "$HG19_ROOT/hg19_random.fa.fai",
      "inputRefFasta": "$HG19_ROOT/hg19_random.fa",
      "inputMutectModules": "gatk/4.1.6.0 hg19/p13 samtools/1.9",
      "combineVariants_modules": "gatk/3.6-0 tabix/0.2.6 hg19/p13",
      "variantEffectPredictor_vep_modules": "vep/105.0 tabix/0.2.6 vep-hg19-cache/105 hg19/p13",
      "variantEffectPredictor_vep_vepCacheDir": "$VEP_HG19_CACHE_ROOT/.vep",
      "variantEffectPredictor_vep_ncbiBuild": "GRCh37",
      "variantEffectPredictor_vcf2maf_modules": "vcf2maf/1.6.21b tabix/0.2.6 hg19/p13 vep-hg19-cache/105",
      "variantEffectPredictor_vcf2maf_vepCacheDir": "$VEP_HG19_CACHE_ROOT/.vep",
      "variantEffectPredictor_vcf2maf_vepPath": "$VEP_ROOT/bin/",
      "variantEffectPredictor_vcf2maf_ncbiBuild": "GRCh37"
      },
    "hg38": {
      "inputRefDict": "$HG38_ROOT/hg38_random.dict",
      "inputRefFai": "$HG38_ROOT/hg38_random.fa.fai",
      "inputRefFasta": "$HG38_ROOT/hg38_random.fa",
      "inputMutectModules": "gatk/4.1.6.0 hg38/p12 samtools/1.9",
      "combineVariants_modules": "gatk/3.6-0 tabix/0.2.6 hg38/p12",
      "variantEffectPredictor_vep_modules": "vep/105.0 tabix/0.2.6 vep-hg38-cache/105 hg38/p12",
      "variantEffectPredictor_vep_vepCacheDir": "$VEP_HG38_CACHE_ROOT/.vep",
      "variantEffectPredictor_vep_ncbiBuild": "GRCh38",
      "variantEffectPredictor_vcf2maf_modules": "vcf2maf/1.6.21b tabix/0.2.6 hg38/p12 vep-hg38-cache/105",
      "variantEffectPredictor_vcf2maf_vepCacheDir": "$VEP_HG38_CACHE_ROOT/.vep",
      "variantEffectPredictor_vcf2maf_vepPath": "$VEP_ROOT/bin/",
      "variantEffectPredictor_vcf2maf_ncbiBuild": "GRCh38"
      }
  }
  

  parameter_meta {
    tumorInputGroup: "partitioned bam files from umiConsensus outputs for tumor sample"
    normalInputGroup: "partitioned bam files from umiConsensus outputs for normal sample"
    outputFileNamePrefix: "Prefix to use for output file"
    intervalFile: "interval file to subset variant calls"
    inputIntervalsToParalellizeBy: "intervals for parallelization"
    tumorName: "Name of the tumor sample"
    normalName: "name of the normal sample"
    reference: "reference version"
  }

  Array[BamAndBamIndex]partitionedBams = [tumorInputGroup.dcsScBamAndIndex, tumorInputGroup.sscsScBamAndIndex, tumorInputGroup.allUniqueBamAndIndex]
  scatter ( bamAndIndex in partitionedBams ) {
    call mutect2.mutect2 {
      input:
        tumorBam = bamAndIndex.bam,
        tumorBai = bamAndIndex.bamIndex,
        filter_refDict = resources[reference].inputRefDict,
        filter_refFai = resources[reference].inputRefFai,
        filter_refFasta = resources[reference].inputRefFasta,
        filter_modules = resources[reference].inputMutectModules,
        mergeVCFs_refFasta = resources[reference].inputRefFasta,
        mergeVCFs_modules = resources[reference].inputMutectModules,
        runMutect2_refDict = resources[reference].inputRefDict,
        runMutect2_refFai = resources[reference].inputRefFasta,
        runMutect2_refFasta = resources[reference].inputRefFasta,
        runMutect2_modules = resources[reference].inputMutectModules,
        intervalFile = intervalFile,
        intervalsToParallelizeBy = inputIntervalsToParalellizeBy
      }
    }

  Array[File] mutect2FilteredVcfFiles = mutect2.filteredVcfFile
  Array[File] mutect2FilteredVcfIndexes = mutect2.filteredVcfIndex

  call getFileName{
    input:
      fileName = mutect2FilteredVcfFiles[0]
  }

  call combineVariants {
    input: 
      inputVcfs = [mutect2FilteredVcfFiles[0],mutect2FilteredVcfFiles[1]],
      inputIndexes = [mutect2FilteredVcfIndexes[0],mutect2FilteredVcfIndexes[1]],
      priority = "mutect2-dcsSc,mutect2-sscsSc",
      outputPrefix = getFileName.outputFileName,
      referenceFasta = resources[reference].inputRefFasta,
      modules = resources[reference].combineVariants_modules
  }

  call annotation {
    input: 
      uniqueVcf = mutect2FilteredVcfFiles[2],
      uniqueVcfIndex = mutect2FilteredVcfIndexes[2],
      mergedVcf = combineVariants.combinedVcf,
      mergedVcfIndex = combineVariants.combinedIndex,
      outputPrefix = getFileName.outputFileName
  }

  call vep.variantEffectPredictor {
    input: 
      vcfFile = annotation.annotatedCombinedVcf,
      vcfIndex = annotation.annotatedCombinedIndex,
      toMAF = true,
      onlyTumor = true,
      tumorOnlyAlign_updateTagValue = true,
      vcf2maf_retainInfoProvided = true,
      vep_referenceFasta = resources[reference].inputRefFasta,
      vcf2maf_referenceFasta = resources[reference].inputRefFasta,
      targetBed = intervalFile,
      tumorName = tumorName,
      vcf2maf_modules = resources[reference].variantEffectPredictor_vcf2maf_modules,
      vcf2maf_ncbiBuild = resources[reference].variantEffectPredictor_vcf2maf_ncbiBuild,
      vcf2maf_vepCacheDir = resources[reference].variantEffectPredictor_vcf2maf_vepCacheDir,
      vcf2maf_vepPath = resources[reference].variantEffectPredictor_vcf2maf_vepPath,
      vep_modules = resources[reference].variantEffectPredictor_vep_modules,
      vep_ncbiBuild = resources[reference].variantEffectPredictor_vep_ncbiBuild,
      vep_vepCacheDir = resources[reference].variantEffectPredictor_vep_vepCacheDir
  }

  File? tumorMaf = variantEffectPredictor.outputMaf

  call filterMaf {
    input:
    mafFile = tumorMaf,
    outputPrefix = outputFileNamePrefix
  }

  meta {
    author: "Alexander Fortuna, Rishi Shah and Gavin Peng"
    email: "alexander.fortuna@oicr.on.ca, rshah@oicr.on.ca, and gpeng@oicr.on.ca"
    description: "The Mutect2Consensus workflow will process umiConsensus outputs for the tumour data through mutect2 in tumour only mode to call variants then use information from the matched normal to identify likely germline variants."
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
      tumorDcsScVcf: "DCS vcf for tumor sample",
      tumorDcsScVcfIndex: "DCS vcf index for tumor sample",
      tumorSscsScVcf: "SSCS vcf for tumor sample",
      tumorSscsScVcfIndex: "SSCS vcf index for tumor sample",
      tumorAllUniqueVcf: "vcf of DCS + singletons for tumor sample",
      tumorAllUniqueVcfIndex: "vcf index for DCS + singletons for tumor sample",
      tumorVepVcf: "vep vcf for tumor sample",
      tumorVepVcfIndex: "vep vcf index for tumor sample",
      tumorMafOutput: "maf output for tumor sample",
      filterredMaf: "maf file after filtering",
    }
  }

  output {
    File tumorDcsScVcf = mutect2FilteredVcfFiles[0]
    File tumorDcsScVcfIndex = mutect2FilteredVcfIndexes[0]
    File tumorSscsScVcf = mutect2FilteredVcfFiles[1]
    File tumorSscsScVcfIndex = mutect2FilteredVcfIndexes[1]
    File tumorAllUniqueVcf = mutect2FilteredVcfFiles[2]
    File tumorAllUniqueVcfIndex = mutect2FilteredVcfIndexes[2]
    File tumorVepVcf = variantEffectPredictor.outputVcf
    File tumorVepVcfIndex = variantEffectPredictor.outputTbi
    File? tumorMafOutput = tumorMaf
    File? filterredMaf = filterMaf.filterredMaf
  }
}

task getFileName {
  input {
    File fileName
    Int jobMemory = 4
    Int timeout = 1
    Int threads = 1
  }

  parameter_meta {
    fileName: "the file to get basename with"
    jobMemory: "memory allocated to preprocessing, in GB"
    timeout: "timeout in hours"
    threads: "number of cpu threads to be used"
    }

  command <<<
    basename ~{fileName} | cut -d. -f1 
  >>>

  output {
    String outputFileName = read_string(stdout())
  }

  runtime {
    memory:  "~{jobMemory} GB"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
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
    String genesToKeep = "$MAF_FILTERING_ROOT/genes_to_keep.txt"
    String outputPrefix 
    String modules = "python/3.9 pandas/1.4.2 maf-filtering/2023-10-06"
    Int jobMemory = 8
    Int timeout = 1
    Int threads = 1
  }

  parameter_meta {
    mafFile: "input maf file for tumor sample"
    mafNormalFile: "input file for normal sample"
    freqList: "frequency list used in maf annotation"
    genesToKeep: "gene list in maf filtering"
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
  >>>

  runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  cpu:     "~{threads}"
  timeout: "~{timeout}"
  }

  output {
    File filterredMaf = "~{outputPrefix}_filtered_maf_for_tar.maf.gz"
  }
}