version 1.0
workflow run_full {
  input {
    File fullin_TEs
    Array[String] file_id
    Array[File] allcs
    Array[File] allcs_idx
    File chrom_size
    Array[String] fullin_sample_name_column
    Int memory_GB
    Int nCPUs
    Float threshold_QC
  }
  String fullin_sample_name = fullin_sample_name_column[0]
    parameter_meta {
      fullin_TEs: "Path to bed file containing TEs"
      file_id: "terra table column containing file IDs"
      allcs: "terra table column containing location of 'allc_*.tsv.gz'"
      allcs_idx: "terra table column containing location of 'allc.tsv.idx'"
      chrom_size: "Path to chrom.sizes obtained by UCSC fetchChromSizes.sh"
      fullin_sample_name_column: "name of sample"
      memory_GB: "memory, in gigabytes"
      nCPUs: "CPUs for parallel execution"
      threshold_QC: "Threshold for discarding methylation fraction"
  }
  call mangle_bed {
    input:
      bed = fullin_TEs,
      mem = memory_GB
  }
  call generate_dataset {
    input:
        fileIDs = file_id,
        allc_list = allcs,
        SampleName = fullin_sample_name,
        nCPU = nCPUs,
        mangledTEs = mangle_bed.bed_mangled,
        chromSize = chrom_size,
        mem = memory_GB
  }
  call calculate_fractions {
    input:
        tempzarr = generate_dataset.zarrTar,
        SampleName = fullin_sample_name,
        mem = memory_GB,
        threshold = threshold_QC
  }
}

task mangle_bed {
    input {
        File bed
        Int mem
    }
    command <<<
      run_steamer mangle-bed-file-ids ~{bed}
    >>>
    output {
        File bed_mangled = sub(bed, basename(bed, ".bed"), basename(bed, ".bed") + "_mangled")
    }
    runtime {
    docker: "ghcr.io/welch-lab/steamer:latest"
    memory: mem + "GB"
  }
}

task generate_dataset {
    input {
        Array[String] fileIDs
        Array[File]   allc_list
        String SampleName
        Int nCPU
        File mangledTEs
        File chromSize
        Int mem
    }
    Int nCPUscale = floor(nCPU*0.75)
    Array[Array[String]] tsvPair = [fileIDs, allc_list]
    File allc_table = write_tsv(tsvPair)
    command <<<
    allcools generate-dataset --allc_table ~{allc_table} --output_path=~{SampleName}.mcds --obs_dim cell \
    --cpu ~{nCPUscale} --chunk 50 --regions TEs ~{mangledTEs} --chrom_size_path ~{chromSize} \
    --quantifiers TEs count CHN; tar -cf tempzarr.tar ~{SampleName}.mcds
    >>>
    output {
#this should be a Directory but cromwell doesn't support WDL 1.2
        File zarrTar = "tempzarr.tar"
    }
    runtime {
      docker: "ghcr.io/welch-lab/steamer:latest"
      memory: mem + "GB"
      cpu: nCPUs
      disks: "local-disk 375 LOCAL"
      
    }
}

task calculate_fractions {
    input {
        File tempzarr
        String SampleName
        Int mem
        Float threshold
    }
    command <<<
    tar -xf ~{tempzarr}; \
    run_steamer mc-fractions ~{SampleName}.mcds ~{threshold}
    >>>
    output {
        File count_mat = SampleName + ".mcds.mtx"
    }
    runtime {
        docker: "ghcr.io/welch-lab/steamer:latest"
        memory: mem + "GB"
  }
}
