version 1.0
workflow run_full {
  input {
    File fullin_TEs
    Array[String] file_id
    Array[File] allcs
    Array[File] allcs_idx
    File chrom_size
    String fullin_sample_name
    Int memory_GB
  }
    parameter_meta {
      fullin_TEs: "Path to bed file containing TEs"
      file_id: "terra table column containing file IDs"
      allcs: "terra table column containing location of 'allc_*.tsv.gz'"
      allcs_idx: "terra table column containing location of 'allc.tsv.idx'"
      chrom_size: "Path to chrom.sizes obtained by UCSC fetchChromSizes.sh"
      fullin_sample_name: "name of sample"
      memory_GB: "memory, in gigabytes"
  }
  call mangle_bed {
    input:
      bed = fullin_TEs,
      mem = memory_GB
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
        File bed_mangled = basename(bed, ".bed") + "_mangled.bed"
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
    Array[Array[String]] tsvPair = [fileIDs, allc_list]
    File allc_table = write_tsv(tsvPair)
    command <<<
    allcools generate-dataset --allc_table ~{allc_table} --output_path=~{SampleName}.mcds --obs_dim cell \
    --cpu ~{nCPU} --chunk 50 --regions TEs ~{mangledTEs} --chrom_size_path ~{chromSize} \
    --quantifiers TEs count CHN; tar -cf tempzarr.tar ~{SampleName}.mcds
    >>>
    output {
#this should be a Directory but cromwell doesn't support WDL 1.2
        File zarrTar = "/tempzarr.tar"
    }
    runtime {
      docker: "ghcr.io/welch-lab/steamer:latest"
      memory: mem + "GB"
      #cores?
    }
}

task calculate_fractions {
    input {
        File tempzarr
        String SampleName
        Int mem
        Int threshold
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
