version 1.0
workflow run_full {
  input {
    File fullin_TEs
    Array[String] file_id
    Array[File] allcs
    File chrom_size
    Array[String] fullin_sample_name_column
    Int memory_GB
    Int nCPUs
    Int threshold_QC
  }
  String fullin_sample_name = fullin_sample_name_column[0]
    parameter_meta {
      fullin_TEs: "Path to bed file containing TEs"
      file_id: "terra table column containing file IDs"
      allcs: "terra table column containing location of 'allc_*.tsv.gz'"
      chrom_size: "Path to chrom.sizes obtained by UCSC fetchChromSizes.sh"
      fullin_sample_name_column: "name of sample"
      memory_GB: "memory, in gigabytes"
      nCPUs: "CPUs for parallel execution"
      threshold_QC: "Threshold for discarding methylation value"
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
  output {
    File mtx = calculate_fractions.count_mat
  }
}

task mangle_bed {
    input {
        File bed
        Int mem
    }
    command <<<
      run_steamer mangle-bed-file-ids ~{bed} TEs_mangled.bed
    >>>
    output {
        File bed_mangled = "TEs_mangled.bed"
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
    parameter_meta {
        allc_list: {
            description: "erra table column containing location of 'allc_*.tsv.gz",
            localization_optional: true
        }
    }
    Int nCPUscale = ceil(nCPU*0.75)
    Int disk = ceil(size(allc_list, "G")) + 375
    String disk_string = "local-disk " + disk + " LOCAL"
    Array[Array[String]] initial_paired = [fileIDs, allc_list]
    Array[Array[String]] tsvPaired = transpose(initial_paired)
    File allc_table = write_tsv(tsvPaired)
    command <<<
    CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt \
    GCS_REQUESTER_PAYS_PROJECT=$(curl -s "http://metadata.google.internal/computeMetadata/v1/project/project-id" -H "Metadata-Flavor: Google") \
    GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token) \
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
      cpu: nCPU
      disks: disk_string
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
