version 1.0

workflow run_full {
  input {
    File full_index
    File full_t2gs
    Array[File] full_fastqs
    String fullin_sample_name
    Int memory_GB
  }
    parameter_meta {
      full_t2gs: "Path to t2gs file"
      full_index:"Path to index file"
      full_fastqs: "List of FASTQs to process"
      fullin_sample_name: "name of sample"
      memory_GB: "memory, in gigabytes"
  }
    call run_rna {
      input:
        index = full_index,
        t2gs = full_t2gs,
        SampleName = fullin_sample_name,
        mem = memory_GB,
        fastqs = full_fastqs
    }
    output {
    File Mat = run_rna.mtx
    File DF =  run_rna.genes
    File Bar = run_rna.barcodes
  }
}
task run_rna {
    input {
      File index
      File t2gs
      String SampleName
      Int mem
      Array[File] fastqs
    }
    command <<<
      kb count -i ~{index} -g ~{t2gs} --h5ad -x 10XV3 -o ~{SampleName} --verbose ~{sep=' ' fastqs}
    >>>
    output {
      File TenxWhiteList = SampleName + "/10x_version3_whitelist.txt"
      File inspect = SampleName + "/inspect.json"
      File kb_info = SampleName + "/kb_info.json"
      File matrixEC = SampleName + "/matrix.ec"
      File outputBus = SampleName + "/output.bus"
      File outputBusUnFiltered = SampleName + "/output.unfiltered.bus"
      File run_info = SampleName + "/run_info.json"
      File transcripts = SampleName + "/transcripts.txt"
      File adata = SampleName + "/counts_unfiltered/adata.h5ad"
      File barcodes = SampleName + "/counts_unfiltered/cells_x_genes.barcodes.txt"
      File genesNames = SampleName + "/counts_unfiltered/cells_x_genes.genes.names.txt"
      File genes = SampleName + "/counts_unfiltered/cells_x_genes.genes.txt"
      File mtx = SampleName + "/counts_unfiltered/cells_x_genes.mtx"
    }
  runtime {
    docker: "ghcr.io/welch-lab/steamer:latest"
    memory: mem + "GB"
    disks: "local-disk 50 SSD"
  }
}

