version 1.0
workflow steamerprep {
    input {
      File? in_TEs
      File? in_Frags
      File? in_QCbar
      Int in_mem
    }
    parameter_meta {
      in_TEs: "Path to .gz or .csv file containing TEs"
      in_Frags: "Path to .gz or .csv file containing Fragments"
      in_QCbar: "Path to .gz or .csv file containing QC'd barcodes"
      in_mem: "Memory in GB"
    }
    call make_bedfiles {
      input:
        TEs = in_TEs,
        Frags = in_Frags,
        QCbar = in_QCbar,
        mem = in_mem
    }
    call intersect_bedfiles {
      input:
        Frags = make_bedfiles.outFrags,
        TEs = make_bedfiles.outTEs,
        mem = in_mem
    }
    output {
      File prepped = intersect_bedfiles.intersected
    }
}

task make_bedfiles {
    input {
        File? TEs
        File? Frags
        File? QCbar
        Int mem
    }
    command <<<
      run_steamer create-bed-for-fragments ~{Frags} ~{QCbar} && run_steamer create-bed-for-tes ~{TEs}
    >>>
    output {
      File outFrags = "Frag.bed"
      File outTEs = "TEs.bed"
    }
  runtime {
    docker: "us-central1-docker.pkg.dev/cobalt-entropy-358220/welch-lab/steamer:latest"
    memory: mem + "GB"
  }
}

task intersect_bedfiles {
  input {
    File? Frags
    File? TEs
    Int mem
  }
  command <<<
    bedtools intersect -a ~{TEs} -b ~{Frags} -wb -sorted > bed_intersect.bed
  >>>
  output {
    File intersected = "bed_intersect.bed"
  }
    runtime {
    docker: "us-central1-docker.pkg.dev/cobalt-entropy-358220/welch-lab/steamer:latest"
    memory: mem + "GB"
  }
}
