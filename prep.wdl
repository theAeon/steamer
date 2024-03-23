version 1.0
workflow steamerprep {
    input {
      File in_TEs
      File in_Frags
      File? in_QCbar
    }
    parameter_meta {
      in_TEs: "Path to .gz or .csv file containing TEs"
      in_Frags: "Path to .gz or .csv file containing Fragments"
      in_QCbar: "Path to .gz or .csv file containing QC'd barcodes"
    }
    call make_bedfiles {
      input:
        TEs = in_TEs,
        Frags = in_Frags,
        QCbar = in_QCbar
    }
    call intersect_bedfiles {
      input:
        Frags = make_bedfiles.outFrags,
        TEs = make_bedfiles.outTEs
    }
    output {
      File prepped = intersect_bedfiles.intersected
    }
}

task make_bedfiles {
    input {
        File TEs
        File Frags
        File? QCbar
    }
    command <<<
      run_steamer create-bed-for-fragments ~{Frags} ~{QCbar} && run_steamer create-bed-for-tes ~{TEs}
    >>>
    output {
      File outFrags = "Frag.bed"
      File outTEs = "TEs.bed"
    }
  runtime {
    docker: "ghcr.io/welch-lab/steamer:latest"
  }
}

task intersect_bedfiles {
  input {
    File Frags
    File TEs
  }
  command <<<
    bedtools intersect -a ~{TEs} -b ~{Frags} -wb -sorted > bed_intersect.bed
  >>>
  output {
    File intersected = "bed_intersect.bed"
  }
    runtime {
    docker: "ghcr.io/welch-lab/steamer:latest"
  }
}
