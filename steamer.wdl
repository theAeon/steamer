version 1.0
workflow steamer {
    input {
      File in_TEs
      File in_Frags
    }
    parameter_meta {
      in_TEs: "Path to .gz or .csv file containing TEs"
      in_Frags: "Path to .gz or .csv file containing Fragments"
    }
    call make_bedfiles {
      input:
        TEs = in_TEs,
        Frags = in_Frags
    }
}

task make_bedfiles {
    input {
        File TEs
        File Frags
    }
    command <<<
      run_steamer create-bed-for-fragments ~{Frags} && run_steamer create-bed-for-tes ~{TEs}
    >>>
    output {
      File outFrag = "Frag.bed"
      File outTEs = "TEs.bed"
    }
  runtime {
    docker: "ghcr.io/welch-lab/steamer:latest"
  }
}
