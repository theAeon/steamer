version 1.0
import "prep.wdl" as prep
import "analysis.wdl" as analysis


workflow run_full {
  input {
    File? fullin_TEs
    File? fullin_Frags
    File? fullin_Intersected
    File? fullin_barcode_list
    String fullin_sample_name
    String Mode = "ATAC"
    Int memory_GB
  }
  parameter_meta {
      fullin_TEs: "Path to bed file containing TEs"
      fullin_Frags: "Path to bed file containing Fragments"
      fullin_Intersected: "Pre-intersected bed file"
      fullin_barcode_list: "Path to .gz or .csv file containing QC'd barcodes"
      fullin_sample_name: "name of sample"
      memory_GB: "memory, in gigabytes"
  }
  if (!defined(fullin_Intersected)) {
      call prep.intersect_bedfiles {
          input:
            TEs = fullin_TEs,
            Frags = fullin_Frags,
            mem = memory_GB
        }
      call analysis.runsteamer as prepped {
          input:
            in_bedfile = intersect_bedfiles.intersected,
            in_barcode_list = fullin_barcode_list,
            in_sample_name = fullin_sample_name,
            in_mem = memory_GB
        }
    }
  if (defined(fullin_Intersected)) {
      call analysis.runsteamer as premade {
          input:
            in_bedfile = fullin_Intersected,
            in_barcode_list = fullin_barcode_list,
            in_sample_name = fullin_sample_name,
            in_mem = memory_GB
        }
  }
  output {
    Array[File] Mat = [select_first([prepped.FamMat, premade.FamMat]),
                      select_first([prepped.UniqueMat, premade.UniqueMat])]
    Array[File] DF = [select_first([prepped.UniqueDF, premade.UniqueDF]),
                      select_first([prepped.FamDF, premade.FamDF])]
    Array[File] Bar = [select_first([prepped.UniqueBar, premade.UniqueBar]),
                       select_first([prepped.FamBar, premade.FamBar])]
  }
}
  task run_rna {
    input {
      File index
      File t2gs
      String SampleName
      Int mem
      Array[File]+ fastqs
    }
    command <<<
      kb count -i ~{index} -g ~{t2gs} --h5ad -x 10XV3 -o ~{SampleName} -m ~{mem}G --verbose ${sep=' ' fastqs}
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
  }
  }
