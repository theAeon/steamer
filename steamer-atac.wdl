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
    Array[File] Mat = select_all([prepped.FamMat, premade.FamMat, prepped.UniqueMat, premade.UniqueMat])
    Array[File] DF =  select_all([prepped.UniqueDF, premade.UniqueDF, prepped.FamDF, premade.FamDF])
    Array[File] Bar = select_all([prepped.UniqueBar, premade.UniqueBar, prepped.FamBar, premade.FamBar])
  }
}
