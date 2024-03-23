version 1.0
import "prep.wdl" as prep
import "analysis.wdl" as analysis


workflow run_full {
  input {
    File fullin_TEs
    File fullin_Frags
    File fullin_barcode_list
    String fullin_sample_name
  }
  parameter_meta {
      fullin_TEs: "Path to .gz or .csv file containing TEs"
      fullin_Frags: "Path to .gz or .csv file containing Fragments"
      fullin_barcode_list: "Path to .gz or .csv file containing QC'd barcodes"
      fullin_sample_name: "name of sample"
  }
  call prep.steamerprep {
      input:
        in_TEs = fullin_TEs,
        in_Frags = fullin_Frags,
        in_QCbar = fullin_barcode_list
    }
  call analysis.runsteamer {
    input:
      in_bedfile = steamerprep.prepped,
      in_barcode_list = fullin_barcode_list,
      in_sample_name = fullin_sample_name
  }
  output {
    File FamMat = runsteamer.FamMat
    File UniqueMat = runsteamer.UniqueMat
    File UniqueDF = runsteamer.UniqueDF
    File UniqueBar = runsteamer.UniqueBar
    File FamDF = runsteamer.FamDF
    File FamBar = runsteamer.FamBar
  }
}
