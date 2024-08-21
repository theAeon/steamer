version 1.0
workflow run_full {
  input {
    File fullin_TEs
    Array[File] allcs
    Array[File] allcs_idx
    File chrom_size
    String fullin_sample_name
    Int memory_GB
  }
    parameter_meta {
      fullin_TEs: "Path to bed file containing TEs"
      allcs: "terra table column containing location of 'allc_*.tsv.gz'"
      allcs_idx: "terra table column containing location of 'allc.tsv.idx'"
      chrom_size: "Path to chrom.sizes obtained by UCSC fetchChromSizes.sh"
      fullin_sample_name: "name of sample"
      memory_GB: "memory, in gigabytes"
  }
}
