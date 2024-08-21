version 1.0
workflow run_full {
  input {
    File fullin_TEs
    Array[File] allcs
    File chrom_size
    String fullin_sample_name
    Int memory_GB
  }
}
