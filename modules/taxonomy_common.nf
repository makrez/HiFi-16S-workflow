process taxonomy_add_md5 {
    conda (params.enable_conda ? "$projectDir/env/jq.yml" : null)
    container "makrezdocker/alpine-jq:1.0"

    publishDir "${params.outdir}/final", mode: params.publish_dir_mode

    input:
    path tax_table
    val method

    output:
    path "${method}_tax_merged_freq_tax.tsv", emit: tax_with_md5

    script:
    """
    set -euo pipefail

    awk 'BEGIN{FS=OFS="\\t"}
    NR==1{
      print "id", \$0
      next
    }
    {
      cmd = "printf \\"%s\\" \\"" \$1 "\\" | md5sum"
      cmd | getline hashline
      close(cmd)
      split(hashline, a, " ")
      print a[1], \$0
    }' ${tax_table} > ${method}_tax_merged_freq_tax.tsv
    """
}

process taxonomy_summary {

    conda (params.enable_conda ? "$projectDir/env/Rdata_table.yml" : null)
    container "quay.io/biocontainers/r-data.table:1.12.2"
    tag "${prefix}"

    publishDir "${params.outdir}/final",
    mode: 'copy',
    overwrite: true

    input:
    path taxonomy_table
    val prefix

    output:
    path "${prefix}_taxonomy_assignment_summary.tsv",
        emit: assignment_summary

    path "${prefix}_taxonomy_deepest_rank_summary.tsv",
        emit: deepest_rank_summary

    script:
    """
    summarize_taxonomy.R \
        "${taxonomy_table}" \
        "${prefix}_taxonomy_assignment_summary.tsv" \
        "${prefix}_taxonomy_deepest_rank_summary.tsv"
    """
}
