process taxonomy_nb_assign {
    label 'highcpu'

    conda (params.enable_conda ? "$projectDir/env/dada2.yml" : null)
    container "quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"

    publishDir "${params.outdir}/nb_tax", mode: params.publish_dir_mode

    input:
    tuple path(asv_fasta), val(db_name), path(db_fasta), path(assign_script)

    output:
    tuple val(db_name), path("${db_name}_nb.tsv"), emit: nb_tax

    script:
    """
    Rscript --vanilla ${assign_script} \\
      ${asv_fasta} \\
      ${db_fasta} \\
      ${db_name} \\
      ${task.cpus}

    harmonise_taxonomy.R \\
      ${db_name}_nb.tsv \\
      ${db_name} \\
      ${db_name}_nb.harmonised.tmp.tsv

    mv ${db_name}_nb.harmonised.tmp.tsv ${db_name}_nb.tsv
    """
}

process taxonomy_vsearch_assign {
    label 'highcpu'

    conda (params.enable_conda ? "$projectDir/env/vsearch.yml" : null)
    container "quay.io/biocontainers/vsearch:2.8.0--hfc679d8_0"

    publishDir "${params.outdir}/vsearch_tax", mode: params.publish_dir_mode

    input:
    tuple path(asv_fasta), val(db_name), path(vsearch_fasta), path(vsearch_taxonomy), path(assign_script)

    output:
    tuple val(db_name), path("${db_name}_vsearch.tsv"), emit: vsearch_tax

    script:
    """
    bash ${assign_script} \\
      ${asv_fasta} \\
      ${vsearch_fasta} \\
      ${vsearch_taxonomy} \\
      ${db_name} \\
      ${task.cpus} \\
      ${params.maxreject} \\
      ${params.maxaccept} \\
      ${params.vsearch_identity} \\
      ${db_name}_vsearch.tsv
    """
}

process taxonomy_best {
    conda (params.enable_conda ? "$projectDir/env/Rdata_table.yml" : null)
    container "quay.io/biocontainers/r-data.table:1.12.2"
    
    publishDir "${params.outdir}/final", mode: params.publish_dir_mode

    input:
    path nb_tax_files
    val db_priority
    path best_script

    output:
    path "best_taxonomy.tsv", emit: best_tax
    path "best_taxonomy_withDB.tsv", emit: best_tax_with_db

    script:
    """
    Rscript --vanilla ${best_script} \\
      ${nb_tax_files.join(' ')} \\
      ${db_priority}
    """
}

process merge_taxonomy_with_table {
    conda (params.enable_conda ? "$projectDir/env/Rdata_table.yml" : null)
    container "quay.io/biocontainers/r-data.table:1.12.2"

    publishDir "${params.outdir}/final", mode: params.publish_dir_mode

    input:
    path taxonomy
    path asv_table
    path merge_script

    output:
    path "best_tax_merged_no_id.tsv", emit: merged_no_id

    script:
    """
    Rscript --vanilla ${merge_script} \\
      ${taxonomy} \\
      ${asv_table} \\
      best_tax_merged_no_id.tsv
    """
}

process add_md5_to_taxonomy_table {
    conda (params.enable_conda ? "$projectDir/env/jq.yml" : null)
    container "makrezdocker/alpine-jq:1.0"

    publishDir "${params.outdir}/final", mode: params.publish_dir_mode

    input:
    path merged_tax_table

    output:
    path "best_tax_merged_freq_tax.tsv", emit: best_tax_merged_freq_tax

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
    }' ${merged_tax_table} > best_tax_merged_freq_tax.tsv
    """
}

