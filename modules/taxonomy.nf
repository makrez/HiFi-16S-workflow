process taxonomy_nb_assign {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"

    publishDir "${params.outdir}/nb_tax", mode: params.publish_dir_mode

    memory { db_name == 'gtdb' ? 24.GB : 8.GB }


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
    """
}

process taxonomy_best {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    
    publishDir "${params.outdir}/results", mode: params.publish_dir_mode

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
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"

    publishDir "${params.outdir}/results", mode: params.publish_dir_mode

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
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "makrezdocker/alpine-jq:1.0"

    publishDir "${params.outdir}/results", mode: params.publish_dir_mode

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

