process download_gtdb_db {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "makrezdocker/alpine-jq:1.0"

    publishDir "${params.gtdb_dir}/nb", pattern: "*.fa.gz", mode: "copy"
    publishDir "${params.gtdb_dir}/vsearch", pattern: "sequences.fasta", mode: "copy"
    publishDir "${params.gtdb_dir}/vsearch", pattern: "taxonomy.tsv", mode: "copy"

    label 'cpu_def'

    input:
    val nb_url
    val nb_filename
    val vsearch_seq_url

    output:
    path "${nb_filename}"
    path "sequences.fasta"
    path "taxonomy.tsv"

    script:
    """
    wget -O ${nb_filename} "${nb_url}"
    wget -O ssu_all_r220.fna.gz "${vsearch_seq_url}"

    gunzip -c ssu_all_r220.fna.gz > sequences.fasta

    zgrep "^>" ssu_all_r220.fna.gz | \\
      awk '{
        id=\$1
        sub(/^>/,"",id)
        tax=\$0
        sub(/^>[^ ]+ /,"",tax)
        sub(/ \\[.*/, "", tax)
        print id "\\t" tax
      }' > taxonomy.tsv
    """
}

process download_silva_db {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "makrezdocker/alpine-jq:1.0"

    publishDir "${params.silva_dir}/nb", pattern: "*.fa.gz", mode: "copy"
    publishDir "${params.silva_dir}/vsearch", pattern: "sequences.fasta", mode: "copy"
    publishDir "${params.silva_dir}/vsearch", pattern: "taxonomy.tsv", mode: "copy"

    label 'cpu_def'

    input:
    val nb_url
    val nb_filename
    val vsearch_seq_url
    val vsearch_tax_url

    output:
    path "${nb_filename}"
    path "sequences.fasta"
    path "taxonomy.tsv"

    script:
    """
    set -euo pipefail

    echo "Downloading SILVA NB database..."
    wget -O ${nb_filename} "${nb_url}"

    echo "Downloading SILVA VSEARCH sequences..."
    wget -O silva_sequences.fasta.gz "${vsearch_seq_url}"
    gunzip -c silva_sequences.fasta.gz > sequences.fasta

    echo "Downloading SILVA taxonomy mapping..."
    wget -O silva_taxmap.txt.gz "${vsearch_tax_url}"
    gunzip -c silva_taxmap.txt.gz > silva_taxmap.txt

    # Convert SILVA taxmap to simple TSV: sequence_id <tab> taxonomy
    # Expected SILVA taxmap format: accession, ..., taxonomy in column 3
    awk -F '\\t' 'BEGIN{OFS="\\t"} NF >= 3 {print \$1, \$3}' silva_taxmap.txt > taxonomy.tsv

    # Clean up intermediates
    rm -f silva_sequences.fasta.gz silva_taxmap.txt.gz silva_taxmap.txt
    """
}

process download_gg2_db {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "makrezdocker/alpine-jq:1.0"

    publishDir "${params.gg2_dir}/nb", pattern: "*.fa.gz", mode: "copy"
    publishDir "${params.gg2_dir}/vsearch", pattern: "sequences.fasta", mode: "copy"
    publishDir "${params.gg2_dir}/vsearch", pattern: "taxonomy.tsv", mode: "copy"

    label 'cpu_def'

    input:
    val nb_url
    val nb_filename
    val vsearch_seq_url
    val vsearch_tax_url

    output:
    path "${nb_filename}"
    path "sequences.fasta"
    path "taxonomy.tsv"

    script:
    """
    wget -O ${nb_filename} "${nb_url}"
    wget -O gg2_sequences.fna.gz "${vsearch_seq_url}"
    wget -O gg2_taxonomy.tsv.gz "${vsearch_tax_url}"

    gunzip -c gg2_sequences.fna.gz > sequences.fasta
    gunzip -c gg2_taxonomy.tsv.gz > taxonomy.tsv
    """
}
