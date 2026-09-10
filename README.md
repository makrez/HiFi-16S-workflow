This repository is a rewritten fork of the PacBio HiFi-16S workflow, redesigned for improved performance, modularity, and flexibility on HPC systems.

# Overview

This pipeline processes PacBio HiFi 16S amplicon sequencing data. It includes quality control, filtering, denoising with DADA2, and taxonomic assignment. The main output of the pipeline is a count table with taxonomic assignment.

The repository includes database download utilities for Naive Bayes classifier from dada2 and for vsearch.

The refactor focuses on:

- Improved parallelisation

- Separation of pipeline stages

- Compatibility wit PacBio Revio / Vega data (binned quality scores)

- Database management


# Pipeline Structure

```
Input FASTQ
   ↓
QC + Filtering
   ↓
Primer trimming
   ↓
Error model learning
   ↓
Denoising (independent)
   ↓
Sequence table construction
   ↓
Chimera removal
   ↓
ASV filtering
   ↓
Taxonomic assignment (per database)
```

# Installation and Usage

Clone the repository and adjust the `nextflow.config`. Further options are located in `conf/base.config`.

Before adjusting the config or downloading any data, you can run a test run with a tiny DB and some test data that are part of the repository:


```
# 16s
nextflow run main.nf -profile test_16s,docker # or signularity, or conda

# its
nextflow run main.nf -profile test_its,docker # or signularity, or conda

# euk_long 
nextflow run main.nf -profile test_euk_long,docker # or signularity, or conda
```

## Database download

All supported databases can be installed through the workflow:

```bash
nextflow run main.nf \
    --download_db true \
    --download_targets silva,gtdb,gg2,euk_ssu,euk_lsu,euk_long,euk_its
```

The location of the installed databases is controlled by `--db_base_dir`, also avaliable as a parameter in the `nextflow.config`.

The database manifest in:

```text
conf/databases.yml
```

defines the database versions, filenames, marker types, and official download
locations.

For example, an installation may have the following structure:

```text
databases/
├── silva/
│   └── nb/
├── gtdb/
│   └── nb/
├── gg2/
│   └── nb/
├── euk_ssu/
│   ├── nb/
│   └── vsearch/
├── euk_lsu/
│   ├── nb/
│   └── vsearch/
├── euk_its/
│   ├── nb/
│   └── vsearch/
└── euk_long/
    ├── nb/
    └── vsearch/
```

## Running the Pipeline

The basic command is like this:

```
nextflow run main.nf \
  --input samplesheet.tsv \
  --meta_data meta_data.tsv
  --outdir results \
  -profile docker

```


You can execute the pipeline on a cluster by using either a custom nextflow profile, or 
by adding a custom config. See `custom_slurm.config` for an example.

```
nextflow run main.nf \
  --input samplesheet.tsv \
  --metadata meta_data.tsv
  --outdir results \
  -profile docker
  -c custom_slurm.config
```


## Configuration


Important options to change before running the pipeline:

```
// QC / preprocessing
filterQ = 20
downsample = 0
skip_primer_trim = false
max_ee = 2
min_len = 500
max_len = 6000

// Primer sequences (5' -> 3')
// Can be either a single string or a list of strings.

forward_p = 'AGRGTTYGATYMTGGCTCAG'
reverse_p = 'AAGTCGTAACAAGGTARCY'

// Error model learning
error_model_reads_per_sample = 20 // 
learn_nbases = null
error_model = null // If provided, this error model is used.
binned_quality_scores = '3,10,17,22,27,35,40'

// Denoising parameters
band_size = 16
homopolymer_gap_penalty = 1
omegac = 1e-40

// Chimera removal parameters
chimera_method = 'consensus'
min_parent_fold = 1.0

// ASV filtering
min_asv_total_freq = 10
min_asv_sample = 0

// Databases
db_base_dir = ""
databases_yaml = "${projectDir}/conf/databases.yml"
db_to_prioritize = 'GG2'

download_db = false
download_targets = null

// VSEARCH LCA
vsearch_lca_id         = 0.90
vsearch_lca_query_cov  = 0.90
vsearch_lca_maxaccepts = 0
vsearch_lca_maxrejects = 0
vsearch_lca_cutoff     = 0.90

// Output
outdir = "results"
publish_dir_mode = "copy"

```

### Parameter Documentation

`--filterQ`: Minimum read quality threshold used during preprocessing. Reads with quality values below this threshold are removed. The default value is `20`.

`--downsample`: Maximum number of reads retained per sample during preprocessing. Setting this to `0` disables downsampling and keeps all reads that pass the filtering criteria. The default value is `0`.

`--skip_primer_trim`: If set to `true`, primer trimming is skipped. By default, primer sequences are removed before downstream processing. The default value is `false`.

`--max_ee`: Maximum number of expected errors allowed per read. Reads with more expected errors than this threshold are removed. The default value is `2`.

`--min_len`: Minimum read length retained after preprocessing and primer trimming. Reads shorter than this value are removed. The default value is `500` bp.

`--max_len`: Maximum read length retained after preprocessing and primer trimming. Reads longer than this value are removed. The default value is `6000` bp.

`--forward_p` and `--reverse_p`: Forward and reverse primer sequences in 5′ to 3′ orientation. These are used for primer trimming. Each parameter can contain either a single primer sequence or a list of primer sequences.

`--error_model_reads_per_sample`: Number of reads sampled from each sample for error-model learning. Reads from the individual samples are combined to construct the dataset used to estimate the DADA2 error model. Reducing this value can speed up error-model learning for datasets containing many samples. The default value is `20`.

`--learn_nbases`: Number of bases used by DADA2 for error-model learning. If set to `null`, the workflow determines the number of bases from the reads selected using `--error_model_reads_per_sample`. The default value is `null`.

`--error_model`: Path to a previously calculated DADA2 error model stored as an RDS file, for example `error.rds`. If provided, this error model is used and error-model learning is skipped. The default value is `null`.

`--binned_quality_scores`: Comma-separated quality-score values used when constructing the error model for sequencing data with binned quality scores, such as PacBio Revio and Vega data. The default value is `3,10,17,22,27,35,40`.

`--band_size`: Band size used by DADA2 during pairwise sequence alignment in the denoising step. Restricting the alignment to a band around the diagonal reduces computational cost while allowing small insertions and deletions. The default value is `16`.

`--homopolymer_gap_penalty`: Gap penalty used by DADA2 for insertions and deletions within homopolymer regions. A reduced penalty accommodates the higher frequency of indel errors in homopolymers in long-read sequencing data. The default value is `1`.

`--omegac`: DADA2 abundance p-value threshold used when determining whether a sequence is sufficiently abundant to be inferred as a new sequence variant rather than explained by sequencing errors. The default value is `1e-40`.

`--chimera_method`: Method used by DADA2 for chimera detection and removal. The default value is `consensus`, in which chimeras are identified independently in each sample and the evidence across samples is combined to determine whether an ASV should be removed.

`--min_parent_fold`: Minimum abundance of each potential parent sequence relative to a candidate chimera during chimera detection. A value of `1.0` allows parent sequences with at least the same abundance as the candidate chimera. The default value is `1.0`.

`--min_asv_total_freq`: Minimum total abundance required for an ASV across the complete dataset. ASVs with fewer reads than this threshold across all samples are removed. The default value is `10`.

`--min_asv_sample`: Minimum number of samples in which an ASV must be observed to be retained. Setting this to `0` disables filtering based on the number of samples containing the ASV. The default value is `0`.

`--db_base_dir`: Base directory where the reference databases are stored. This should point to the directory containing the databases defined in `--databases_yaml`.

`--databases_yaml`: YAML configuration file describing the available reference databases and their associated files. By default, the workflow uses `conf/databases.yml` included with the pipeline.

`--db_to_prioritize`: Reference database to prioritize when resolving ties between taxonomic assignments produced by the Naive Bayes classifier. If multiple databases result in equally supported assignments, the assignment from the specified database is selected. The default value is `GG2`.

`--download_db`: If set to `true`, reference databases are downloaded instead of using databases already present under `--db_base_dir`. The databases to download can be selected using `--download_targets`. The default value is `false`.

`--download_targets`: Specifies which reference databases should be downloaded when `--download_db` is enabled. If set to `null`, no explicit subset of download targets is specified. The default value is `null`.

`--vsearch_lca_id`: Minimum sequence identity required for a reference sequence to be considered a match during VSEARCH LCA classification. The default value of `0.90` requires at least 90% sequence identity.

`--vsearch_lca_query_cov`: Minimum fraction of the query sequence that must be covered by the alignment. The default value of `0.90` requires at least 90% of the query sequence to align to the reference sequence.

`--vsearch_lca_maxaccepts`: Maximum number of matching reference sequences accepted for each query before VSEARCH stops searching. Setting this to `0` allows an unlimited number of accepted matches. Together with `--vsearch_lca_maxrejects 0`, this causes VSEARCH to search the complete reference database.

`--vsearch_lca_maxrejects`: Maximum number of non-matching candidate reference sequences considered before VSEARCH stops searching for a query. Setting this to `0` disables this limit. Together with `--vsearch_lca_maxaccepts 0`, the complete reference database is searched.

`--vsearch_lca_cutoff`: Fraction of matching reference sequences required to support a taxonomic assignment when determining the lowest common ancestor (LCA). The default value of `0.90` requires 90% of the accepted matches to support the reported taxonomic lineage.


### Supported databases

The appropriate reference database depends on the marker being sequenced.

| Database       | Marker         | Typical target                 | Naive Bayes | VSEARCH-LCA |
| -------------- | -------------- | ------------------------------ | ----------: | ----------: |
| SILVA          | 16S rRNA       | Bacteria and Archaea           |           ✓ |           — |
| GTDB           | 16S rRNA       | Bacteria and Archaea           |           ✓ |           — |
| Greengenes2    | 16S rRNA       | Bacteria and Archaea           |           ✓ |           — |
| EUKARYOME SSU  | 18S / SSU rRNA | Eukaryotes                     |           ✓ |           ✓ |
| EUKARYOME LSU  | 28S / LSU rRNA | Eukaryotes                     |           ✓ |           ✓ |
| EUKARYOME ITS  | ITS            | Eukaryotes, particularly fungi |           ✓ |           ✓ |
| EUKARYOME long | SSU–ITS–LSU    | Long eukaryotic amplicons      |           ✓ |           ✓ |

The database identifiers used by the workflow are:

```text
silva
gtdb
gg2
euk_ssu
euk_lsu
euk_its
euk_long
```

### 16S rRNA data

For bacterial and archaeal 16S rRNA amplicons, the workflow supports three
reference databases:

* **SILVA**
* **GTDB**
* **Greengenes2**

These databases are currently used with the Naive Bayes taxonomy
classification workflow.

#### SILVA

SILVA provides a curated collection of aligned small- and large-subunit
ribosomal RNA sequences across all domains of life. In this workflow, the
SILVA SSU reference is used for taxonomic classification of bacterial and
archaeal 16S rRNA sequences.

The currently configured classifier is based on **SILVA 138.2**.

#### GTDB

The Genome Taxonomy Database (GTDB) provides a standardized genome-based
taxonomy for Bacteria and Archaea.

The workflow uses an officially published DADA2-compatible SSU reference
derived from **GTDB release R220**. GTDB is particularly useful when a
genome-based bacterial and archaeal taxonomy is desired.

#### Greengenes2

Greengenes2 provides an updated reference taxonomy integrating microbial
genome and marker-gene information.

The workflow currently uses the officially published
**Greengenes2 2024.09** DADA2-compatible reference.

#### Why VSEARCH-LCA is not used for these databases

The workflow does not construct its own VSEARCH/SINTAX versions of SILVA,
GTDB, or Greengenes2.

Although it is technically possible to convert sequence and taxonomy files
into a SINTAX-compatible FASTA, doing so would require the workflow to define
its own rules for taxonomy parsing, rank normalization, identifier matching,
and reference construction.

Instead, the workflow uses the officially published classifier-ready
references directly. Consequently, SILVA, GTDB, and Greengenes2 are currently
available through the Naive Bayes classification branch only.

### Eukaryotic marker data

Eukaryotic amplicons are classified using the **EUKARYOME** reference
database.

EUKARYOME provides marker-specific reference sets for several commonly used
eukaryotic amplicons. Importantly, EUKARYOME publishes databases specifically
prepared for different taxonomic classification algorithms.

The workflow therefore uses the official EUKARYOME DADA2 references for the
Naive Bayes classifier and the official EUKARYOME SINTAX references for
VSEARCH.

No conversion between these formats is performed by the workflow.

#### 18S / SSU

For 18S rRNA amplicon sequencing, use:

```text
euk_ssu
```

This corresponds to the EUKARYOME small-subunit (SSU) reference database.

Typical input consists of eukaryotic 18S rRNA amplicons. Both classification
methods are supported:

```text
euk_ssu/
├── nb/
│   └── DADA2_EUK_SSU_v2.1.fa.gz
└── vsearch/
    └── SINTAX_EUK_SSU_v2.1.fa.gz
```

#### 28S / LSU

For 28S rRNA amplicon sequencing, use:

```text
euk_lsu
```

This corresponds to the EUKARYOME large-subunit (LSU) reference database.

Both Naive Bayes and VSEARCH-LCA classification are supported:

```text
euk_lsu/
├── nb/
│   └── DADA2_EUK_LSU_v2.1.fa.gz
└── vsearch/
    └── SINTAX_EUK_LSU_v2.1.fa.gz
```




This repository is a rewritten fork of the official PacBio HiFi-16S workflow, redesigned for improved performance, modularity, and flexibility on HPC systems.

