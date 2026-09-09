This repository is a rewritten fork of the official PacBio HiFi-16S workflow, redesigned for improved performance, modularity, and flexibility on HPC systems.

# Overview

This pipeline processes PacBio HiFi 16S amplicon sequencing data using a modular Nextflow DSL2 workflow. It includes quality control, filtering, denoising with DADA2, and taxonomic assignment. The main output of the pipeline is a count table with taxonomic assignment.

Please create an issue for bugs and feature requests.

The refactor focuses on:

- Improved parallelisation

- Separation of pipeline stages

- Compatibility wit PacBio Revio data (binned quality scores)

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

# Usage

## Test pipeline

```
nextflow run main.nf -profile test_16s,docker # or signularity, or conda
```
## Basic Command

```
nextflow run main.nf \
  --input samplesheet.tsv \
  --meta_data meta_data.tsv
  --outdir results \
  -profile conda

```

## Cluster execution

You can execute the pipeline on a cluster by using either a custom nextflow profile, or 
by adding a custom config. See `custom_slurm.config` for an example.

```
nextflow run main.nf \
  --input samplesheet.tsv \
  --metadata meta_data.tsv
  --outdir results \
  -profile conda
  -c custom_slurm.config
```

## Database download

All supported databases can be installed through the workflow:

```bash
nextflow run main.nf \
    --download_db true \
    --download_targets silva,gtdb,gg2,euk_ssu,euk_lsu,euk_long,euk_its
```

The location of the installed databases is controlled by `--db_base_dir`.

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

This workflow supports taxonomic classification of PacBio HiFi amplicon
sequence variants (ASVs) using reference databases appropriate for the
sequenced marker.

A deliberate design principle of the workflow is that **only officially
published reference databases or officially published classifier-specific
derivatives are used**. This keeps the taxonomic references reproducible and avoids introducing
pipeline-specific decisions into the construction of reference databases.

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



