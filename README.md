# Nextclade dataset template

A Snakemake workflow for building a [Nextclade](https://github.com/nextstrain/nextclade)
dataset for a virus of your choice. It downloads sequences from NCBI, curates
the metadata, builds a reference tree, and packages everything Nextclade needs
into a single `dataset.zip`.

The template ships with picornavirus defaults, because that is what it grew out
of, but nothing in it is enterovirus-specific. Every value you need to change is
marked `CHANGE ME` in the `Snakefile`.

> [!NOTE]
> If you're not familiar with Nextclade, please read its [documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/) first.

---

## Contents

- [What you get](#what-you-get)
- [Prerequisites](#prerequisites)
- [Setup, step by step](#setup-step-by-step)
- [Running the workflow](#running-the-workflow)
- [Looking at the result](#looking-at-the-result)
- [Sub-workflows](#sub-workflows)
- [Configuration reference](#configuration-reference)
- [Troubleshooting](#troubleshooting)
- [Releasing a dataset](#releasing-a-dataset)
- [Authors and citation](#authors-and-citation)

---

## What you get

Running the workflow produces:

| Output | What it is |
|---|---|
| `dataset.zip` | The packaged Nextclade dataset.  |
| `out-dataset/` | The same files unpacked (used to serve locally to test in Nextclade Web) |
| `results/auspice.json` | The reference tree. Drag it onto [auspice.us](https://auspice.us) to inspect the build. |
| `results/example_sequences.fasta` | A handful of sequences shipped with the dataset as the "Load example" button. |

A Nextclade dataset is a reference sequence, its annotation, a reference tree,
clade definitions and a set of QC rules. Nextclade places new sequences onto
that tree, calls their mutations and assigns them a clade. The
[Nextclade documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/user/datasets.html)
explains the format in detail.

---

## Prerequisites

### Software

Everything comes from conda-forge and bioconda. `environment.yml` in the project
root lists all of it — create the environment from that file:

```bash
conda env create -f environment.yml
conda activate nextclade
```

If the environment already exists and the file has changed since, update it with
`conda env update -f environment.yml --prune`. Micromamba works as a replacement (`micromamba env create -f environment.yml`).

### An email address for NCBI

NCBI asks automated downloads to identify themselves, and throttles those that
do not. Put your email address in a `.env` file in the project root:

```bash
echo "EMAIL=you@example.org" > .env
```

`.gitignore` already keeps `.env` out of version control.

### Things to know before you start

You will need, for your virus:

- the **NCBI taxon id** — look it up in the [NCBI Taxonomy browser](https://www.ncbi.nlm.nih.gov/taxonomy);
- a **reference accession**, ideally the RefSeq entry, with its version suffix (e.g. `AY426531.1`);
- a rough idea of your **genome length**, to set a sensible minimum sequence length;
- optionally, a **clade nomenclature** and the mutations that define each clade. You can build a
  dataset without this and add clades later.

---

## Setup, step by step

### 1. Generate the reference files

Nextclade needs the reference as a FASTA file and its annotation as GFF3. This
script fetches both from GenBank and helps you tidy up the CDS names:

```bash
python3 scripts/generate_from_genbank.py --reference "<accession>" --output-dir dataset/
```

The script asks a few questions while it runs. When a genome is annotated both as a whole CDS and as mature peptides, the mature peptides are usually the more
useful choice (`mature_protein [2]`); otherwise pick `on a case-by-case basis [0]` and then `[product]` or `[leave empty for manual choice]` to choose the
features yourself. The right answers depend on how the reference is annotated.

CDS names end up as labels on amino acid mutations, so keep them short and recognisable — `VP1`, not `VP1_capsid_protein_precursor`. The script strips
boilerplate like a trailing "protein" and replaces characters that augur and Nextclade cannot handle, warning you when it does.

This writes three files:

- `dataset/reference.fasta` — the reference sequence
- `dataset/genome_annotation.gff3` — the annotation
- `dataset/reference.gbk` — the annotation in GenBank format

Copy the last one into `resources/`, where `augur ancestral` looks for it:

```bash
cp dataset/reference.gbk resources/reference.gbk
```

> [!IMPORTANT]
> The CDS names must match **exactly** in three places: the `Name` attributes in
> `dataset/genome_annotation.gff3`, the `product` qualifiers in
> `resources/reference.gbk`, and the `GENES` list in the `Snakefile`. A mismatch
> makes `augur ancestral` fail, because it cannot line the features up.
>
> Check them with:
> ```bash
> grep -o 'Name=[^;]*' dataset/genome_annotation.gff3 | sort -u
> grep '/product=' resources/reference.gbk | sort -u
> ```

### 2. Tell the ingest workflow which virus to download

Edit [`ingest/config/config.yaml`](ingest/config/config.yaml):

```yaml
ncbi_taxon_id: 42789          # your virus, from the NCBI Taxonomy browser
entrez_search_term: "your virus name"
```

Only `ncbi_taxon_id` matters for the default download route.

### 3. Describe your dataset in `pathogen.json`

Edit [`dataset/pathogen.json`](dataset/pathogen.json). At a minimum, fill in the
`attributes` block:

```json
"attributes": {
  "name": "Enterovirus D68",
  "reference name": "Fermon",
  "reference accession": "AY426531.1"
}
```

The rest of the file controls alignment and quality control. The defaults are
reasonable for a ~7.5 kb picornavirus; for a genome of a different size, scale
`qc.missingData.missingDataThreshold` and `qc.privateMutations.typical`/`cutoff`
accordingly. The
[pathogen config reference](https://docs.nextstrain.org/projects/nextclade/en/latest/user/input-files/05-pathogen-config.html)
documents every field.

> [!WARNING]
> QC only runs for the checks you enable here. If you drop the `qc` block,
> Nextclade will not check your sequences for quality issues.

### 4. Set the workflow parameters

Open the [`Snakefile`](Snakefile) and work through the blocks at the top.

```python
REFERENCE_ACCESSION = "AY426531.1"
```

Then sanity-check the filtering thresholds, especially `MIN_LENGTH` (the default
`6000` assumes a picornavirus-sized genome and will discard everything for a
smaller virus). See the [configuration reference](#configuration-reference) below.

### 5. Decide how to root the tree

For many Enteroviruses the official reference sequence was deposited decades ago and sits a long way from anything circulating today. Mutations called against it are numerous and not very informative. As an alternative, the workflow can root the
tree on a **Static Inferred Ancestor**: a sequence reconstructed from the most
recent common ancestor of your own sequences, obtained by rooting a tree on an
outgroup of related viruses. Nextclade can then report mutations relative to
that ancestor instead of the distant reference.

The template does not ship such an ancestor, so you have to choose:

**Option A — use the inferred ancestor (recommended).** Requires a bit of setup:
you need outgroup sequences. Follow
[`inferred-root/README.md`](inferred-root/README.md), then set in the `Snakefile`:

```python
STATIC_ANCESTRAL_INFERRENCE = True
INFERRENCE_RERUN = True          # only for the first build
ROOTING = "ancestral_sequence"
```

**Option B — skip it for now.** Gets you a working dataset immediately, and you
can switch to Option A later:

```python
STATIC_ANCESTRAL_INFERRENCE = False
ROOTING = "mid_point"
```

---

## Running the workflow

```bash
snakemake --cores 9 all
```

Adjust `--cores` to your machine. Tree building dominates the runtime; expect
tens of minutes to a few hours depending on `MAX_SEQS`.

Before committing to a full run, check what Snakemake intends to do:

```bash
snakemake --cores 9 all --dry-run
```

The workflow validates its configuration while parsing, so mistakes surface in
seconds rather than after an hour of tree building.

### The first run, with root inference

If you chose Option A above, the first build has to compute the ancestor. That
is destructive and slow, so it asks for confirmation:

```bash
snakemake --cores 9 all --config static_inference_confirmed=true
```

Afterwards set `INFERRENCE_RERUN = False`, so later runs reuse the ancestor in
`resources/inferred-root.fasta`. Commit that file: it is what lets someone else
reproduce your dataset without rerunning the inference.

### Adding clade labels

Clades are defined in [`resources/clades.tsv`](resources/clades.tsv): one row per
clade, listing the mutations that are unique to it. The template ships it empty,
so fill it in once you have done an initial run and can see the tree. Until it
holds real definitions, every mutation Nextclade finds in a query sequence counts
as private; afterwards, Nextclade can tell clade-defining mutations apart from
genuinely private ones.

After building the dataset, run:

```bash
snakemake --cores 9 mutLabels
```

The rule writes the labels into the dataset's `pathogen.json`. Use
`clades_high_threshold` in the `mutLabels` rule to control how specific a mutation
must be to a clade before it is labelled — raise it to label only your main
clades, lower it to include finer ones. Mutations that fall just short of the
threshold are collected in `results/newly_relevant_mutations.tsv` for you to
review.

### Starting over

```bash
snakemake --cores 1 clean
```

This removes `results/`, `out-dataset/`, `dataset.zip` and the downloaded data.
It deliberately leaves `resources/inferred-root.fasta` alone, since recomputing
the ancestor is a separate decision.

---

## Looking at the result

### The reference tree

Drag `results/auspice.json` onto [auspice.us](https://auspice.us). This is the
quickest way to spot problems: sequences on implausibly long branches are
usually misannotated or poor quality. Add their accessions to
`resources/exclude.txt` and rerun.

Alternatively, view the tree locally:

```bash
auspice view --datasetDir results
```
Or run:
```bash
snakemake -c 9 viz
```

What Auspice displays — colourings and filters — is controlled by [auspice_config.json](resources/auspice_config.json). Further configuration options are described in the [Nextclade reference-tree documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/user/input-files/04-reference-tree.html).

### The dataset in Nextclade Web

Serve `out-dataset/` locally and point Nextclade at it:

```bash
serve --cors out-dataset -l 3000
```

Or run:
```bash
snakemake -c 9 serve
```

Then open:

```text
https://master.clades.nextstrain.org/?dataset-url=http://localhost:3000
```

Click **Load example**, then **Run**. With a large tree, reducing
*Settings → Sequence view → Max. nucleotide markers* to 500 keeps the interface
responsive.

### The dataset in the Nextclade CLI

```bash
nextclade3 run \
    --input-dataset dataset.zip \
    --output-all test_out/ \
    results/example_sequences.fasta
```

`test_out/nextclade.tsv` then holds one row per sequence, with clade
assignments, QC verdicts and mutation calls. This is the check worth running
before any release: see [`Nextclade_Dataset_Checklist.md`](Nextclade_Dataset_Checklist.md)
for what else to test, and [`nextclade_cli.md`](nextclade_cli.md) for the full
CLI reference.

---

## Sub-workflows

### `ingest/` — getting the data

Downloads sequences and metadata for your taxon from NCBI, normalises dates and
place names, and writes `data/sequences.fasta` and `data/metadata.tsv`. The main
workflow calls it for you while `FETCH_SEQUENCES = True`.

To supply your own data instead, set `FETCH_SEQUENCES = False` and put the two
files in `data/` yourself. The metadata needs at least `accession`, `strain`,
`date`, `country` and `region` columns.

See [`ingest/README.md`](ingest/README.md).

### `inferred-root/` — reconstructing the ancestor

Builds a tree of your sequences plus an outgroup, roots it on the outgroup, and
reconstructs the sequence at the MRCA of your sequences. Run by the main
workflow when `INFERRENCE_RERUN = True`.

See [`inferred-root/README.md`](inferred-root/README.md).

---

## Configuration reference

All of these live at the top of the [`Snakefile`](Snakefile).

### Required

| Parameter | Description |
|---|---|
| `REFERENCE_ACCESSION` | GenBank accession of your reference, with version suffix. |
| `GENES` | CDS names to translate. Must match your annotation and GenBank file exactly. |

### Filtering and subsampling

| Parameter | Default | Description |
|---|---|---|
| `MIN_DATE` | `"1960-01-01"` | Drop sequences collected earlier. `""` disables. |
| `MIN_LENGTH` | `"6000"` | Drop shorter sequences. Scale to your genome; `""` disables. |
| `MAX_SEQS` | `"10000"` | Upper bound on tree size. Lower it for a faster build. |
| `ALLOWED_DIVERGENCE` | `"3000"` | Substitution count above which a sequence is treated as an outlier and excluded automatically. Lower it once you know what normal looks like for your virus. |
| `ROOTING` | `"ancestral_sequence"` | Tree root: `"ancestral_sequence"`, `"mid_point"`, `"outgroup"`, or an accession. |
| `ID_FIELD` | `"accession"` | Metadata column holding sequence IDs: `"accession"` or `"strain"`. |

### Switches

| Parameter | Default | Description |
|---|---|---|
| `FETCH_SEQUENCES` | `True` | Download data via `ingest/`. Set `False` to supply your own. |
| `STATIC_ANCESTRAL_INFERRENCE` | `True` | Include the inferred ancestor in the tree. |
| `INFERRENCE_RERUN` | `False` | Recompute the ancestor. Needed for the first build; slow and destructive. |

### Files you will edit

| File | Purpose |
|---|---|
| `dataset/pathogen.json` | Dataset attributes, alignment parameters, QC thresholds. |
| `dataset/README.md` | Description shown to users in Nextclade. |
| `dataset/CHANGELOG.md` | Release history. |
| `resources/clades.tsv` | Clade-defining mutations, one per row. |
| `resources/exclude.txt` | Accessions to keep out of the tree. |
| `resources/include.txt` | Accessions that must survive subsampling. |
| `resources/include_examples.txt` | Accessions to force into the example set. |
| `resources/accession_strain.tsv` | Corrected strain names, where NCBI's are wrong. |
| `resources/auspice_config.json` | Tree display: colourings, filters, reference nodes. |

---

## Troubleshooting

**`CONFIGURATION ERROR: REFERENCE_ACCESSION is still the template placeholder`**
Work through [Setup](#setup-step-by-step). The workflow checks its configuration
up front rather than failing deep into a run; the message names what to fix.

**`Missing reference FASTA: dataset/reference.fasta`**
Run `scripts/generate_from_genbank.py` — step 1 of the setup.

**`Missing inferred ancestral root: resources/inferred-root.fasta`**
Nothing has computed the ancestor yet. Either set `INFERRENCE_RERUN = True` and
confirm the run, or take Option B in [step 5](#5-decide-how-to-root-the-tree).

**`augur ancestral` fails on a missing feature**
The CDS names disagree between `GENES`, the GFF3 and the GenBank file. Compare
them with the `grep` commands in [step 1](#1-generate-the-reference-files).

**`augur filter` returns no sequences**
Usually `MIN_LENGTH` is longer than your genome, or `MIN_DATE` postdates your
sequences. Check what you actually downloaded:

```bash
seqkit stats data/sequences.fasta
cut -f1,4 data/metadata.tsv | head
```

**The tree is dominated by one long branch**
One or a few bad sequences. Find them in `results/outliers.txt`, or open
`results/auspice.json` on auspice.us, and add the accessions to
`resources/exclude.txt`.

**`nextclade3: command not found`**
See [Prerequisites](#prerequisites) for the symlink.

**Snakemake says the directory is locked**
A previous run was interrupted. `snakemake --unlock` clears it.

Still stuck? [Open an issue](https://github.com/enterovirus-phylo/dataset-template-inferred-root/issues),
or check the [Nextclade documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/index.html).

---

## Releasing a dataset

[`Nextclade_Dataset_Checklist.md`](Nextclade_Dataset_Checklist.md) walks through
what to verify before publishing, and how to submit a dataset to
[`nextclade_data`](https://github.com/nextstrain/nextclade_data) so it appears in
the Nextclade dataset list.

---

## Authors and citation

Maintainers: Nadia Neuner-Jehle, Alejandra González-Sánchez and Emma B. Hodcroft
([eve-lab.org](https://eve-lab.org/)).

Questions and suggestions are welcome as
[issues](https://github.com/enterovirus-phylo/dataset-template-inferred-root/issues),
or by email to <eve-group@swisstph.ch>.

If you use this template, we would appreciate a citation:

> Neuner-Jehle, N., González Sánchez, A., Hodcroft, E. B., & European Non-Polio
> Enterovirus Network (ENPEN). (2025). *enterovirus-phylo/nextclade_d68:
> Enterovirus D68 Nextclade Dataset v1.0.0* (v1.0.0--2025-11-18). Zenodo.
> https://doi.org/10.5281/zenodo.17642338

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17642338.svg)](https://doi.org/10.5281/zenodo.17642338)
