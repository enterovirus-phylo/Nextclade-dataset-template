# Inferred ancestral root

This sub-workflow reconstructs a **Static Inferred Ancestor**: a sequence at the
most recent common ancestor (MRCA) of your sequences, obtained by rooting a tree
on an outgroup of related viruses.

## Why bother

For many viruses the official reference was deposited decades ago and has since
drifted far from anything circulating. Mutations called against it are numerous,
largely shared by every modern sequence, and not very informative.

An inferred ancestor sits at the base of your own diversity instead. Mutations
called against it are the ones that actually distinguish your sequences from
each other. The reference still defines the coordinate system — positions do not
change — but Nextclade can report mutations relative to the ancestor as well,
which is what the `ref_nodes` block in
[`../resources/auspice_config.json`](../resources/auspice_config.json) configures.

The sequence is *static*: computed once, committed to the repository, and reused
by every later build. That is deliberate. A root that shifted each time new data
arrived would silently change every mutation call in the dataset.

---

## When to run it

- **First build.** The template ships no ancestor, so you have to compute one.
- **After adding a lot of new data**, especially sequences from a lineage that
  was not represented before and might pull the root.
- **After changing the outgroup.**

Not otherwise. Between those occasions, leave `INFERRENCE_RERUN = False` in the
top-level Snakefile and let the workflow reuse `resources/inferred-root.fasta`.

---

## Setup

### 1. Choose an outgroup

You need two or three viruses that are related to yours but clearly outside the
diversity you are describing — for enteroviruses, typically other serotypes in
the same species, or a sister species.

The choice matters. Too distant, and the alignment degrades to the point where
the reconstruction is meaningless. Too close, and the root lands inside your own
diversity rather than below it, which defeats the purpose.

A reasonable way to find candidates: BLAST your reference, and take hits that
are clearly homologous but fall outside the group you are building the dataset
for.

### 2. Download the outgroup sequences

One FASTA file per accession, in `resources/outgroup/`, named after the
accession:

```text
resources/outgroup/
├── KC631740.1.fasta
├── MH341887.1.fasta
└── MT703464.1.fasta
```

```bash
cd resources/outgroup
datasets download virus genome accession KC631740.1 --filename tmp.zip
unzip -jp tmp.zip ncbi_dataset/data/genomic.fna > KC631740.1.fasta && rm tmp.zip
```

The FASTA header inside each file should also be the accession, because
`pick_ancestral_sequence.py` identifies the outgroup tips by name. See
[`resources/outgroup/README.md`](resources/outgroup/README.md).

### 3. Configure the sub-workflow

Edit the settings block at the top of [`snakefile`](snakefile):

| Parameter | Default | Description |
|---|---|---|
| `REFERENCE_ACCESSION` | placeholder | Must match the top-level Snakefile. |
| `OUTGROUP` | placeholders | Your outgroup accessions; each needs a matching FASTA. |
| `MIN_DATE` | `"1990-01-01"` | Earliest collection date to include. |
| `MIN_LENGTH` | `"6000"` | Minimum sequence length in bp. Scale to your genome. |
| `MAX_SEQS` | `"1000"` | Sequences to subsample for the tree. A thousand is usually plenty; the root position is not sensitive to more. |
| `ID_FIELD` | `"accession"` | Metadata column holding sequence IDs. |

The workflow checks all of this before starting and names anything missing,
including every outgroup FASTA it cannot find.

---

## Running it

Normally you do not run this directly. Set `INFERRENCE_RERUN = True` in the
top-level Snakefile and run the main workflow from the project root:

```bash
snakemake --cores 9 all --config static_inference_confirmed=true
```

The confirmation flag is required because the run overwrites
`resources/inferred-root.fasta`.

To run it on its own — useful when experimenting with different outgroups — the
sub-workflow needs files the main workflow prepares, so generate those first:

```bash
# from the project root
snakemake --cores 9 results/metadata.tsv results/include.txt

cd inferred-root
snakemake --cores 9 all_sub
```

---

## What it does

| Rule | Step |
|---|---|
| `index_sequences`, `filter` | Subsample the ingroup by date, length and `MAX_SEQS`. |
| `join_fastas` | Add the outgroup, dropping duplicate headers. |
| `align_fastas` | Align everything with MAFFT. Nextclade's reference-based alignment is no good here, because the outgroup is too divergent from the reference. |
| `tree` | Build a maximum-likelihood tree with `augur tree` (IQ-TREE). |
| `pick_ancestral_sequence` | Reroot on the outgroup, reconstruct the sequence at the MRCA of the ingroup with TreeTime, and fill gaps from the reference. |
| `all_sub` | Copy the result to `../resources/inferred-root.fasta`. |

The gap-filling step matters: a reconstruction with gaps in it cannot serve as a
Nextclade reference. Columns where the ancestor is ambiguous take the reference
base, and insertions relative to the reference are dropped, so the output keeps
the reference coordinate system.

Details are in [`../scripts/pick_ancestral_sequence.py`](../scripts/pick_ancestral_sequence.py).

---

## Checking the result

### Is the tree rooted where you think?

Open `results/nwk_tree_outgroup.png`. You want to see the outgroup sequences at
the base, and your own sequences forming a single clade above them. If the
outgroup is scattered through the tree, it is too close to your ingroup and the
reconstruction is not trustworthy.

### Is the sequence sane?

```bash
seqkit stats ../resources/inferred-root.fasta
```

Expect a length close to your reference, and no gaps. A much shorter sequence
usually means the alignment went wrong.

### Does it look like an ancestor?

The point of the exercise is that the ancestor should be closer to your
sequences than the reference is. Worth confirming after the main build: in the
Auspice tree, the ancestor tip should sit near the root of your diversity, not
out on a branch of its own.

---

## Requirements

Beyond the main workflow's dependencies: **MAFFT**, **TreeTime**, and the Python
packages **biopython**, **typer** and **matplotlib**. All are in the conda
environment described in the [main README](../README.md#prerequisites).

---

## Troubleshooting

**`CONFIGURATION ERROR: Missing outgroup sequences`**
The filenames in `resources/outgroup/` do not match the `OUTGROUP` list. The
message lists both what is missing and what is present.

**`These outgroup sequences are not in the tree`**
The files exist, but the FASTA *headers* inside them differ from the accessions
in `OUTGROUP`. Check with `grep '>' resources/outgroup/*.fasta`.

**`CONFIGURATION ERROR: Missing ... ../results/metadata.tsv`**
You are running the sub-workflow standalone without the files the main workflow
produces. Generate them first, as shown above.

**`augur filter` produces no sequences**
`MIN_LENGTH` or `MIN_DATE` is too strict for your data. Check what you have with
`seqkit stats ../data/sequences.fasta`.

**MAFFT fails, or the alignment is very long**
The outgroup is too divergent, or covers a different genomic region than your
sequences. Use closer relatives, and make sure the outgroup FASTAs are whole
genomes if your ingroup is.

**`No ingroup tips left after removing the outgroup and the reference`**
The tree contains only outgroup and reference sequences, so `augur filter`
returned nothing usable. See the filtering point above.

**The inferred sequence has many ambiguous sites**
Some is expected in divergent regions, and those positions take the reference
base. A lot of it suggests poor-quality or misaligned sequences in the ingroup;
add the offenders to `../resources/exclude.txt` and rerun.

---

## Good practice

- **Commit `resources/inferred-root.fasta`.** It is what lets someone else
  reproduce your dataset without repeating the inference.
- **Record why you chose those outgroups**, as a comment next to the `OUTGROUP`
  list. Future you will want to know.
- **Test the root's stability** by rerunning with a different `MAX_SEQS`. If the
  ancestor changes noticeably, the outgroup is probably too close.
- **Keep `resources/static_inferred_metadata.tsv` in step** with your metadata
  columns. It supplies the one metadata row describing the ancestor, and its
  `accession` must stay `ancestral_sequence` to match the FASTA header.

---

## Further reading

- **Augur documentation:** https://docs.nextstrain.org/projects/augur/
- **TreeTime documentation:** https://treetime.readthedocs.io/

---

## Authors & Contact

- Maintainers: Nadia Neuner-Jehle, Alejandra González Sánchez and Emma B. Hodcroft ([eve-lab.org](https://eve-lab.org/))
- For questions or suggestions, please [open an issue](https://github.com/enterovirus-phylo/dataset-template-inferred-root/issues)