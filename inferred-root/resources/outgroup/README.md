# Outgroup sequences

Put one FASTA file here per outgroup accession, named after the accession:

```bash
resources/outgroup/
├── KC631740.1.fasta
├── MH341887.1.fasta
└── MT703464.1.fasta
```

The filenames must match the `OUTGROUP` list in [`../../snakefile`](../../snakefile)
exactly, and the FASTA header inside each file should be the accession too, since
`pick_ancestral_sequence.py` identifies the outgroup tips by name.

A quick way to fetch one:

```bash
datasets download virus genome accession KC631740.1 --filename tmp.zip
unzip -jp tmp.zip ncbi_dataset/data/genomic.fna > KC631740.1.fasta && rm tmp.zip
```

## Choosing an outgroup

Pick two or three viruses that are related to yours but clearly outside the
diversity you want to describe — for enteroviruses, typically other serotypes
within the same species, or a sister species.

Too distant, and the alignment falls apart; too close, and the inferred root
lands inside your own diversity instead of below it. After a run, open
`../../results/nwk_tree_outgroup.png` and check that the outgroup really sits at
the base of the tree.
