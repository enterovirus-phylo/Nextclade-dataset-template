# Ingest

Downloads sequences and metadata for one taxon from NCBI and curates them into
the two files the dataset build needs:

- `data/sequences.fasta`
- `data/metadata.tsv`

Curation normalises dates to ISO format, resolves place names to a consistent
region/country/division hierarchy, tidies author lists, and applies any manual
corrections you record. Both files are then copied up to the project's top-level
`data/` directory.

The main workflow runs all of this for you while `FETCH_SEQUENCES = True`, so
you normally do not need to invoke it directly. Read on if you want to run it on
its own or change what it collects.

---

## Layout

| Path | Contents |
|---|---|
| `snakefile` | The workflow itself. |
| `config/config.yaml` | Which taxon to fetch, which fields to keep, how to curate them. |
| `workflow/snakemake_rules/` | The rules, split by stage: fetch, curate, sort. |
| `source-data/` | Manual corrections: per-accession annotations and local geolocation rules. |
| `bin/` | Standalone helper scripts, kept for reference; the default route does not use them. |
| `data/` | Downloaded and intermediate files. Not version-controlled. |

---

## Configuration

Everything lives in [`config/config.yaml`](config/config.yaml).

### The one setting you must change

```yaml
ncbi_taxon_id: 42789
```

Find the numeric id for your virus in the
[NCBI Taxonomy browser](https://www.ncbi.nlm.nih.gov/taxonomy). The workflow
refuses to start while this is still a placeholder.

Pick a taxon at the level you actually want. A species-level id pulls in every
serotype below it; if you want one serotype, use its own id.

### Metadata fields

`ncbi_datasets_fields` lists what to request from NCBI. The
[NCBI Datasets field reference](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/dataformat/tsv/dataformat_tsv_virus-genome/)
documents the available names. `curate.field_map` then renames them to the
column names augur and the dataset build expect — so if you add a field, add it
to both.

The main workflow's tree colourings rely on `date`, `country` and `region`
being present, so keep those.

### Manual corrections

NCBI records contain mistakes, and the fix belongs in version control rather
than in a hand-edited TSV that the next download overwrites.

[`source-data/annotations.tsv`](source-data/annotations.tsv) overrides individual
fields: one correction per line, tab-separated, keyed on the **unversioned**
accession.

```
AY741551	country	Sierra Leone
AY603973	date	1961-XX-XX
```

[`source-data/geolocation-rules.tsv`](source-data/geolocation-rules.tsv) maps raw
NCBI location strings onto a consistent hierarchy, with `*` matching anything at
that level:

```
raw_region/raw_country/raw_division/raw_location	region/country/division/location
Africa/Cote d'Ivoire/*/*	Africa/Côte d'Ivoire/*/*
```

These are applied on top of the
[shared Nextstrain rules](https://github.com/nextstrain/ncov-ingest/blob/master/source-data/gisaid_geoLocationRules.tsv),
which the workflow downloads, so you only need entries for the cases those get
wrong.

> [!NOTE]
> Both files arrive pre-populated with corrections from the pathogen this
> template grew out of. They are keyed on specific accessions, so they do
> nothing for an unrelated virus, but you may want to clear them out to keep
> your own corrections easy to find.

---

## Running it on its own

From inside this directory:

```bash
snakemake --cores 9 all
```

Override a setting for one run without editing the file:

```bash
snakemake --cores 9 all --config ncbi_taxon_id=42789
```

Some scripts may need the executable bit, depending on how you checked the
repository out:

```bash
chmod +x ./bin/*
```

---

## How the download works

Two routes exist, and this template is set up for the first:

1. **NCBI Datasets** (default) — takes a taxon id, returns sequences and a
   metadata table in one package. Fast, and needs no custom parsing.
2. **Entrez** — takes a free-text search term and returns full GenBank records,
   which carry fields NCBI Datasets does not expose.

The Entrez route is left in place as a starting point but is **not wired up**:
it depends on a `vendored/fetch-from-ncbi-entrez` script that this template does
not ship, and on a parser you would have to write yourself. Use route 1 unless
you specifically need a field only GenBank records carry, in which case see
[`nextstrain/ingest`](https://github.com/nextstrain/ingest) for the vendored
scripts and the `git subrepo` workflow for keeping them current.

---

## Stages

| Rule | What it does |
|---|---|
| `fetch_ncbi_dataset_package` | Downloads the NCBI Datasets zip for the taxon. Retries on network errors. |
| `extract_ncbi_dataset_sequences` | Unpacks the genomic FASTA. |
| `format_ncbi_dataset_report` | Turns the metadata into TSV and splits the versioned accession into `accession` and `accession_version`. |
| `format_ncbi_datasets_ndjson` | Joins sequences and metadata into one NDJSON stream. |
| `curate` | Renames fields, normalises dates and place names, abbreviates authors, applies manual corrections. |
| `sort` | Removes duplicate sequences. |
| `metadata` | Keeps the metadata rows whose sequences survived deduplication. |

The curate step is a single pipeline of `augur curate` commands streaming NDJSON
between them. To add a step of your own, read NDJSON records from stdin and
write them to stdout, and slot the command into the pipeline in
[`workflow/snakemake_rules/curate.smk`](workflow/snakemake_rules/curate.smk).

---

## Troubleshooting

**`CONFIGURATION ERROR: ncbi_taxon_id ... is not set`**
Put your taxon id in `config/config.yaml`.

**`datasets: command not found`**
The NCBI Datasets CLI is missing. It is part of `environment.yml`, so activate
the `nextclade` environment, or install it on its own with
`conda install -c conda-forge ncbi-datasets-cli`.

**The download returns far fewer sequences than expected**
Check the taxon id resolves to what you think, and at the level you want:

```bash
datasets summary virus genome taxon 42789
```

**`N sequences have no metadata row`**
The curate step dropped records the download had provided. `logs/curate.txt`
says why — usually an unparseable date or a strain name that failed the regex.

**Dates come out as `XXXX-XX-XX`**
That is the marker for an unknown or ambiguous date, and is expected for older
records. If it affects sequences that do have a date in GenBank, add the format
to `curate.expected_date_formats`.

**A country or region is wrong or empty**
Add a rule to `source-data/geolocation-rules.tsv`, or correct the individual
record in `source-data/annotations.tsv`.
