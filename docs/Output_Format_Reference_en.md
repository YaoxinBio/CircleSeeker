# CircleSeeker Output Format Reference

Step 16 writes the following results directly into the directory supplied with `-o`. `<prefix>` is the sample prefix supplied with `-p`.

```text
<output>/
├── <prefix>_eccDNA_summary.csv
├── <prefix>_eccDNA_regions.csv
├── <prefix>_eccDNA_reads.csv
├── <prefix>_eccDNA_all.fasta
├── UeccDNA/   (<prefix>_uecc.bed, <prefix>_uecc.fasta)
├── MeccDNA/   (<prefix>_mecc_sites.bed, <prefix>_mecc.fasta)
├── CeccDNA/   (<prefix>_cecc_segments.bed, <prefix>_cecc_junctions.bedpe, <prefix>_cecc.fasta)
├── <prefix>_report.html
└── <prefix>_summary.txt
```

The packager creates all three type directories. It writes BED/BEDPE files when the corresponding regions are present, writes each type's FASTA when sequences are available, and copies the reports produced by the preceding step.

## Packaged CSV files

`<prefix>_eccDNA_summary.csv` has one row per eccDNA. Its header is

```csv
eccDNA_id,type,state,chr,start,end,strand,length,location,segment_count,read_count,copy_number,confidence,inferred_reads
```

| Fields | Interpretation |
| --- | --- |
| `eccDNA_id`, `type`, `state` | IDs start at 1 within each type using at least four digits, e.g. `UeccDNA0001`. `type` is `Uecc`, `Mecc` or `Cecc`. `state` is `Confirmed` or `Inferred`. Use `eccDNA_id` to join files from the same run. |
| `chr`, `start`, `end`, `strand` | Coordinates are zero-based, half-open. Uecc uses its source locus. Mecc shows a representative candidate. A Cecc row has `chr=multi` and `.` for the other three fields. |
| `length`, `location`, `segment_count` | `length` is in bp. Uecc has one locus. For Mecc, `\|` separates alternative source loci for **one** sequence and `segment_count` is 1. For Cecc, `;` separates the constituent segments in traversal order and `segment_count` counts these segments. |
| `read_count`, `inferred_reads` | `read_count` is total supporting reads. A Confirmed call can receive additional split-read support when overlapping calls are merged. `inferred_reads` records the portion added to that call. For a standalone Inferred call, use `read_count` for its split-read support. |
| `copy_number` | Confirmed calls use the repeat/copy estimate and can receive the estimate of an overlapping Inferred call. Standalone Inferred calls use rounded coverage as a proxy, with a fallback of 1 for a missing or nonpositive estimate. |
| `confidence` | Confirmed calls carry the upstream score when available. The built-in SplitReads-Core adapter supplies fixed compatibility values of `prob_present=0.95` and `prob_artifact=0.0`. Inferred `confidence` inherits `prob_present`. A value of 0.95 provides no candidate-specific probability estimate or ranking. |

`<prefix>_eccDNA_regions.csv` lists one row per locus or segment.

```csv
eccDNA_id,region_idx,chr,start,end,strand,length,role
```

`region_idx` starts at 1. Uecc has role `source`. Mecc has `primary` and `candidate` loci. The origin remains ambiguous across the candidate loci. Cecc uses `head`, `middle` and `tail` for the ordered segments.

`<prefix>_eccDNA_reads.csv` lists distinct `(eccDNA_id, read_name)` pairs from the Confirmed evidence files. Its header is

```csv
eccDNA_id,read_name,copy_number,eccDNA_copy_number,match_degree,identity,mapq
```

The two copy-number columns refer to the read and the Confirmed eccDNA, respectively. Missing per-read estimates can be blank. Split-read support added by unification is counted in the summary, so counting read-detail rows can give a smaller total. With no named Confirmed reads this CSV can be empty and have no header.

## Packaged BED, BEDPE and FASTA

The BED files contain six tab-separated fields in the order `chrom, start, end, name, score, strand`, without a header. The score is fixed at 1. The name is the ID for Uecc, `ID|role` for Mecc and `ID|segI/N|role` for Cecc. The Cecc BEDPE contains `chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2`. Each row links a strand-aware segment exit to the next entry using one-base intervals, including the closing junction from the last segment to the first.

The combined and per-type FASTA files use `>ID|location|length=N|state=STATE`. A Mecc record can have additional `|` characters inside its location. Read the ID up to the first `|` and the two metadata fields from the right, or look up the ID in the regions CSV.

```text
>MeccDNA0001|chr1:100-200(+)|chr2:300-400(-)|length=100|state=Confirmed
```

## Intermediate files and fallback output

The following sections describe per-caller files. Use `--keep-tmp` to retain intermediate results. When the unified CSV is missing, Step 16 can copy available caller files into `<prefix>_Confirmed_UeccDNA/`, `<prefix>_Confirmed_MeccDNA/`, `<prefix>_Confirmed_CeccDNA/` and `<prefix>_Inferred_eccDNA/`. In this fallback, `<prefix>_merged_output.csv` is a copy of the available **Confirmed** summary.

The intermediate `reads_count` records read support. `Hit_count` in the Confirmed Mecc table counts candidate sites and defaults to 1 for Confirmed Uecc/Cecc. Inferred Uecc/Cecc rows use the split-read count for both fields.

## 2. Confirmed intermediate files

### 2.1 UeccDNA

#### `*_UeccDNA.core.csv`

| Field | Description |
|-------|-------------|
| `eccDNA_id` | Unique ID |
| `chr` / `start0` / `end0` / `strand` | 0-based coordinates and strand |
| `length` | Length (bp) |
| `match_degree` | Match degree (0-100) |
| `confidence_score` / `query_cov_best` / `query_cov_2nd` | Evidence scores (may be empty) |
| `mapq_best` / `identity_best` | Best alignment quality (may be empty) |
| `low_mapq` / `low_identity` | Low-quality flags (may be empty) |
| `copy_number` / `repeat_number` | Copy/repeat estimates |
| `eccdna_type` | `Uecc` |
| `num_merged` / `merged_from_ids` | Merge bookkeeping |
| `reads_count` / `read_name` | Read count and list |

#### `*_UeccDNA.bed`

Columns: `chrom, chromStart, chromEnd, name, score, strand, length, eccdna_type, repeat_number, match_degree`

#### `*_UeccDNA_C.fasta`

Header format:
```
>UeccDNA1|chr:start-end(strand)|length=...|repeats=...|reads=...
```

---

### 2.2 MeccDNA

#### `*_MeccSites.core.csv`

Adds the following fields on top of Uecc core fields:

| Field | Description |
|-------|-------------|
| `hit_index` | Site index (1-based) |
| `hit_count` | Total sites for the eccDNA |

#### `*_MeccSites.bed`

Columns: `chrom, chromStart, chromEnd, name, score, strand, length, eccdna_type, copy_number`

#### `*_MeccBestSite.bed`

Columns (legacy names): `chrom, chromStart, chromEnd, name, score, strand, eLength, eClass, copyNum`

#### `*_MeccDNA_C.fasta`

Header format:
```
>MeccDNA1|multi_loci:{sites}_sites|length=...|copies=...|reads=...
```

---

### 2.3 CeccDNA

#### `*_CeccSegments.core.csv`

Adds the following fields on top of Uecc core fields:

| Field | Description |
|-------|-------------|
| `seg_index` / `seg_total` | Segment index and total |
| `junction_role` | `head` / `middle` / `tail` |

#### `*_CeccSegments.bed`

Columns: `chrom, chromStart, chromEnd, name, score, strand, length, eccdna_type, copy_number, match_degree`

#### `*_CeccJunctions.bedpe`

In 1.6.1, N segments produce N directed junctions, including the last-to-first closure. Endpoints use one-base, zero-based half-open intervals at each segment's strand-aware entry and exit. LAST-derived MAPQ, low-MAPQ and confidence fields are missing because LAST does not provide MAPQ; identity is measured from alignment columns. Per-read RCA support remains distinct from molecule-level totals.

Columns: `chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2`

#### `*_CeccDNA_C.fasta`

Header format:
```
>CeccDNA1|segments:{seg_total}|junctions:{chr-chain}|length=...|copies=...|reads=...
```

---

## 3. Inferred intermediate files

### 3.1 `*_simple.csv`

Fallback packaging renames this file to `*_UeccDNA_I.csv`.

Fields:
`eccDNA_id, chr, start0, end0, strand, length, eccdna_type, state, num_split_reads, prob_present, prob_artifact, hifi_abundance, copy_number`

### 3.2 `*_chimeric.csv`

Fields:
`eccDNA_id, chr, start0, end0, strand, length, eccdna_type, state, seg_index, seg_total, junction_role, read_count, num_split_reads, prob_present, prob_artifact, hifi_abundance, copy_number`

### 3.3 FASTA

- `*_UeccDNA_I.fasta`: `>ID|chr:start-end|length=...|type=simple`
- `*_CeccDNA_I.fasta`: `>ID|segments=N|length=...|type=chimeric|seg1:chr:start-end;...`

---

## 4. Reports

- `*_summary.txt`: plain-text summary
- `*_report.html`: HTML report
