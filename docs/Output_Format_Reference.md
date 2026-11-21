# CircleSeeker 输出格式参考手册

第 16 步直接把结果写入 `-o` 指定的目录。`<prefix>` 是用 `-p` 指定的样本前缀。

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

打包时会创建三个类型目录。有相应区域时才写 BED/BEDPE，有序列时才写各类型的 FASTA。HTML 和文本报告由前一步生成，打包时复制到输出目录。

## 最终 CSV

`<prefix>_eccDNA_summary.csv` 中每条 eccDNA 占一行，表头如下。

```csv
eccDNA_id,type,state,chr,start,end,strand,length,location,segment_count,read_count,copy_number,confidence,inferred_reads
```

| 字段 | 含义 |
| --- | --- |
| `eccDNA_id`、`type`、`state` | 每种类型从 1 独立编号，至少四位，如 `UeccDNA0001`。`type` 取 `Uecc`、`Mecc` 或 `Cecc`，`state` 取 `Confirmed` 或 `Inferred`。同一次运行的结果表可用 `eccDNA_id` 关联。 |
| `chr`、`start`、`end`、`strand` | 坐标采用零基半开区间。Uecc 对应单一来源，Mecc 显示代表性候选位点。Cecc 的 `chr` 为 `multi`，其余三列为 `.`。 |
| `length`、`location`、`segment_count` | `length` 单位为 bp。Uecc 对应一个位点。Mecc 是一条环序列，`\|` 分隔它的多个候选来源位点，`segment_count` 为 1。Cecc 由多个片段组成，`;` 按环上顺序分隔组成片段，`segment_count` 计片段数。 |
| `read_count`、`inferred_reads` | `read_count` 计支持读段总数。合并重叠结果时，Confirmed 记录可以获得额外的 split-read 支持，`inferred_reads` 记录追加的数量。独立 Inferred 记录的 split-read 支持数看 `read_count`。 |
| `copy_number` | Confirmed 记录使用重复拷贝估计，合并时可计入重叠 Inferred 记录的估计值。独立 Inferred 记录使用取整后的覆盖深度作为代理值，缺失或非正数时取 1。 |
| `confidence` | Confirmed 有上游评分时沿用该值。内置 SplitReads-Core 适配器使用固定兼容值 `prob_present=0.95` 和 `prob_artifact=0.0`，Inferred 的 `confidence` 继承 `prob_present`。0.95 不能用来判断单条候选的存在概率或给候选排序。 |

`<prefix>_eccDNA_regions.csv` 每行对应一个位点或片段。

```csv
eccDNA_id,region_idx,chr,start,end,strand,length,role
```

`region_idx` 从 1 开始。Uecc 的 `role` 为 `source`。Mecc 的 `primary` 是代表性位点，`candidate` 是其他候选位点，真实来源仍有多种可能。Cecc 使用 `head`、`middle`、`tail` 表示各片段在环上的顺序。

`<prefix>_eccDNA_reads.csv` 每行对应 Confirmed 分支中一个去重后的 `(eccDNA_id, read_name)`。表头如下。

```csv
eccDNA_id,read_name,copy_number,eccDNA_copy_number,match_degree,identity,mapq
```

两个拷贝数字段分别对应读段和 Confirmed eccDNA，逐读段估计缺失时可以为空。整合时追加的 split-read 支持计入主表，逐读段明细行数因此可能小于主表的 `read_count`。没有具名 Confirmed 读段时，明细 CSV 可能为空且没有表头。

## 最终 BED、BEDPE 和 FASTA

BED 含六列，用制表符分隔，无表头，顺序为 `chrom, start, end, name, score, strand`。`score` 固定为 1。Uecc 的 `name` 是 ID，Mecc 使用 `ID|role`，Cecc 使用 `ID|segI/N|role`。Cecc 的 BEDPE 列为 `chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2`。每行用单碱基区间连接上一片段的出口与下一片段的入口，保留链方向，并包含末段到首段的闭环连接。

总 FASTA 与各类型 FASTA 的表头均为 `>ID|location|length=N|state=STATE`。Mecc 的 `location` 本身可以含多个 `|`。解析时从第一个 `|` 左边取 ID，从右边取最后两个字段，也可以直接用 ID 查询区域表。

```text
>MeccDNA0001|chr1:100-200(+)|chr2:300-400(-)|length=100|state=Confirmed
```

## 中间文件与回退输出

下面列出各证据分支的中间文件。使用 `--keep-tmp` 可保留中间结果。缺少统一 CSV 时，第 16 步可以把已有文件复制到 `<prefix>_Confirmed_UeccDNA/`、`<prefix>_Confirmed_MeccDNA/`、`<prefix>_Confirmed_CeccDNA/` 和 `<prefix>_Inferred_eccDNA/`。此时 `<prefix>_merged_output.csv` 来自已有的 **Confirmed** 汇总表。

中间表的 `reads_count` 计支持读段数。Confirmed Mecc 的 `Hit_count` 计候选位点数，Confirmed Uecc/Cecc 的 `Hit_count` 默认为 1。Inferred Uecc/Cecc 的两个字段均取 split-read 支持数。

## 2. Confirmed 中间文件

### 2.1 UeccDNA

#### `*_UeccDNA.core.csv`

| 字段 | 说明 |
|------|------|
| `eccDNA_id` | 唯一 ID |
| `chr` / `start0` / `end0` / `strand` | 0-based 坐标与方向 |
| `length` | 序列长度（bp） |
| `match_degree` | 匹配度（0-100） |
| `confidence_score` / `query_cov_best` / `query_cov_2nd` | 证据打分与覆盖度（可能为空） |
| `mapq_best` / `identity_best` | 最佳比对质量（可能为空） |
| `low_mapq` / `low_identity` | 低质量标记（可能为空） |
| `copy_number` / `repeat_number` | 拷贝数/重复数估计 |
| `eccdna_type` | `Uecc` |
| `num_merged` / `merged_from_ids` | 合并信息 |
| `reads_count` / `read_name` | 支持 reads 数与列表 |

#### `*_UeccDNA.bed`

列顺序：`chrom, chromStart, chromEnd, name, score, strand, length, eccdna_type, repeat_number, match_degree`

#### `*_UeccDNA_C.fasta`

FASTA header 格式：
```
>UeccDNA1|chr:start-end(strand)|length=...|repeats=...|reads=...
```

---

### 2.2 MeccDNA

#### `*_MeccSites.core.csv`

在 Uecc 核心字段基础上追加：

| 字段 | 说明 |
|------|------|
| `hit_index` | 当前位点序号（1-based） |
| `hit_count` | 同一 eccDNA 命中总数 |

#### `*_MeccSites.bed`

列顺序：`chrom, chromStart, chromEnd, name, score, strand, length, eccdna_type, copy_number`

#### `*_MeccBestSite.bed`

列顺序（保留 legacy 字段名）：`chrom, chromStart, chromEnd, name, score, strand, eLength, eClass, copyNum`

#### `*_MeccDNA_C.fasta`

FASTA header 格式：
```
>MeccDNA1|multi_loci:{sites}_sites|length=...|copies=...|reads=...
```

---

### 2.3 CeccDNA

#### `*_CeccSegments.core.csv`

在 Uecc 核心字段基础上追加：

| 字段 | 说明 |
|------|------|
| `seg_index` / `seg_total` | 片段序号与总数 |
| `junction_role` | `head` / `middle` / `tail` |

#### `*_CeccSegments.bed`

列顺序：`chrom, chromStart, chromEnd, name, score, strand, length, eccdna_type, copy_number, match_degree`

#### `*_CeccJunctions.bedpe`

1.6.1 中，N 个片段输出 N 条有向连接，包含末段到首段的闭环连接；端点使用按链方向确定的入口/出口，表示为单碱基、零基半开区间。LAST 不提供 MAPQ，因此对应 MAPQ、低 MAPQ 标记和置信度字段保留缺失；identity 根据实际比对列计算。逐读段 RCA 支持量与分子总量分别记录。

列顺序：`chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2`

#### `*_CeccDNA_C.fasta`

FASTA header 格式：
```
>CeccDNA1|segments:{seg_total}|junctions:{chr-chain}|length=...|copies=...|reads=...
```

---

## 3. Inferred 中间文件

### 3.1 `*_simple.csv`

备用打包路径将此文件重命名为 `*_UeccDNA_I.csv`。

字段：
`eccDNA_id, chr, start0, end0, strand, length, eccdna_type, state, num_split_reads, prob_present, prob_artifact, hifi_abundance, copy_number`

### 3.2 `*_chimeric.csv`

字段：
`eccDNA_id, chr, start0, end0, strand, length, eccdna_type, state, seg_index, seg_total, junction_role, read_count, num_split_reads, prob_present, prob_artifact, hifi_abundance, copy_number`

### 3.3 FASTA

- `*_UeccDNA_I.fasta`：`>ID|chr:start-end|length=...|type=simple`
- `*_CeccDNA_I.fasta`：`>ID|segments=N|length=...|type=chimeric|seg1:chr:start-end;...`

---

## 4. 统计与报告

- `*_summary.txt`：文本统计摘要
- `*_report.html`：HTML 报告
