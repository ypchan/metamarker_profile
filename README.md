# meta_marker_count

`meta_marker_count` 是一个 paired-end 宏基因组/扩增子 clean reads 的 marker reads 计数流程。主脚本会按样本统计 clean reads 数量，用 BBDuk 从 reads 中提取 16S/18S/ITS 候选 reads，用 minimap2 对候选 reads 比对参考库，最后按 taxonomy 表生成各分类单元的 marker RPM。

默认 marker 是 `16S,ITS`。如果需要真核 18S，运行时显式使用 `--markers 16S,18S,ITS`。

## 功能概览

- `reads_count`: 用 `seqkit stats` 统计 R1 read pairs，并计算 `clean_reads_total = R1 read_pairs * 2`。
- `extract`: 用 `bbduk.sh` 和 marker FASTA 从 R1/R2 中提取候选 reads。
- `align`: 用 `minimap2 -x sr` 将候选 reads 比对到 marker index，输出 PAF。
- `abundance`: 用 taxonomy TSV 解析 PAF，按 `domain/phylum/class/order/family/genus/species` 生成 long table、matrix、domain total 和 assignment stats。

## 安装

### 1. 下载 repo

```bash
git clone <repo-url>
cd meta_marker_count
```

如果已经在本目录内，可以直接从依赖安装开始。

### 2. 安装运行依赖

推荐使用 conda/mamba：

```bash
mamba env create -f environment.yml
mamba activate meta-marker-count
```

也可以使用 conda：

```bash
conda env create -f environment.yml
conda activate meta-marker-count
```

必须能在 `PATH` 中找到这些程序：

- `seqkit`
- `bbduk.sh`
- `minimap2`
- `awk`
- `sort`
- `find`
- `xargs`
- `python3`

检查依赖：

```bash
./meta_marker_count.sh --check-deps --sample-id test --r1 R1.fq.gz --r2 R2.fq.gz
```

上面命令还需要真实 FASTQ 才会通过完整参数检查；更常用的是安装后用实际运行命令加 `--dry-run` 检查。

### 3. 安装命令到 PATH

```bash
make install PREFIX="$HOME/.local"
export PATH="$HOME/.local/bin:$PATH"
```

安装后会得到两个命令：

- `meta_marker_count`: 主流程。
- `meta_marker_build_refs`: 从 SILVA/UNITE 原始 FASTA 生成 short-id FASTA 和 taxonomy TSV。

验证：

```bash
meta_marker_count --help
meta_marker_build_refs --help
```

### 4. 更新

```bash
cd meta_marker_count
git pull --ff-only
make install PREFIX="$HOME/.local"
```

## 输入样本表

多样本输入使用 TSV，必须有表头。列按列名识别，不按列顺序识别。

必需列：

- `sample_id`
- `r1_path`
- `r2_path`

可选列：

- `year`
- `month`
- `depth`
- 任意其他样本 metadata，例如 `site_id`, `site_type`, `group`, `habitat`, `replicate`

示例：

```tsv
sample_id	year	month	depth	site_type	r1_path	r2_path
202311_MF1_00-10	2023	11	00-10	MF	/data/202311_MF1_00-10_R1.fq.gz	/data/202311_MF1_00-10_R2.fq.gz
202311_MF1_10-20	2023	11	10-20	MF	/data/202311_MF1_10-20_R1.fq.gz	/data/202311_MF1_10-20_R2.fq.gz
```

单样本也可以不用 TSV：

```bash
meta_marker_count \
  --sample-id 202311_MF1_00-10 \
  --r1 /data/202311_MF1_00-10_R1.fq.gz \
  --r2 /data/202311_MF1_00-10_R2.fq.gz \
  --ref-16s ref/SILVA_138.2_SSURef_NR99_tax_silva.dna.arc_bac.shortid.fasta \
  --index-16s ref/SILVA_138.2_SSURef_NR99_tax_silva.dna.arc_bac.shortid.fasta.mmi \
  --ref-its ref/UNITE_public_19.02.2025.shortid.fasta \
  --index-its ref/UNITE_public_19.02.2025.shortid.fasta.mmi \
  --taxonomy ref/ref_taxonomy.tsv
```

## 参考库文件生成

主流程需要三类参考文件：

1. BBDuk 使用的 marker FASTA，也就是 short-id FASTA。
2. minimap2 使用的 `.mmi` index。
3. abundance 步骤使用的合并 taxonomy TSV。

本 repo 提供 `scripts/build_ref_files.py` / `meta_marker_build_refs` 从原始 SILVA/UNITE FASTA 生成这些 taxonomy 相关文件。

### 原始文件

准备两个原始文件：

```text
UNITE_public_19.02.2025.fasta.gz
SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
```

### 一次性生成 SILVA + UNITE 参考文件

在 repo 中运行：

```bash
mkdir -p ref

python3 scripts/build_ref_files.py \
  --silva SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
  --unite UNITE_public_19.02.2025.fasta.gz \
  --outdir ref \
  --force
```

如果已经 `make install`，也可以使用：

```bash
meta_marker_build_refs \
  --silva SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
  --unite UNITE_public_19.02.2025.fasta.gz \
  --outdir ref \
  --force
```

生成结果：

```text
ref/UNITE_public_19.02.2025.shortid.fasta
ref/UNITE_public_19.02.2025.taxonomy.tsv

ref/SILVA_138.2_SSURef_NR99_tax_silva.dna.arc_bac.shortid.fasta
ref/SILVA_138.2_SSURef_NR99_tax_silva.dna.arc_bac.taxonomy.tsv

ref/SILVA_138.2_SSURef_NR99_tax_silva.dna.euk.shortid.fasta
ref/SILVA_138.2_SSURef_NR99_tax_silva.dna.euk.taxonomy.tsv

ref/ref_taxonomy.tsv
```

含义：

- `UNITE_public_19.02.2025.fasta.gz -> UNITE_public_19.02.2025.shortid.fasta + UNITE_public_19.02.2025.taxonomy.tsv`
- `SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz -> SILVA_138.2_SSURef_NR99_tax_silva.dna.arc_bac.shortid.fasta + SILVA_138.2_SSURef_NR99_tax_silva.dna.arc_bac.taxonomy.tsv`
- `SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz -> SILVA_138.2_SSURef_NR99_tax_silva.dna.euk.shortid.fasta + SILVA_138.2_SSURef_NR99_tax_silva.dna.euk.taxonomy.tsv`
- `SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz + UNITE_public_19.02.2025.fasta.gz -> ref_taxonomy.tsv`

`ref_taxonomy.tsv` 格式：

```tsv
ref_id	marker	domain	phylum	class	order	family	genus	species
AB000393.1.1510	16S	Bacteria	Pseudomonadota	Gammaproteobacteria	Enterobacterales	Vibrionaceae	Vibrio	Vibrio_halioticoli
```

说明：

- SILVA 中 `Bacteria` 和 `Archaea` 被写入 `.dna.arc_bac.*`，marker 标记为 `16S`。
- SILVA 中 `Eukaryota` 被写入 `.dna.euk.*`，marker 标记为 `18S`。
- UNITE 被写入 `.shortid.*`，marker 标记为 `ITS`。
- FASTA header 会被改成短 ID，并与 taxonomy 表的 `ref_id` 保持一致。
- 如果不同记录产生重复短 ID，后续重复项会自动加数字后缀，避免 `ref_taxonomy.tsv` 中重复 `ref_id` 覆盖。

### 只生成单个数据库

只生成 UNITE：

```bash
meta_marker_build_refs \
  --unite UNITE_public_19.02.2025.fasta.gz \
  --outdir ref \
  --force
```

只生成 SILVA：

```bash
meta_marker_build_refs \
  --silva SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
  --outdir ref \
  --force
```

如果不想生成合并 taxonomy：

```bash
meta_marker_build_refs \
  --silva SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
  --unite UNITE_public_19.02.2025.fasta.gz \
  --outdir ref \
  --skip-combined \
  --force
```

### 生成 minimap2 index

主流程的 `--index-16s/--index-18s/--index-its` 需要 minimap2 index。用生成的 short-id FASTA 建索引：

```bash
minimap2 -d ref/SILVA_138.2_SSURef_NR99_tax_silva.dna.arc_bac.shortid.fasta.mmi \
  ref/SILVA_138.2_SSURef_NR99_tax_silva.dna.arc_bac.shortid.fasta

minimap2 -d ref/SILVA_138.2_SSURef_NR99_tax_silva.dna.euk.shortid.fasta.mmi \
  ref/SILVA_138.2_SSURef_NR99_tax_silva.dna.euk.shortid.fasta

minimap2 -d ref/UNITE_public_19.02.2025.shortid.fasta.mmi \
  ref/UNITE_public_19.02.2025.shortid.fasta
```

如果只跑 `16S,ITS`，可以不生成 18S index。

## 运行主流程

### 16S + ITS 默认流程

```bash
meta_marker_count \
  --input data_path.tsv \
  --outdir marker_count_out \
  --markers 16S,ITS \
  --ref-16s ref/SILVA_138.2_SSURef_NR99_tax_silva.dna.arc_bac.shortid.fasta \
  --index-16s ref/SILVA_138.2_SSURef_NR99_tax_silva.dna.arc_bac.shortid.fasta.mmi \
  --ref-its ref/UNITE_public_19.02.2025.shortid.fasta \
  --index-its ref/UNITE_public_19.02.2025.shortid.fasta.mmi \
  --taxonomy ref/ref_taxonomy.tsv \
  --rank genus \
  --jobs 8 \
  --threads-per-sample 4
```

### 16S + 18S + ITS

```bash
meta_marker_count \
  --input data_path.tsv \
  --outdir marker_count_out \
  --markers 16S,18S,ITS \
  --ref-16s ref/SILVA_138.2_SSURef_NR99_tax_silva.dna.arc_bac.shortid.fasta \
  --index-16s ref/SILVA_138.2_SSURef_NR99_tax_silva.dna.arc_bac.shortid.fasta.mmi \
  --ref-18s ref/SILVA_138.2_SSURef_NR99_tax_silva.dna.euk.shortid.fasta \
  --index-18s ref/SILVA_138.2_SSURef_NR99_tax_silva.dna.euk.shortid.fasta.mmi \
  --ref-its ref/UNITE_public_19.02.2025.shortid.fasta \
  --index-its ref/UNITE_public_19.02.2025.shortid.fasta.mmi \
  --taxonomy ref/ref_taxonomy.tsv \
  --rank genus \
  --jobs 8 \
  --threads-per-sample 4
```

### 试运行

不实际执行分析，只检查参数、依赖、输入路径并打印 run plan：

```bash
meta_marker_count \
  --input data_path.tsv \
  --outdir marker_count_out \
  --markers 16S,ITS \
  --ref-16s ref/SILVA_138.2_SSURef_NR99_tax_silva.dna.arc_bac.shortid.fasta \
  --index-16s ref/SILVA_138.2_SSURef_NR99_tax_silva.dna.arc_bac.shortid.fasta.mmi \
  --ref-its ref/UNITE_public_19.02.2025.shortid.fasta \
  --index-its ref/UNITE_public_19.02.2025.shortid.fasta.mmi \
  --taxonomy ref/ref_taxonomy.tsv \
  --dry-run
```

## 常用参数

- `--steps all`: 默认运行全部步骤。也可以用 `--steps reads_count,extract,align,abundance` 的子集。
- `--markers 16S,ITS`: 默认 marker。可用 `16S`, `18S`, `ITS` 任意组合。
- `--rank genus`: 输出分类层级。可选 `domain,phylum,class,order,family,genus,species,all`。
- `--jobs 8`: 并行样本数。
- `--threads-per-sample 4`: 每个样本内 BBDuk/minimap2 使用线程数。
- `--bbmap-mem 64G`: BBDuk Java heap。
- `--force`: 覆盖并重跑已有输出。
- `--no-resume`: 忽略 checkpoint。
- `--clean-tmp`: 成功后删除临时目录。

## 输出目录

默认输出到 `marker_count_out/`：

```text
marker_count_out/
  00_manifest/samples.tsv
  00_logs/marker_count.YYYYmmdd_HHMMSS.log
  00_logs/run_config.tsv
  01_reads_count/clean_reads.tsv
  02_marker_reads/
    16S/
    18S/
    ITS/
    stats/
    logs/
  03_align/
    16S/
    18S/
    ITS/
    logs/
  04_abundance/
    marker_rpm.<rank>.long.tsv
    marker_rpm.<rank>.matrix.tsv
    marker_rpm.<rank>.domain_total.tsv
    marker_rpm.<rank>.assignment_stats.tsv
  99_checkpoints/
```

关键结果：

- `04_abundance/marker_rpm.genus.long.tsv`: 每个样本、marker、domain、taxon 的 reads 和 RPM。
- `04_abundance/marker_rpm.genus.matrix.tsv`: 样本 x feature 矩阵。
- `04_abundance/marker_rpm.genus.domain_total.tsv`: 每个样本每个 domain/marker 的总 marker reads/RPM。
- `04_abundance/marker_rpm.genus.assignment_stats.tsv`: PAF hit 过滤和 ambiguous 分配统计。

## 质量过滤默认值

主脚本默认值适合 genus-level 的生态比较：

- 16S identity >= `0.97`
- 18S identity >= `0.97`
- ITS identity >= `0.95`
- alignment length >= `80`
- query coverage >= `0.60`
- mapQ >= `0`
- near-top identity diff <= `0.010`
- near-top alignment length diff <= `10`

如果一个 read 的 near-best hits 对应多个不同 taxon，会标记为 ambiguous，不进入 taxon abundance。

查看完整默认参数：

```bash
meta_marker_count --print-defaults
```

## 基础检查

开发或更新后运行：

```bash
make check
```

它会执行：

- `bash -n meta_marker_count.sh`
- `python3 -m py_compile scripts/build_ref_files.py`
