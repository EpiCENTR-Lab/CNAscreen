# CNAscreen — Copy Number Alteration QC with QDNAseq

## Overview

This pipeline performs a genome integrity screen using low-pass whole-genome sequencing (WGS) BAM files.

It is designed for cultured cell quality control, including iPSCs, ESCs, primary epithelial cells, fibroblasts, and engineered cell lines.

The script produces genome-wide copy number plots and summary tables that help detect:

-   Whole chromosome gains or losses
-   Chromosome arm-level changes
-   Large structural abnormalities

This tool is not intended for full CNV discovery.

## Recommended Sequencing Parameters

-   Coverage: 0.1×–5×
-   Read length: 75–150 bp
-   Paired-end sequencing preferred
-   Reference genome: hg38 (required)

## Installation

CNAscreen needs **R 4.5** (its packages come from Bioconductor 3.22, which requires R 4.5; they will not install on R 4.6 or later). We use **rig** to install R 4.5.3 alongside any other R version you have.

1.  **Install RStudio Desktop** from <https://posit.co/download/rstudio-desktop/>.

2.  **Install rig** by following the instructions for your system at <https://github.com/r-lib/rig#installation>.

3.  **Install R 4.5.3 in rig and make it the default.** Open a terminal (Mac: *Terminal* app; Windows: *Command Prompt*) and run:

    ```         
    rig add 4.5.3
    rig default 4.5.3
    ```

    *Windows only:* also run `rig rtools add`. *Mac only:* also run `sudo rig system make-orthogonal`

4.  **Download CNAscreen** and unzip it.

5.  **Set the working directory and open the project.** In *Terminal* (Mac) or *Command Prompt* (Windows).

    ```         
    cd /path/to/CNAscreen
    rig rstudio CNAscreen.Rproj
    ```

    On Windows, if the folder is on another drive, use `cd /d D:\path\to\CNAscreen`. RStudio opens. Check that the window title says `CNAscreen` and the Console starts with `R version 4.5.3`.

6.  **Install the packages** (one-off, needs internet, \~10–30 min). In the **RStudio Console**, type:

    ``` r
    renv::restore()
    ```

    If prompted about project activation, select `1`. Answer `y` when asked if you would like to proceed with updating packages.

N.B. To switch back to another R version later, run `rig list` to see what's installed and `rig default <version>` to choose one.

## Running CNAscreen

1.  Put your BAM files (and their paired `.bai` index files) in the `bam` folder, and fill in `sample_sex.csv`.

2.  In RStudio, click the **Terminal** tab (next to Console) and run:

    ```         
    Rscript run_CNAscreen.R --bam-dir bam --out-dir CNAscreen_output --bin-size 100 --genome hg38 --sex-map sample_sex.csv --map-cutoff 80 --bases-cutoff 99.9 --seed 170826
    ```

3.  Results appear in the `CNAscreen_output` folder.

## Input Data

The script expects BAM files and BAM index files in the `bam` folder inside CNAscreen. BAM files must be indexed (`.bam.bai` files), and the BAM filename becomes the sample name.

Example:

```         
bam/
    H9_p116.bam
    H9_p116.bam.bai
```

### Sex Annotation File

Modify the file called `sample_sex.csv` to replace the placeholder example text with sample names and sexes.

Example:

| sample      | sex |
|-------------|-----|
| iPSC_A      | XX  |
| iPSC_B      | XY  |
| fibroblast1 | XX  |

Sex must be `XX` or `XY`.

Sample names must match BAM filenames without `.bam`.

## Pipeline Steps

1.  **Read binning** – reads are counted in genomic bins using QDNAseq
2.  **Filtering** – problematic regions are removed (low mappability, blacklist regions)
3.  **Bias correction** – GC and mappability bias correction
4.  **Normalization** – counts converted to relative copy number
5.  **Segmentation** – contiguous regions with similar copy number are identified
6.  **Visualization** – genome-wide plots and a summary PDF are generated

## Options

| Option | Default | Description |
|------------------------|------------------------|------------------------|
| `--bam-dir` | `bam` | Folder containing `.bam` files and their `.bam.bai` indexes |
| `--out-dir` | `CNAscreen_output` | Output folder (created if it doesn't exist) |
| `--sex-map` | *required* | CSV with columns `sample,sex` (`XX` or `XY`) covering every BAM |
| `--genome` | `hg38` | Genome build of the QDNAseq bin annotation (`hg38`) |
| `--bin-size` | `100` | Bin size in kb |
| `--min-mapq` | `37` | Minimum read mapping quality counted into bins |
| `--map-cutoff` | `80` | Bins with mappability below this percentage are filtered |
| `--bases-cutoff` | `99.9` | Bins with characterised bases below this percentage are filtered |
| `--exclude-regions` | none | BED file of additional regions to mask (see below) |
| `--no-sex-ploidy-correction` | off | Switch: disable the chrX/chrY correction for XY samples (see below) |
| `--segment-alpha` | `1e-10` | Significance level for segment change-points (see below) |
| `--segment-undo-sd` | `1` | Merge adjacent segments differing by less than this many SDs (see below) |
| `--segment-tol` | `0.01` | Merge adjacent segments with near-identical log2 levels when drawing (see below) |
| `--seed` | `170826` | Random seed for segmentation (see Reproducibility) |
| `--threads` | `1` | Number of samples processed in parallel (Mac/Linux only; see below) |
| `--plot-filtered-bins` | `TRUE` | Draw filtered bins as light grey points (`TRUE`/`FALSE`) |
| `--ymin`, `--ymax` | `-2.2`, `2.2` | y-axis limits of the plots (log2 ratio) |
| `--width-mm`, `--height-mm` | `95`, `60` | Size of each per-sample PNG (600 dpi). Does not affect the summary PDF |
| `--verbose` | `TRUE` | Print progress messages (`TRUE`/`FALSE`) |

Switches take no value. `--plot-filtered-bins` and `--verbose` need an explicit `TRUE` or `FALSE`.

### Bin size

```         
--bin-size 100
```

Recommended:

-   50 kb – higher resolution but noisier
-   100 kb – recommended default
-   500 kb – smoother but lower resolution

### QC thresholds

```         
--map-cutoff 80
--bases-cutoff 99.9
```

### Sex-chromosome ploidy correction

```         
--no-sex-ploidy-correction
```

XY samples carry one copy of chrX and one copy of chrY, versus two copies of each autosome. Left uncorrected, this \~50% dosage difference reads on the genome-wide plot as an apparent whole-chromosome loss on X and Y, even in a chromosomally normal sample.

By default, the pipeline rescales chrX/chrY copy number by 2x for XY samples before computing the log2 ratio, so a normal XY sample sits at \~0 on X and Y, just like the autosomes. A genuine loss or gain on X/Y still shows up as a deviation from 0 after this correction. Correction is applied automatically to any sample classified as XY (via `--sex-map`). XX samples are unaffected.

To see the raw, uncorrected QDNAseq output instead, pass `--no-sex-ploidy-correction`.

Whether the correction was applied is recorded per sample in the `sex_ploidy_corrected` column of `sample_qc_summary.tsv`.

### Segmentation sensitivity

```         
--segment-alpha 1e-10
--segment-undo-sd 1
--segment-tol 0.01
```

Segmentation uses DNAcopy circular binary segmentation via QDNAseq's `segmentBins()` (`undo.splits = "sdundo"`). Two parameters control how readily the genome is split into segments:

-   `--segment-alpha` — significance level for accepting a change-point. Smaller values give fewer, more conservative segments.
-   `--segment-undo-sd` — adjacent segments whose means differ by less than this many standard deviations of the data are merged back together. Because it is measured in SD units, it scales with each sample's noise.

The pipeline uses QDNAseq's defaults (`1e-10` and `1`).

**Trade-off.** Stricter settings are more robust across sequencing depths but less sensitive to small or low-level subclonal events — for example, the recurrent \~1 Mb 20q11.21 gain seen in hPSC cultures spans only a handful of 100 kb bins. The defaults suit the broad genome-integrity QC this pipeline is designed for. If you are looking for focal or mosaic events, loosen `--segment-alpha` in steps, inspect the segment table, and confirm any finding orthogonally (see Limitations).

**Checking for over-segmentation.** `sample_qc_summary.tsv` reports `total_reads`, `noise_mapd` and `n_segments` per sample. A sample with far more segments than others in the run, but no large events visible in its plot, is a candidate for over-segmentation; compare it against its `noise_mapd` and `total_reads`.

`--segment-tol` (default `0.01`) is separate: it only merges adjacent, already-called segments whose log2 levels are essentially identical so they draw as a single line.

### Excluding additional problem regions

```         
--exclude-regions problem_regions.bed
```

Beyond QDNAseq's built-in blacklist, mappability, and bases filters, you can mask additional genomic regions with a BED file (`chrom`, `start`, `end`, tab-delimited, 0-based half-open per standard BED convention; extra columns are ignored; a `chr` prefix on chromosome names is stripped automatically). This is for regions you've identified empirically — for example, a bin position that recurs as an outlier across many of your own samples regardless of biology, or a known problem region such as a segmental duplication exported from the UCSC Table Browser or the ENCODE blacklist — that you want masked consistently across every sample in a run without changing the global `--map-cutoff`/`--bases-cutoff` thresholds (which apply genome-wide).

Bins overlapping any listed region are flagged `Excluded region` (or `Excluded + QC fail` if they also fail mappability/bases) in the `status` column, and excluded from the plot and QC summary the same way as any other failing bin. The per-sample `frac_excluded_region` column in `sample_qc_summary.tsv` reports how much of that sample's genome was masked this way.

### Parallelising across samples

```         
--threads 4
```

By default, samples are processed one at a time. With `--threads` \> 1, that many samples are processed concurrently (Mac/Linux only; on Windows the script falls back to serial, so use WSL for parallel runs).

Each worker holds one sample's read data in memory, so N workers need roughly N times the RAM of a serial run. Choose the value based on BAM size and available memory, not sample count. The value is automatically capped at the number of samples and at one fewer than the number of CPU cores.

Output is identical to a serial run. Progress messages from different workers may interleave.

### Reproducibility

DNAcopy segmentation uses permutation testing, so breakpoints can vary between runs unless a seed is fixed. Each sample is segmented with its own seed, derived from `--seed` (default `170826`) and the sample name. A sample's result therefore depends only on its name and `--seed`, not on which other samples are in the run, their order, or the `--threads` setting.

`--seed` is recorded in `SESSION_INFO.txt`. To reproduce a previous run, use the same `--seed` and sample names.

## Output Files

Example output directory:

```         
CNAscreen_output/
    sample1.plot.png
    sample1.qdnaseq.rds
    sample1.bins.tsv.gz
    sample1.segments.tsv
    sample2.plot.png
    sample2.qdnaseq.rds
    sample2.bins.tsv.gz
    sample2.segments.tsv
    sample_qc_summary.tsv
    CNAscreen_summary.pdf
    SESSION_INFO.txt
```

## Output Descriptions

### Genome plot — `sample.plot.png`

-   Black points = high-quality bins
-   Purple lines = segmented copy number
-   Dashed line = diploid baseline

Interpretation:

-   Flat near 0 → diploid genome
-   Whole chromosome above 0 → gain
-   Whole chromosome below 0 → loss

### Bin-level data — `sample.bins.tsv.gz`

Contains one row per genomic bin including genomic coordinates, QC status, raw copy number, and segmented copy number. The `status` column is `OK` for bins used in the plot/segmentation, or one of `Bases < cutoff`, `Low mappability`, `Both low`, `Excluded region`, or `Excluded + QC fail` for bins filtered out (the last two only apply when `--exclude-regions` is used).

### Segment table — `sample.segments.tsv`

Contains segmented copy number regions. One row per segment, with hg38 coordinates and copy ratio; see Locating and Reporting Abnormalities.

### QDNAseq object — `sample.qdnaseq.rds`

Stores the full QDNAseq analysis object for future analysis.

### QC summary — `sample_qc_summary.tsv`

One row per sample. Columns:

| Column | Description |
|------------------------------------|------------------------------------|
| `sample` | Sample name (BAM filename without `.bam`) |
| `total_reads` | Reads with MAPQ ≥ `--min-mapq` counted into genomic bins, before bin filtering. A direct measure of sequencing depth |
| `noise_mapd` | Bin-to-bin noise: median absolute difference in log2 ratio between consecutive passing bins on the same chromosome (MAPD-style, unscaled). Lower is cleaner. Comparable between samples run at the same `--bin-size` |
| `n_segments` | Number of segments across the genome, after merging adjacent segments within `--segment-tol` |
| `n_bins` | Total bins on the chromosomes analysed for this sample (chrY is included only for XY samples) |
| `n_finite` | Bins with a usable copy-number value |
| `frac_negative` | Fraction of usable bins with log2 ratio \< 0 |
| `median_log2` | Median log2 ratio across usable bins; expected to be close to 0 |
| `median_bases_pct` | Median percentage of characterised (non-N) reference bases per bin |
| `median_mappability_pct` | Median mappability per bin (%) |
| `frac_low_bases` | Fraction of bins below `--bases-cutoff` |
| `frac_low_map` | Fraction of bins below `--map-cutoff` |
| `frac_excluded_region` | Fraction of bins overlapping `--exclude-regions` (0 if not used) |
| `sex_mode` | `XX` or `XY`, as assigned in `--sex-map` |
| `sex_ploidy_corrected` | Whether the chrX/chrY ploidy correction was applied |
| `plot_file`, `qdnaseq_rds`, `bins_tsv_gz`, `segments_tsv` | Names of this sample's output files |

`median_bases_pct`, `median_mappability_pct`, `frac_low_bases` and `frac_low_map` come from the reference bin annotation, not from the sequencing data. They are therefore the same for all samples of the same sex at a given `--bin-size` and genome build. The sample-specific quality measures are `total_reads`, `noise_mapd` and `n_segments`.

### Session information — `SESSION_INFO.txt`

Records software versions and command used to run the pipeline.

## Locating and Reporting Abnormalities

The plot shows whether a sample is abnormal; `sample.segments.tsv` shows where. Each row is one segment: `chromosome`, `start`, `end` (1-based hg38 coordinates), `size_mb`, `n_bins` (bins supporting the segment), `log2_ratio` and `copy_ratio`.

To find an abnormality, look for segments that differ clearly from 0, especially long ones with many bins. Then:

1.  **Compare with expected values.** A one-copy change present in every cell of a diploid sample gives:

    | Change           | copy_ratio | log2_ratio |
    |------------------|------------|------------|
    | Loss of one copy | 0.5        | −1.0       |
    | Gain of one copy | 1.5        | +0.58      |

    Values between 0 and these suggest the change is in only some cells. For a single-copy change, the fraction of cells affected is roughly `2 × |copy_ratio − 1|`; for example, 1.25 suggests \~50% of cells. Treat this as approximate at low coverage.

2.  **Find the cytogenetic location.** Paste `chromosome:start-end` into the [UCSC Genome Browser](https://genome.ucsc.edu/) (hg38) to see the cytobands it spans, and report it in standard form, e.g. "+12" for a whole-chromosome gain or "gain 20q11.21".

3.  **Check for neighbouring segments.** A whole chromosome or arm may be split into two or more segments either side of the centromere, or where noise breaks it up. Adjacent segments at similar levels usually describe one event.

4.  **Confirm what matters.** Treat events under 1 Mb, and any finding that will affect whether a line is used, as needing orthogonal confirmation (karyotype, SNP array, or ddPCR/qPCR for known hotspots such as 20q11.21). See Limitations.

In the demo sample (H9 p116), H9_p116.segments.tsv shows a gain from 12p to distal 12q (chr12:100,001–110,700,000, copy ratio 1.50, one extra copy in all cells) and a gain at the end of 17q (chr17:77,100,001–83,200,000, copy ratio 1.36).

## Limitations

This pipeline cannot reliably detect small CNVs.

Detection reliability:

| Size    | Reliability                  |
|---------|------------------------------|
| \>10 Mb | Very reliable                |
| 1–10 Mb | Usually detectable           |
| \<1 Mb  | Unreliable with low-pass WGS |

Given 100 kb bins and observed genome-wide noise, fully clonal copy-number changes are statistically detectable below 1 Mb; however, we adopt 1 Mb as a conservative calling threshold to guard against locus-specific mapping artifacts and treat sub-Mb findings as requiring orthogonal investigation.

## Troubleshooting

**Missing BAM index**

Create indexes with:

``` bash
samtools index sample.bam
```

OR in R Studio console:

``` r
Rsamtools::indexBam("bam/sample1.bam")
```

**X/Y chromosome issues**

A sample showing a clean whole-chromosome-scale offset (\~±1.0 log2) confined to chrX/chrY, with autosomes flat, most likely indicates a sex-map/sample mismatch rather than a true copy-number event. Verify `sample_sex.csv` against the sample's known genotype.

## Demo Data

1.  **Download** `H9_p116.bam` and `H9_p116.bam.bai` from Zenodo
    ([doi:10.5281/zenodo.24480952](https://doi.org/10.5281/zenodo.24480952))
    into `demo_data/bam/`.
    
2.  **Check the download** (optional but recommended). In the RStudio Console:

    ``` r
    tools::md5sum(c("demo_data/bam/H9_p116.bam", "demo_data/bam/H9_p116.bam.bai"))
    ```

    The values should match those in `demo_data/md5sums.txt`. If they don't, download the files again.

3.  **Run the demo** in the Terminal tab:

    ```         
    Rscript run_CNAscreen.R --bam-dir demo_data/bam --sex-map demo_data/sample_sex.csv --out-dir demo_output --seed 170826
    ```

4.  **Compare with the expected result** in the Console:

    ``` r
    all.equal(read.delim("demo_output/H9_p116.segments.tsv"),
              read.delim("demo_data/expected_output/H9_p116.segments.tsv"))
    ```

    `TRUE` means your installation reproduces the published result.

**Source.** H9 human embryonic stem cells, passage 116 (karyotypically abnormal subline), SRA run [SRR14354893](https://www.ncbi.nlm.nih.gov/sra/SRR14354893) (BioProject PRJNA726033, BioSample SAMN18918172), from <doi:10.1186/s12864-021-07686-z>. The reads are open access in SRA.

**Processing.** Reads were aligned to GRCh38 (no-alt analysis set; GCA_000001405.15) with BWA-MEM (Galaxy Version 0.7.19+galaxy2) on Galaxy. Coverage before read removal is \~0.11×. Unmapped, secondary, supplementary and MAPQ \< 37 reads were removed.

**Expected result.** A gain of chr12 from 12p to distal 12q (chr12:0.1–110.7 Mb, copy ratio \~1.5, consistent with one extra copy in all cells) and a smaller gain at the end of 17q (chr17:77.1–83.2 Mb, copy ratio \~1.36, consistent with a gain in most cells). Both are recurrent in hPSCs.

## License and citation

CNAscreen is released under the MIT License (see `LICENSE`). The R packages it uses are installed separately by renv and are covered by their own licences.

If you use CNAscreen, please cite: \> Hynds, R. E. (2026). CNAscreen (version \<x.y.z\>) [Software]. Zenodo. doi:<Zenodo DOI>
