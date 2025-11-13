# Spatial single-molecule DNA methylation and transcriptome co‑profiling

## 1. Overview
This repository contains reproducible preprocessing, downstream analysis, and visualization scripts used in the study **"Spatial single-molecule DNA methylation and transcriptome co‑profiling in mammalian tissues"**. The workflow covers raw sequencing to pixel-/cell-resolved DNA methylation (5mC), hydroxymethylation (5hmC) and transcriptome co‑analysis.

**Design principles:**  
- Global workflow: FAST5/FASTQ → BAM/matrices → visualizations.  
- Explicit parameters: software versions, thresholds, filtering rules.  
- Modular and reusable components

---

## 2. Quick start
1. Configure dependencies and confirm key tool versions (see Section 3).  
2. Run RNA and DNA preprocessing scripts in order.  
3. Run downstream integration and visualization scripts.

---

## 3. Dependencies & recommended versions
- Dorado (Nanopore basecalling)
- samtools
- Skewer
- UMI‑tools (recommended v1.0.1)
- seqkit
- bamCoverage / mosdepth
- Sinto
- bedtools
- modkit
- minimap2
- Python3 and R with: Seurat, Signac, rGREAT, and other commonly used bioinformatics packages

> We recommend using Conda environments or containerized deployments for reproducibility.

---

## 4. RNA (spatial transcriptomics) preprocessing  
**Goal:** Extract spatial barcodes and UMIs from raw reads, produce pixel × gene expression matrices, and retain intronic reads (nascent transcripts) to increase sensitivity.

```bash
# RNA preprocessing (generate expression matrices)
bash STpipe_clipping.sh
```

**Key steps:**  
1. Use `UMI‑tools v1.0.1` to extract Barcode B, Barcode A and UMI from Read 2.    
2. Reformatted FASTQ are processed through the ST pipeline for adapter trimming, alignment (human → GRCh38; mouse → GRCm39), barcode demultiplexing, read annotation and UMI counting.  
3. To retain nuclear nascent transcripts, CDS annotations are extended to full-transcript regions and exonic + intronic counts are combined to form final pixel × gene matrices.  
4. Tissue boundaries and cavities are defined by manual alignment of adjacent H&E images to the spatial pixel map, and by per-pixel coverage (reads/UMIs) to exclude pixels outside tissue or inside obvious lumina. Only pixels within the defined tissue regions are carried forward.

**Key scripts:** `1_extract_barcode.sh`, `3_st_pipeline.soft_clipping.sh`.

---

## 5. RNA downstream analysis (Seurat-centric)
**Workflow highlights:**
- Normalization: `SCTransform()` (retain raw counts for differential testing).  
- Sample merging: `merge()` followed by standard Seurat steps (scaling, PCA, neighbor graph construction, UMAP, Louvain clustering).  
- Second‑round subclustering: applied to complex compartments (e.g., spermatogonial stem cell regions) to resolve rare or finely partitioned populations.  
- Sparsity mitigation: use MAGIC to impute and denoise expression profiles for visualization and lineage inference. (Do not use MAGIC‑imputed data for statistical testing.)  

**Practical notes:**
- Log clustering resolution parameters and random seeds for reproducibility.  
- When presenting gene expression figures, always state which data type is used (raw counts / SCTransform / MAGIC‑inferred).  

**Key scripts:** `merge_data.R`, `filter_pixel.R`.

---

## 6. Third-generation (Nanopore) spatial DNA prepocessing pipeline
**Goal:** Call modified bases (5mC, 5hmC) from Nanopore reads, recover spatial barcodes, and produce pixel × CpG or pixel × interval modification matrices.

```bash
# DNA preprocessing (spatial long-read methylation)
bash pipeline.sh config
```

Key steps and notes:  
1. **Basecalling & modified‑base calling:** Dorado with the `hac,5mC_5hmC` model to predict 5mC and 5hmC probabilities; outputs processed BAM files.  
2. **Alignment:** `minimap2` to align reads to GRCh38 / GRCm39.  
3. **Spatial barcode recovery:** convert BAM to FASTQ, trim adapters with Skewer, and extract Barcode B, Barcode A and I7 index using UMI‑tools. The I7 index is used to demultiplex pooled libraries.  
4. **Error correction:** compare each read barcode to a whitelist; accept a whitelist barcode as the corrected barcode if its edit distance to the read is less than half the minimum pairwise edit distance among whitelist barcodes — this substantially improves barcode recovery and read utilization.  
5. **BAM tagging:** write sample identity and spatial pixel coordinates into BAM as `SP` and `CB` tags respectively using Sinto, enabling downstream per‑pixel aggregation.  

**Unified modification-calling thresholds (after inspecting per-sample modkit probability distributions):** C = 0.7, m = 0.7, h = 0.6.

---

## 7. Pixel- and molecule-level metrics & multi-scale aggregation
Due to extreme sparsity of genome-wide coverage, we aggregate CpG‑level calls across multiple genomic scales:  
- 2 kb genomic windows  
- Gene promoters and gene bodies  
- Repeat‑family intervals  
- CpG islands and cell‑type‑specific enhancers  

**Implementation:** `modkit pileup` split by `CB` tag produces pixel-level single-CpG modification proportions. Use `intersectBed` plus custom scripts to aggregate CpG calls into the intervals above, producing pixel × interval matrices.

**Parallel path:** `modkit extract calls` recovers single‑molecule calls; per‑read methylation and hydroxymethylation levels are computed for analyses requiring molecule resolution.

**Key scripts:** `merge_pixel_or_region.GW6.R` , `modkit_work.sh&modkit_density_plot.R`.

---

## 8. Calling DMRs and high‑5hmC regions
**Reproducible calling rules:**  
- **DMR (2 kb windows):** absolute difference in methylation between cell type of interest and background ≥ 0.2; at least 5 observed CpG calls in the window; Fisher’s exact test BH‑adjusted p‑value < 0.05.  
- **High‑5hmC region (2 kb windows):** target cell type 5hmC level > 90th percentile of global distribution and > 10 observed CpG calls.  

**Downstream:** Aggregate cell‑type DMRs back to pixels to produce a pixel × DMR‑group matrix; reduce dimensionality with Signac to create DMR‑based UMAPs. Apply the same pipeline to high‑5hmC regions.

**Key scripts:** `DMR.GW6.R`.

---

## 9. Functional annotation & enrichment
- Map DMRs / high‑5hmC regions to nearby genes using rGREAT (basal‑plus‑extension rule).  
- Use GREAT’s binomial test for enrichment and apply Benjamini‑Hochberg correction.  
- Retain GO terms with binomial fold enrichment > 2 and BH‑adjusted p < 0.05.  

**Key scripts:** `DMR.GREAT.R`,`5hmc_high_bins.R`.

---

## 10. Integrative RNA–DNA analysis strategy
**Aggregation principles:** summarize gene-body methylation and 5hmC as mean modification level per gene within each cell type. Obtain gene expression by `Seurat::AggregateExpression()` and quantile‑normalize aggregated expression for comparability.  

**Analysis flow:**  
1. Keep genes with coverage in >70% of cell types.  
2. Compute pairwise Pearson correlations for each retained gene across cell types:
   - gene‑body methylation vs expression  
   - gene‑body hydroxymethylation vs expression  
   - methylation vs hydroxymethylation  
3. Report genes with absolute correlation > 0.4 for downstream interpretation.  

**Key scripts:** `joint_RNA_hmC_analysis.signac.zscore0_as_NA.dmr_deg_feature.R`,`joint_RNA_mC_analysis.signac.zscore0_as_NA.dmr_deg_feature.R`.

---

## 11. Example analysis: TF enrichment in high‑5hmC regions
- Identify enriched sequence motifs within high‑5hmC intervals using HOMER (`findMotifsGenome.pl`) and record enrichment p‑values.  
- Assign motifs to genes via automated name extraction followed by manual literature curation to resolve ambiguous matches.  
- Quantify TF expression by aggregating counts to cell types and converting to TPM. Combine motif enrichment p‑values with TF TPM for integrated visualizations that prioritize TFs likely regulated via 5hmC.

**Key scripts:** `correlation.body_mC_hmC_level.cell_type.GW6.use_3_cor.R`. 

---

## 12. Contact & contributing
To contribute, report issues, or request help, please open an issue in the repository or contact the corresponding author listed in the paper.

---
