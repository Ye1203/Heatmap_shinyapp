# Heatmap Helper package

Heatmap Helper package
This code uses shinyapp to generate heatmaps. It can use 09d data (from Waxman's Lab Data), segex data (from Waxman's Lab Data), and other suitable data. This will generate a series of data, including heatmaps (original dimensions and print-fit dimensions), processed data, and more. For further information, please contact Bingtian Ye (btye@bu.edu), including how to build the environment with renv and how to run the content.
To this: This code uses shinyapp to generate heatmaps using pheatmap(1.0.12), and runs in three different modes: The 09d mode is a customized R-based analytical pipeline to visualize differentially expressed genes (DEGs) across multiple experimental conditions. The Segex mode is a dedicated R-based visualization module to generate heatmaps from transcriptomic signal intensity data without applying differential expression filtering. The ChIP-seq modevisualizes signal intensities from ChIP-seq or similar genomic region-based data using hierarchical clustering heatmaps. Input matrices are derived from one, or multiple, sheets of an Excel file, where rows represent genomic regions and columns correspond to distinct experimental samples and conditions. Optionally, log2-to-linear or linear-to-log2 transformations are performed prior to normalization. For further information, please contact Bingtian Ye (btye@bu.edu), including how to build the environment with renv and how to run the content. For other questions, contact David Waxman (djw@bu.edu).

---

## Acknowledgements

Created by: Bingtian Ye, laboratory of David Waxman, Dept. of Biology, Boston University

Grant support: Development of this script was carried out in the laboratory of David J Waxman, Boston University, supported in part by NIH grants R01-DK121998 (Growth Hormone Regulation of Sex Differences in Liver Metabolism) and R01-ES024421 (Epigenetic Actions of Environmental Chemicals) to DJW.

---
