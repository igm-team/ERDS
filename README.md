# ERDS

ERDS is designed for detection of copy number variants (CNVs) on human genomes from next generation sequence data, utilizing information from read depth of short reads and SNV heterozygosity.

The original version can be found from subdirectory named "original".

Brett Trost improved memory usage and performance of ERDS and his team published a paper in AJHG describing a workflow for detecting CNVs from WGS data. The modifed code version can be found from subdirectory named "tcag". A simple compilation step is required - just run "make" within the "src" directory.

### Citation

If you use ERDS for publication, please cite:

Zhu, Mingfu, et al. "Using ERDS to infer copy-number variants in high-coverage genomes." The American Journal of Human Genetics 91.3 (2012): 408-421.
