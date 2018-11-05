# ERDS

ERDS is designed for detection of copy number variants (CNVs) on human genomes from next generation sequence data, utilizing information from read depth of short reads and SNV heterozygosity.

The original version can be found from subdirectory named "erds1.1".

Brett Trost and Joe Whitney (The Centre for Applied Genomics, Hospital for Sick Children, Toronto) improved memory usage and performance of ERDS and their team published a paper describing a workflow for detecting CNVs from WGS data (The American Journal of Human Genetics 102.1 (2018): 142-155). The modifed code version can be found from subdirectory named "tcag". A simple compilation step is required - just run "make" within the "src" directory.

### Citation

If you use ERDS for publication, please cite:

Zhu M, Need AC, Han Y, Ge D, Maia JM, Zhu Q, Heinzen EL, Cirulli ET, Pelak K, He M, Ruzzo EK, Gumbs CE, Singh A, Feng S, Shianna KV and Goldstein DB. Using ERDS to Infer Copy-Number Variants in High-Coverage Genomes. The American Journal of Human Genetics. Volume 91, Issue 3, 7 September 2012, Pages 408-421.
