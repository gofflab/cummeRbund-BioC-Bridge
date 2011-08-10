\name{PINK1}
\alias{PINK1}
\docType{data}
\title{
PINK1
}
\description{
A sample 'CuffGene' dataset
}
\usage{data(sampleData)}
\format{
  The format is:
Formal class 'CuffGene' [package "cummeRbund"] with 7 slots
  ..@ id        : chr "PINK1"
  ..@ isoforms  :Formal class 'CuffFeature' [package "cummeRbund"] with 3 slots
  .. .. ..@ annotation:'data.frame':	2 obs. of  16 variables:
  .. .. .. ..$ isoform_id     : chr [1:2] "TCONS_00000480" "TCONS_00000481"
  .. .. .. ..$ gene_id        : chr [1:2] "XLOC_000172" "XLOC_000172"
  .. .. .. ..$ CDS_id         : chr [1:2] "NA" "NA"
  .. .. .. ..$ TSS_group_id   : chr [1:2] "TSS264" "TSS265"
  .. .. .. ..$ class_code     : chr [1:2] "=" "="
  .. .. .. ..$ nearest_ref_id : chr [1:2] "uc001bdm.2" "uc001bdn.2"
  .. .. .. ..$ locus          : chr [1:2] "chr1:20959947-20978003" "chr1:20959947-20978003"
  .. .. .. ..$ length         : int [1:2] 2659 1658
  .. .. .. ..$ coverage       : num [1:2] NA NA
  .. .. .. ..$ gene_id        : chr [1:2] "XLOC_000172" "XLOC_000172"
  .. .. .. ..$ class_code     : chr [1:2] NA NA
  .. .. .. ..$ nearest_ref_id : chr [1:2] NA NA
  .. .. .. ..$ gene_short_name: chr [1:2] "PINK1" "PINK1"
  .. .. .. ..$ locus          : chr [1:2] "chr1:20959947-20978003" "chr1:20959947-20978003"
  .. .. .. ..$ length         : int [1:2] NA NA
  .. .. .. ..$ coverage       : num [1:2] NA NA
  .. .. ..@ fpkm      :'data.frame':	6 obs. of  6 variables:
  .. .. .. ..$ isoform_id  : chr [1:6] "TCONS_00000480" "TCONS_00000481" "TCONS_00000480" "TCONS_00000481" ...
  .. .. .. ..$ sample_name : Factor w/ 3 levels "hESC","Fibroblasts",..: 2 2 1 1 3 3
  .. .. .. ..$ fpkm        : num [1:6] 2102 818 574 120 1598 ...
  .. .. .. ..$ conf_hi     : num [1:6] 3111 1392 669 153 2282 ...
  .. .. .. ..$ conf_lo     : num [1:6] 1091.9 243.7 478.3 87.2 913.7 ...
  .. .. .. ..$ quant_status: chr [1:6] "OK" "OK" "OK" "OK" ...
  .. .. ..@ diff      :'data.frame':	6 obs. of  11 variables:
  .. .. .. ..$ isoform_id    : chr [1:6] "TCONS_00000480" "TCONS_00000481" "TCONS_00000480" "TCONS_00000481" ...
  .. .. .. ..$ sample_1      : Factor w/ 3 levels "hESC","Fibroblasts",..: 1 1 1 1 2 2
  .. .. .. ..$ sample_2      : Factor w/ 3 levels "hESC","Fibroblasts",..: 2 2 3 3 3 3
  .. .. .. ..$ status        : chr [1:6] "OK" "OK" "OK" "OK" ...
  .. .. .. ..$ value_1       : num [1:6] 574 120 574 120 2102 ...
  .. .. .. ..$ value_2       : num [1:6] 2102 818 1598 0 1598 ...
  .. .. .. ..$ ln_fold_change: num [1:6] 1.87 2.77 1.48 -1.80e+308 -3.95e-01 ...
  .. .. .. ..$ test_stat     : num [1:6] -5.11 -5.10 -4.46 -1.80e+308 8.51e-01 ...
  .. .. .. ..$ p_value       : num [1:6] 3.22e-07 3.44e-07 8.10e-06 1.14e-13 3.95e-01 ...
  .. .. .. ..$ q_value       : num [1:6] 4.77e-06 5.04e-06 9.20e-05 4.63e-12 6.76e-01 ...
  .. .. .. ..$ significant   : chr [1:6] "yes" "yes" "yes" "yes" ...
  ..@ TSS       :Formal class 'CuffFeature' [package "cummeRbund"] with 3 slots
  .. .. ..@ annotation:'data.frame':	2 obs. of  14 variables:
  .. .. .. ..$ TSS_group_id   : chr [1:2] "TSS264" "TSS265"
  .. .. .. ..$ class_code     : chr [1:2] NA NA
  .. .. .. ..$ nearest_ref_id : chr [1:2] NA NA
  .. .. .. ..$ gene_id        : chr [1:2] "XLOC_000172" "XLOC_000172"
  .. .. .. ..$ locus          : chr [1:2] "chr1:20959947-20978003" "chr1:20959947-20978003"
  .. .. .. ..$ length         : int [1:2] NA NA
  .. .. .. ..$ coverage       : num [1:2] NA NA
  .. .. .. ..$ gene_id        : chr [1:2] "XLOC_000172" "XLOC_000172"
  .. .. .. ..$ class_code     : chr [1:2] NA NA
  .. .. .. ..$ nearest_ref_id : chr [1:2] NA NA
  .. .. .. ..$ gene_short_name: chr [1:2] "PINK1" "PINK1"
  .. .. .. ..$ locus          : chr [1:2] "chr1:20959947-20978003" "chr1:20959947-20978003"
  .. .. .. ..$ length         : int [1:2] NA NA
  .. .. .. ..$ coverage       : num [1:2] NA NA
  .. .. ..@ fpkm      :'data.frame':	6 obs. of  6 variables:
  .. .. .. ..$ TSS_group_id: chr [1:6] "TSS264" "TSS265" "TSS264" "TSS265" ...
  .. .. .. ..$ sample_name : Factor w/ 3 levels "hESC","Fibroblasts",..: 2 2 1 1 3 3
  .. .. .. ..$ fpkm        : num [1:6] 2102 818 574 120 1598 ...
  .. .. .. ..$ conf_hi     : num [1:6] 3080 1237 663 149 2282 ...
  .. .. .. ..$ conf_lo     : num [1:6] 1123.4 398.5 483.7 90.8 913.7 ...
  .. .. .. ..$ quant_status: chr [1:6] "OK" "OK" "OK" "OK" ...
  .. .. ..@ diff      :'data.frame':	6 obs. of  11 variables:
  .. .. .. ..$ TSS_group_id  : chr [1:6] "TSS264" "TSS265" "TSS264" "TSS265" ...
  .. .. .. ..$ sample_1      : Factor w/ 3 levels "hESC","Fibroblasts",..: 1 1 1 1 2 2
  .. .. .. ..$ sample_2      : Factor w/ 3 levels "hESC","Fibroblasts",..: 2 2 3 3 3 3
  .. .. .. ..$ status        : chr [1:6] "OK" "OK" "OK" "OK" ...
  .. .. .. ..$ value_1       : num [1:6] 574 120 574 120 2102 ...
  .. .. .. ..$ value_2       : num [1:6] 2102 818 1598 0 1598 ...
  .. .. .. ..$ ln_fold_change: num [1:6] 1.87 2.77 1.48 -1.80e+308 -3.95e-01 ...
  .. .. .. ..$ test_stat     : num [1:6] -5.29 -6.76 -4.49 -1.80e+308 8.66e-01 ...
  .. .. .. ..$ p_value       : num [1:6] 1.23e-07 1.34e-11 6.96e-06 1.04e-16 3.86e-01 ...
  .. .. .. ..$ q_value       : num [1:6] 1.32e-06 2.43e-10 4.87e-05 4.10e-15 5.70e-01 ...
  .. .. .. ..$ significant   : chr [1:6] "yes" "yes" "yes" "yes" ...
  ..@ CDS       :Formal class 'CuffFeature' [package "cummeRbund"] with 3 slots
  .. .. ..@ annotation:'data.frame':	2 obs. of  15 variables:
  .. .. .. ..$ CDS_id         : chr [1:2] "P363" "P364"
  .. .. .. ..$ class_code     : chr [1:2] NA NA
  .. .. .. ..$ nearest_ref_id : chr [1:2] NA NA
  .. .. .. ..$ gene_id        : chr [1:2] "XLOC_000172" "XLOC_000172"
  .. .. .. ..$ TSS_group_id   : chr [1:2] "TSS265" "TSS264"
  .. .. .. ..$ locus          : chr [1:2] "chr1:20959947-20978003" "chr1:20959947-20978003"
  .. .. .. ..$ length         : int [1:2] NA NA
  .. .. .. ..$ coverage       : num [1:2] NA NA
  .. .. .. ..$ gene_id        : chr [1:2] "XLOC_000172" "XLOC_000172"
  .. .. .. ..$ class_code     : chr [1:2] NA NA
  .. .. .. ..$ nearest_ref_id : chr [1:2] NA NA
  .. .. .. ..$ gene_short_name: chr [1:2] "PINK1" "PINK1"
  .. .. .. ..$ locus          : chr [1:2] "chr1:20959947-20978003" "chr1:20959947-20978003"
  .. .. .. ..$ length         : int [1:2] NA NA
  .. .. .. ..$ coverage       : num [1:2] NA NA
  .. .. ..@ fpkm      :'data.frame':	6 obs. of  6 variables:
  .. .. .. ..$ CDS_id      : chr [1:6] "P363" "P364" "P363" "P364" ...
  .. .. .. ..$ sample_name : Factor w/ 3 levels "hESC","Fibroblasts",..: 2 2 1 1 3 3
  .. .. .. ..$ fpkm        : num [1:6] 818 2102 120 574 0 ...
  .. .. .. ..$ conf_hi     : num [1:6] 1237 3080 149 663 0 ...
  .. .. .. ..$ conf_lo     : num [1:6] 398.5 1123.4 90.8 483.7 0 ...
  .. .. .. ..$ quant_status: chr [1:6] "OK" "OK" "OK" "OK" ...
  .. .. ..@ diff      :'data.frame':	6 obs. of  11 variables:
  .. .. .. ..$ CDS_id        : chr [1:6] "P363" "P364" "P363" "P364" ...
  .. .. .. ..$ sample_1      : Factor w/ 3 levels "hESC","Fibroblasts",..: 1 1 1 1 2 2
  .. .. .. ..$ sample_2      : Factor w/ 3 levels "hESC","Fibroblasts",..: 2 2 3 3 3 3
  .. .. .. ..$ status        : chr [1:6] "OK" "OK" "OK" "OK" ...
  .. .. .. ..$ value_1       : num [1:6] 120 574 120 574 818 ...
  .. .. .. ..$ value_2       : num [1:6] 818 2102 0 1598 0 ...
  .. .. .. ..$ ln_fold_change: num [1:6] 2.77 1.87 -1.80e+308 1.48 -1.80e+308 ...
  .. .. .. ..$ test_stat     : num [1:6] -6.76 -5.29 -1.80e+308 -4.49 -1.80e+308 ...
  .. .. .. ..$ p_value       : num [1:6] 1.34e-11 1.23e-07 1.04e-16 6.96e-06 4.79e-05 ...
  .. .. .. ..$ q_value       : num [1:6] 2.63e-10 1.40e-06 4.57e-15 5.52e-05 3.10e-04 ...
  .. .. .. ..$ significant   : chr [1:6] "yes" "yes" "yes" "yes" ...
  ..@ annotation:'data.frame':	1 obs. of  7 variables:
  .. ..$ gene_id        : chr "XLOC_000172"
  .. ..$ class_code     : chr NA
  .. ..$ nearest_ref_id : chr NA
  .. ..$ gene_short_name: chr "PINK1"
  .. ..$ locus          : chr "chr1:20959947-20978003"
  .. ..$ length         : int NA
  .. ..$ coverage       : num NA
  ..@ fpkm      :'data.frame':	3 obs. of  6 variables:
  .. ..$ gene_id     : chr [1:3] "XLOC_000172" "XLOC_000172" "XLOC_000172"
  .. ..$ sample_name : Factor w/ 3 levels "hESC","Fibroblasts",..: 2 1 3
  .. ..$ fpkm        : num [1:3] 2919 693 1598
  .. ..$ conf_hi     : num [1:3] 4003 814 2282
  .. ..$ conf_lo     : num [1:3] 1836 573 914
  .. ..$ quant_status: chr [1:3] "OK" "OK" "OK"
  ..@ diff      :'data.frame':	3 obs. of  11 variables:
  .. ..$ gene_id       : chr [1:3] "XLOC_000172" "XLOC_000172" "XLOC_000172"
  .. ..$ sample_1      : Factor w/ 3 levels "hESC","Fibroblasts",..: 1 1 2
  .. ..$ sample_2      : Factor w/ 3 levels "hESC","Fibroblasts",..: 2 3 3
  .. ..$ status        : chr [1:3] "OK" "OK" "OK"
  .. ..$ value_1       : num [1:3] 693 693 2919
  .. ..$ value_2       : num [1:3] 2919 1598 1598
  .. ..$ ln_fold_change: num [1:3] 2.07 1.2 -0.87
  .. ..$ test_stat     : num [1:3] -7.02 -3.61 2.13
  .. ..$ p_value       : num [1:3] 2.29e-12 3.02e-04 3.35e-02
  .. ..$ q_value       : num [1:3] 4.04e-11 1.22e-03 7.65e-02
  .. ..$ significant   : chr [1:3] "yes" "yes" "no"
}
\details{
Sample CuffGene data for gene 'PINK1'
}
\source{
None
}
\references{
None
}
\examples{
data(sampleData)
PINK1
}
\keyword{datasets}
