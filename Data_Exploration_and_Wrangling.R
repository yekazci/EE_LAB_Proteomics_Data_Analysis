setwd("D:/2_POST_DOC/2_PROTEOMICS/OXIDATIVE_STRESS_ANALYSIS_1")


data <- read.delim("EE_LAB_data.txt")
head(data)
library(tidyverse)

data_selected <- data %>% select(-c("X.3","X.4"))## we remove the columns that we dont want.
head(data_selected)
data_abundances <- data_selected %>% select(colnames(data_selected)[c(1,8:ncol(data_selected))])## we remove the colums that we dont want.
data_abundances %>% head
colnames(data_abundances) <- c("IDs", "mito_n", "NES_n", "NLS_n", "mito_r", "NES_r", "NLS_r")
data_abundances <- data_abundances[-(1:2),] ## remove first two rows.

head(data_abundances)
data_abundances <- data.frame(data_abundances[,-1], row.names = data_abundances$IDs)
head(data_abundances)

#                         mito_n       NES_n       NLS_n      mito_r       NES_r       NLS_r
# Q9NZL9               1498.287199 1431.553645 1866.794285 1559.735039 1431.553645 1429.003029
# A1A5D9               2173.752159  3331.36228 9937.810657 2262.902207  3331.36228 7607.245021
# Q9NR45               3499.173026 3185.488001 3361.262315  3642.68131 3185.488001 2572.995893
# Q01085;C9JTN7;F8WE16 328.9601445 290.5244742 195.2229413 342.4514767 290.5244742 149.4402338
# Q9UHV9               4186.473921 3692.417489  4469.70679 4358.169828 3692.417489 3421.493515
# Q6TFL3;H0Y5M5;H0Y701 1142.721589 935.6942709 1038.638704 1189.586952 935.6942709 795.0623517

data_abundances[,1:6] <- apply(data_abundances, 2, function(x) as.numeric(as.character(x))) ## convert charaters to numbers.

data_abundances %>% head()

data_abundances %>% class()

table(apply(data_abundances,2,function(x) is.numeric(x))) ### check if all values are numeric.

# TRUE 
# 6

data_raw <- data_abundances[,4:6]
head(data_raw)



nrow(data_raw)
table(is.na(data_raw))
data_log2 <- log2(data_raw+1) ### log2 normalisation. we add one to avoid log2 of zero.

annotation_col <- data.frame(names=colnames(data_log2), row.names = colnames(data_log2))

pheatmap::pheatmap(mat = as.matrix(data_log2), show_rownames=FALSE, annotation_col = annotation_col, dendrogram=c("row"))

data_log2_mat <- as.matrix(data_log2)

top100 <- order(matrixStats::rowVars(data_log2_mat), decreasing = TRUE)[1:100] ### proteins with top variance across samples:
pheatmap::pheatmap(mat = data_log2_mat[top100,], show_rownames=FALSE, annotation_col = annotation_col)



### prepare gene ids for top 100 genes:

top_100_genes <- data.frame(protein_IDs = row.names(data_log2_mat[top100,]))
top_100_genes %>% head()


data_selected_id_desc <- select(data_selected, c("X", "X.8")) ### selct cos with ids and descriptions.
data_selected_id_desc %>% head()
colnames(data_selected_id_desc) <- c("IDs", "Descriptions")
data_selected_id_desc %>% head()
data_selected_id_desc <- data_selected_id_desc[-(1:2),] ### remove unnecessary rows.
data_selected_id_desc %>% head()

row.names(data_selected_id_desc) <- data_selected_id_desc$IDs ### add row.names as ids in order to subset decsriptions in top100 list order.

data_selected_id_desc %>% head()
data_selected_id_desc[top_100_genes$protein_IDs,]
data_top100_table <- data_selected_id_desc[top_100_genes$protein_IDs,]

### short the id only selecting first ones:

data_top100_ids_short <- strsplit(data_top100_table$IDs, split = ";")


head(data_top100_ids_short)

data_top100_ids_short[1]

data_top100_table$short_IDs <- unlist(lapply(data_top100_ids_short, function(x) as.vector(x)[[1]]))

library(org.Hs.eg.db)

keytypes(org.Hs.eg.db)  ### look at the options.



data_top100_table$genes <- mapIds(org.Hs.eg.db, keys=data_top100_table$short_IDs, column='SYMBOL', keytype='UNIPROT', multiVals='first')



data_top100_table %>% dplyr::select(c("genes", "short_IDs")) ### select function was masked by annotationdb package. So specify the select fun from dplyr.


which(is.na(data_top100_table$genes))  ### indices with no gene name. I manually added this searching google.
# [1]  1  2  3  8 20 22 23 33 38 42 59 61 71 74 79 82 83 88 93 

data_top100_table$genes[which(is.na(data_top100_table$genes))] <- c("EIF2AK2","GOLT1B","GATD3","LRRC28","SDCCAG8","RPS15","CACNA1A","MYO15B","HSD17B12","MRM2","MYO15B","DYNLRB1","CCDC38","RPA3","NME1-NME2","FDFT1","SEC11A","RER1","40S_ribosomal_pro")

which(is.na(data_top100_table$genes))  ### check if any NA remained:
# integer(0)

mat_top100 <- data_log2_mat[top100,]
row.names(mat_top100) <- data_top100_table$genes


my_clustered_data <- pheatmap::pheatmap(mat = mat_top100, show_rownames=TRUE, annotation_col = annotation_col, fontsize_row = 5);

### Take the ordered rows:

my_clustered_data$tree_row$order

# [1]  41  15  89  71  90  20  26  37  18  19  46  58  45  53  50  54  99  77
# [19]  81  88   9  49  80  32  27  31  40  55  25  39  43  83  76  64  72   3
# [37]   6  11   1   2   4   5   7   8  16  14  61  38  74  63  97  57  82  62
# [55]  70  91  94  69  84  86  96  60  67  92  17  22  21  28  30  23  48  79
# [73]  42  56  59  36  98  51  29  93  85 100  12  10  13  75  66  68  65  95
# [91]  47  33  35  78  87  52  73  24  34  44

mat_top100[my_clustered_data$tree_row$order,]

write.csv(mat_top100[my_clustered_data$tree_row$order,], file = "top100_genes.txt")

sessionInfo()

# R version 4.2.1 (2022-06-23 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19043)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=Turkish_Turkey.utf8  LC_CTYPE=Turkish_Turkey.utf8   
# [3] LC_MONETARY=Turkish_Turkey.utf8 LC_NUMERIC=C                   
# [5] LC_TIME=Turkish_Turkey.utf8    
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods  
# [8] base     
# 
# other attached packages:
#   [1] org.Hs.eg.db_3.15.0  AnnotationDbi_1.58.0 IRanges_2.30.1      
# [4] S4Vectors_0.34.0     Biobase_2.56.0       BiocGenerics_0.42.0 
# [7] forcats_0.5.2        stringr_1.4.1        dplyr_1.0.10        
# [10] purrr_0.3.4          readr_2.1.2          tidyr_1.2.1         
# [13] tibble_3.1.8         ggplot2_3.3.6        tidyverse_1.3.2     
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.9             lubridate_1.8.0        Biostrings_2.64.1     
# [4] png_0.1-7              assertthat_0.2.1       utf8_1.2.2            
# [7] GenomeInfoDb_1.32.4    R6_2.5.1               cellranger_1.1.0      
# [10] backports_1.4.1        reprex_2.0.2           RSQLite_2.2.17        
# [13] httr_1.4.4             pillar_1.8.1           zlibbioc_1.42.0       
# [16] rlang_1.0.5            googlesheets4_1.0.1    readxl_1.4.1          
# [19] rstudioapi_0.14        blob_1.2.3             googledrive_2.0.0     
# [22] RCurl_1.98-1.8         pheatmap_1.0.12        bit_4.0.4             
# [25] munsell_0.5.0          broom_1.0.1            compiler_4.2.1        
# [28] modelr_0.1.9           pkgconfig_2.0.3        tidyselect_1.1.2      
# [31] KEGGREST_1.36.3        GenomeInfoDbData_1.2.8 matrixStats_0.62.0    
# [34] fansi_1.0.3            crayon_1.5.1           tzdb_0.3.0            
# [37] dbplyr_2.2.1           withr_2.5.0            bitops_1.0-7          
# [40] grid_4.2.1             jsonlite_1.8.0         gtable_0.3.1          
# [43] lifecycle_1.0.2        DBI_1.1.3              magrittr_2.0.3        
# [46] scales_1.2.1           cli_3.4.0              stringi_1.7.8         
# [49] cachem_1.0.6           farver_2.1.1           XVector_0.36.0        
# [52] fs_1.5.2               xml2_1.3.3             ellipsis_0.3.2        
# [55] generics_0.1.3         vctrs_0.4.1            RColorBrewer_1.1-3    
# [58] tools_4.2.1            bit64_4.0.5            glue_1.6.2            
# [61] hms_1.1.2              fastmap_1.1.0          colorspace_2.0-3      
# [64] gargle_1.2.1           rvest_1.0.3            memoise_2.0.1         
# [67] haven_2.5.1   