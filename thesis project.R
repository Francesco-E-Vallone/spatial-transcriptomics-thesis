library(spatialLIBD)
library(SingleCellExperiment)
library(SpatialExperiment)
library(nnSVG)
library(PlackettLuce)
library(dplyr)
library(ggplot2)
library(rmarkdown)
library(knitr)
library(kableExtra)

#data management ####
#upload Visium data
get(load("Human_DLPFC_Visium_processedData_sce_scran_spatialLIBD.Rdata"))
get(load("Human_DLPFC_Visium_processedData_sce_scran_sce_layer_spatialLIBD.Rdata"))
get(load("Human_DLPFC_Visium_modeling_results.Rdata" ))
sce <- spatialLIBD:::.update_sce(sce) #solving the check_sce issue (solution from github)
sce_layer <- spatialLIBD:::.update_sce_layer(sce_layer)

#exploring data
check_sce(sce)
sce_layer
modeling_results

#convert to SpatialExperiment obj
spe <- sce_to_spe(sce)

#converting the modeling statistics into a long format (more manageable)
system.time(
  sig_genes <-
    sig_genes_extract_all(
      n = nrow(sce_layer),
      modeling_results = modeling_results,
      sce_layer = sce_layer
    )
)

head(sig_genes)

#visualisation of spatial infos & genes####
#visualisation of spatial info for 1 sample
vis_clus(
  spe = spe,
  clustervar = "layer_guess_reordered",
  sampleid = "151673",
  colors = libd_layer_colors,
  ... = " LIBD Layers"
)

#setting marker genes for the white matter
white_matter_genes <- c("GFAP", "AQP4", "MBP", "PLP1", "NCDN", "ELOVL1")
white_matter_genes <- rowData(spe)$gene_search[
  rowData(spe)$gene_name %in% white_matter_genes
]

white_matter_genes

#plotting one gene (e.g. GFAP)
vis_gene(
  spe,
  geneid = white_matter_genes[2],
  point_size = 1.5
) #we can visually check where the gene in most enriched

#checking for its ranking 
i_gfap <- subset(sig_genes, gene == "ELOVL1" &
                   test == "WM")$top
i_gfap #position in the rank

#checking quantitatively in which layer the gene is most enriched
layer_boxplot(
  i = i_gfap,
  sig_genes = sig_genes,
  sce_layer = sce_layer
)

#multi-gene visualisation
vis_gene(
  spe,
  geneid = white_matter_genes,
  multi_gene_method = "pca",
  point_size = 1.5
)  #multi_gene_method = z-score / pca
gc()

#spatial registration ####
#loading single nucleus RNA-seq data
bfc <- BiocFileCache::BiocFileCache()
url <- paste0(
  "https://libd-snrnaseq-pilot.s3.us-east-2.amazonaws.com/",
  "SCE_DLPFC-n3_tran-etal.rda"
)
local_data <- BiocFileCache::bfcrpath(url, x = bfc)

get(load(local_data, verbose = T))

#inspect obj (sce.dlpfc.tran from "local_data")
table(sce.dlpfc.tran$cellType)

#performing spatial registration (computes pseudo-bulking, Bayesian correlation, and t-statistics)
sce_modeling_results <- registration_wrapper(
  sce = sce.dlpfc.tran,
  var_registration = "cellType",
  var_sample_id = "donor",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

#extracting t-statistics
registration_t_stats <- sce_modeling_results$enrichment[, grep("^t_stat", colnames(sce_modeling_results$enrichment))]
colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats)) #renaming colnames

#visualising cell type x gene
dim(registration_t_stats)
gc()

#correlating t-statistics with layer info (layer reference)
cor_layer <- layer_stat_cor(
  stats = registration_t_stats,
  modeling_results = modeling_results,
  model_type = "enrichment",
  top_n = 100
)
cor_layer

layer_stat_cor_plot(cor_layer, max = max(cor_layer)) #visualisation

#annotating cell types based on correlation values (confidence about a cell type)
ann <- annotate_registered_clusters(
  cor_stats_layer = cor_layer,
  confidence_threshold = 0.25,
  cutoff_merge_ratio = 0.25
)

ann

#writing table
rmd_text <- "
---
title: 'Table Output'
output: html_document
---

```{r, echo=FALSE}
library(kableExtra)
# Use the existing dataframe (e.g., ann)
kable(ann, format = 'html', booktabs = TRUE) %>%
  kable_styling('striped', full_width = F)"
writeLines(rmd_text, con = "output_table.Rmd")
render("output_table.Rmd", output_format = "html_document")

#gene ranking analyses ####
#keeping spots over tissue only
spe2 <- spe[, colData(spe)$in_tissue == T]
dim(spe)
dim(spe2)

##downsampling spots
available_spots <- colnames(spe2)

#setting a seed and sampling the spots
set.seed(100) #for reproducibility (to ensure that the same spots are taken if the analysis is re-run)
sampled_spots <- sample(available_spots, size = 20000)

#subsetting using the sampled spots 
spe2 <- spe2[, sampled_spots]
dim(spe2)

spe2 <- filter_genes(spe2, filter_genes_ncounts = 3, filter_genes_pcspots = 0.7)
dim(spe2)

##dealing with duplicates (there was a problem of spot duplicates that didn't allow for the nnSVG computation)
coords <- spatialCoords(spe2)
duplicates <- duplicated(coords) #keeping only the first instance of a duplicated coord

#keeping only unique coordinates
spe2 <- spe2[, !duplicates]
dim(spe2)

##running nnSVG
#since the spe (now spe2 object already containsthe logcounts,I'll skip the calc of logcounts using "scran")
spe2 <- nnSVG(spe2, 
             n_neighbors = 10,  
             n_threads = 1,
             verbose = TRUE)

#number of significant SVGs
table(rowData(spe2)$padj<=0.05)

#checking results
rowData(spe2)[order(rowData(spe2)$rank)[1:15],]

#plotting spatial expression of top-ranked SVG
ix <- which(rowData(spe2)$gene_name == "SSC1")
ix_name <- rowData(spe2)$gene_name[ix]
ix_name

df <- as.data.frame(cbind(spatialCoords(spe2), 
        expr = counts(spe2)[ix, ]))

ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = expr)) + 
  geom_point(size = 0.8) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_gradient(low = "gray90", high = "blue", 
                       trans = "sqrt", breaks = range(df$expr), 
                       name = "counts") + 
  ggtitle(ix_name) + 
  theme_bw() + 
  theme(plot.title = element_text(face = "italic"), 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

##running an nnSVG analysis for the white matter
#filtering out NA rows 
valid_spots <- !is.na(colData(spe)$layer_guess_reordered)

#subsetting the data for valid spots
spe_valid <- spe[, valid_spots]

#taking only spots located in the white matter
white_matter_spots <- spe_valid[, colData(spe_valid)$layer_guess_reordered == "WM"]
rownames(white_matter_spots) <- rowData(white_matter_spots)$gene_name
gc()

#filtering genes
white_matter_spots <- filter_genes(white_matter_spots,
                                   filter_genes_ncounts = 5,
                                   filter_genes_pcspots = 0.7)
dim(white_matter_spots)
gc()

#running nnSVG analysis for the white matter spots
white_matter_spots <- nnSVG(white_matter_spots, 
                            n_neighbors = 10,  
                            n_threads = 1,
                            verbose = TRUE)
gc()

#checking results
table(rowData(white_matter_spots)$padj <= 0.05)

#checking top nnSVG results
top_svg_wm <- rowData(white_matter_spots)[order(rowData(white_matter_spots)$rank)[1:15],]
print(top_svg_wm)

#plotting data
ix_wm <- which(rowData(white_matter_spots)$gene_name == "NCDN")
ix_name_wm <- rowData(white_matter_spots)$gene_name[ix_wm]

df_wm <- as.data.frame(cbind(spatialCoords(white_matter_spots), 
                             expr = counts(white_matter_spots)[ix_wm, ]))

ggplot(df_wm, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = expr)) + 
  geom_point(size = 0.8) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_gradient(low = "gray90", high = "blue", 
                       trans = "sqrt", breaks = range(df_wm$expr), 
                       name = "counts") + 
  ggtitle(ix_name_wm) + 
  theme_bw() + 
  theme(plot.title = element_text(face = "italic"), 
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())


##running nnSVG analysis of Layer 2
#subsetting layer2 data
layer_index <- colData(spe)$layer_guess_reordered %in% "Layer2"
layer_spots_spe <- spe[, layer_index]
rownames(layer_spots_spe) <- rowData(layer_spots_spe)$gene_name
dim(layer_spots_spe)
gc()

#filtering genes
layer_spots_spe <- filter_genes(layer_spots_spe,
                                filter_genes_ncounts = 5,
                                filter_genes_pcspots = 0.7)
dim(layer_spots_spe)
gc()

#running nnSVG 
layer_spots_spe <- nnSVG(layer_spots_spe, 
                         n_neighbors = 10,  
                         n_threads = 1,
                         verbose = TRUE)
gc()

##building the Plackett-Luce model for two Layers (WM and Layer 2)
#select significant genes from the WM
significant_genes_wm <- rowData(white_matter_spots)[rowData(white_matter_spots)$padj <= 0.05, ]

#extracting logcounts for significant genes
logcounts_wm <- assay(white_matter_spots, "logcounts")[rownames(significant_genes_wm), ] 
logcounts_layer <- assay(layer_spots_spe, "logcounts") #I'll compare the genes that appear as the most significant in the WM

#creating a ranking for White Matter
rankings_wm <- apply(logcounts_wm, 1, function(x) rank(-x, ties.method = "first"))
rankings_wm <- rankings_wm[, 1:15]

rankings_wm_pl <- as.rankings(rankings_wm)
pl_model_wm <- PlackettLuce(rankings_wm_pl) #PL model

coef_wm <- coef(pl_model_wm) #extracting coefficients (log-worths)

#creating a ranking for Layer 2
rankings_layer <- apply(logcounts_layer, 1, function(x) rank(-x, ties.method = "first"))
rankings_layer <- rankings_layer[, colnames(rankings_layer) %in% colnames(rankings_wm)]

rankings_layer_pl <- as.rankings(rankings_layer)
pl_model_layer <- PlackettLuce(rankings_layer_pl) #PL model

coef_layer <- coef(pl_model_layer) #extracting coefficients (log-worths)

##single plots
#layer 2
df_layer <- data.frame(Gene = names(coef_layer), Coefficient = coef_layer)
df_layer <- df_layer[!df_layer$Gene %in% c("tie2", "tie3"), ]

#plotting layer 2
ggplot(df_layer, aes(x = reorder(Gene, Coefficient), y = Coefficient)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Coefficients for Layer 2",
       x = "Genes",
       y = "Coefficient Value") +
  theme_minimal()

#MW
df_wm <- data.frame(Gene = names(coef_wm), Coefficient = coef_wm)
df_wm <- df_wm[!df_wm$Gene %in% c("tie2", "tie3", "tie4"), ]

#plotting WM
ggplot(df_wm, aes(x = reorder(Gene, Coefficient), y = Coefficient)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  coord_flip() +
  labs(title = "Coefficients for White Matter",
       x = "Genes",
       y = "Coefficient Value") +
  theme_minimal()
  