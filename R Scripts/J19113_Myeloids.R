#Jayalaxmi Suresh Babu/ Won Jin Ho
#Last updated 9/27/2021
#R version 3.6.2 (2019-12-12)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 19042)

rm(list = ls())
library(reshape2)
library(randomcoloR)
library(pals)
library(ggplot2)
library(Hmisc)
library(edgeR)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(corrplot)
library(readxl)
library(ggridges)
library(ggpubr)
memory.limit(size=56000) #may need this to run locally

####READ and CLUSTER FUNCTIONS####
returnfcs <- function(FDR_cutoff=.05,
                      metaDataFile='~/',
                      panelDataFile='~/',
                      dataDirectory='~/',
                      shape_timepoint=NULL,
                      color_timepoint=NULL){
  #This function generates an fcs file, subtype_markers, colors and shapes for clustering 
  require(scales);require(readxl);require(dplyr);require(flowCore)
  ##directory and metadatafile checking
  if(!dir.exists(dataDirectory)) {stop('ERR: cannot find data directory')}
  if(!file.exists(metaDataFile)) {stop('ERR: cannot find metadata.xlsx or .csv file')}
  ##readin metadata and clean
  ifelse(grepl(metaDataFile,pattern='.xls'),md <- read_excel(metaDataFile),md <- read.csv(metaDataFile,header = TRUE))#must be in xl format or csv
  md$benefit <- factor(md$benefit)
  md$batch <- factor(md$batch)
  md$timepoint <- factor(md$timepoint)
  md$OS <- factor(md$OS_6mo)
  md$patientID <- factor(md$patient_id)
  md$recist <- factor(md$recist)
  rownames(md) = md$sample_id;md$sample_id <- md$sample_id
  #Make sure all files in metadata present in datadirectory
  if(!all(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])){
    print(paste('ERR: not all filenames in metadata present in data folder - missing',md$file_name[!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])],'Subsetting...'))
    md <- md[-c(!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])),]
  }
  ##Define shapes for timepoint
  if(is.null(shape_timepoint)){shape_timepoint <- c(0:25)[1:length(levels(md$timepoint))]}#can specify as long as number is same
  if(length(shape_timepoint)!=length(levels(md$timepoint))){stop(paste0('ERR no. shapes specified is less than no. of timepoint (',length(levels(md$timepoint)),')'))}
  names(shape_timepoint) <- levels(md$timepoint)
  ## Define colors for the timepoint
  if(is.null(color_timepoint)){color_timepoint <- hue_pal()(length(levels(md$timepoint)))}#can specify as long as number is same
  if(length(color_timepoint)!=length(levels(md$timepoint))){stop(paste0('ERR no. shapes specified is less than no. of timepoint (',length(levels(md$timepoint)),')'))}
  ## read fcs
  fcs_raw <- read.flowSet(md$file_name, path = dataDirectory, transformation = FALSE, truncate_max_range = FALSE)
  sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))
  panel <- read_excel(panelDataFile)
  head(data.frame(panel))
  ## Replace problematic characters
  panel$Metal <- gsub('-', '_', panel$Metal)
  panel_fcs <- pData(parameters(fcs_raw[[1]]))
  panel_fcs$desc <- gsub('-', '_', panel_fcs$desc)
  panel_fcs$desc[is.na(panel_fcs$desc)] <- paste0('NA_',which(is.na(panel_fcs$desc))) #was labelled 'todo'(mapping based on isotope for now just getting rid of NA keeping rownum) 
  # use panel$Antigen to fix description in panel_fcs
  # use metal+isotope as mapping between panel from xlsx and panel from the fcs files
  rownames(panel_fcs) = panel_fcs$name
  panel_fcs[paste0(panel$Metal,panel$Isotope,'Di'),2] <- panel$Antigen
  ## Replace paramater data in flowSet
  pData(parameters(fcs_raw[[1]])) <- panel_fcs
  ## Define variables indicating marker types
  subtype_markers <- panel$Antigen[panel$Subtype == 1]
  functional_markers <- panel$Antigen[panel$Functional == 1]
  if(!all(subtype_markers %in% panel_fcs$desc)){stop('ERR: Not all subtype_markers in panel_fcs$desc (isotopes)')}
  if(!all(functional_markers %in% panel_fcs$desc)){stop('ERR: Not all functional_markers in panel_fcs$desc (isotopes)')}
  ## arcsinh transformation and column subsetting
  fcs <- fsApply(fcs_raw, function(x, cofactor = 5){
    colnames(x) <- panel_fcs$desc
    expr <- exprs(x)
    expr <- asinh(expr[, union(subtype_markers,functional_markers)] / cofactor)
    exprs(x) <- expr
    x
  })
  return(list('fcs'=fcs,
              'subtype_markers'=subtype_markers,
              'functional_markers'=functional_markers,
              'shape_timepoint'=shape_timepoint,
              'color_timepoint'=color_timepoint,
              'sample_ids'=sample_ids,
              'meta_data'=md))
}

clusterfcs <- function(fcs=output$fcs,
                       subtype_markers = output$subtype_markers,
                       seed=1234,plottitle='consensus_plots',
                       numclusters=40){
  ## Cell population identification with FlowSOM and ConsensusClusterPlus
  require(dplyr);require(FlowSOM);require(ConsensusClusterPlus)
  set.seed(seed)
  som <- ReadInput(fcs, transform = FALSE, scale = FALSE) %>% BuildSOM(colsToUse = subtype_markers)
  ## Get the cell clustering into 100 SOM codes
  cell_clustering_som <- som$map$mapping[,1]
  ## Metaclustering into numclusters with ConsensusClusterPlus
  codes <- som$map$codes
  mc <- ConsensusClusterPlus(t(codes), maxK = numclusters, reps = 100,
                             pItem = 0.9, pFeature = 1, title = plottitle, 
                             plot = "png", clusterAlg = "hc", 
                             innerLinkage = "average", finalLinkage = "average",
                             distance = "euclidean", seed = 1234)
  
  ## Get cluster ids for each cell
  code_clustering <- mc[[numclusters]]$consensusClass#metaclusters consensus
  cell_clustering <- code_clustering[cell_clustering_som]#cell clustering from som
  return(list('code_clustering'=code_clustering,'cell_clustering'=cell_clustering,'metaclusters'=mc))
}


####CLUSTER HEATMAP FUNCTIONS ####
plot_clustering_heatmap_wrapper <- function(fcs, cell_clustering, nclusters=40,
                                            color_clusters='colorslusters', cluster_merging = NULL, 
                                            subtype_markers,
                                            clusterMergeFile=NULL,
                                            fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  #if((color_clusters)=='auto'){color_clusters <- hue_pal()(nclusters)}
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  ## Calculate the mean expression##################################################
  pdf(fileName, width=8, height=11) 
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  color_heat <- colorRampPalette(brewer.pal(n = 9, name = "YlOrBr"))(100)
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  labels_row <- paste0(expr01_mean$cell_clustering, " (", clustering_prop ,
                       "%)")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  ## Annotation for the merged clusters
  
  if(!is.null(clusterMergeFile)){
    ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Merged <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Merged <- color_clusters2
  }
  
  p <- pheatmap(expr_heat, color = rev(colorRampPalette(brewer.pal(n = 11, name = "RdYlBu"))(100)), 
                cluster_cols = T,
                cluster_rows = T, 
                labels_row = labels_row,
                #scale="column",
                display_numbers = FALSE, number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 8,
                cellheight = 8
  )
  dev.off() 
  print('Colors:')
  print(color_clusters)
  print(p);return(p)
}

plot_clustering_heatmap_wrapper2 <- function(fcs, cell_clustering, nclusters=40,
                                             color_clusters='clustercolors',
                                             subtype_markers,
                                             fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
 # if((color_clusters)=='auto'){color_clusters <- hue_pal()(nclusters)}
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  ## Calculate the mean expression##################################################
  pdf(fileName, width=8, height=11) 
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  #labels_row <- expr01_mean$cell_clustering
  
  labels_row <- paste0(expr01_mean$cell_clustering, " ")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  p <- pheatmap(expr_heat, 
                color = c(rep(magma(100)[1],25),magma(100)[1:100]), 
                cluster_cols = T,
                cluster_rows = F, 
                labels_row = labels_row,
                #scale="column",
                display_numbers = F, 
                number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 8,
                cellheight = 8,
                border_color = "black",
                annotation_legend = F
  )
  dev.off() 
  print('Colors:')
  print(color_clusters)
  print(p);return(p)
}

####DIAGNOSTICS####
makeDiagnosticPlots = function(exprData, 
                               md = output$meta_data,
                               sample_ids = output$sample_ids,
                               fcs = output$fcs,
                               subtype_markers = output$subtype_markers,
                               color_conditions = clustercolors,
                               shape_conditions = c(1:13),
                               fileName = 'diagnostics.pdf', 
                               tit = '', 
                               fun = mean)
{
  pdf(file = fileName)
  
  # plot 1
  ggdf <- data.frame(sample_id = sample_ids, exprData)
  ggdf <- melt(ggdf, id.var = 'sample_id', value.name = 'expression', 
               variable.name = 'antigen')
  mm <- match(ggdf$sample_id, md$sample_id)
  ggdf$timepoint <- md$timepoint[mm]
  print(ggplot(ggdf, aes(x = expression, color = timepoint, group = sample_id)) + 
          geom_density() +
          facet_wrap(~ antigen, nrow = 4, scales = 'free') + theme_bw() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                strip.text = element_text(size = 7),
                axis.text = element_text(size = 5)) + 
          scale_color_manual(values = color_timepoint) )
  # plot 2
  
  ## Spot check - number of cells per sample
  cell_table <- table(sample_ids)
  ggdf <- data.frame(sample_id = names(cell_table), 
                     cell_counts = as.numeric(cell_table))
  mm <- match(ggdf$sample_id, md$sample_id)
  ggdf$timepoint <- md$timepoint[mm]
  print(ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = timepoint)) + 
          geom_bar(stat = 'identity') + 
          geom_text(aes(label = cell_counts), hjust = 0.5, vjust = -0.5, size = 2.5) + 
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  
          scale_fill_manual(values = color_timepoint, drop = FALSE) + 
          scale_x_discrete(drop = FALSE))
  
  dev.off()
}


####CLUSTER HISTO####

plot_clustering_distr_wrapper <- function(expr = expr, 
                                          cell_clustering){
  # Calculate the median expression
  cell_clustering <- factor(cell_clustering)
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  # Calculate cluster frequencies
  freq_clust <- table(cell_clustering)
  freq_clust <- round(as.numeric(freq_clust)/sum(freq_clust)*100, 2)
  cell_clustering <- factor(cell_clustering,
                            labels = levels(cell_clustering))
  ### Data organized per cluster
  ggd <- melt(data.frame(cluster = cell_clustering, expr),
              id.vars = "cluster", value.name = "expression",
              variable.name = "antigen")
  ggd$antigen <- factor(ggd$antigen, levels = colnames(expr))
  ggd$reference <- "no"
  ### The reference data
  ggd_bg <- ggd
  ggd_bg$cluster <- "reference"
  ggd_bg$reference <- "yes"
  
  ggd_plot <- rbind(ggd, ggd_bg)
  ggd_plot$cluster <- factor(ggd_plot$cluster,
                             levels = c(levels(cell_clustering)[rev(cluster_rows$order)], "reference"))
  
  ggplot() +
    geom_density_ridges(data = ggd_plot, aes(x = expression, y = cluster,
                                             color = reference, fill = reference), alpha = 0.3) +
    facet_wrap( ~ antigen, scales = "free_x", nrow = 2) +
    theme_ridges() +
    theme(axis.text = element_text(size = 7),  
          strip.text = element_text(size = 7), legend.position = "none")
  
}
####UMAP####
do_umap <- function(fcs,subtype_markers,sample_ids,cell_clustering,metadata,
                    clusterMergeFile='~_merging.xlsx',
                    seed = 1234, ncells=2000,sample_subset=NULL){
  require(umap);require(flowCore);require(readxl)
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Create vector to later find and skip duplicates
  dups <- duplicated(expr[, subtype_markers])
  dups <- which(!(dups))## Find and skip duplicates
  ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
  ## New clustering1m
  mm <- match(cell_clustering, cluster_merging$original_cluster)
  cell_clustering1m <- cluster_merging$new_cluster[mm]
  ## Create a data frame of sample_ids and cell_clustering1m
  dtf<-data.frame(ids=sample_ids,type=cell_clustering1m)
  #dtf$B<- dtf$type!="B"#add a column that indicates whether the cell is a B cell or not; TRUE is non-B
  ##Why exclude B CELLS?
  ## WE HAVE NO B CELLS bc dtf$B depends on dtf$type depends on cellclustering1m which is just numbers 1:40 so..?
  #should it be the named parts cluster in merge file corresponding to it like if 30 -> grepl(cluster_merging[30,3],pattern='^B')??
  #Im blocking this out till we know why we have to do this
  ## Create subset columns without B cells (ONLY to generate the correct # of columns in inds2 object)
  #sample_ids2 <- dtf$ids[dtf$type!="B"] #sampleids without B cells
  #cell_clustering1m2 <- dtf$type[dtf$type!="B"] #clusters without B cells
  ## Data subsampling: create indices by sample
  inds <- split(1:length(sample_ids), sample_ids) #to get original indexes belonging to each cluster
  #inds2 <- split(1:length(sample_ids2), sample_ids2) #to fill in the original indexes that do not have B cells
  samplenames <- names(inds) #create a name vector of the files
  #FOR THIS TO WORK MUST BE IN FORMAT PRE/POST Tx
  # for (i in 1:(length(samplenames)/2)){#1:15 was here because ids in dtf was 30 factors and could be divided into 0 and 6wk for each so changed it
  #   templength <- length(inds2[[i]])
  #   inds2[[i]] <- inds[[i]][dtf$B[dtf$ids==samplenames[i]]] #subsets the "B cell or not a B cell" column for each sample by the name
  #   inds2[[i]] <- inds2[[i]][1:templength]
  # }
  custom.settings = umap.defaults
  custom.settings$seed = seed
  custom.settings$n.neighbors = neighbors
  ####umapindex generation####
  #umap ncells = table of sample ids with how many to downsample to by default col = id, row = ncells
  #sample ids = chr [1:2851129] "4927_0wk" "4927_0wk" "4927_0wk" "4927_0wk" ...
  #^ from ## Generate sample IDs corresponding to each cell in the 'expr' matrix sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))
  #can subset sample_ids and rerun umap 
  #can do timepoint or patient number or pre/post if you know corresp sample ids
  #sample_subset = '02_0wk' or c('02_0wk','2973_0wk') for example picking all 0 wk ones or use regex sample_ids[(grepl(sample_ids,pattern = '0wk'))]
  ifelse(is.null(sample_subset),
         umap_ncells <- pmin(table(sample_ids), ncells),
         umap_ncells <- pmin(table(sample_ids), ncells)[sample_subset]
  )
  if(!is.null(sample_subset)){inds <- inds[sample_subset]}
  umap_inds <- lapply(names(inds), function(i){
    s <- sample(inds[[i]], umap_ncells[i], replace = FALSE)
    intersect(s, dups)
  })
  set.seed(seed)
  umap_inds <- unlist(umap_inds)
  umap_out <- umap(expr[umap_inds, subtype_markers], config = custom.settings, method = 'naive')
  umapRes2D = data.frame(umap1 = umap_out$layout[, 1], umap2 = umap_out$layout[, 2], 
                         expr[umap_inds, subtype_markers],
                         sample_id = sample_ids[umap_inds], cell_clustering = factor(cell_clustering1m[umap_inds]), check.names = FALSE)
  
  #exclude any unassigned cluster post umap if needed--this has to be done by looking at the two columns 
  #metadata$sampleid is just a number in this metadatafile so to make unique ones combine samp_id col with timepoint (0wk)
  #to get in format of umapres2d$sample_id which looks like "02_0wk" do:
  return(umapRes2D)
}

plotUmap <- function(umapRes,seed=1234,neighbors=10,midpoint,color_clusters='clustercolors',code_clustering,subtype_markers=NULL)
{require(umap);require(ggplot2);require(viridis);require(ggrepel)
  #if((color_clusters)=='auto'){color_clusters <- hue_pal()(length(unique(code_clustering)))}
  custom.settings = umap.defaults
  custom.settings$seed = seed
  custom.settings$n.neighbors = neighbors
  ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = cell_clustering)) +
    geom_point(size = 1) +
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    
    scale_color_manual(values = color_clusters, name="CLUSTERS") +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  
  print(ggp)
  #other options
  print(ggp + facet_wrap(~ timepoint, ncol = 3)+ggtitle('TIMEPOINTS'))
  print(ggp + facet_wrap(~ batch, ncol = 2)+ggtitle('BATCH'))
  
  
  ggp2 <- ggplot(umapRes,  aes(x = umap1, y = umap2)) +
    geom_point(size = 1) +
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  print(ggp2 + facet_wrap(~ sample_id, ncol = 8)+ggtitle('SAMPLE'))
  
  
  
  
  ggp3 <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = sample_id)) +
    geom_point(size = 1) +
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  
  print(ggp3)
  #can specify which markers to display
  if(!is.null(subtype_markers)){
    for(i in subtype_markers)
    {
      ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = umapRes[,i])) +
        geom_point(size = 1) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        scale_color_gradient2(i, low="dark blue",mid="white",high="dark red", midpoint = mean(unlist(umapRes[,i])))
      print(ggp)
    }
  }
}





#======================
#     RUNNING DATA
#======================


####DATA LOADING AND CLUSTERING####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
workd<-getwd()

####if start with backup outputs
output<-readRDS('backup_output.rds')

#read and cluster

output <- returnfcs(metaDataFile = paste0(workd,"/Config","/metadata_M.xlsx"),
                    panelDataFile = paste0(workd,"/Config","/panel_M.xlsx"),
                    dataDirectory = paste0(workd,'/Data'))
output[8:10] <- clusterfcs(numclusters=40)
names(output)[8:10] <- c('code_clustering','cell_clustering','metaclusters')

saveRDS(output, file = "backup_output.rds")

#annotations

clusterMergeFile = paste0(workd,"/Config","/merged_M.xlsx")

cluster_merging <- read_xlsx(clusterMergeFile)


#set up factor levels

metadata = output$meta_data 
merged = cluster_merging

clusterlevels = unique(merged$new_cluster)

samplevels <- c(metadata$sample_id)

timelevels=unique(metadata$timepoint)

clustercolors <- as.character(c(cols25(n=25),alphabet(n=19)))

mm1 <- match(output$cell_clustering, cluster_merging$original_cluster)

cell_clustering1m <- cluster_merging$new_cluster[mm1]

output$cell_clustering1m <- cell_clustering1m

#metacluster heatmap
plot_clustering_heatmap_wrapper(fcs=output$fcs,
                                color_clusters = c(stepped2(20),stepped3(20),rev(cubicl(20))),
                                cell_clustering = output$cell_clustering, 
                                subtype_markers=output$subtype_markers,
                                clusterMergeFile = clusterMergeFile,
                                fileName = 'M_clusteringheatmap.pdf');dev.off()

plot_clustering_heatmap_wrapper2(fcs=output$fcs,
                                 color_clusters = clustercolors,
                                 cell_clustering = factor(output$cell_clustering1m, levels=clusterlevels), 
                                 subtype_markers=output$subtype_markers,
                                 fileName = 'J19113_clusteringheatmap_merged.pdf');dev.off()

#save output
saveRDS(output, file = "backup_output.rds")


#diagnostic plots

color_timepoint <-c("red","blue","purple")

expr <- fsApply(output$fcs, exprs)
rng <- colQuantiles(expr, probs = c(0.01, 0.99))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1])) #scaling 0-1
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1

gc()
dgplots<-makeDiagnosticPlots(expr01, 
                             md = output$meta_data,
                             sample_ids = output$sample_ids,
                             fcs = output$fcs,
                             subtype_markers = output$subtype_markers,
                             color_conditions = color_timepoint,
                             shape_conditions = c(0:25,0,1),
                             fileName = 'J19113_plot_diagnostics.pdf')

#for mds
expr_mean_sample_tbl <- data.frame(sample_id = output$sample_ids, expr01) %>%
  group_by(sample_id) %>%  summarize_all(funs(mean))

expr_mean_sample <- t(expr_mean_sample_tbl[, -1])
colnames(expr_mean_sample) <- expr_mean_sample_tbl$sample_id

mds <- plotMDS(expr_mean_sample, plot = FALSE)
ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
                   sample_id = colnames(expr_mean_sample))
mm <- match(ggdf$sample_id, output$meta_data$sample_id)
ggdf$timepoint <- output$meta_data$timepoint[mm]
ggdf$batch <- output$meta_data$batch[mm]

mdsggp <- ggplot(ggdf, aes(x = MDS1, y = MDS2, color = timepoint, shape = batch)) +
  geom_point(size = 2, alpha = 0.8) +
  #geom_label_repel(aes(label = sample_id)) +
  theme_bw() +
  theme(axis.text = element_text(color="black")) +
  #ylim(-0.1, 0.08) +
  scale_color_manual(values = c("blue","orange","purple","black")) +
  scale_shape_manual(values = c(2,0,1,3,4,5,6,7,8)) +
  coord_fixed()
pdf("J19113_plot_mds.pdf", width=4, height=5);mdsggp;dev.off()

#for mean heatmap
mm <- match(colnames(expr_mean_sample[,samplevels]), output$meta_data$sample_id)
colnames(expr_mean_sample[,samplevels])
annotation_col <- data.frame(timepoint = factor(output$meta_data$timepoint[mm]), 
                             row.names = colnames(expr_mean_sample[,samplevels]))

annotation_colors <- list(condition = clustercolors[1:length(levels(annotation_col$timepoint))])

color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)

meanht<- pheatmap(expr_mean_sample[,samplevels], color = color, display_numbers = TRUE,
                  number_color = "black", fontsize_number = 5, 
                  cluster_cols = F,
                  annotation_col = annotation_col, main = '',
                  annotation_colors = annotation_colors, clustering_method = "average")
dev.off()
pdf("J19113_plot_meanheatmap.pdf", width=12, height=5);meanht;dev.off()




####DIFFERENTIAL PLOTS####

#set up count and prop matrices
counts_table <- table(output$cell_clustering1m, output$sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

write.csv(counts, file='J19113_counts.csv')
write.csv(props, file='J19113_props.csv')

#set up the data frame for proportional plotting
ggdf <- melt(data.frame(cluster = rownames(props),props, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdf$cluster <- factor(ggdf$cluster, levels=clusterlevels)
ggdf$timepoint <- factor(output$meta_data$timepoint[match(ggdf$sample_id,output$meta_data$sample_id)], levels=timelevels)
ggdf$patient_id <- factor(output$meta_data$patient_id[match(ggdf$sample_id,output$meta_data$sample_id)])
ggdf$batch <- output$meta_data$batch[match(ggdf$sample_id,output$meta_data$sample_id)]
ggdf$recist <- output$meta_data$recist[match(ggdf$sample_id,output$meta_data$sample_id)]
ggdf$benefit <- output$meta_data$benefit[match(ggdf$sample_id,output$meta_data$sample_id)]


#plot box plots
boxplot <- ggplot(ggdf, aes(x=timepoint, y=proportion, fill=timepoint))+
  geom_boxplot(outlier.size=0, lwd=0.25, color="black")+
  #geom_jitter(width=0, aes(shape=sample_id))+
  facet_wrap(~cluster,ncol=7,scales="free")+
  ylab("% of Live Cells")+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=7.5, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white")  
        )
pdf("J19113_plot_abundance_box.pdf", width=10, height=8);boxplot;dev.off()

#line plot for all timepoints
##Abundance plot - Day 1 vs Day 8
ggdf_gvax<-ggdf[ggdf$timepoint %in% c("C1D1", "C1D8"),]
ggdf_gvax$type_tp <- paste(ggdf_gvax$recist, ggdf_gvax$timepoint, sep="_")
ggdf_gvax$type_tp<- factor(ggdf_gvax$type_tp, levels=unique(ggdf_gvax$type_tp))

ggdf_gvax$type_tp<- factor(ggdf_gvax$type_tp, levels=c("NA_C1D1","NA_C1D8","NA_C1D15","PD_C1D1","PD_C1D8","PD_C1D15","SD_C1D1","SD_C1D8","SD_C1D15"))

lineplot_r <- ggplot(ggdf_gvax, aes(x=timepoint, y=proportion, group=patient_id))+
  facet_wrap(~cluster, scales='free',ncol=4)+
  scale_shape_manual(values=c(0:25,0,1))+
  geom_point(aes(shape=patient_id,col = recist))+
  geom_line(show.legend = T, aes(col = recist,linetype=recist))+
  theme(axis.text.x = element_text(size=8, angle = 45, hjust=1, vjust=1, color="black"),
        axis.title.x = element_blank(),
        axis.ticks=element_line(color="black"),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=10),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(1,'lines'),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill="white"))
dev.off()
pdf("M_abundance_C1D1 vs C1D8.pdf", width=11, height=16);lineplot_r;dev.off()

##Abundance plot - Day 8 vs Day 15
ggdf_gvax<-ggdf[ggdf$timepoint %in% c("C1D8", "C1D15"),]
ggdf_gvax$type_tp <- paste(ggdf_gvax$recist, ggdf_gvax$timepoint, sep="_")
ggdf_gvax$type_tp<- factor(ggdf_gvax$type_tp, levels=unique(ggdf_gvax$type_tp))

ggdf_gvax$type_tp<- factor(ggdf_gvax$type_tp, levels=c("NA_C1D1","NA_C1D8","NA_C1D15","PD_C1D1","PD_C1D8","PD_C1D15","SD_C1D1","SD_C1D8","SD_C1D15"))

lineplot_r <- ggplot(ggdf_gvax, aes(x=type_tp, y=proportion, group=patient_id))+
  facet_wrap(~cluster, scales='free',ncol=4)+
  scale_shape_manual(values=c(0:25,0,1))+
  geom_point(aes(shape=patient_id,col = recist))+
  geom_line(show.legend = T, aes(col = recist,linetype=recist))+
  theme(axis.text.x = element_text(size=8, angle = 45, hjust=1, vjust=1, color="black"),
        axis.title.x = element_blank(),
        axis.ticks=element_line(color="black"),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=10),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(1,'lines'),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill="white"))
dev.off()
pdf("M_abundance_C1D8 vs C1D15.pdf", width=11, height=16);lineplot_r;dev.off()

##Abundance plot - Day 1 vs Day 15
ggdf_gvax<-ggdf[ggdf$timepoint %in% c("C1D1", "C1D15"),]
ggdf_gvax$type_tp <- paste(ggdf_gvax$recist, ggdf_gvax$timepoint, sep="_")
ggdf_gvax$type_tp<- factor(ggdf_gvax$type_tp, levels=unique(ggdf_gvax$type_tp))

ggdf_gvax$type_tp<- factor(ggdf_gvax$type_tp, levels=c("NA_C1D1","NA_C1D8","NA_C1D15","PD_C1D1","PD_C1D8","PD_C1D15","SD_C1D1","SD_C1D8","SD_C1D15"))

lineplot_r <- ggplot(ggdf_gvax, aes(x=type_tp, y=proportion, group=patient_id))+
  facet_wrap(~cluster, scales='free',ncol=4)+
  scale_shape_manual(values=c(0:25,0,1))+
  geom_point(aes(shape=patient_id,col = recist))+
  geom_line(show.legend = T, aes(col = recist,linetype=recist))+
  theme(axis.text.x = element_text(size=8, angle = 45, hjust=1, vjust=1, color="black"),
        axis.title.x = element_blank(),
        axis.ticks=element_line(color="black"),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=10),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(1,'lines'),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill="white"))
dev.off()
pdf("M_abundance_C1D1 vs C1D15.pdf", width=11, height=16);lineplot_r;dev.off()

##Abundance plot - Day 1 vs Day 8 vs Day 15
ggdf_gvax<-ggdf[ggdf$timepoint %in% c("C1D1","C1D8", "C1D15"),]
ggdf_gvax$type_tp <- paste(ggdf_gvax$recist, ggdf_gvax$timepoint, sep="_")
ggdf_gvax$type_tp<- factor(ggdf_gvax$type_tp, levels=unique(ggdf_gvax$type_tp))

ggdf_gvax$type_tp<- factor(ggdf_gvax$type_tp, levels=c("NA_C1D1","NA_C1D8","NA_C1D15","PD_C1D1","PD_C1D8","PD_C1D15","SD_C1D1","SD_C1D8","SD_C1D15"))

lineplot_r <- ggplot(ggdf_gvax, aes(x=type_tp, y=proportion, group=patient_id))+
  facet_wrap(~cluster, scales='free',ncol=4)+
  scale_shape_manual(values=c(0:25,0,1))+
  geom_point(aes(shape=patient_id,col = recist))+
  geom_line(show.legend = T, aes(col = recist,linetype=recist))+
  theme(axis.text.x = element_text(size=8, angle = 45, hjust=1, vjust=1, color="black"),
        axis.title.x = element_blank(),
        axis.ticks=element_line(color="black"),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.text.y = element_text(color="black"),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=10),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(1,'lines'),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill="white"))
dev.off()
pdf("M_abundance_C1D1 vs C1D8 vs C1D15.pdf", width=14, height=20);lineplot_r;dev.off()


## Umap
umapRes <- do_umap(fcs=output$fcs,subtype_markers = output$subtype_markers,
                   sample_ids = output$sample_ids,cell_clustering = output$cell_clustering, metadata=output$metadata,
                   clusterMergeFile=clusterMergeFile,
                   seed = 1234, ncells=500,sample_subset=NULL)

mm <- match(as.character(umapRes$sample_id), as.character(output[["meta_data"]]$sample_id))
umapRes$timepoint <- factor(output[["meta_data"]]$timepoint[mm], levels=timelevels)
umapRes$batch <- output[["meta_data"]]$batch[mm]
umapRes$sample_id <- factor(output[["meta_data"]]$sample_id[mm], levels=samplevels)
umapRes$cell_clustering = factor(umapRes$cell_clustering, levels=clusterlevels)

dev.off()
pdf('J19113_plot_umaps.pdf',width=10,height=10)
plotUmap(umapRes = umapRes,
         code_clustering=cell_clustering1m,
         color_clusters = clustercolors,
         subtype_markers = output$subtype_markers)
dev.off()

saveRDS(umapRes, file="backup_umap.rds")


## Functional marker expression line plots

exprmntbl <- data.frame(fsApply(output$fcs,exprs)[, c(output$subtype_markers)],
                        sample_id = output$sample_ids, 
                        cluster = output$cell_clustering1m) %>%
  group_by(sample_id, cluster) %>%
  summarize_all(funs(mean))
exprmntbl$cluster <- factor(exprmntbl$cluster, levels=clusterlevels)
exprmntbl$timepoint <- factor(output$meta_data$timepoint[match(exprmntbl$sample_id,output$meta_data$sample_id)], levels=timelevels)
exprmntbl$patient_id <- factor(output$meta_data$patient_id[match(exprmntbl$sample_id,output$meta_data$sample_id)])
exprmntbl$batch <- output$meta_data$batch[match(exprmntbl$sample_id,output$meta_data$sample_id)]
exprmntbl$recist <- output$meta_data$recist[match(exprmntbl$sample_id,output$meta_data$sample_id)]
exprmntbl$benefit <- output$meta_data$benefit[match(exprmntbl$sample_id,output$meta_data$sample_id)]
ggdf2<-melt(exprmntbl, id.var=c("cluster","sample_id","patient_id","benefit","recist","timepoint"))
ggdf2$value <- as.numeric(ggdf2$value)
fmlistplot <- c(output$subtype_markers)

##Functional Plot - Day 1 vs D8 vs D15
ggdf_gvax<-ggdf[ggdf$timepoint %in% c("C1D1", "C1D8","C1D15"),]
ggdf2_gvax<-ggdf2[ggdf2$patient_id %nin% c("P8"),]  ##Excluded due to missing data point 
ggdf2_gvax$timepoint<- factor(ggdf2_gvax$timepoint, levels=c("C1D1","C1D8","C1D15"))
comps<- list(c("C1D1","C1D8"))
comps2<- list(c("C1D8","C1D15"))

pdf("M_functional_C1D1 vs C1D8 vs C1D15.pdf",width=10,height=16)
for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf2_gvax[ggdf2_gvax$variable==fmlistplot[i],], aes(x=timepoint,y=value, group=patient_id))+
    ggtitle(fmlistplot[i])+
    facet_wrap(~cluster, scales='free',ncol=4)+
    scale_shape_manual(values=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,0,1,2,3,4,5))+
    geom_point(aes(shape=patient_id, col=recist))+
    geom_line(show.legend = T, aes(col=recist, linetype=recist))+
    stat_compare_means(method="t.test",show.legend=FALSE, paired=TRUE, label="p.format", size=2, comparisons=comps)+
    stat_compare_means(method="t.test",show.legend=FALSE, paired=TRUE, label="p.format", size=2, comparisons=comps2)+
    theme(axis.text.x = element_text(size=8, angle = 45, hjust=1, vjust=1, color="black"),
          axis.title.x = element_blank(),
          axis.ticks=element_line(color="black"),
          axis.line.x = element_line(size=0.5),
          axis.line.y = element_line(size=0.5),
          axis.text.y = element_text(color="black"),
          axis.title.y = element_blank(),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=10),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.key.size = unit(1,'lines'),
          legend.text = element_text(size=10),
          legend.key = element_rect(fill="white"))
  print(ggp)
}
dev.off()

##Functional Plot - D8 vs D15
ggdf2_gvax<-ggdf2[ggdf2$timepoint %in% c("C1D8", "C1D15")]
ggdf2_gvax$timepoint<- factor(ggdf2_gvax$timepoint, levels=c("C1D1","C1D8","C1D15"))
                              
comps<- list(c("C1D8", "C1D15"))
pdf("M_functional_C1D8 vs C1D15.pdf",width=14,height=20)
for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf2_gvax[ggdf2_gvax$variable==fmlistplot[i],], aes(x=timepoint,y=value, group=patient_id))+
    ggtitle(fmlistplot[i])+
    facet_wrap(~cluster, scales='free',ncol=4)+
    scale_shape_manual(values=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,0,1,2,3,4,5))+
    geom_point(aes(shape=patient_id, col=recist))+
    stat_compare_means(method="t.test",show.legend=FALSE, paired=TRUE, label="p.format", size=2, comparisons=comps)+
    geom_line(show.legend = T, aes(col=recist, linetype=recist))+
    theme(axis.text.x = element_text(size=8, angle = 45, hjust=1, vjust=1, color="black"),
          axis.title.x = element_blank(),
          axis.ticks=element_line(color="black"),
          axis.line.x = element_line(size=0.5),
          axis.line.y = element_line(size=0.5),
          axis.text.y = element_text(color="black"),
          axis.title.y = element_blank(),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=10),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.key.size = unit(1,'lines'),
          legend.text = element_text(size=10),
          legend.key = element_rect(fill="white"))
  print(ggp)
}
dev.off()

##Funcitonal Plot - D1 vs D15
ggdf2_gvax<-ggdf2[ggdf2$timepoint %in% c("C1D1", "C1D15")]
ggdf2_gvax$timepoint<- factor(ggdf2_gvax$timepoint, levels=c("C1D1","C1D8","C1D15"))

comps<- list(c("C1D1", "C1D15"))
pdf("M_functional_C1D1 vs C1D15.pdf",width=14,height=20)
for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf2_gvax[ggdf2_gvax$variable==fmlistplot[i],], aes(x=timepoint,y=value, group=patient_id))+
    ggtitle(fmlistplot[i])+
    facet_wrap(~cluster, scales='free',ncol=4)+
    scale_shape_manual(values=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,0,1,2,3,4,5))+
    geom_point(aes(shape=patient_id, col=recist))+
    stat_compare_means(method="t.test",show.legend=FALSE, paired=TRUE, label="p.format", size=2, comparisons=comps)+
    geom_line(show.legend = T, aes(col=recist, linetype=recist))+
    theme(axis.text.x = element_text(size=8, angle = 45, hjust=1, vjust=1, color="black"),
          axis.title.x = element_blank(),
          axis.ticks=element_line(color="black"),
          axis.line.x = element_line(size=0.5),
          axis.line.y = element_line(size=0.5),
          axis.text.y = element_text(color="black"),
          axis.title.y = element_blank(),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=10),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.key.size = unit(1,'lines'),
          legend.text = element_text(size=10),
          legend.key = element_rect(fill="white"))
  print(ggp)
}
dev.off()


################################################## edge R

View(counts)
View(ggdf)
View(ggdf2)
View(exprmntbl)

library(edgeR)
library(dplyr)

################single cell level
marker = c()
FDR_adjusted = vector()
clusternames = vector()
comparison = vector()
edgeR_table = data.frame(logFC=NULL,	logCPM=NULL,	PValue=NULL)

for (i in 1:length(clusterlevels)) {
  celltable = exprmntbl[!exprmntbl$timepoint=="C1D1",]
  g1 = celltable[celltable$cluster==clusterlevels[i],]
  g = t(g1)
  ids = unique(g1$sample_id)
  timepoints = c(g1$timepoint)
  g = g[-c(1,2,39:43),]
  class(g) = "numeric"
  d <- DGEList(counts=g,group=factor(timepoints))
  d1 <- estimateCommonDisp(d, verbose=T)
  d1 <- estimateTagwiseDisp(d1)
  
  
  table1cluster = data.frame(logFC=NULL,	logCPM=NULL,	PValue=NULL)
  for (j in 1:(2-1)) {
    for (k in (j+1):2) {
      
      et <- exactTest(d1, pair=c(timepoints[j],timepoints[k]))
      fdr = p.adjust(et$table$PValue,method = "BH")
      FDR_adjusted = c(FDR_adjusted,fdr)
      table1cluster = rbind(table1cluster,et$table)
      edgeR_table = rbind(edgeR_table,et$table)
      
      group = paste(timepoints[j],' and ',timepoints[k])
      marker = c(marker,row.names(g))
      
      for (z in 1:length(fmlistplot)) {
        comparison = c(comparison,group)
    }
  }
  
  for (p in 1:nrow(table1cluster)) {
    clusternames = c(clusternames,clusterlevels[i])
    }
  }
}

final_table = data.frame(cluster = clusternames,
                         comparison = comparison,
                         marker = marker,
                         logFC = edgeR_table$logFC,
                         logCPM = edgeR_table$logCPM,
                         PValue = edgeR_table$PValue,
                         FDR = FDR_adjusted)

write.csv(final_table,paste0("singleCell edgeR for C1D8vsC1D15.csv"))



############abundance level

abundance_table = counts
ids = colnames(counts)
tp = data.frame(timepoint = NULL)
for (i in 1:ncol(counts)) {
  tp = rbind(tp,as.character(unique(ggdf[which(ggdf$sample_id==ids[i]),4])))
}
tp = t(tp)
colnames(tp) = colnames(abundance_table)
abundance_table = rbind(abundance_table,tp)
colnames(abundance_table) = tp

one8_table = abundance_table[,!colnames(abundance_table)=="C1D15"]
eight15_table = abundance_table[,!colnames(abundance_table)=="C1D1"]
groups = c(as.character(eight15_table[21,]))
timepointlevels = unique(groups)

one8_table = one8_table[-21,]
one8_table[,] = apply(one8_table[,c(1:ncol(one8_table))],2,function(x) as.numeric(as.character(x)))
one8_table1 = data.matrix(one8_table)

eight15_table = eight15_table[-21,]
eight15_table[,] = apply(eight15_table[,c(1:ncol(eight15_table))],2,function(x) as.numeric(as.character(x)))
eight15_table1 = data.matrix(eight15_table)


d = DGEList(counts=eight15_table1,group=groups)
d0 <- estimateCommonDisp(d, verbose=T)
d0 <- estimateTagwiseDisp(d0)

finaltable = data.frame(Cluster=NULL,Comparison=NULL,logFC=NULL,logCPM=NULL,PValue=NULL,FDR=NULL)
mergedEDGER = data.frame(logFC=NULL,	logCPM=NULL,	PValue=NULL)
FDR = c()
edgeR_table = data.frame(logFC=NULL,	logCPM=NULL,	PValue=NULL)
fourCluster = data.frame(logFC=NULL,	logCPM=NULL,	PValue=NULL)
clusternames = c()
compare = c()
funcmarker = c()
FDR_adjusted = c()

for (x in 1:(length(timepointlevels)-1)) {
  for (y in (x+1):length(timepointlevels)) {
    et <- exactTest(d0, pair=c(timepointlevels[x],timepointlevels[y]))
    fdr = p.adjust(et$table$PValue,method = "BH")
    FDR_adjusted = c(FDR_adjusted,fdr)
    fourCluster = rbind(fourCluster,et$table)
    mergedEDGER = rbind(mergedEDGER,et$table)
    
    clusternames = c(clusternames,rownames(eight15_table1))
    
    for (z in 1:nrow(et$table)) {
      compa = paste(timepointlevels[x],'and',timepointlevels[y])
      compare = c(compare,compa)
    }
    
  }
  
}

finaltable = data.frame(Cluster=clusternames,
                        Comparison=compare,
                        logFC=mergedEDGER$logFC,
                        logCPM=mergedEDGER$logCPM,
                        PValue=mergedEDGER$PValue,
                        FDR=FDR_adjusted)

write.csv(finaltable,paste0("edgeR_counts_C1D8vsC1D15.csv"))
