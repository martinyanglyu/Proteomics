
library(dplyr)
library(ggplot2)
library(yaml)
library(fgsea)
library(dplyr)
library(tibble)
#library(rJava) # need to insall for pathvisio
#library("XMLRPC") # need to insall for pathvisio

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



pheatmap_pathway_NES<- function (expression_data,pathway_name) {
  Color = colorRampPalette(c("green", "white", "red"))(50)
  #  paletteLength <- 50
  # NES_breaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
  #            seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
  #NES_breaks <- c(seq(-3, 0, length.out=ceiling(paletteLength/2) + 1),
  # seq(3/paletteLength, 3, length.out=floor(paletteLength/2)))
  NES_breaks = seq(-3, 3, length.out =50)
  
  pheatmap::pheatmap(expression_data, cellwidth=20, cellheight = 10, scale = "none",color = Color,
                     treeheight_row = 5,main=pathway_name,
                     show_rownames = T,show_colnames = T,
                     #annotation_col= Annotation,
                     # annotation_row=Annotation,
                     #annotation_legend=Label_def,
                     cluster_rows = T,  cluster_cols = F,clustering_distance_rows = "euclidean",breaks=NES_breaks)
  
}

matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[,1]
  m
}


load_config <- function(config_path = "config.yaml") {
  if (!file.exists(config_path)) {
    stop("Configuration file not found at: ", config_path)
  }
  
  config <- yaml::read_yaml(config_path)
  
  # Create directories if they don't exist
  dirs_to_create <- c(
    config$default$data$raw_data,
    config$default$data$processed_data,
    config$default$output$figures,
    config$default$output$tables,
    config$default$output$results
  )
  
  lapply(dirs_to_create, function(dir) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  })
  
  return(config$default)
}

UMAP_legend_4=function(data,pathway_name,x_name,fill_name,color){
  my_hist <-ggplot(data,aes_string(
    x =x_name,
    y = pathway_name,
    fill=fill_name
  ) )+  geom_violin( )+scale_fill_manual(name = fill_name,values=color)
  
  legend <- cowplot::get_legend(my_hist)
  
}

Valcano_plot_NES=function(Data,Name){
   
  
  library(ggrepel)
  # add a column of NAs
  de=Data %>%mutate(diffexpressed=ifelse(pvalue>0.05,"NO",ifelse(log2FoldChange>0,"UP","DOWN")))

  de$delabel <- NA
 
  de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]
   # plot adding up all layers we have seen so far
 
  ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
    geom_point() + ggtitle(Name)+
    geom_text_repel(size=6) +  theme_minimal() +
    scale_color_manual(values=c("DOWN"="blue", "NO"= "black","UP"= "red")) +
    geom_vline(xintercept=c(-0.14, 0.14), col=c("blue" ,"red")) +
    geom_hline(yintercept=-log10(0.05), col="green")+ theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=4, fill=NA))
  
  
  
  
}

pheatmap_singscore<- function (pathways,data,Annotation) {
  Gene_select_anno= data[,colnames(data) %in% pathways] %>%t()%>%na.omit()%>%t()%>%as.data.frame()%>%.[,rownames(Annotation)]

  pheatmap::pheatmap(Gene_select_anno, cellwigermline=5, cellheight = 10,cellwidth = 10, scale = "row",
                                    treeheight_row = 5,
                                    show_rownames = T,show_colnames = F,
                                    annotation_col= Annotation,
                                    # annotation_row=Annotation,
                                    #annotation_legend=Label_def,
                                    cluster_rows = T,  cluster_cols = F,clustering_distance_rows = "euclidean")
  
}

Valcano_plot_all_proteins=function(Data,Title){
   
  
  library(ggrepel)
  # add a column of NAs
  de=Data %>%mutate(diffexpressed=ifelse(pvalue>0.05,"NO",ifelse(log2FoldChange>0,"UP","DOWN")))

  de$delabel <- NA
 
  de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]
   # plot adding up all layers we have seen so far
 
  ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
    geom_point() + ggtitle(Title)+
    theme_minimal() +
    geom_text_repel() +
    scale_color_manual(values=c("DOWN"="blue", "NO"= "black","UP"= "red")) +
    geom_vline(xintercept=c(-0.26, 0.26), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")
  
  
  
  
}

Singscore_get=function(Normalized_data,pathway_list) {
    rankdata<-singscore::rankGenes(Normalized_data)
  # return(rankdata)
  # rownames Should be the names of genes
  out=list()
  for (i in 1: length(pathway_list)) {
     upgenelist=pathway_list[[i]]
     pathway_name=names(pathway_list)[i]
     scoredf <- singscore::simpleScore(rankdata, upSet=upgenelist)
     if( nrow(scoredf)==ncol(Normalized_data)&!is.null(scoredf)) {
       Output= scoredf
       colnames(Output)=paste(pathway_name,colnames(Output),sep="_")
       out[[pathway_name]]=Output[,1]  }
   
  }
  combine=do.call("cbind",out)
}


Wilcox_test_singscore<- function (Data,Group1_ID,Group2_ID) {
  plist=list()
  plist2=list()
  for (Gene in rownames(Data)) {
    Data_select=Data[rownames(Data)==Gene,]
    a=Data_select[,Group1_ID]%>%as.numeric()
    b=Data_select[,Group2_ID]%>%as.numeric()
    
    Test=wilcox.test(b,a, correct=FALSE)
    pvalue=Test$p.value
    plist[[Gene]]=pvalue
    Mean_a=mean(a)
    Mean_b=mean(b)
    if (Mean_a >0 & Mean_b >0) {
      log2FoldChange=log2(mean(b)/mean(a))
      
    }
    
    if (Mean_a <0 & Mean_b <0) {
      log2FoldChange=log2( 1/(mean(b)/mean(a)))
      
    }
    
    if (Mean_a <0 & Mean_b >0) {
      log2FoldChange=log2( (mean(b)-mean(a)/abs(mean(a))))
      
    }
    
    if (Mean_b <0 & Mean_a >0) {
      log2FoldChange=log2(1/ (mean(a)-mean(b)/abs(mean(b))))
      
    }
     plist2[[Gene]]=data.frame(Mean_a=mean(a),Mean_b=mean(b),log2FoldChange=log2FoldChange)
  }
  # return(plist2)
  test=plist%>%as.data.frame%>%t()%>%as.data.frame()%>%setNames(c("pvalue")) %>%mutate (padj=p.adjust(.$pvalue, method = "fdr"))
  test2=do.call(rbind,plist2)
  output=cbind(test,test2)%>%tibble::rownames_to_column(.,var="SYMBOL")

}

 
Valcano_plot_new_2=function(Data,Pathway_list=NULL,pathway_name=NULL){
  
  if (!is.null(Pathway_list)) {
    pathway_genes=Pathway_list [[pathway_name]]%>%as.character()
    
    library(ggrepel)
    Data=Data%>%.[.$gene_symbol%in%pathway_genes,]
    
  }  
  # if colnames (Data) is not log2Foldchange, then change it
  if (!"log2Foldchange" %in% colnames(Data)) {
    colnames(Data)[colnames(Data)=="log2FoldChange"] <- "log2Foldchange"
  }
  # add a column of NAs
  de=Data %>%mutate(diffexpressed=ifelse(pvalue>0.05,"NO",ifelse(log2Foldchange>0,"UP","DOWN")))
  de$delabel <- NA
  de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]
  
  # plot adding up all layers we have seen so far
  ggplot(data=de, aes(x=log2Foldchange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
    geom_point() + ggtitle(pathway_name)+
    theme_minimal() +
    ggrepel:: geom_text_repel(max.overlaps = Inf ) +
    scale_color_manual(values=c("DOWN"="blue", "NO"= "black","UP"= "red")) +
    geom_vline(xintercept=c(-0.6, 0.6), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")
  
  
  
  
}

Valcano_plot_FC_cut=function(Data,Pathway_list=NULL,pathway_name=NULL,FC){
  
  if (!is.null(Pathway_list)) {
    pathway_genes=Pathway_list [[pathway_name]]%>%as.character()
    
    library(ggrepel)
    Data=Data%>%.[.$gene_symbol%in%pathway_genes,]
    
  }  
  
  # add a column of NAs
  FC=log2(FC)
  de=Data %>%mutate(diffexpressed=ifelse(pvalue>0.05|abs(log2Foldchange)<FC,"NO",ifelse(log2Foldchange>0,"UP","DOWN")))
  de$delabel <- NA
  de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]
  top_genes_p <- de %>%
    dplyr::filter(diffexpressed != "NO") %>%
    arrange(pvalue) %>%
    head(10)
  
  top_upgenes <- de %>%
    dplyr::filter(diffexpressed != "NO") %>%dplyr::filter(log2Foldchange>0)%>%
    arrange(desc(log2Foldchange)) %>%
    head(10)
  
  top_dngenes <- de %>%
    dplyr::filter(diffexpressed != "NO") %>%dplyr::filter(log2Foldchange<0)%>%
     arrange(log2Foldchange) %>%
    head(10)
  
  top_genes <- dplyr::bind_rows(top_genes_p, top_upgenes, top_dngenes) %>%
    distinct() 
  
  # plot adding up all layers we have seen so far
  ggplot(data=de, aes(x=log2Foldchange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
    geom_point() + ggtitle(pathway_name)+
    theme_minimal() +
    ggrepel::geom_text_repel(data = top_genes, aes(label = gene_symbol), 
                             max.overlaps = Inf, box.padding = 0.35, 
                             point.padding = 0.3, segment.color = "grey50") +
    scale_color_manual(values=c("DOWN"="blue", "NO"= "black","UP"= "red")) +
    geom_vline(xintercept=c(-FC, FC), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")
  
  
  
  
}

UMAP_output=function(Data_select ) {
  set.seed(1234)
  # kmeans_1_data=read.table(file.path( Rawdata_path,"Proteomics","54_proteomics_singscore_kmeans_1.txt"))
  
  
  Data.umap = umap::umap(Data_select)
  umap_plot_df <- data.frame(Data.umap$layout)%>%setNames(c("UMAP1","UMAP2"))
 
  umap_plot_df_2=umap_plot_df%>%tibble::rownames_to_column(.,var="Sample_ID")
  Data_select2=Data_select%>%mutate(Sample_ID=rownames(.))
  umap_plot_df_2=umap_plot_df%>%tibble::rownames_to_column(.,var="Sample_ID")%>%left_join(.,Data_select2,by="Sample_ID")
  
}



UMAP_plot=function(Data_umap,shape_label,color_label,color_defination) {
  set.seed(1234)
  if(shape_label=="") {
    ggplot(
      Data_umap,
      aes_string(
        x = "UMAP1",
        y = "UMAP2",
        color=color_label
      )
    ) +  geom_point( size = 6)+scale_colour_manual(values=color_defination) +theme(
    panel.background = element_blank(),   # Remove the panel background
    panel.grid.major = element_blank(),   # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    axis.line = element_line(color = "black") # Optionally keep axis lines
  )
    
    
  } else {
    ggplot(
      Data_umap,
      aes_string(
        x = "UMAP1",
        y = "UMAP2",
        color=color_label,
        shape=shape_label# label points with different colors for each `subgroup`
      )
    ) +  geom_point( size = 6)+scale_colour_manual(values=color_defination) +theme(
    panel.background = element_blank(),   # Remove the panel background
    panel.grid.major = element_blank(),   # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    axis.line = element_line(color = "black") # Optionally keep axis lines
  )
  }
  
  
}


get_diff_expressed_genes <- function(data, FC_cutoff = 1.2, pvalue_cutoff = 0.05) {
  FC_log2_cutoff <- log2(FC_cutoff)
  
  # Apply T-test and filter genes based on FC and p-value cutoffs
  diff_genes <- data %>%
    dplyr::filter(pvalue < pvalue_cutoff & abs(log2Foldchange) >= FC_log2_cutoff) %>%
    mutate(direction = ifelse(log2Foldchange > 0, "UP", "DOWN"))
  
  # Split the results into upregulated and downregulated genes
  upregulated_genes <- diff_genes %>%
    dplyr::filter(direction == "UP") %>%
    pull(gene_symbol)  # Get the gene names
  
  downregulated_genes <- diff_genes %>%
    dplyr::filter(direction == "DOWN") %>%
    pull(gene_symbol)  # Get the gene names
  
  return(list(upregulated_genes = upregulated_genes, downregulated_genes = downregulated_genes))
}



perform_fgsea_analysis <- function(Ttest_list, GMT_files, stat_value="pvalue",min_size = 5, max_size = 500, n_perms = 10000) {
  
  # Get GMT files from the specified directory
  # GMT_files <- list.files(file.path(GMT_directory), full.names = TRUE) %>%
  #   .[grepl("gmt.txt", .)] %>%
  #   .[grepl("c2|h.", .)]
 

  
  # Create an empty list to store FGSEA output
  Fgsea_output <- list()
  
  # Loop through each T-test result in the list
  for (i in names(Ttest_list)) {
    
    # Filter the T-test results for significant genes (p-value < 0.05)
    res_p005 <- Ttest_list[[i]] %>%
      dplyr::filter(.data[[stat_value]] < 0.05)
  
    # Prepare the ranks based on log2Foldchange, and remove NA values
    # if the log2FoldChange, correct it to log2Foldchange
    if (!"log2Foldchange" %in% colnames(res_p005)) {
      colnames(res_p005)[colnames(res_p005) == "log2FoldChange"] <- "log2Foldchange"
    }
    
    res_p005_rank <- res_p005 %>%
      arrange(log2Foldchange) %>%
      dplyr::select(log2Foldchange) %>%
      tibble::rownames_to_column(var = "gene_symbol") %>%
      na.omit()
    
    # Convert the ranked data frame into a named vector (required by FGSEA)
    ranks <-tibble::deframe(res_p005_rank)
    
    # Create a list to store FGSEA results for each GMT file
    fgsea_list <- list()
    
    # Loop through GMT files and perform FGSEA for each gene set
    for (k in GMT_files) {
      name <- gsub(".*/", "", k)  # Extract the file name
      name2 <- gsub(".symbols.gmt.txt", "", name)  # Remove the extension
   
      # Load the pathways from the GMT file
      pathway <- fgsea::gmtPathways(k)
      
      # Perform FGSEA multilevel analysis
      fgseaRes <- fgsea::fgseaMultilevel(pathways = pathway, 
                                         stats = ranks,
                                         minSize = min_size,
                                         maxSize = max_size,
                                         nPermSimple = n_perms,
                                         eps = 0)
      
      # Store the FGSEA result for the current gene set
      fgsea_list[[name2]] <- fgseaRes
    }
    
    # Combine the results from all GMT files and store in the final output list
    Fgsea_output[[i]] <- do.call("rbind", fgsea_list)
  }
  
  # Return the FGSEA output list
  return(Fgsea_output)
}

# Example usage:
# Fgsea_output <- perform_fgsea_analysis(Ttest_list, Reference_path)

 

# Function to get up and downregulated pathways and find common pathways
find_common_pathways <- function(Fgsea_output) {
  
  up_pathways_list <- list()
  down_pathways_list <- list()
  
  # Loop through each sample's FGSEA results
  for (sample in names(Fgsea_output)) {
    # Get the FGSEA results for the current sample
    sample_data <- Fgsea_output[[sample]]
    
    # Filter for upregulated pathways (padj < 0.05 and NES > 0)
    up_pathways <- sample_data %>%
      dplyr::filter(padj < 0.05 & NES > 0) %>%
      pull(pathway)  # Extract only the pathway names
    
    # Filter for downregulated pathways (padj < 0.05 and NES < 0)
    down_pathways <- sample_data %>%
      dplyr::filter(padj < 0.05 & NES < 0) %>%
      pull(pathway)  # Extract only the pathway names
    
    # Store the up and down pathways in respective lists
    up_pathways_list[[sample]] <- up_pathways
    down_pathways_list[[sample]] <- down_pathways
  }
  
  # Find the common upregulated pathways across all samples
  common_up_pathways <- Reduce(intersect, up_pathways_list)
  
  # Find the common downregulated pathways across all samples
  common_down_pathways <- Reduce(intersect, down_pathways_list)
  
  # Return the common up and downregulated pathways
  return(list(common_up_pathways = common_up_pathways, 
              common_down_pathways = common_down_pathways))
}

# Example usage:
# Assuming Fgsea_output is the result from the perform_fgsea_analysis function
#common_pathways <- find_common_pathways(Fgsea_output)


combine_selected_pathways <- function(Fgsea_output, selected_pathways) {
  
  # Create an empty list to store the selected pathway data
  combined_data <- list()
  
  # Loop through each sample in the Fgsea_output
  for (sample in names(Fgsea_output)) {
    
    # Get FGSEA results for the current sample
    sample_data <- Fgsea_output[[sample]]
    
    # Filter the FGSEA results to only include the selected pathways
    selected_data <- sample_data %>%
      dplyr::filter(pathway %in% selected_pathways) %>%
      mutate(Group = sample)  # Add the group (sample) column
    
    # Append the selected data to the list
    combined_data[[sample]] <- selected_data
    
  }  
  combined_data_df <- bind_rows(combined_data)
  
  return(combined_data_df)
}

pheatmap_singscore<- function (pathways,data,Annotation,annotation_colors) {
  Gene_select_anno= data%>%.[,pathways] %>%t()%>%na.omit()%>%as.data.frame()

  pheatmap::pheatmap(Gene_select_anno, cellwigermline=5, cellheight = 10,cellwidth = 10, scale = "row",
                     treeheight_row = 5,
                     show_rownames = T,show_colnames = F,
                     annotation_col= Annotation, annotation_colors = annotation_colors,
                     # annotation_row=Annotation,
                     #annotation_legend=Label_def,
                     cluster_rows = T,  cluster_cols = F,clustering_distance_rows = "euclidean")
  
}



Singscore_Bar_plot=function(data,pathway_name,x_name,fill_name,compare_list,color) {
  my_comparisons <- compare_list
  ggplot(data,aes_string(
    x =x_name,
    y = pathway_name,
    fill=fill_name
  ) )+  geom_violin( ) +ggtitle(gsub("_Normalized_singscore","", pathway_name))+ylab("Singscore")+xlab("")+    geom_boxplot(width=0.1)+ theme_bw() +
    theme(plot.title=element_text(size=14),axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1,size=0),axis.text.y = element_text(size = 8),legend.position="none") +
    scale_fill_manual(values=color)+ggpubr::stat_compare_means(comparisons =  my_comparisons)
}


Batch_singscore_boxplot=function(pathway_list,data,x_name,fillby,nCol,compare_list,color){
  library(gridExtra)
  
  plist=list()
  for(i in 1:length(pathway_list)) {
    Pathway=pathway_list[i]
    # return(Pathway)
    plist[[i]]=Singscore_Bar_plot(data,Pathway,x_name,fillby,compare_list,color)
    
  }
  #plist[["legend"]]=UMAP_legend_2(data, pathway_list[1] ,fillby,fillby)
  plist[["legend"]]=UMAP_legend_4(data,pathway_list[1] ,fillby,fillby,color)
  do.call("grid.arrange", c(plist, ncol=nCol))
  
}



Enrich_dotplot=function(FGSEA_out){
  library(ggbreak)
  FGSEA_out$NES=round(FGSEA_out$NES,2)
  #mycolor=RColorBrewer::brewer.pal(3,"YlOrRd")
  # data=FGSEA_out%>%mutate(Norm_padj=-log10(padj))%>%mutate(State=ifelse(NES>0&padj<0.05,"Activated",ifelse(NES<0&padj<0.05,"Inhibited", ifelse(NES>0&padj<0.07,"P_Activated",ifelse(NES<0&padj<0.07,"P_Inhibited","no_change")))))%>%arrange(desc(NES))
  data=FGSEA_out%>%mutate(Norm_padj=-log10(padj))%>%mutate(State=ifelse(NES>0&padj<0.05,"Activated",ifelse(NES<0&padj<0.05,"Inhibited","no_change")))%>%arrange(desc(NES))
  # Up_data=data%>%dplyr::filter(NES>0)%>%arrange(NES)
  #   DN_data=data%>%dplyr::filter(NES<0)
  #   data=rbind(DN_data,Up_data)
  # data$pathway=factor(data$pathway,levels =data$pathway )
  MAX_NES=max(data$NES)%>%ceiling()
  min_NES=min(data$NES)%>%floor()
  # ggplot(data) +geom_point(aes(x = NES, y = pathway, color = -log10(padj),size=gene_counts))+scale_size("gene_counts") +
  # color_set= setNames(c('red','green',"wheat","turquoise" ,'grey'),c("Activated", "Inhibited","P_Activated","P_Inhibited","no_change"))
  color_set= setNames(c('red','green','grey'),c("Activated", "Inhibited","no_change"))
  P= ggplot(data) +geom_point(aes(x = NES, y = pathway, color =State ,size=Norm_padj)) +
    theme_bw() +labs(x="Normalized enrichment score",y="")+
    theme(axis.text.x = element_text(size=rel(1.15)),
          axis.title = element_text(size=rel(1.15)),plot.title = element_text(hjust=0.5, 
                                                                              face = "bold"),legend.title = element_text(size=rel(1.15),hjust=0.5, face="bold"),axis.text.y = element_text(size =15))+
    scale_size_continuous (name="-log10(padj)",range=c(5,12))+scale_color_manual(name="State",values=color_set)+facet_wrap(~Group,nrow=1)
  # if (min_NES>0|MAX_NES<0) {
  #   return(P)
  #   } else {
  #   Breaks=c(min_NES,-1,1,MAX_NES)
  #  return(P+scale_x_break(c(-1, 1), scales = 1.5))
  # xlab("Normalized enrichment score") +ylab("")+scale_x_break(c(-1, 1), scales = 1.5)
  #   }
}


combined_self_define_merge <- function(Pathway_list_select, GSEA_pathway_list, Ttest_list) {
  
  # Internal function to define the term and process the data
  Self_define_term <- function(Pathway_names, Pathway_Data, Gene_input) {
    out <- list()
    colnames(Gene_input)[grepl("log", colnames(Gene_input))] <- "logFC"
    colnames(Gene_input)[grepl("pv", colnames(Gene_input))] <- "pvalue"
    colnames(Gene_input)[grepl("symbol", colnames(Gene_input))] <- "Gene.symbol"
    
    for (i in Pathway_names) {
      Gene_in_pathway <- Pathway_Data[[i]] %>% as.character()
      Gene_list <- Gene_input %>%
        dplyr::filter(Gene.symbol %in% Gene_in_pathway) %>%
        dplyr::filter(pvalue < 0.05 & abs(logFC) > 0.27) %>%
        mutate(Change = ifelse(logFC > 0.27, "Up_regulated", "Down_regulated"))
      
      Up_Genes <- Gene_list %>% dplyr::filter(Change == "Up_regulated")
      DN_Genes <- Gene_list %>% dplyr::filter(Change != "Up_regulated")
      
      Gene_anno <- data.frame(
        Term_Description = i,
        lowest_p = NA,
        Up_regulated = toString(Up_Genes$Gene.symbol, collapse = ", "),
        Down_regulated = toString(DN_Genes$Gene.symbol, collapse = ", ")
      )
      out[[i]] <- Gene_anno
    }
    
    out_all <- do.call("rbind", out)
    return(out_all)
  }
  
  # Internal function to merge and rename genes for multiple groups
  merge_and_rename_genes_multi <- function(all_pathways_data, group_names) {
    Out_all <- data.frame()
    
    for (pathway_name in unique(all_pathways_data[[1]]$Term_Description)) {
      merged_list <- list()
      all_genes <- c()
      
      # Collect data from all groups
      for (i in seq_along(all_pathways_data)) {
        Data <- all_pathways_data[[i]] %>% 
          dplyr::filter(Term_Description == pathway_name) %>% 
          mutate(group = group_names[i])
        
        Up_regulated <- stringr::str_split(Data$Up_regulated, ", ") %>% unlist() %>% as.character()
        Down_regulated <- stringr::str_split(Data$Down_regulated, ", ") %>% unlist() %>% as.character()
        
        merged_genes <- c(Up_regulated, Down_regulated)
        all_genes <- unique(c(all_genes, merged_genes))
        
        merged_list[[i]] <- list(data = Data, merged_genes = merged_genes)
      }
      
      # Identify missing genes for each group
      for (i in seq_along(merged_list)) {
        missing_genes <- setdiff(all_genes, merged_list[[i]]$merged_genes)
        merged_list[[i]]$data$no_change <- toString(missing_genes, collapse = ", ")
        Out_all <- rbind(Out_all, merged_list[[i]]$data)
      }
    }
    
    Out_all$Term_Description <- toupper(Out_all$Term_Description)
    return(Out_all)
  }
  
  # Get the group names automatically from the names of the Ttest_list
  group_names <- names(Ttest_list)
  
  # Process each dataset in Ttest_list using Self_define_term
  all_pathways_data <- lapply(seq_along(Ttest_list), function(i) {
    Self_define_term(Pathway_list_select, GSEA_pathway_list, Ttest_list[[i]])
  })
  
  
  # Merge the results for all groups using merge_and_rename_genes_multi
  Self_pathway_merge <- merge_and_rename_genes_multi(all_pathways_data, group_names)
  
  return(Self_pathway_merge)
}

Term_gene_graph_3_update_2 <- function(result_df, Ttest_list, num_terms = 10, layout = "stress", 
                                       use_description = FALSE, node_size = "Fold change", 
                                       node_colors = c("#E5D7BF", "#E7211A","#18499E",  "grey"), 
                                       facet_variable = NULL) {
  # Check input validity
  if (!is.numeric(num_terms) & !is.null(num_terms)) stop("`num_terms` must either be numeric or NULL!")
  if (!is.logical(use_description)) stop("`use_description` must either be TRUE or FALSE!")
  
  ID_column <- ifelse(use_description, "Term_Description", "ID")
  val_node_size <- c("num_genes", "logFC")
  #  if (!node_size %in% val_node_size) stop("`node_size` should be one of ", paste(dQuote(val_node_size), collapse = ", "))
  if (!is.data.frame(result_df)) stop("`result_df` should be a data frame")
  
  necessary_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated", "no_change", "group")
  if (!all(necessary_cols %in% colnames(result_df))) stop("`result_df` must have all necessary columns!")
  
  if (length(node_colors) != 4) stop("`node_colors` must contain exactly 4 colors")
  
  if (!is.null(num_terms)) {
    if (nrow(result_df) < num_terms) num_terms <- NULL
  }
  
  Group <- unique(result_df$group)
  glist <- list()
  input <- result_df
  
  for (G in Group) {
    result_df <- input[order(input$Term_Description, decreasing = FALSE), ] %>% dplyr::filter(group == G)
    
    if (!is.null(num_terms)) result_df <- result_df[1:num_terms, ]
    
    graph_df <- data.frame()
    ###############################################
    for (i in base::seq_len(nrow(result_df))) {
      up_genes <- unlist(strsplit(result_df$Up_regulated[i], ", "))
      down_genes <- unlist(strsplit(result_df$Down_regulated[i], ", "))
      no_change <- unlist(strsplit(result_df$no_change[i], ", "))
      
      for (gene in c(up_genes, down_genes, no_change)) {
        graph_df <- rbind(graph_df, data.frame(Term = result_df[i, ID_column], Gene = gene)) %>% arrange(Gene)
      }
    }
    
    up_genes <- unlist(lapply(result_df$Up_regulated, function(x) strsplit(x, ", ")))
    DN_genes <- unlist(lapply(result_df$Down_regulated, function(x) strsplit(x, ", ")))
    no_change_genes <- unlist(lapply(result_df$no_change, function(x) strsplit(x, ", ")))
    
    # Corrected logFC_genes filtering for the current group G
    logFC_genes <- Ttest_list[[G]] %>%
      dplyr::filter(gene_symbol %in% c(up_genes, DN_genes, no_change_genes)) %>%
      dplyr::select(gene_symbol, log2Foldchange, pvalue)%>%mutate(FC=2^log2Foldchange)
    
    # Create the graph
    g <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
    cond_term <- names(igraph::V(g)) %in% result_df[, ID_column]
    cond_up_gene <- names(igraph::V(g)) %in% up_genes
    cond_down_gene <- names(igraph::V(g)) %in% DN_genes
    cond_no_change_gene <- names(igraph::V(g)) %in% no_change_genes
    
    node_type <- ifelse(cond_term, "term", 
                        ifelse(cond_up_gene, "up", 
                               ifelse(cond_down_gene, "down", 
                                      ifelse(cond_no_change_gene, "no_change", NA))))
    node_type <- factor(node_type, levels = c("term", "up", "down", "no_change"))
    
    igraph::V(g)$type <- node_type
    
    type_descriptions <- c(term = "enriched term", up = "up-regulated gene", 
                           down = "down-regulated gene", no_change = "no_change gene")
    names(node_colors) <- c("term", "up", "down", "no_change")
    
    node_colors <- node_colors[levels(node_type)]
    
    if (node_size == "num_genes") {
      sizes <- igraph::degree(g)
      sizes <- ifelse(igraph::V(g)$type == "term", sizes, 2)
      size_label <- "# genes"
    } else {
      # Replace sizes with logFC (absolute)
      logFC_idx <- match(names(igraph::V(g)), logFC_genes$gene_symbol)
      sizes <- abs(logFC_genes$FC[logFC_idx])
      sizes[is.na(sizes)] <- 5
      size_label <- "Fold change"
    }
    
    igraph::V(g)$size <- sizes
    igraph::V(g)$label.cex <- 0.5
    igraph::V(g)$frame.color <- "gray"
    glist[[G]] <- g
  }
  
  plist <- list()
  ######################################################################3
  for (s in Group) {
    g <- glist[[s]]
    p <- ggraph::ggraph(g, layout = layout) +
      ggraph::geom_edge_link(alpha = 0.8, colour = "darkgrey") +
      ggraph::geom_node_point(ggplot2::aes(color = type, size = size)) +
      ggplot2::scale_size(range = c(5, 10), breaks = c(0, 0.25, 0.5, 1, 2), name = size_label)+
      ggplot2::scale_color_manual(values = node_colors, name = NULL, labels = type_descriptions) +
      ggplot2::theme_void() +
      suppressWarnings(ggraph::geom_node_text(ggplot2::aes(label = .data$name), nudge_y = 0.05, repel = F, max.overlaps = 20)) +
      #suppressWarnings(ggraph::geom_node_text(ggrepel::geom_text_repel(aes(label = igraph::V(g)$name), max.overlaps = 20))) +  # Use ggrepel for label repelling
      ggplot2::ggtitle(s) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), 
                     plot.subtitle = ggplot2::element_text(hjust = 0.5))
    plist[[s]] <- p
  }
  
  return(plist)
}

run_pathview_analysis <- function(test_data, pathway_name, pathway_id, output_dir,outname) {
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Map gene symbols to Entrez IDs
  test_data$entrez_id <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db, 
    keys = test_data$gene_symbol, 
    column = "ENTREZID", 
    keytype = "SYMBOL", 
    multiVals = "first"
  )
  
  # Remove rows with missing Entrez IDs
  test_data <- test_data %>% dplyr::filter(!is.na(entrez_id))
  
  # if log2FoldChange is not in the column names, change it to log2FoldChange
  if (!"log2FoldChange" %in% colnames(test_data)) {
    colnames(test_data)[colnames(test_data) == "log2Foldchange"] <- "log2FoldChange"
  }
  
  # Prepare gene data for pathview
  gene_data <- test_data$log2FoldChange
  names(gene_data) <- test_data$entrez_id
  
  # Set output filename based on pathway name
  
  # Run pathview for the given pathway
  pathview_result <- pathview::pathview(
    gene.data = gene_data, 
    pathway.id = pathway_id,  # Provide pathway ID (e.g., "04110")
    species = "hsa", 
    out.dir = output_dir, 
    bins = list(gene = 5, cpd = 5),  # Adjust bins if needed
    kegg.native = TRUE
  )
  
  outfie_name=file.path(".",paste0(pathway_id,".pathview.png"))
  newname=file.path(output_dir,outname)
  file.rename(outfie_name,newname)
  png=file.path(".",paste0(pathway_id,".png"))
  xml=file.path(".",paste0(pathway_id,".xml"))
  file.remove(png)
  file.remove(xml)
  # Return the pathview result object
  
}


################################
# BiocManager::install("BridgeDbR")
# install.packages("XMLRPC", repos="http://R-Forge.R-project.org")
# library(BridgeDbR), derby <- getDatabase("Rattus norvegicus",location=getwd()) in tutorial.Rnw does not work. need change with 
# derby <- ("C:/R/Packages/RPathVisio/RPathVisio/vignettes/Rn_Derby_Ensembl_91.bridge")
# BiocManager::install("rWikiPathways")
###################
# wikipathway_list
# TNF_pathways <- rWikiPathways::findPathwayIdsByXref("TNF",'H')
# hs.pathways <- rWikiPathways::listPathways('Homo sapiens') 
# hs.pathwyas_unlist=lapply(hs.pathways,function(x) unlist(x) %>% as.data.frame()) %>% do.call("cbind",.) %>% t()
# 
# write.table(hs.pathwyas_unlist,file= "C:/R/project/8_IgA_proteome/Raw_data/wikipathways_list.txt",row.names = F,sep="\t")
#################################
##########################################################################################################################################


create_gex_file_EnID <- function(test_data, output_file, row.names = TRUE) {
  # Map gene symbols to Entrez IDs
  test_data$entrez_id <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db, 
    keys = test_data$gene_symbol, 
    column = "ENTREZID", 
    keytype = "SYMBOL", 
    multiVals = "first"
  )
  
  # Remove rows with missing Entrez IDs
  test_data <- test_data %>% dplyr::filter(!is.na(entrez_id))
  
  # Select and rename columns for GEX format
  test_data2 <- test_data %>%
    dplyr::select(entrez_id, log2Foldchange) %>%
    dplyr::rename(ID = entrez_id, `Expression Value` = log2Foldchange) %>%
    mutate(`System Code` = "L") %>%
    dplyr::select(ID, `System Code`, `Expression Value`)  # Ensure column order
  
  # Write the data to a tab-delimited file
  write.table(test_data2, file = file.path("./Pathvisio",output_file), sep = "\t", row.names = row.names, 
              col.names = TRUE, quote = FALSE)
  
  message("GEX file created at: ", output_file)
}



data_reorg= function (data_input) {
  
  data_input2<- data_input %>% mutate (Sig=ifelse (padj<0.05,1,0)) %>%
    mutate(FC_range=ifelse(log2FoldChange>=1,"4",
                           ifelse(log2FoldChange>=0.67&log2FoldChange<1,"3",
                                  ifelse(log2FoldChange>= 0.27&log2FoldChange<0.67,"2",
                                         ifelse(log2FoldChange>0, "1",
                                                ifelse(log2FoldChange <= -1,"-4",
                                                       ifelse(log2FoldChange<= -0.67  &log2FoldChange> -1,"-3",
                                                              ifelse(log2FoldChange<= -0.27&log2FoldChange> -0.67,"-2",
                                                                     ifelse(log2FoldChange<0,"-1", "0")))))))))
  
  data_input3=data_input2 %>% mutate(FC_label=as.integer (FC_range) * Sig)
  
  # (Sig_fold_change=ifelse(log2FC>=1 & pvalue<0.05,"3",
  #                                                             ifelse(log2FC>=0.67&pvalue<0.05$log2FC<1,"2",
  #                                                                    ifelse(log2FC>=0.27&pvalue<0.05$log2FC<0.67,"1",
  #                                                                           ifelse(log2FC<=-1 & pvalue<0.05,"-3", 
  #                                                                                  ifelse(log2FC> -1$log2FC<= -0.67 & pvalue<0.05,"-2",
  #                                                                                         ifelse(log2FC >-0.67 &log2FC<= -0.27 &pvalue<0.05,"-1",
  #                                                                                                ifelse(log2FC > 0 &pvalue<0.05,"-1","" 
  #                                                                                                )))))))
  
}  

pathvisioServer <- function(port=9000) {
  ##work around, otherwise the driver can't be found.
  .jnew("org/apache/derby/jdbc/EmbeddedDriver")
  
  server = .jnew("org/pathvisio/xmlrpc/JavaServer")
  start = .jcall(server, "V", "main", c(port,""))
  paste("http://localhost:", port, "/", sep="") 
  port
}



getParticipants <- function(pathway, type, host="localhost", port=9000, path=NA, outputdir=NA) {
  if (missing(pathway)) stop("You must provide a pathway name.");
  if (missing(type)) stop("You must provide which participants you wants; i.e GeneProduct, Protein, Metabolite or Reactions/Interactions.");
  if (is.na(path)) path = getwd();
  if (is.na(outputdir)) outputdir = tempdir();
  
  pwyPath = paste(path, "/", pathway, ".gpml", sep="")
  hostUrl = paste("http://", host, ":", port, "/", sep="")
  xml.rpc(hostUrl, "PathVisio.getPathwayParticipants", pwyPath, type, outputdir)
}



createVisualization <- function(gexname, gsample="", gcolors="", gvalues="", rsample="", rcolors="", rexprs="", host="localhost", port="9000", gexpath=NA) {
  if (missing(gexname)) stop("You must provide the name of the gexfile to use.");
  if ((gsample=="" | gcolors=="" | gvalues=="")
      && (rsample=="" | rcolors=="" | rexprs=="")) 
    stop("You must either provide a complete set of gradient values or a complete set of color rules.");
  if (is.na(gexpath)) gexpath = getwd();
  
  exts = c("txt",".pgex")
  hostUrl = paste("http://", host, ":", port, "/", sep="")
  # for (ext in exts) {
  #   g = paste0(gexpath,"/",gexname,".",ext)
  #   if (file.exists(g)) gex = g
  # }
  gex=paste0(gexpath,"/",gexname,".txt",".pgex")
  rexprs = paste("[",unlist(strsplit(rsample,";")),"]",unlist(strsplit(rexprs,";")),sep="",collapse=";")
  xml.rpc(hostUrl,"PathVisio.createVisualization", gex, gsample, gcolors, gvalues, rsample, rcolors, rexprs)
}

importData <- function(name, dataframe, dbname, host="localhost", port=9000, filepath=NA, dbpath=NA, outputdir=NA, row.names=TRUE, source=NA) {
  if (missing(name)) stop("You must provide a name for the pgex file")
  if (missing(dataframe)) stop("You must provide a table.");
  if (missing(dbname)) stop("You must provide the name of the database to use for mapping the data.");
  if (is.na(filepath)) filepath = tempdir();
  if (row.names) {
    dataframe <- cbind(rownames(dataframe),dataframe)
    colnames(dataframe)[1] <- "ID"
  }
  hostUrl=paste("http://",host,":",port,sep="")
  # if there's no source given, have bridgedb guess the source
  if (is.na(source)) {
    firstid = dataframe[1,1]
    list <- BridgeDbR::getMatchingSources(as.character(firstid))
    if (is.na(list[1])) stop("identifier pattern not recognized")
    # if there's just a single system code possible, use that system code
    # but if there are more, stop and return the possible source's
    # when the list is empty, report an error
    if (length(list)==1) {
      source = BridgeDbR::getSystemCode(list[1])
    } else if (!class(list) == "try-error") {
      strlist = paste(list,collapse=(", "))
      stop(paste("Unable to detect source, possible source's are:",strlist))
    }
    else stop("Incorrect data frame, unable to match identifiers to a source");
  }
  # else try if source contains the full name 
  else {
    tryname = BridgeDbR::getSystemCode(source)
    if (!is.na(tryname)) source = tryname;
  }
  # check if the system code is set correctly
  res = BridgeDbR::getFullName(source)
  if (is.na(res)) stop ("Invalid source");
  l = ncol(dataframe)
  scnum = l + 1
  dataframe["System Code"] <- source
  dataframe = dataframe[,c(1,scnum,2:l)]
  file = paste(filepath,"/",name,".txt",sep="")
  write.table(dataframe,file,sep="\t",row.names=FALSE,quote=FALSE)
  importDataByFile(name,dbname,host=host,port=port,filepath=filepath, dbpath=dbpath, outputdir=outputdir)
}






importDataByFile <- function(filename, dbname, host="localhost", port=9000, filepath=NA, dbpath=NA, outputdir=NA) {
  if (missing(filename)) stop("You must provide the name of a tab delimited data file.");
  if (missing(dbname)) stop("You must provide the name of the database to use for mapping the data.");
  if (is.na(filepath)) filepath = getwd();
  if (is.na(dbpath)) dbpath = getwd();
  if (is.na(outputdir)) outputdir = getwd();
  
  fexts = c(".txt","",".csv")
  dbexts = c(".bridge","",".pgdb")
  hostUrl = paste("http://", host, ":", port, "/", sep="")
  for (fext in fexts) {
    f = paste(filepath,"/",filename,fext,sep="")
    if (file.exists(f)) file = f
  }
  for (dbext in dbexts) {
    d = paste(dbpath,"/",dbname,dbext,sep="")
    if (file.exists(d)) db = d
  }
  if(is.na(file)) stop("expression file not found")
  if(is.na(db)) stop("database file not found")
  xml.rpc(hostUrl,"PathVisio.importData",file,db,outputdir)
}



visualizeDataByURI <- function(uri, gexname, dbname, host="localhost", port=9000, pwypath=NA, gexpath=NA, dbpath=NA, outputdir=NA) {
  if (missing(uri)) stop("You must provide WikiPathway ID for the pathway");
  if (missing(gexname)) stop("You must provide the name of the gexfile to use.");
  if (missing(dbname)) stop("You must provide the name of the database to use for mapping the data.");
  if (is.na(gexpath)) gexpath = getwd();
  if (is.na(dbpath)) dbpath = getwd();
  if (is.na(outputdir)) outputdir = getwd();
  
  gexexts = c(".pgex","")
  dbexts = c(".bridge","",".pgdb")
  for (ext in gexexts) {
    g = paste(gexpath,"/",gexname,".txt",ext,sep="")
    if (file.exists(g)) gex = g
  }
  for (ext in dbexts) {
    d = paste(dbpath,"/",dbname,ext,sep="")
    if (file.exists(d)) db = d
  } 
  hostUrl = paste("http://", host, ":", port, "/", sep="")
  xml.rpc(hostUrl, "PathVisio.visualizeDataByURI", uri, gex, db, outputdir)
}

visualizeData <- function(pathway, gexname, dbname, host="localhost", port=9000, pwypath=NA, gexpath=NA, dbpath=NA, outputdir=NA) {
  if (missing(pathway)) stop("You must provide a pathway name.");
  if (missing(gexname)) stop("You must provide the name of the gexfile to use.");
  if (missing(dbname)) stop("You must provide the name of the database to use for mapping the data.");
  if (is.na(pwypath)) pwypath = getwd();
  if (is.na(gexpath)) gexpath = getwd();
  if (is.na(dbpath)) dbpath = getwd();
  if (is.na(outputdir)) outputdir = getwd();
  
  gexexts = c(".pgex","")
  dbexts = c(".bridge","",".pgdb")
  for (ext in gexexts) {
    g = paste(gexpath,"/",gexname,".txt",ext,sep="")
    if (file.exists(g)) gex = g
  }
  for (ext in dbexts) {
    d = paste(dbpath,"/",dbname,ext,sep="")
    if (file.exists(d)) db = d
  } 
  pwy = paste(pwypath,"/",pathway,".gpml",sep="")
  hostUrl = paste("http://", host, ":", port, "/", sep="")
  xml.rpc(hostUrl, "PathVisio.visualizeData", pwy, gex, db, outputdir)
}

visualizeData_2 <- function(pathway_uri, gexname, dbname, host="localhost", port=9000, Wikipath,pwypath=NA, gexpath=NA, dbpath=NA, outputdir=NA) {
  # if (missing(pathway)) stop("You must provide a pathway name.");
  if (missing(gexname)) stop("You must provide the name of the gexfile to use.");
  if (missing(dbname)) stop("You must provide the name of the database to use for mapping the data.");
  if (is.na(pwypath)) pwypath = getwd();
  if (is.na(gexpath)) gexpath = getwd();
  if (is.na(dbpath)) dbpath = getwd();
  if (is.na(outputdir)) outputdir = getwd();
  gmpl=list.files(wikipath,full.names = T)
  pwy_select= gmpl[grep (pathway_uri,gmpl)]
  # return(pwy)
  # return(pwy)
  gexexts = c(".pgex","")
  dbexts = c(".bridge","",".pgdb")
  for (ext in gexexts) {
    g = paste(gexpath,"/",gexname,ext,sep="")
    if (file.exists(g)) gex = g
  }
  
  for (ext in dbexts) {
    d = paste(dbpath,"/",dbname,ext,sep="")
    if (file.exists(d)) db = d
  }
  # db=file.path(Pathvisio_path,dbname)
  # return(db)
  # pwy = paste(pwypath,"/",pathway,".gpml",sep="")
  hostUrl = paste("http://", host, ":", port, "/", sep="")
  for (i in 1:length(pwy_select)) {
    pwy=pwy_select[i]
    xml.rpc(hostUrl, "PathVisio.visualizeData", pwy, gex, db, outputdir)
  }
}

pathway_visual =function (outputdir,gexdata_name,pathway_index) {
  visualizeData_2 (pathway=pathway_index, gexname=gexdata_name, dbname="Hs_Derby_Ensembl_91.bridge", host="localhost",  port=9000,Wikipath=wikipath, 
                   gexpath=outputdir,dbpath=Pathvisio_install_path, outputdir=outputdir)
  
}



Enrich_dotplot_pval=function(FGSEA_out){
  library(ggbreak)
  FGSEA_out$NES=round(FGSEA_out$NES,2)
  #mycolor=RColorBrewer::brewer.pal(3,"YlOrRd")
  data=FGSEA_out%>%mutate(Norm_pval=-log10(pval))%>%mutate(State=ifelse(NES>0&pval<0.05,"Activated",ifelse(NES<0&pval<0.05,"Inhibited", ifelse(NES>0&pval<0.07,"P_Activated",ifelse(NES<0&pval<0.07,"P_Inhibited","no_change")))))%>%arrange(desc(NES))
  
  Up_data=data%>%dplyr::filter(NES>0)%>%arrange(NES)
  DN_data=data%>%dplyr::filter(NES<0)
  data=rbind(DN_data,Up_data)
  data$pathway=factor(data$pathway,levels =data$pathway )
  MAX_NES=max(data$NES)%>%ceiling()
  min_NES=min(data$NES)%>%floor()
  # ggplot(data) +geom_point(aes(x = NES, y = pathway, color = -log10(pval),size=gene_counts))+scale_size("gene_counts") +
  color_set= setNames(c('red','green',"wheat","turquoise" ,'grey'),c("Activated", "Inhibited","P_Activated","P_Inhibited","no_change"))
  return(data)
  P= ggplot(data) +geom_point(aes(x = NES, y = pathway, color =State ,size=Norm_pval)) +
    theme_bw() +labs(x="Normalized enrichment score",y="")+
    theme(axis.text.x = element_text(size=rel(1.15)),
          axis.title = element_text(size=rel(1.15)),plot.title = element_text(hjust=0.5, 
                                                                              face = "bold"),legend.title = element_text(size=rel(1.15),hjust=0.5, face="bold"),axis.text.y = element_text(size =15))+
    scale_size_continuous (name="-log10(pval)",range=c(2,12))+scale_color_manual(name="State",values=color_set) 
  P+facet_wrap(~Group,nrow=1)
  # if (min_NES>0|MAX_NES<0) {
  #   return(P)
  #   } else {
  #   Breaks=c(min_NES,-1,1,MAX_NES)
  #  return(P+scale_x_break(c(-1, 1), scales = 1.5))
  # xlab("Normalized enrichment score") +ylab("")+scale_x_break(c(-1, 1), scales = 1.5)
  #   }
}


Enrich_dotplot=function(FGSEA_out){
  library(ggbreak)
  FGSEA_out$NES=round(FGSEA_out$NES,2)
  #mycolor=RColorBrewer::brewer.pal(3,"YlOrRd")
  # data=FGSEA_out%>%mutate(Norm_padj=-log10(padj))%>%mutate(State=ifelse(NES>0&padj<0.05,"Activated",ifelse(NES<0&padj<0.05,"Inhibited", ifelse(NES>0&padj<0.07,"P_Activated",ifelse(NES<0&padj<0.07,"P_Inhibited","no_change")))))%>%arrange(desc(NES))
  data=FGSEA_out%>%mutate(Norm_padj=-log10(padj))%>%mutate(State=ifelse(NES>0&padj<0.05,"Activated",ifelse(NES<0&padj<0.05,"Inhibited","no_change")))%>%arrange(desc(NES))
  # Up_data=data%>%dplyr::filter(NES>0)%>%arrange(NES)
  #   DN_data=data%>%dplyr::filter(NES<0)
  #   data=rbind(DN_data,Up_data)
  # data$pathway=factor(data$pathway,levels =data$pathway )
  MAX_NES=max(data$NES)%>%ceiling()
  min_NES=min(data$NES)%>%floor()
  # ggplot(data) +geom_point(aes(x = NES, y = pathway, color = -log10(padj),size=gene_counts))+scale_size("gene_counts") +
  # color_set= setNames(c('red','green',"wheat","turquoise" ,'grey'),c("Activated", "Inhibited","P_Activated","P_Inhibited","no_change"))
  color_set= setNames(c('red','green','grey'),c("Activated", "Inhibited","no_change"))
  P= ggplot(data) +geom_point(aes(x = NES, y = pathway, color =State ,size=Norm_padj)) +
    theme_bw() +labs(x="Normalized enrichment score",y="")+
    theme(axis.text.x = element_text(size=rel(1.15)),
          axis.title = element_text(size=rel(1.15)),plot.title = element_text(hjust=0.5, 
                                                                              face = "bold"),legend.title = element_text(size=rel(1.15),hjust=0.5, face="bold"),axis.text.y = element_text(size =15))+
    scale_size_continuous (name="-log10(padj)",range=c(5,12))+scale_color_manual(name="State",values=color_set)+facet_wrap(~Group,nrow=1)
  # if (min_NES>0|MAX_NES<0) {
  #   return(P)
  #   } else {
  #   Breaks=c(min_NES,-1,1,MAX_NES)
  #  return(P+scale_x_break(c(-1, 1), scales = 1.5))
  # xlab("Normalized enrichment score") +ylab("")+scale_x_break(c(-1, 1), scales = 1.5)
  #   }
}

Enrich_dotplot_bydataset = function(FGSEA_out){
  library(ggplot2)
  library(dplyr)
  library(ggbreak)
  
  # Round NES
  FGSEA_out$NES = round(FGSEA_out$NES, 2)
  
  # Create data with extra columns
  data = FGSEA_out %>%
    mutate(Norm_padj = -log10(padj)) %>%
    mutate(State = ifelse(NES > 0 & padj < 0.05, "Activated",
                          ifelse(NES < 0 & padj < 0.05, "Inhibited", "no_change"))) %>%
    arrange(Group, Database, desc(NES))  # Sort by Group, Database, NES
  
  # Insert a gap between Databases inside each Group
  data = data %>%
    group_by(Group) %>%
    group_modify(~ {
      db_ordered = .x %>%
        arrange(Database, desc(NES))
      
      if (length(unique(db_ordered$Database)) > 1) {
        gap_row = db_ordered[1, ]
        gap_row$pathway = " "  # Insert empty line
        gap_row$NES = NA
        gap_row$Norm_padj = NA
        gap_row$State = NA
        gap_row$Database = NA
        gap_row$Group = unique(db_ordered$Group)  # IMPORTANT: Keep Group
        db_ordered = bind_rows(
          db_ordered[db_ordered$Database == unique(db_ordered$Database)[1], ],
          gap_row,
          db_ordered[db_ordered$Database == unique(db_ordered$Database)[2], ]
        )
      }
      return(db_ordered)
    }) %>%
    ungroup()
  
  # Set pathway as factor with custom levels
  data$pathway = factor(data$pathway, levels = rev(unique(data$pathway)))
  
  # Color settings
  color_set = setNames(c('red', 'green', 'grey'), c("Activated", "Inhibited", "no_change"))
  
  
  gap_y = which(levels(data$pathway) == " ")
  
  # Make ggplot
  p = ggplot(data, aes(x = NES, y = pathway, color = State, size = Norm_padj)) +
    geom_point(na.rm = TRUE)+ geom_hline(yintercept = gap_y,size=1 ,linetype = "solid", color = "black") +
    theme_bw() +
    labs(x = "Normalized enrichment score", y = "") +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),
      axis.text.x = element_text(size = rel(1.15)),
      axis.title = element_text(size = rel(1.15)),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(size = rel(1.15), hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 15)
    ) +
    scale_size_continuous(name = "-log10(padj)", range = c(5, 12)) +
    scale_color_manual(name = "State", values = color_set) +
    facet_wrap(~ Group, nrow = 1)  # <-- Keep facet by Group
  
  # Optional: Add x-axis break at (-1,1)
  MAX_NES = ceiling(max(data$NES, na.rm = TRUE))
  min_NES = floor(min(data$NES, na.rm = TRUE))
  
  # if (min_NES <= -1 && MAX_NES >= 1) {
  #   p = p + scale_x_break(c(-1, 1), scales = 1.5)
  # }
  
  return(p)
}


Self_define_term=function(Pathway_names,Pathway_Data,Gene_input){
  necessary_cols <- c("Term_Description", "lowest_p", "Up_regulated", 
                      "Down_regulated", "no_change","group")
  out=list()
  colnames(Gene_input)[grepl("log",colnames(Gene_input))]="logFC"
  colnames(Gene_input)[grepl("pv",colnames(Gene_input))]="pvalue"
  colnames(Gene_input)[grepl("symbol",colnames(Gene_input))]="Gene.symbol"
  for(i in Pathway_names){
    Gene_in_pathway=Pathway_Data[[i]]%>%as.character()
    #Gene_list=Gene_input%>%dplyr::filter(Gene.symbol%in%Gene_in_pathway)%>%mutate(Change=ifelse(logFC>0.27&pvalue<0.05,"Up_regulated",ifelse(abs(logFC)>0.27&pvalue<0.05,"Down_regulated","no_change")))
    Gene_list=Gene_input%>%dplyr::filter(Gene.symbol%in%Gene_in_pathway)%>%dplyr::filter(pvalue<0.05&abs(logFC)>0.27)%>%mutate(Change=ifelse(logFC>0.27,"Up_regulated","Down_regulated"))
    Up_Genes=Gene_list%>%dplyr::filter(Change=="Up_regulated")
    DN_Genes=Gene_list%>%dplyr::filter(Change!="Up_regulated")
    Gene_anno=data.frame(Term_Description=i,lowest_p=NA,Up_regulated=toString (Up_Genes$Gene.symbol, collapse = ", "),
                         Down_regulated=toString  (DN_Genes$Gene.symbol, collapse = ", "))
    out[[i]]=Gene_anno
    
  }
  out_all=do.call("rbind",out)
}



merge_and_rename_genes <- function(pathway_out_1, pathway_out_2,group1_name="Group1",group2_name="Group2") {
  Out_2 <- data.frame()
  
  # Add missing genes column and rename for pathway_out_1
  for (pathway_name in pathway_out_1$Term_Description) {
    # Extract genes for the current pathway
    Data1 <- pathway_out_1 %>% 
      dplyr::filter(Term_Description == pathway_name) %>% 
      mutate(group = group1_name)
    
    Data2 <- pathway_out_2 %>% 
      dplyr::filter(Term_Description == pathway_name) %>% 
      mutate(group = group2_name)
    
    # Split Up_regulated and Down_regulated by ","
    Data1_Up_regulated <- stringr::str_split(Data1$Up_regulated, ", ")%>%unlist() %>%as.character()
    Data1_Down_regulated <- stringr::str_split(Data1$Down_regulated, ", ")%>%unlist() %>%as.character()
    Data2_Up_regulated <- stringr::str_split(Data2$Up_regulated, ", ")%>%unlist() %>%as.character()
    Data2_Down_regulated <- stringr::str_split(Data2$Down_regulated, ", ")%>%unlist() %>%as.character()
    
    # Merge genes for pathway_out_1 and pathway_out_2 for the current pathway
    merged_genes_1 <- c(Data1_Up_regulated, Data1_Down_regulated)
    merged_genes_2 <- c(Data2_Up_regulated, Data2_Down_regulated)
    all_genes=unique(c(merged_genes_1,merged_genes_2))
    
    missing_1 <- setdiff( all_genes, merged_genes_1)
    
    Data1$no_change <- toString  (missing_1, collapse = ", ")
    
    missing_genes_2 <- setdiff( all_genes, merged_genes_2)
    Data2$no_change <- toString(missing_genes_2,collapse = ", ")
    
    out <- rbind(Data1, Data2)
    Out_2 <- rbind(Out_2, out)
  }
  Out_2$Term_Description=toupper(Out_2$Term_Description)
  return(Out_2)
}


combined_self_define_merge <- function(Pathway_list_select, GSEA_pathway_list, Ttest_list) {
  
  # Internal function to define the term and process the data
  Self_define_term <- function(Pathway_names, Pathway_Data, Gene_input) {
    out <- list()
    colnames(Gene_input)[grepl("log", colnames(Gene_input))] <- "logFC"
    colnames(Gene_input)[grepl("pv", colnames(Gene_input))] <- "pvalue"
    colnames(Gene_input)[grepl("symbol", colnames(Gene_input))] <- "Gene.symbol"
    
    for (i in Pathway_names) {
      Gene_in_pathway <- Pathway_Data[[i]] %>% as.character()
      Gene_list <- Gene_input %>%
        dplyr::filter(Gene.symbol %in% Gene_in_pathway) %>%
        dplyr::filter(pvalue < 0.05 & abs(logFC) > 0.27) %>%
        mutate(Change = ifelse(logFC > 0.27, "Up_regulated", "Down_regulated"))
      
      Up_Genes <- Gene_list %>% dplyr::filter(Change == "Up_regulated")
      DN_Genes <- Gene_list %>% dplyr::filter(Change != "Up_regulated")
      
      Gene_anno <- data.frame(
        Term_Description = i,
        lowest_p = NA,
        Up_regulated = toString(Up_Genes$Gene.symbol, collapse = ", "),
        Down_regulated = toString(DN_Genes$Gene.symbol, collapse = ", ")
      )
      out[[i]] <- Gene_anno
    }
    
    out_all <- do.call("rbind", out)
    return(out_all)
  }
  
  # Internal function to merge and rename genes for multiple groups
  merge_and_rename_genes_multi <- function(all_pathways_data, group_names) {
    Out_all <- data.frame()
    
    for (pathway_name in unique(all_pathways_data[[1]]$Term_Description)) {
      merged_list <- list()
      all_genes <- c()
      
      # Collect data from all groups
      for (i in seq_along(all_pathways_data)) {
        Data <- all_pathways_data[[i]] %>% 
          dplyr::filter(Term_Description == pathway_name) %>% 
          mutate(group = group_names[i])
        
        Up_regulated <- stringr::str_split(Data$Up_regulated, ", ") %>% unlist() %>% as.character()
        Down_regulated <- stringr::str_split(Data$Down_regulated, ", ") %>% unlist() %>% as.character()
        
        merged_genes <- c(Up_regulated, Down_regulated)
        all_genes <- unique(c(all_genes, merged_genes))
        
        merged_list[[i]] <- list(data = Data, merged_genes = merged_genes)
      }
      
      # Identify missing genes for each group
      for (i in seq_along(merged_list)) {
        missing_genes <- setdiff(all_genes, merged_list[[i]]$merged_genes)
        merged_list[[i]]$data$no_change <- toString(missing_genes, collapse = ", ")
        Out_all <- rbind(Out_all, merged_list[[i]]$data)
      }
    }
    
    Out_all$Term_Description <- toupper(Out_all$Term_Description)
    return(Out_all)
  }
  
  # Get the group names automatically from the names of the Ttest_list
  group_names <- names(Ttest_list)
  
  # Process each dataset in Ttest_list using Self_define_term
  all_pathways_data <- lapply(seq_along(Ttest_list), function(i) {
    Self_define_term(Pathway_list_select, GSEA_pathway_list, Ttest_list[[i]])
  })
  
  
  # Merge the results for all groups using merge_and_rename_genes_multi
  Self_pathway_merge <- merge_and_rename_genes_multi(all_pathways_data, group_names)
  
  return(Self_pathway_merge)
}

# Example usage:
# result <- combined_self_define_merge(Pathway_list_select, GSEA_pathway_list, Ttest_list)


Term_gene_graph_3_update=function (result_df, num_terms = 10, layout = "stress", use_description = FALSE, 
                                   node_size = "num_genes", node_colors = c("#E5D7BF","#18499E","#E7211A","grey"),facet_variable = NULL) {
  # For paired plot
  
  if (!is.numeric(num_terms) & !is.null(num_terms)) {
    stop("`num_terms` must either be numeric or NULL!")
  }
  if (!is.logical(use_description)) {
    stop("`use_description` must either be TRUE or FALSE!")
  }
  ID_column <- ifelse(use_description, "Term_Description", 
                      "ID")
  val_node_size <- c("num_genes", "p_val")
  if (!node_size %in% val_node_size) {
    stop("`node_size` should be one of ", paste(dQuote(val_node_size), 
                                                collapse = ", "))
  }
  if (!is.data.frame(result_df)) {
    stop("`result_df` should be a data frame")
  }
  necessary_cols <- c(ID_column, "lowest_p", "Up_regulated", 
                      "Down_regulated", "no_change","group")  # Include "no_change"
  if (!all(necessary_cols %in% colnames(result_df))) {
    stop(paste(c("All of", paste(necessary_cols, collapse = ", "), 
                 "must be present in `results_df`!"), collapse = " "))
  }
  if (length(node_colors) != 4) {  # Corrected length check
    stop("`node_colors` must contain exactly 4 colors")  # Updated error message
  }
  if (!is.null(num_terms)) {
    if (nrow(result_df) < num_terms) {
      num_terms <- NULL
    }
  }
  Group=unique(result_df$group)
  
  ####################
  glist=list()
  input=result_df
  
  #####################################
  for (G in Group)    {
    result_df<- input[order(input$Term_Description, decreasing = FALSE),  ]%>%dplyr::filter(group==G)
    if (!is.null(num_terms)) {
      result_df <- result_df[1:num_terms, ]
    }
    graph_df <- data.frame()
    for (i in base::seq_len(nrow(result_df))) {
      up_genes <- unlist(strsplit(result_df$Up_regulated[i], 
                                  ", "))
      down_genes <- unlist(strsplit(result_df$Down_regulated[i], 
                                    ", "))
      
      no_change= unlist(strsplit(result_df$no_change[i], 
                                 ", "))
      
      for (gene in c(up_genes, down_genes,no_change)) {
        graph_df <- rbind(graph_df, data.frame(Term = result_df[i, 
                                                                ID_column], Gene = gene))%>%arrange(Gene)
      }
    }
    
    up_genes <- unlist(lapply(result_df$Up_regulated, function(x) strsplit(x, ", ")))
    DN_genes=  unlist(lapply(result_df$Down_regulated, function(x) strsplit(x, ", ")))
    no_change_genes <- unlist(lapply(result_df$no_change, function(x) strsplit(x, ", ")))  # Fix here
    
    # Assign node type correctly, including "no_change" condition
    
    
    g <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
    cond_term <- names(igraph::V(g)) %in% result_df[, ID_column]
    cond_up_gene <- names(igraph::V(g)) %in% up_genes
    cond_down_gene=  names(igraph::V(g)) %in% DN_genes
    cond_no_change_gene <- names(igraph::V(g)) %in% no_change_genes # fix here
    # node_type <- ifelse(cond_term, "term", 
    #                     ifelse(cond_up_gene, "up", 
    #                            ifelse(cond_down_gene, "down", "no_change")))
    
    node_type <- ifelse(cond_term, "term", 
                        ifelse(cond_up_gene, "up", 
                               ifelse(cond_down_gene, "down", 
                                      ifelse(cond_no_change_gene, "no_change", NA)))) # fix here
    
    node_type <- factor(node_type, levels = c("term", "up", "down", "no_change"))
    node_type <- factor(node_type, levels = c("term", "up", "down","no_change"))
    igraph::V(g)$type <- node_type
    type_descriptions <- c(term = "enriched term", up = "up-regulated gene", 
                           down = "down-regulated gene", no_change = "no_change gene")  # Include "no_change" description
    type_descriptions <- type_descriptions[levels(node_type)]
    names(node_colors) <- c("term", "up", "down", "no_change")  # Added "no_change"
    node_colors <- node_colors[levels(node_type)]
    
    if (node_size == "num_genes") {
      sizes <- igraph::degree(g)
      sizes <- ifelse(igraph::V(g)$type == "term", sizes, 2)
      size_label <- "# genes"
    }
    else {
      idx <- match(names(igraph::V(g)), result_df[, ID_column])
      sizes <- -log10(result_df$lowest_p[idx])
      sizes[is.na(sizes)] <- 2
      size_label <- "-log10(p)"
    }
    
    igraph::V(g)$size <- sizes
    igraph::V(g)$label.cex <- 0.5
    igraph::V(g)$frame.color <- "gray"
    glist[[G]]=g
  }    
  
  
  ##################################################
  
  plist=list()
  
  for (s in Group) {
    g=glist[[s]]
    #xy <- graphlayouts::layout_with_stress(glist[[1]])
    p <- ggraph::ggraph(g, layout = layout)
    
    p <- p + ggraph::geom_edge_link(alpha = 0.8, colour = "darkgrey")
    p <- p + ggraph::geom_node_point(ggplot2::aes(color =type, size =size))
    p <- p + ggplot2::scale_size(range = c(5, 10), breaks = round(seq(round(min(igraph::V(g)$size)), 
                                                                      round(max(igraph::V(g)$size)), length.out = 4)), name = size_label)
    p <- p + ggplot2::scale_color_manual(values = node_colors, 
                                         name = NULL, labels = type_descriptions)
    
    breaks = round(seq(round(min(igraph::V(g)$size))))
    return(list(igraph::V(g)$size ,breaks))
    p <- p + ggplot2::theme_void()
    p <- p + suppressWarnings(ggraph::geom_node_text(ggplot2::aes(label = .data$name), 
                                                     nudge_y = 0.2, repel = TRUE, max.overlaps = 20))
    
    
    if (is.null(num_terms)) {
      p <- p + ggplot2::ggtitle("Term-Gene Graph")
    }
    else {
      p <- p + ggplot2::ggtitle("Term-Gene Graph", subtitle = paste(c("Top", 
                                                                      num_terms, "terms"), collapse = " "))
    }
    p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), 
                            plot.subtitle = ggplot2::element_text(hjust = 0.5))+ggtitle(s)
    plist[[s]]=p
  }
  
  return(plist)
}

Term_gene_graph_3_update_2 <- function(result_df, Ttest_list, num_terms = 10, layout = "stress", 
                                       use_description = FALSE, node_size = "Fold change", 
                                       node_colors = c("#E5D7BF", "#E7211A","#18499E",  "grey"), 
                                       facet_variable = NULL) {
  # Check input validity
  if (!is.numeric(num_terms) & !is.null(num_terms)) stop("`num_terms` must either be numeric or NULL!")
  if (!is.logical(use_description)) stop("`use_description` must either be TRUE or FALSE!")
  
  ID_column <- ifelse(use_description, "Term_Description", "ID")
  val_node_size <- c("num_genes", "logFC")
  #  if (!node_size %in% val_node_size) stop("`node_size` should be one of ", paste(dQuote(val_node_size), collapse = ", "))
  if (!is.data.frame(result_df)) stop("`result_df` should be a data frame")
  
  necessary_cols <- c(ID_column, "lowest_p", "Up_regulated", "Down_regulated", "no_change", "group")
  if (!all(necessary_cols %in% colnames(result_df))) stop("`result_df` must have all necessary columns!")
  
  if (length(node_colors) != 4) stop("`node_colors` must contain exactly 4 colors")
  
  if (!is.null(num_terms)) {
    if (nrow(result_df) < num_terms) num_terms <- NULL
  }
  
  Group <- unique(result_df$group)
  glist <- list()
  input <- result_df
  
  for (G in Group) {
    result_df <- input[order(input$Term_Description, decreasing = FALSE), ] %>% dplyr::filter(group == G)
    
    if (!is.null(num_terms)) result_df <- result_df[1:num_terms, ]
    
    graph_df <- data.frame()
    ###############################################
    for (i in base::seq_len(nrow(result_df))) {
      up_genes <- unlist(strsplit(result_df$Up_regulated[i], ", "))
      down_genes <- unlist(strsplit(result_df$Down_regulated[i], ", "))
      no_change <- unlist(strsplit(result_df$no_change[i], ", "))
      
      for (gene in c(up_genes, down_genes, no_change)) {
        graph_df <- rbind(graph_df, data.frame(Term = result_df[i, ID_column], Gene = gene)) %>% arrange(Gene)
      }
    }
    
    up_genes <- unlist(lapply(result_df$Up_regulated, function(x) strsplit(x, ", ")))
    DN_genes <- unlist(lapply(result_df$Down_regulated, function(x) strsplit(x, ", ")))
    no_change_genes <- unlist(lapply(result_df$no_change, function(x) strsplit(x, ", ")))
    
    # Corrected logFC_genes filtering for the current group G
    logFC_genes <- Ttest_list[[G]] %>%
      dplyr::filter(gene_symbol %in% c(up_genes, DN_genes, no_change_genes)) %>%
      dplyr::select(gene_symbol, log2Foldchange, pvalue)%>%mutate(FC=2^log2Foldchange)
    
    # Create the graph
    g <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
    cond_term <- names(igraph::V(g)) %in% result_df[, ID_column]
    cond_up_gene <- names(igraph::V(g)) %in% up_genes
    cond_down_gene <- names(igraph::V(g)) %in% DN_genes
    cond_no_change_gene <- names(igraph::V(g)) %in% no_change_genes
    
    node_type <- ifelse(cond_term, "term", 
                        ifelse(cond_up_gene, "up", 
                               ifelse(cond_down_gene, "down", 
                                      ifelse(cond_no_change_gene, "no_change", NA))))
    node_type <- factor(node_type, levels = c("term", "up", "down", "no_change"))
    
    igraph::V(g)$type <- node_type
    
    type_descriptions <- c(term = "enriched term", up = "up-regulated gene", 
                           down = "down-regulated gene", no_change = "no_change gene")
    names(node_colors) <- c("term", "up", "down", "no_change")
    
    node_colors <- node_colors[levels(node_type)]
    
    if (node_size == "num_genes") {
      sizes <- igraph::degree(g)
      sizes <- ifelse(igraph::V(g)$type == "term", sizes, 2)
      size_label <- "# genes"
    } else {
      # Replace sizes with logFC (absolute)
      logFC_idx <- match(names(igraph::V(g)), logFC_genes$gene_symbol)
      sizes <- abs(logFC_genes$FC[logFC_idx])
      sizes[is.na(sizes)] <- 5
      size_label <- "Fold change"
    }
    
    igraph::V(g)$size <- sizes
    igraph::V(g)$label.cex <- 0.5
    igraph::V(g)$frame.color <- "gray"
    glist[[G]] <- g
  }
  
  plist <- list()
  ######################################################################3
  for (s in Group) {
    g <- glist[[s]]
    p <- ggraph::ggraph(g, layout = layout) +
      ggraph::geom_edge_link(alpha = 0.8, colour = "darkgrey") +
      ggraph::geom_node_point(ggplot2::aes(color = type, size = size)) +
      ggplot2::scale_size(range = c(5, 10), breaks = c(0, 0.25, 0.5, 1, 2), name = size_label)+
      ggplot2::scale_color_manual(values = node_colors, name = NULL, labels = type_descriptions) +
      ggplot2::theme_void() +
      suppressWarnings(ggraph::geom_node_text(ggplot2::aes(label = .data$name), nudge_y = 0.05, repel = F, max.overlaps = 20)) +
      #suppressWarnings(ggraph::geom_node_text(ggrepel::geom_text_repel(aes(label = igraph::V(g)$name), max.overlaps = 20))) +  # Use ggrepel for label repelling
      ggplot2::ggtitle(s) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), 
                     plot.subtitle = ggplot2::element_text(hjust = 0.5))
    plist[[s]] <- p
  }
  
  return(plist)
}

pheatmap_singscore<- function (pathways,data,Annotation,annotation_colors) {
  Gene_select_anno= data[,colnames(data) %in% pathways] %>%t()%>%na.omit()%>%as.data.frame()%>%.[,rownames(Annotation)]
  pheatmap::pheatmap(Gene_select_anno, cellwigermline=5, cellheight = 10,cellwidth = 10, scale = "row",
                     treeheight_row = 5,
                     show_rownames = T,show_colnames = F,
                     annotation_col= Annotation,annotation_colors = annotation_colors,
                     # annotation_row=Annotation,
                     #annotation_legend=Label_def,
                     cluster_rows = T,  cluster_cols = F,clustering_distance_rows = "euclidean")
  
}


pheatmap_singscore_2<- function (pathways,data,Annotation,annotation_colors) {
  Gene_select_anno= data[,colnames(data) %in% pathways] %>%t()%>%na.omit()%>%as.data.frame()%>%.[,rownames(Annotation)]
  quantile_breaks <- function(xs, n = 6) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
  }
  
  mat_breaks <- quantile_breaks(as.matrix(Gene_select_anno), n = 10)
  Color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 6, name = "RdYlBu")))(length(mat_breaks))
  pheatmap::pheatmap(Gene_select_anno, cellwigermline=5, cellheight = 25,cellwidth = 10, scale = "row",
                     treeheight_row = 5,color= Color,fontsize_row = 20,
                     show_rownames = T,show_colnames = F,
                     breaks = mat_breaks,
                     annotation_col= Annotation,annotation_colors = annotation_colors,fontsize = 20,
                     # annotation_row=Annotation,
                     #annotation_legend=Label_def,
                     cluster_rows = T,  cluster_cols = F,clustering_distance_rows = "euclidean")
  
}
###################################
