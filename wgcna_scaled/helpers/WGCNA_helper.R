pick_beta <- function(GSE_WGCNA_matrix, n_samples){
    powers = c(c(1:10), seq(from = 12, to=30, by=2))
    
    sft = pickSoftThreshold(GSE_WGCNA_matrix, powerVector = powers, verbose = 5)
    sft$fitIndices$plot_param <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
    if(max(sft$fitIndices$plot_param)<0.85){
      if(n_samples>20 & n_samples <30)
        beta <- 16
      else if(n_samples >30 & n_samples <40)
        beta <- 14
      else 
        beta <- 12
    }
    else
      beta <- sft$powerEstimate
    return(beta)
}

get_modules <- function(GSE_WGCNA_matrix, beta){
    #calclute the adjacency matrix
    adj= adjacency(GSE_WGCNA_matrix,type = "signed", power = beta);

    #turn adjacency matrix into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations
    TOM=TOMsimilarityFromExpr(GSE_WGCNA_matrix,networkType = "signed", 
                          TOMType = "signed", power = beta)
    
    # get dissimilarity matrix from TOM matrix
    colnames(TOM) =rownames(TOM) 
    dissTOM=1-TOM

    #hierarchical clustering of the genes based on the TOM dissimilarity measure
    geneTree = hclust(as.dist(dissTOM),method="average");

    # Set the minimum module size
    minModuleSize = 20;

    # cut the tree 
    dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);

    #the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
    dynamicMods <- dynamicMods + 1

    diag(dissTOM) = NA;
    dynamicColors = labels2colors(dynamicMods)
    print("Module Lengths")
    print(table(dynamicColors))
    module_colors= setdiff(unique(dynamicColors), "grey")
    
    module_list = list()
    for (color in module_colors){
        module= colnames(GSE_WGCNA_matrix)[which(dynamicColors==color)]
        
        # expand list if slashes exist in entrez ids  
        module <- unlist(sapply(module, function(x) sapply(strsplit(x, split = "///"), trimws), 
                                USE.NAMES = F))

        module <- module[grep("^[0-9]+$", module)]
	module <- module[!is.na(module)]
	module <- module[module!=""]
	module <- module[module!=" "]
	print(module)
	if(length(module) >300)
		next
        module_list[[color]] <- paste(module, collapse = ',')
    }
    cat(module_list[[1]][1:50])
    module_list <- module_list[sapply(module_list, nchar)>0]
    return(module_list)
}

write_modules <- function(module_list, modules_dir_path, gse_id, platform_name){
  
    module_file_path = paste(modules_dir_path,gse_id,".gmt", sep = "")  
    for(index in 1:length(module_list)){
        dataset_platform <- paste(gse_id, platform_name, sep = '_')
        platform_cluster <- paste(dataset_platform, index, sep = '#')

        
        
        line_to_write <- paste(platform_cluster, module_list[[index]], sep = '\t')

        write(line_to_write, file = module_file_path, append=TRUE)
        print(paste(platform_cluster, "module_added"))
        }
    }

write_gse_info <- function(gse_info_file_path, gse_id, abstract){
  
    info = file.info(gse_info_file_path)
    empty <- (info$size==0)
  
    if(!empty){
      gse_existing <- read.table(gse_info_file_path, sep = "\t", stringsAsFactors = F)
      if(!(gse_id %in% gse_existing[, 1]))
        write(paste(gse_id, abstract, sep = '\t'),
              file = gse_info_file_path, append = T)
    }
    else
      write(paste(gse_id, abstract, sep = '\t'),
            file = gse_info_file_path, append = T)
    }

