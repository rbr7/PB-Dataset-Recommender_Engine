run_module_pipeline_RNASeq <- function(GSE_id, module_file_path, gse_info_file_path, logs_path){
    
    source("helpers/GEO_helper.R")
    source("helpers/WGCNA_helper.R")
    
    stop_messages <- list("low_samples" = "less than 12 samples present.",
                          "high_samples" = "more than 400 samples present(maybe a single cell dataset)",
                          "low_genes" = "No rows in the matrix,GEOQuery did not get data properly",
                          "NA_genes" = "All rows contain columns with NAs",
                          "platform_annotation" = "Not able to annotate platform",
                          "success" = "WGCNA successfully completed!!!", 
                          "failure" = "failed", 
                          "platform_not_MA" = "platform not a Microrarray platform or not in annotations file",
                          "annot_not_defined" = "Annotation is not defined in platform annotation file")
    
    log_GSE_sink(GSE_id = GSE_id, logs_path = logs_path)
    
    
    tryCatch({
    GSE_object <- download_geo(GSE_id)
    message(cat("downloading data for ---------", GSE_id))
    
    fail <- 0
   
    #platform_annot_table <- read.csv("helpers/MA_platform_annotation.csv", stringsAsFactors = F)

    for(platform_index in 1:length(GSE_object)){
        
        module_list <- list()
        
        message("checking number of samples")
        GSE_platform <- get_GSE_platform(GSE_object, platform_index)
        
        GSE_n_samples <- ncol(GSE_object[[1]]@assayData$exprs)
        message(paste("no. of samples", GSE_n_samples))
        if(GSE_n_samples < 12){
            logerror(stop_messages[['low_samples']], logger=paste(GSE_id, GSE_platform, sep = "#"))
            fail <- 1
            next
        }
        
        if(GSE_n_samples > 400){
            logerror(stop_messages[['high_samples']], logger=paste(GSE_id, GSE_platform, sep = "#"))
            fail <- 1
            next
        }
        
        
        # download RAW counts data
        download_command <- paste("wget -r -nH -nc --cut-dirs=3 ftp://ftp.ncbi.nlm.nih.gov/geo/series/",
                                  paste(substr(GSE_id, 1,5), "nnn/", GSE_id, sep = ""), sep = "")
        print(download_command)
        
        system(download_command)
        
        message("Getting data for platform----- ", platform_index)
        
        sapply(GSE_object[[1]]@phenoData@data$supplementary_file_1, function(x) 
            system(paste("wget", as.character(x), "-P", GSE_id)))
        
        counts_file_names <- list.files(path = GSE_id)
        print(counts)
        
        GSE_combined_matrix <- edgeR::readDGE(counts_file_names, path = GSE_id, columns = c(1,3))
        dds <- DESeqDataSetFromMatrix(GSE_combined_matrix$counts, GSE_combined_matrix$samples, design = ~1)
        dds <- DESeq(dds)
        GSE_matrix <- assay(rlog(dds))
        
        print(dim(GSE_matrix))
        if(nrow(GSE_matrix)<100){
            logerror(stop_messages[['low_genes']], logger=paste(GSE_id, GSE_platform, sep = "#"))
            fail <- 1
            next
        }
        
        message(paste("GSE_platform ------", GSE_platform))
        GSE_tax_id <- get_GSE_taxonomy(GSE_object, platform_index)
        GSE_title <- get_GSE_title(GSE_object, platform_index)
        
        # annotate with entrez ID's GSE_matrix
        
        ## see which column to use from platform file        
        GSE_matrix <- GSE_matrix[!apply(GSE_matrix,1, anyNA),]
        print(dim(GSE_matrix))
        if(nrow(GSE_matrix)==0){
                logerror(stop_messages[['NA_genes']], logger=paste(GSE_id, GSE_platform, sep = "#"))
                fail <- 1
                next
            }
        
        message("After NA removal")
        message(paste(dim(GSE_matrix)[1], dim(GSE_matrix)[2]))

        # get top MAD genes
        GSE_WGCNA_matrix <- filter_MAD_genes(GSE_matrix)
        print(paste(dim(GSE_WGCNA_matrix)[1], dim(GSE_WGCNA_matrix)[2]))
        
        # choose beta
        message("Finding Soft thresholding power...")
        beta <- pick_beta(GSE_WGCNA_matrix, GSE_n_samples)
        message(message("power chosen= ", beta))
        
        #get modules
        message("finding modules...")
        module_list <- get_modules(GSE_WGCNA_matrix, beta = beta)
        message("WGCNA run successful")
        loginfo(stop_messages[["success"]], logger=paste(GSE_id, GSE_platform, sep = "#"))
        
        # write modules to file
        message("Writing modules to file")
        write_modules(module_list, modules_dir_path = "modules/", gse_id = GSE_id, platform_name = GSE_platform)
        message("modules written")
        write_gse_info(gse_info_file_path = "module-info.txt", gse_id = GSE_id, abstract = GSE_title)
        message("module info written")
    }
    }, finally = {
        if(length(module_list) == 0 & fail ==0){
            logerror(stop_messages[["failure"]], logger=paste(GSE_id, GSE_platform, sep = "#"))
        }
    })
}

