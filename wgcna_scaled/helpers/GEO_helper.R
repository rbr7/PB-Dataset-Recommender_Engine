library(WGCNA)
library(GEOquery)
allowWGCNAThreads(8)

log_GSE_sink <- function(GSE_id, logs_path){  
    logs_path <- paste(logs_path, "/",sep = "")
    GSE_log_file_path <- paste(logs_path, GSE_id, '.log', sep = '')
    if(!dir.exists(logs_path))
        dir.create(logs_path)
    
    print(GSE_log_file_path)
    GSE_sink_file <- file(GSE_log_file_path, open = "w")
    sink(GSE_sink_file, type = "message")
    
    message(paste0("Executing: ", date(), "\n"),GSE_log_file_path)
}

download_geo <- function(geo_file_path, matrix_dir_path="matrix/"){
    geo_object <- getGEO(filename = paste(matrix_dir_path, geo_file_path, sep = ""), getGPL = F)
    return(geo_object)
}

get_GSE_matrix <- function(geo_object){
    return(geo_object@assayData$exprs)
}

get_GSE_metadata <- function(geo_object){
    return(GSE_object@phenoData@data)
}

get_GSE_platform <- function(geo_object){
    return(geo_object@annotation)
}

get_GSE_taxonomy <- function(geo_object){
    return(geo_object@experimentData@other$sample_taxid)
}

get_GSE_n_samples <- function(geo_object){
    return(length(strsplit(geo_object@experimentData@other$sample_id," ")[[1]]))
}

get_GSE_title <- function(geo_object){
    return(geo_object@experimentData@title)
}

download_GPL <- function(platform_name){
    gseurl <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    myurl <- paste(gseurl,'?targ=self&acc=',platform_name,'&form=text&view=full', sep='')
    message("Downloading data for platform ", platform_name)
    download_command <- paste("mkdir soft/; wget -O soft/", platform_name, ".txt ", "\"", myurl, "\"", sep = "")
    system(download_command, intern = T)
    message("Platform data downloaded")
}



annotate_GSE_matrix <- function(GSE_matrix, platform_name, annot_column, annot_type){
    
    source("helpers/GEO_helpers/ensembl_to_entrez.R")
    source("helpers/GEO_helpers/gb_acc_to_entrez.R")
    source("helpers/GEO_helpers/gb_list_to_entrez.R")
    source("helpers/GEO_helpers/symbol_to_entrez.R")
    source("helpers/GEO_helpers/gene_assign_to_entrez.R")
    source("helpers/GEO_helpers/assoc_gene_to_entrez.R")
    
    #download platform annotation soft file using wget
    download_GPL(platform_name)
    
    platform_object <- getGEO(filename = paste("soft/", platform_name, ".txt", sep = ""))
    
    message(platform_object@dataTable@columns)
    # identify entrez gene column
    
    if(annot_type == "entrez"){
        match_positions <- match(rownames(GSE_matrix), platform_object@dataTable@table$ID)
        
        if(sum(!is.na(as.numeric(platform_object@dataTable@table[,annot_column][match_positions])))<100){
            message("entrez column does not have numbers")
            return(NA)
        }
        
        rownames(GSE_matrix) <- platform_object@dataTable@table[,annot_column][match_positions]
        
        GSE_matrix <- GSE_matrix[!is.na(rownames(GSE_matrix)), ]
        GSE_matrix <- GSE_matrix[rownames(GSE_matrix)!="", ]
        
        GSE_matrix_df <- data.frame(GSE_matrix,stringsAsFactors =F) 
        GSE_matrix_df$mad_score <- apply(GSE_matrix_df, 1, function(x) mad(x, na.rm = T))
        GSE_matrix_df$entrez_id <- rownames(GSE_matrix)
        GSE_matrix_df <- GSE_matrix_df[order(GSE_matrix_df$mad_score, decreasing = T),]
        GSE_matrix_df <- GSE_matrix_df[!duplicated(GSE_matrix_df$entrez_id), ]
        rownames(GSE_matrix_df) <- GSE_matrix_df$entrez_id
        GSE_matrix_df <- GSE_matrix_df[, !grepl("entrez_id|mad_score", colnames(GSE_matrix_df))]
        return(as.matrix(GSE_matrix_df))
    }
    
    if(annot_type == "symbol")
        return(symbol_to_entrez(GSE_matrix, platform_object, annot_column))
    
    if(annot_type == "ensembl")
        return(ensembl_to_entrez(GSE_matrix, platform_object, annot_column))
    
    if(annot_type == "gene_ass")
        return(gene_assign_to_entrez(GSE_matrix, platform_object, annot_column))
    
    if(annot_type == "assoc_gene")
        return(assoc_gene_to_entrez(GSE_matrix, platform_object, annot_column))
    
    if(annot_type == "genbank")
        return(gb_acc_to_entrez(GSE_matrix, platform_object, annot_column))
    
}
    
check_log2 <- function(GSE_matrix){
    
    if(abs(range(GSE_matrix)[2]) - abs(range(GSE_matrix)[1]) >100)
        return(FALSE)
    else
        return(TRUE)
}

tranform_log2 <- function(GSE_matrix){
    if(sum(GSE_matrix <= 0) >0){
        message("The dataset has negative values which is weird")
        message("Make negative values and zeros as 1")
        GSE_matrix[GSE_matrix <=0 ] <- 1
    }
    return(log2(GSE_matrix))
}

filter_MAD_genes <- function(GSE_matrix){
    GSE_WGCNA_matrix = t(GSE_matrix[order(apply(GSE_matrix,1,mad), decreasing = T)[1:5000],])
    return(GSE_WGCNA_matrix)
}
