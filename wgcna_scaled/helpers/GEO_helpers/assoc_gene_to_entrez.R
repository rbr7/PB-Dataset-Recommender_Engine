#6th element of 'Associated Gene' column is the entrez id
assoc_gene_to_entrez <- function(GSE_matrix, platform_object, annot_column){
  match_positions <- match(rownames(GSE_matrix), platform_object@dataTable@table$ID)
  
  GSE_matrix_df <- data.frame(GSE_matrix,stringsAsFactors =F) 
  GSE_matrix_df$associated_gene <- platform_object@dataTable@table[, annot_column][match_positions]
  
  GSE_matrix_df <- GSE_matrix_df[!is.na(GSE_matrix_df$associated_gene), ]
  GSE_matrix_df <- GSE_matrix_df[GSE_matrix_df$associated_gene!="", ] 
  
  GSE_matrix_df$entrez_id <- sapply(GSE_matrix_df$associated_gene, function(x){
    id_vec <- unlist(sapply(strsplit(x, split = " /// ")[[1]], function(x) strsplit(x, split = ' // ')[[1]][6]))
    return(paste(unique(id_vec[id_vec!='---']), collapse = ','))
  })
  GSE_matrix_df <- GSE_matrix_df[!is.na(GSE_matrix_df$entrez_id), ]
  print(dim(GSE_matrix_df))                                      
  
  GSE_matrix_df <- GSE_matrix_df[GSE_matrix_df$entrez_id != "", ]
  
  GSE_matrix_df$mad_score <- apply(GSE_matrix_df[, !grepl("entrez_id|associated_gene", colnames(GSE_matrix_df))], 
                                   1, function(x) mad(x, na.rm = T))
  
  GSE_matrix_df <- GSE_matrix_df[order(GSE_matrix_df$mad_score, decreasing = T),]
  GSE_matrix_df <- GSE_matrix_df[!duplicated(GSE_matrix_df$entrez_id), ]
  
  rownames(GSE_matrix_df) <- GSE_matrix_df$entrez_id
  GSE_matrix_df <- GSE_matrix_df[, !grepl("entrez_id|mad_score|associated_gene", colnames(GSE_matrix_df))]
  return(as.matrix(GSE_matrix_df))
}