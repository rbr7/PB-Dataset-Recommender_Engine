symbol_to_entrez <- function(GSE_matrix, platform_object, annot_column){
  match_positions <- match(rownames(GSE_matrix), platform_object@dataTable@table$ID)
  rownames(GSE_matrix) <- platform_object@dataTable@table[, annot_column][match_positions]
  GSE_matrix <- GSE_matrix[!is.na(rownames(GSE_matrix)), ]
  GSE_matrix <- GSE_matrix[rownames(GSE_matrix)!="", ] 
  
  ensembl_to_entrez_df <- read.table("helpers/annotation_files/symbol-to-entrez.tsv",
                                     sep = '\t', stringsAsFactors = F)
  
  m <- match(rownames(GSE_matrix), ensembl_to_entrez_df[,3])
  
  GSE_matrix_df <- data.frame(GSE_matrix,stringsAsFactors =F) 
  GSE_matrix_df$entrez_id <- ensembl_to_entrez_df[m, 2]
  GSE_matrix_df <- GSE_matrix_df[!is.na(GSE_matrix_df$entrez_id),]
  GSE_matrix_df$mad_score <- apply(GSE_matrix_df[, !grepl("entrez_id", colnames(GSE_matrix_df))], 
                                   1, function(x) mad(x, na.rm = T))
  
  GSE_matrix_df <- GSE_matrix_df[order(GSE_matrix_df$mad_score, decreasing = T),]
  GSE_matrix_df <- GSE_matrix_df[!duplicated(GSE_matrix_df$entrez_id), ]
  
  
  rownames(GSE_matrix_df) <- GSE_matrix_df$entrez_id
  GSE_matrix_df <- GSE_matrix_df[, !grepl("entrez_id|mad_score", colnames(GSE_matrix_df))]
  return(as.matrix(GSE_matrix_df))
}
