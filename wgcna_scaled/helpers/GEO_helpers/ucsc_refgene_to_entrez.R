ucsc_refgene_to_entrez <- function(GSE_matrix, platform_object){
  match_positions <- match(rownames(GSE_matrix), platform_object@dataTable@table$ID)
  rownames(GSE_matrix) <- platform_object@dataTable@table[,"UCSC_RefGene_Accession"][match_positions]
  GSE_matrix <- GSE_matrix[!is.na(rownames(GSE_matrix)), ]
  GSE_matrix <- GSE_matrix[rownames(GSE_matrix)!="", ] 
  
  refseq_to_entrez_df <- read.table('~/genequery_rnaseq_integration/gqcmd-internal_new/refseq-to-entrez.tsv',
                                    sep = '\t', stringsAsFactors = F)
  
  GSE_matrix_df <- data.frame(GSE_matrix,stringsAsFactors =F)
  
  print(dim(GSE_matrix))
  GSE_matrix_df$entrez_id <- sapply(rownames(GSE_matrix), function(gc_acc_comma){
    gb_acc_vec <- strsplit(gc_acc_comma, split = ";")[[1]]
    convert_df <- refseq_to_entrez_df[refseq_to_entrez_df[,3] %in% gb_acc_vec,]
    
    if(nrow(convert_df)==0)
      return(NA)
    else
      return(paste(convert_df[,2], collpase = ","))
  })
  
  GSE_matrix_df <- GSE_matrix[!is.na(GSE_matrix_df$entrez_id), ]
  GSE_matrix_df$mad_score <- apply(GSE_matrix_df[, !grepl("entrez_id", colnames(GSE_matrix_df))], 
                                   1, function(x) mad(x, na.rm = T))
  GSE_matrix_df <- GSE_matrix_df[order(GSE_matrix_df$mad_score, decreasing = T),]
  GSE_matrix_df <- GSE_matrix_df[!duplicated(GSE_matrix_df$entrez_id), ]
  rownames(GSE_matrix_df) <- GSE_matrix_df$entrez_id
  GSE_matrix_df <- GSE_matrix_df[, !grepl("entrez_id|mad_score", colnames(GSE_matrix_df))]
  return(as.matrix(GSE_matrix_df))
}