run_module_pipeline_MA <- function(GSE_id, module_file_path, gse_info_file_path, logs_path) {
  source("helpers/GEO_helper.R")
  source("helpers/WGCNA_helper.R")

  stop_messages <- list(
    "low_samples" =
      "less than 12 samples present.",
    "low_genes" = "No rows in the matrix,GEOQuery did not get data properly",
    "NA_genes" = "All rows contain columns with NAs",
    "platform_annotation" = "Not able to annotate platform",
    "success" = "WGCNA successfully completed!!!",
    "failure" = "failed",
    "platform_not_MA" = "platform not a Microrarray platform or not in annotations file",
    "annot_not_defined" = "Annotation is not defined in platform annotation file",
    "entrez_col_no_digit" = "Entrez column does not contain integers in platform annotation file",
    "low_genes_WGCNA" = "Very less number of genes in matrix for WGCNA to run"
  )

  log_GSE_sink(GSE_id = GSE_id, logs_path = logs_path)
  platform_annot_table <- read.csv("helpers/MA_platform_annotation.csv", stringsAsFactors = F)[, c(1, 3, 4)]

  for (platform_index in 1:length(list.files("matrix/"))) {
    tryCatch({
      fail <- 0

      GSE_object <- download_geo(list.files("matrix/")[[platform_index]])
      module_list <- list()

      message("checking number of samples")
      GSE_platform <- get_GSE_platform(GSE_object)
      if (!(GSE_platform %in% platform_annot_table[, 1])) {
        logerror(stop_messages[["platform_not_MA"]], logger = paste(GSE_id, GSE_platform, sep = "#"))
        fail <- 1
        next
      }

      GSE_n_samples <- get_GSE_n_samples(GSE_object)
      message(paste("no. of samples", GSE_n_samples))
      if (GSE_n_samples < 12) {
        logerror(stop_messages[["low_samples"]], logger = paste(GSE_id, GSE_platform, sep = "#"))
        fail <- 1
        next
      }

      message("Getting data for platform----- ", platform_index)
      GSE_matrix <- get_GSE_matrix(GSE_object)
      print(dim(GSE_matrix))
      if (nrow(GSE_matrix) < 100) {
        logerror(stop_messages[["low_genes"]], logger = paste(GSE_id, GSE_platform, sep = "#"))
        fail <- 1
        next
      }

	GSE_n_samples <- ncol(GSE_object)
        message(paste("actual no. of samples", GSE_n_samples))
            if (GSE_n_samples < 12) {
		  logerror(stop_messages[["low_samples"]], logger = paste(GSE_id, GSE_platform, sep = "#"))
	          fail <- 1
		  next
	    }





      message(paste("GSE_platform ------", GSE_platform))
      GSE_tax_id <- get_GSE_taxonomy(GSE_object)
      GSE_title <- get_GSE_title(GSE_object)

      # annotate with entrez ID's GSE_matrix



      ## see which column to use from platform file
      annot_column <- platform_annot_table[platform_annot_table[, 1] == GSE_platform, 2]
      annot_type <- platform_annot_table[platform_annot_table[, 1] == GSE_platform, 3]

      if (nchar(annot_type) == 0) {
        message("annot not defined")
        fail <- 1
        logerror(stop_messages[["annot_not_defined"]], logger = paste(GSE_id, GSE_platform, sep = "#"))
        next
      }

      message(paste("column for annotation------", annot_column))
      message(paste("annotation type------", annot_type))

      GSE_matrix <- annotate_GSE_matrix(GSE_matrix, GSE_platform, annot_column, annot_type)

      if (nrow(GSE_matrix) == 0) {
        fail <- 1
        message("annotation not done due to entrez column error")
        logerror(stop_messages[["entrez_col_no_digit"]], logger = paste(GSE_id, GSE_platform, sep = "#"))
        next
      }

      message("GSE matrix dims after annotation")
      message(paste(dim(GSE_matrix)[1], dim(GSE_matrix)[2]))

      GSE_matrix <- GSE_matrix[!apply(GSE_matrix, 1, anyNA), ]
      print(dim(GSE_matrix))
      if (nrow(GSE_matrix) == 0) {
        logerror(stop_messages[["NA_genes"]], logger = paste(GSE_id, GSE_platform, sep = "#"))
        fail <- 1
        next
      }

      message("After NA removal")
      message(paste(dim(GSE_matrix)[1], dim(GSE_matrix)[2]))


      log2_tranformed <- check_log2(GSE_matrix)

      if (!log2_tranformed) {
        message("matrix is not log2 tranformed..so transforming it")
        GSE_matrix <- tranform_log2(GSE_matrix)
      }
      # get top MAD genes
      GSE_WGCNA_matrix <- filter_MAD_genes(GSE_matrix)
      print(paste(dim(GSE_WGCNA_matrix)[1], dim(GSE_WGCNA_matrix)[2]))


      GSE_WGCNA_matrix <- GSE_WGCNA_matrix[, !apply(GSE_WGCNA_matrix, 2, anyNA)]


      if (ncol(GSE_WGCNA_matrix) < 2000) {
        logerror(stop_messages[["low_genes_WGCNA"]], logger = paste(GSE_id, GSE_platform, sep = "#"))
        fail <- 1
        next
      }


      # choose beta
      message("Finding Soft thresholding power...")
      beta <- pick_beta(GSE_WGCNA_matrix, GSE_n_samples)
      message(message("power chosen= ", beta))

      # get modules
      message("finding modules...")
      module_list <- get_modules(GSE_WGCNA_matrix, beta = beta)
      message("WGCNA run successful")
      loginfo(stop_messages[["success"]], logger = paste(GSE_id, GSE_platform, sep = "#"))

      # write modules to file
      message("Writing modules to file")
      write_modules(module_list, modules_dir_path = "modules/", gse_id = GSE_id, platform_name = GSE_platform)
      message("modules written")
      write_gse_info(gse_info_file_path = "module-info.txt", gse_id = GSE_id, abstract = GSE_title)
      message("module info written")
    },
    finally = {
      if (length(module_list) == 0 & fail == 0) {
        logerror(stop_messages[["failure"]], logger = paste(GSE_id, GSE_platform, sep = "#"))
        next
      }
    }
    )
  }
  gc()
}
