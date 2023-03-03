library(SingleCellExperiment);library(Matrix)

#' Read 10X hdf5 file and turn it into an SCE object
#' 
#' @param filename Path to H5 file
#' @param sample_name Name of the sample, which will be used to uniquify the 
#' colnames
#' @return SingleCellExperimentObject with counts, colData, rowData
#' 
Read10X_h5_2SCE <- function(filename, sample_name){
    
    ## read the counts in------------------------------
    mat <- Read10X_h5(filename, use.names=FALSE, unique.features=FALSE)
    ncells <- dim(mat)[2]
    
    ## colData-----------------------------------------
    cd <- DataFrame(Barcodes = colnames(mat), 
        cell = paste(sample_name, 1:ncells, sep = "_"),
        Sample = sample_name)
    rownames(cd) <- cd$cell
    colnames(mat) <- cd$cell
    
    ## rowData-----------------------------------------
    rd <- Read10X_h5_rowData(filename)
    
    ## SCE---------------------------------------------
    outsce <- SingleCellExperiment(
        assays = list(counts = mat), 
        colData = cd[colnames(mat),],
        rowData = rd[rownames(mat),])
    return(outsce)
}


#' Read 10X hdf5 file
#'
#' Read count matrix from 10X CellRanger hdf5 file.
#' This can be used to read both scATAC-seq and scRNA-seq matrices.
#'
#' @param filename Path to h5 file
#' @param use.names Label row names with feature names rather than ID numbers.
#' @param unique.features Make feature names unique (default TRUE)
#' 
#' @details Originally from Seurat.
#'
#' @return Returns a sparse matrix with rows and columns labeled. If multiple
#' genomes are present, returns a list of sparse matrices (one per genome).
#'
#' @examples \dontrun{
#'  sce <- Read10X_h5_2SCE("JonesLab_IL15J_10Xdata/SampleX/filtered_feature_bc_matrix.h5", "Our_sample_x")
#' }
#' @export
#' @concept preprocessing
#'
Read10X_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
    if (!requireNamespace('hdf5r', quietly = TRUE)) {
        stop("Please install hdf5r to read HDF5 files")
    }
    if (!file.exists(filename)) {
        stop("File not found")
    }
    infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
    genomes <- names(x = infile)
    output <- list()
    if (hdf5r::existsGroup(infile, 'matrix')) {
        # cellranger version 3
        if (use.names) {
            feature_slot <- 'features/name'
        } else {
            feature_slot <- 'features/id'
        }
    } else {
        if (use.names) {
            feature_slot <- 'gene_names'
        } else {
            feature_slot <- 'genes'
        }
    }
    for (genome in genomes) {
        counts <- infile[[paste0(genome, '/data')]]
        indices <- infile[[paste0(genome, '/indices')]]
        indptr <- infile[[paste0(genome, '/indptr')]]
        shp <- infile[[paste0(genome, '/shape')]]
        features <- infile[[paste0(genome, '/', feature_slot)]][]
        barcodes <- infile[[paste0(genome, '/barcodes')]]
        sparse.mat <- sparseMatrix(
            i = indices[] + 1,
            p = indptr[],
            x = as.numeric(x = counts[]),
            dims = shp[],
            repr = "T"
        )
        if (unique.features) {
            features <- make.unique(names = features)
        }
        rownames(x = sparse.mat) <- features
        colnames(x = sparse.mat) <- barcodes[]
        sparse.mat <- as(object = sparse.mat, Class = 'dgCMatrix')
        # Split v3 multimodal
        if (infile$exists(name = paste0(genome, '/features'))) {
            types <- infile[[paste0(genome, '/features/feature_type')]][]
            types.unique <- unique(x = types)
            if (length(x = types.unique) > 1) {
                message(
                    "Genome ",
                    genome,
                    " has multiple modalities, returning a list of matrices for this genome"
                )
                sparse.mat <- sapply(
                    X = types.unique,
                    FUN = function(x) {
                        return(sparse.mat[which(x = types == x), ])
                    },
                    simplify = FALSE,
                    USE.NAMES = TRUE
                )
            }
        }
        output[[genome]] <- sparse.mat
    }
    infile$close_all()
    if (length(x = output) == 1) {
        return(output[[genome]])
    } else{
        return(output)
    }
}


Read10X_h5_rowData <- function(filename) {
    if (!requireNamespace('hdf5r', quietly = TRUE)) {
        stop("Please install hdf5r to read HDF5 files")
    }
    if (!file.exists(filename)) {
        stop("File not found")
    }
    infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
    genomes <- names(x = infile)
    output <- list()
    if (hdf5r::existsGroup(infile, 'matrix')) {
        # cellranger version 3
            feature_slot_name <- 'features/name'
            feature_slot_id <- 'features/id'
            } else {
                feature_slot_name <- 'gene_names'
                feature_slot_id <- 'genes'
            }
    
    for (genome in genomes) {
        features_name <- infile[[paste0(genome, '/', feature_slot_name)]][]
        features_id <- infile[[paste0(genome, '/', feature_slot_id)]][]
        
        rowd <- DataFrame(ID = features_id, Symbol = features_name)
        rownames(rowd) <- features_id ## predicated on  use.names = FALSE for Read10X_hf5!
        
        output[[genome]] <- rowd
        #rownames(x = sparse.mat) <- features
        #colnames(x = sparse.mat) <- barcodes[]
        # Split v3 multimodal
    }
    infile$close_all()
    if (length(x = output) == 1) {
        return(output[[genome]])
    } else{
        return(output)
    }
}
