% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{save_monocle_objects}
\alias{save_monocle_objects}
\title{Save a Monocle3 full cell_data_set.}
\usage{
save_monocle_objects(
  cds,
  directory_path,
  hdf5_assays = FALSE,
  comment = "",
  verbose = TRUE
)
}
\arguments{
\item{cds}{a cell_data_set to save.}

\item{directory_path}{a string giving the name of the directory
in which to write the object files.}

\item{hdf5_assays}{a boolean determining whether the
non-HDF5Array assay objects are saved as HDF5 files. At this
time cell_data_set HDF5Array assay objects are stored as
HDF5Assay files regardless of the hdf5_assays parameter value.}

\item{comment}{a string with optional notes that is saved with
the objects.}

\item{verbose}{a boolean determining whether to print information
about the saved files.}
}
\value{
none.
}
\description{
Save a Monocle3 full cell_data_set to a specified directory
by writing the R objects to RDS files and the nearest
neighbor indexes to index files. The assays
objects are saved as HDF5Array files when hdf5_assays=TRUE
or when the cell_data_set assays are HDF5Array objects. If
any assay in the cell_data set is an HDF5 object, all assays
must be. When save_monocle_objects is run with hdf5_assays=TRUE,
the load_monocle_objects function loads the saved assays into
HDF5Array objects in the resulting cell_data_set. Note:
operations such as preprocess_cds that are run on assays stored
as HDF5Arrays are much, much slower than the same operations
run on assays stored as in-memory matrices. You may want to
investigate parameters related to the Bioconductor DelayedArray
and BiocParallel packages in this case.
}
\examples{
  \dontrun{
    cds <- load_a549()
    save_monocle_objects(cds, 'mo')
  }

}
