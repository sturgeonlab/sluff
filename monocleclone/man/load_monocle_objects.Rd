% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{load_monocle_objects}
\alias{load_monocle_objects}
\title{Load a full Monocle3 cell_data_set.}
\usage{
load_monocle_objects(directory_path)
}
\arguments{
\item{directory_path}{a string giving the name of the directory
from which to read the saved cell_data_set files.}
}
\value{
a cell_data_set.
}
\description{
Load a full Monocle3 cell_data_set, which was saved using
save_monocle_objects. For more information read the help
information for save_monocle_objects.
}
\examples{
  \dontrun{
    cds <- load_a549()
    save_monocle_objects(cds, 'mo')
    cds1 <- load_monocle_objects('mo')
  }

}
