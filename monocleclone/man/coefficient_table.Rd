% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expr_models.R
\name{coefficient_table}
\alias{coefficient_table}
\title{Extract coefficient table from a fit_models result.}
\usage{
coefficient_table(model_tbl)
}
\arguments{
\item{model_tbl}{A tibble of model objects, generally the output of
\code{\link{fit_models}}.}
}
\value{
A table of coefficient data for each gene.
}
\description{
Extracts a table of coefficients from a tibble
containing model objects. It tests whether each coefficient
differs significantly from zero under the Wald test and
adjusts the p-values for multiple hypothesis testing using
the method of Benjamini and Hochberg, placing these
adjusted values in the q-value column.
}
\examples{
  \donttest{
    cell_metadata <- readRDS(system.file('extdata',
                                         'worm_embryo/worm_embryo_coldata.rds',
                                         package='monocle3'))
    gene_metadata <- readRDS(system.file('extdata',
                                         'worm_embryo/worm_embryo_rowdata.rds',
                                         package='monocle3'))
    expression_matrix <- readRDS(system.file('extdata',
                                             'worm_embryo/worm_embryo_expression_matrix.rds',
                                             package='monocle3'))
   
    cds <- new_cell_data_set(expression_data=expression_matrix,
                             cell_metadata=cell_metadata,
                             gene_metadata=gene_metadata)

    cds <- preprocess_cds(cds, num_dim=50)
    cds <- align_cds(cds, alignment_group = "batch",
                     residual_model_formula_str = "~ bg.300.loading + bg.400.loading +
                     bg.500.1.loading + bg.500.2.loading + bg.r17.loading +
                     bg.b01.loading + bg.b02.loading")
    cds <- reduce_dimension(cds)
    ciliated_genes <- c("che-1", "hlh-17", "nhr-6", "dmd-6", "ceh-36", "ham-1")
    cds_subset <- cds[rowData(cds)$gene_short_name \%in\% ciliated_genes,]
    gene_fits <- fit_models(cds_subset, model_formula_str = "~embryo.time")
    fit_coefs <- coefficient_table(gene_fits)
  }

}
