#' Ligand-Receptor Table from OmniPath
#' 
#' A data frame of ligand-receptor pairs. Each monomeric ligand/receptor is represented by its corresponding gene symbol. 
#' Each Multimeric ligand/receptor is represented by concatenating the gene symbols of all constituent subunits using underscores (`_`).
#'
#' @format ## `omnipathr`
#' a data frame with 1946 rows and 2 columns: `ligand` and `receptor`.
#' 
#' @source the dataset was obtained using the \pkg{OmnipathR} package (version 3.14.0).
#' The package is available on Bioconductor: \url{https://www.bioconductor.org/packages/release/bioc/html/OmnipathR.html}.
#' 
#' @references D Turei, A Valdeolivas, L Gul, N Palacio-Escat, M Klein, O Ivanova, M Olbei, A Gabor, F
#' Theis, D Modos, T Korcsmaros and J Saez-Rodriguez (2021) Integrated intra- and
#' intercellular signaling knowledge for multicellular omics analysis. Molecular Systems
#' Biology 17: e9923; DOI: \url{https://doi.org/10.15252/msb.20209923}
#' 
#' D Turei, T Korcsmaros and J Saez-Rodriguez (2016) OmniPath: guidelines and gateway for
#' literature-curated signaling pathway resources. Nature Methods 13 (12); PMID: 27898060
#' 
#' A Valdeolivas, D Turei, A Gabor (2019) “OmnipathR: client for the OmniPath web service.”
#' Bioconductor Package
#' 
"omnipathr"