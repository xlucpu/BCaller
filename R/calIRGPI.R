#' @name calIRGPI
#' @title Calculate immune-related genes pair index (IRGPI) for bladder cancer
#' @description This function calculates an IRGPI risk score individually for bladder cancer only.
#' @param expr A numerical expression matrix or data frame with row for gene symbol name and column for sample ID. Note: In principle, the expression profile does not need any normalization, but since the amount of gene expression is affected by the gene length, the original count or the normalized count data may not be suitable for this analysis. It is recommended to provide FPKM or TPM value and log transformation is not necessary.
#' @param res.path A string value to indicate the output path for storing the Estimated IRGPI.txt file.
#'
#' @return A data frame stored sample name and estimated IRGPI.
#' @export
#' @references Prognosis stratification and personalized treatment in bladder cancer through a robust immune gene-pair based signature
#' @author Xiaofan Lu
#' @examples
#' library(BCaller)
#' load(system.file("extdata", "demo.RData", package = "BCaller", mustWork = TRUE)) # load example data
#' IRGPI  <- calIRGPI(expr = demo)
#' head(IRGPI)
calIRGPI <- function(expr     = NULL,
                     res.path = getwd()) {

  # genes in IRGPI1
  IRGP1 <- c("AMMECR1L",
             "BIRC5",
             "G6PD",
             "BTK",
             "C5",
             "EGR2",
             "LY9",
             "CXCL1",
             "HLA-DPB1",
             "CLEC4A",
             "EGR1",
             "ITGA5",
             "CCL11",
             "ITGA1",
             "AXL",
             "NUP107",
             "SYT17",
             "SAA1",
             "AMMECR1L",
             "CARD11",
             "CD96",
             "CD3G",
             "CCL17",
             "AMMECR1L",
             "JAK2",
             "TLR4",
             "INPP5D",
             "MEFV",
             "ICAM2")

  # genes in IRGPI2
  IRGP2 <- c("MST1R",
             "ST6GAL1",
             "SYT17",
             "CD247",
             "CD96",
             "IL15",
             "MAGEC2",
             "IL1A",
             "RRAD",
             "EOMES",
             "NFKB1",
             "MR1",
             "MARCO",
             "TLR4",
             "CEACAM1",
             "THY1",
             "THY1",
             "TGFB1",
             "MAP2K1",
             "IL1R1",
             "CNOT4",
             "DPP4",
             "PBK",
             "IL4R",
             "PDGFC",
             "VEGFC",
             "PVR",
             "SERPINB2",
             "THY1")

  # coefficient for each IRGP
  coeff <- c(-0.278125761,
             -0.240534714,
             -0.231481299,
             -0.219619661,
             -0.14952529,
             -0.134399494,
             -0.123209559,
             -0.109416714,
             -0.107499057,
             -0.103013301,
             -0.088255513,
             -0.078881388,
             -0.040480343,
             -0.036568213,
             -0.026672769,
             0.007069121,
             0.030336853,
             0.066360407,
             0.132744554,
             0.136897427,
             0.159178132,
             0.18447765,
             0.211337792,
             0.233469355,
             0.261049599,
             0.327501543,
             0.33704137,
             0.397375749,
             0.46981137)

  # initial check
  if(max(expr) >= 25){
    message("--please make sure a properly normalized expression data has been provided (e.g., FPKM or TPM); count data is not suitable because it does not consider gene length.")
  }

  if(!all(is.element(unique(c(IRGP1,IRGP2)), rownames(expr)))) {
    missgene <- setdiff(unique(c(IRGP1,IRGP2)), rownames(expr))
    stop(paste0(length(missgene)," genes cannot be mapped in your data, please check for the following missing feature(s):\n",
         paste(missgene, collapse = "\n")))
  }

  # extract expression
  IRGP1.expr <- expr[IRGP1,]
  IRGP2.expr <- expr[IRGP2,]

  # calculate difference and convert to binary matrix
  IRGP.diff <- IRGP1.expr - IRGP2.expr
  IRGP.binary <- ifelse(IRGP.diff < 0, 1, 0)

  # estimate IRGPI
  IRGPI <- apply(t(IRGP.binary),1,function(x) {x %*% coeff})

  # output
  outTab <- data.frame(SampleID = names(IRGPI),
                       IRGPI = as.numeric(IRGPI),
                       stringsAsFactors = F)
  write.table(outTab,
              file = file.path(res.path,"Estimated IRGPI.txt"),
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE)

  return(outTab)
}
