#' @importFrom abind abind
#' @importFrom devtools install_github
#' @importFrom dplyr  %>%  filter group_by mutate select summarize
#' @importFrom fmsb radarchart
#' @importFrom ggplot2 aes element_blank expand_limits facet_wrap geom_boxplot ggplot ggplotGrob geom_rect geom_point labs theme theme_classic xlim ylim xlab ylab   
#' @importFrom ggrepel geom_text_repel
#' @importFrom graphics abline arrows axis axTicks barplot boxplot contour hist identify layout legend
#'  lines matplot mtext par plot plot.new points polygon segments text title text
#' @importFrom grDevices col2rgb colorRampPalette dev.off jpeg rainbow rgb xy.coords
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid unit.c unit grid.newpage grid.draw
#' @importFrom kableExtra cell_spec kable_styling column_spec add_header_above 
#' @importFrom knitr kable
#' @importFrom MASS mvrnorm kde2d
#' @importFrom methods getClassDef getSlots .hasSlot new show slot slot<- slotNames
#' @importFrom mvtnorm rmvnorm
#' @importFrom parallel detectCores 
#' @importFrom snowfall sfClusterEval sfInit sfExportAll sfIsRunning sfExport sfSapply sfLibrary
#' @importFrom stats acf approx coef dbeta density dnorm dlnorm lm loess loess.smooth
#' median nlm optim optimise optimize plogis pnorm predict qlnorm quantile rbeta
#' rlnorm rmultinom rnorm runif sd
#' @importFrom utils  browseURL capture.output combn flush.console packageVersion ls.str lsf.str read.csv read.csv2
#' 
NULL


.onUnload <- function (libpath) {
  library.dynam.unload("DLMtool", libpath)
}

