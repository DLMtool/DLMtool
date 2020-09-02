#' @importFrom abind abind
#' @importFrom dplyr %>%  arrange filter group_by left_join mutate select summarize
#' @importFrom ggplot2 aes element_blank expand_limits facet_wrap geom_boxplot ggplot ggplotGrob geom_rect geom_point labs theme theme_classic xlim ylim xlab ylab   
#' @importFrom graphics abline arrows axis axTicks barplot boxplot contour hist identify layout legend
#'  lines matplot mtext par plot plot.new points polygon segments text title text
#' @importFrom grDevices col2rgb colorRampPalette dev.off jpeg rainbow rgb xy.coords
#' @importFrom grid unit.c unit grid.newpage grid.draw
#' @importFrom methods getClassDef getSlots .hasSlot new show slot slot<- slotNames
#' @importFrom parallel detectCores 
#' @importFrom snowfall sfClusterEval sfInit sfExportAll sfIsRunning sfExport sfSapply sfLibrary
#' @importFrom stats acf approx coef cor dbeta density dnorm dlnorm lm loess loess.smooth nls setNames SSasympOff
#' median nlm optim optimise optimize plogis pnorm predict qlnorm quantile rbeta
#' rlnorm rmultinom rnorm runif sd
#' @importFrom utils  browseURL capture.output combn flush.console packageVersion ls.str lsf.str read.csv read.csv2
#' 
NULL


.onUnload <- function (libpath) {
  library.dynam.unload("DLMtool", libpath)
}

# global variable names
Names <- c("maxage", "R0", "Mexp", "Msd", "dep", "D", "Mgrad", "SRrel", "hs", "procsd",
           "L50", "L95", "L50_95", "CAL_binsmid", "Len_age", "maxlen", "Linf", 
           "M_at_Length", "Frac_area_1", "Prob_staying", "M_ageArray", "Mat_age",
           "Wt_age", "V", "Spat_targ", "procmu", "recMulti", "Linfrand", "Krand",
           "Abias Aerr", "Brefbias", "CAA_ESS", "CAA_nsamp", "CAL_ESS", "CAL_bins", "CAL_nsamp",
           "Cbias", "Crefbias", "Csd", "Dbias", "Derr", "TAEFrac", "TAESD", "EffLower",
           "EffUpper", "EffYears", "FMSY_Mbias", "Frac_area_1", "Irefbias", "Isd", "K", "Kbias", "Kgrad",
           "Krand", "Ksd", "L5", "L5s", "LFCbias", "LFS", "LFSbias", "LFSs", "LatASD", "Linfbias", "Linfgrad",
           "Linfrand", "Linfsd", "M", "M_ageArray", "Mat_age", "Mbias", "Mrand", "Prob_staying", "Recsd",
           "SLarray", "SizeLimFrac", "SizeLimSD", "Spat_targ", "TACFrac", "TACSD", 
           "Vmaxlen", "Vmaxlens", "Wt_age", "ageM", "betas", "lenMbias", "nCALbins", "procmu", "qcv", "qinc",
           "recMulti",  "t0", "t0bias", "Abias", "Aerr", "Perr", "Esd", "qvar", "Marray",
           "Linfarray", "Karray", "t0array", "mov",  "nareas", "AC", "LenCV", "a", "b", "FinF", 
           "Fdisc", "R50", "Rslope", "retA", "retL", "LR5", "LFR", "Rmaxlen",
           "V2", "SLarray2", "DR", "Asize", "Size_area_1", "L50array", "L95array",
           "Fdisc_array", "Fdisc_array2", "Pinitdist", "DataOut",
           'Perr_y', "Cobs", "Iobs", "Dobs", "Btbiascv", 'Btobs', "h", 'Index',
           '.', 'MP', 'Data', 'DataClass', "Type", "Recs", "DominatedMPs"
)

if(getRversion() >= "2.15.1") utils::globalVariables(Names)


# change messages to blue text instead of default red
message <- function(...) {
  if (requireNamespace("crayon", quietly = TRUE)) {
    x <- paste(...)
    return(base::message(crayon::blue(paste(base::strwrap(x), collapes="\n"))))
  } else {
    return(base::message(...))
  }
}