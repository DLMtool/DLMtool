#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the MP recommendation(s)
#' @param plot Logical. Show the plot?
#' 
#' @section Required Data:
#' See \linkS4class{Data} for information on the \code{Data} object \cr 
#' 
#' @return
#'  <%= if(MPtype(mp)[2] == 'Output') paste0("An object of class \\code{\\link[=Rec-class]{Rec}} with the \\code{TAC} slot populated with a numeric vector of length \\code{reps}") %>
#'  <%= if(MPtype(mp)[2] == 'Input') paste0("An object of class \\code{\\link[=Rec-class]{Rec}} with the ", paste(MPtype(mp)[3], collapse=','), " slot(s) populated") %>
#'  <%= if(MPtype(mp)[2] == 'Mixed') paste0("An object of class \\code{\\link[=Rec-class]{Rec}} with the ", paste(MPtype(mp)[3], collapse=','), " slot(s) populated") %>
#' 




