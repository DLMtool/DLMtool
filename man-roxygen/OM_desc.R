#' @slot OM A table of sampled parameter of the operating model. Table object of nsim rows. Real numbers\cr
#'   \itemize{
#'   <% for (x in 1:nrow(OM_desc)) { %>
#'   <%   cat(paste0("\\item ", OM_desc$Variable[x], ": ", OM_desc$Description[x])) %>
#'   <% } %>
#'   }
#'   
#'   

