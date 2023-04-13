
.onLoad <- function(...){
  requireNamespace("ggplot2")
  requireNamespace("compiler")
  compiler::enableJIT(3)
  return(NULL)
}
