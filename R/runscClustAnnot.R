#' Launch Shiny App for scClustAnnot
#'
#' A function that launches the Shiny app for scClustAnnot.
#' The code has been placed in \code{./inst/shiny-scripts}.
#'
#' @return No return value but open up a Shiny page.
#'
#' @examples
#' \dontrun{
#'
#' scClustAnnot::runscClustAnnot()
#' }
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials. \href{https://shiny.rstudio.com/tutorial/}{Link}
#'
#' @export
#' @import shiny

runscClustAnnot <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "scClustAnnot")
  actionShiny <- shiny::runApp(appDir, display.mode = "normal")

  return(actionShiny)
}
# [END]
