# app.R — entry point, run this file to launch the dashboard
# setwd to project root first, then: shiny::runApp("app")

library(shiny)
source("ui.R")
source("server.R")
shinyApp(ui = ui, server = server)
