library(shiny)
library(ggplot2)


fluidPage(
  
  titlePanel("MVN curvature"),
  
  sidebarPanel(
    
    sliderInput(
      'variance', label = 'variance:',
      max = 5, 
      min = 0.5,
      value = 1, 
      step = .1
    )
    
  ),
  
  mainPanel(
    plotOutput('plot')
  )
)