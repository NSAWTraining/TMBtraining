library(shiny)
library(ggplot2)
library(mvtnorm)

function(input, output) {
  
  output$plot <- renderPlot({
    
    cor2cov <- function(R,sd2){
      S <- c(sd2, sd2)
      diag(S) %*% R %*% diag(S)
    }
    
    C <- cbind(c(1, .2), c(.2, 1))
    Sigma <- cor2cov(C, input$variance)
    x <- y <- seq(-3, 3,.5)
    z <- matrix(0,length(x),length(y))
    for(i in 1:length(x)){
      for(j in 1:length(y)){
        z[i,j] <- dmvnorm(c(x[i],y[j]), sigma = Sigma)
      }
    }
    colnames(Sigma) <- c("Sigma", " ")
    H <- -solve(Sigma)
    colnames(H) <- c("H", " ")
    
    persp(x, y, -log(z), zlim = range(0,12))
    output$table1 <- renderTable(round(Sigma, 3))
    output$table2 <- renderTable(round(H,3))
    
  })
}