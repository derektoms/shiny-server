########################################
#$#$#$#$#$#    HEADER      #$$#$#$#$#$#$
########################################
## last update: 29-09-2022
## Initial start for a bioreactor data visualization app

## Libraries

require(shiny)
require(ggplot2)
require(dplyr)
require(tidyr)
require(magrittr)
require(RColorBrewer)

########################################
#$#$#$#$#$#    SERVER      #$$#$#$#$#$#$
########################################

## Define server logic
shinyServer(function(input, output, session) {

data <- read.csv('data/toy_data.csv')
meta <- read.csv('data/toy_meta.csv')

data$duration <- as.numeric(data$duration)
datacols <- colnames(data[,c(-1:-2)])

data2 <- reactive({
    data %>% filter(data$cell.culture %in% input$cc3d)
})

mycols <- reactive({
    brewer.pal(length(input$cc3d),"Set1")
})

output$scatterplot <- renderPlot({
       ggplot(data=data2(), aes_string(x=input$xaxis, y=input$yaxis)) + geom_point(aes(colour=cell.culture), size=4) + geom_line(aes(colour=cell.culture)) + theme_minimal() + scale_colour_manual(values=mycols())
})

output$data2 <- renderDataTable  ({data2()})

## Kill shinyApp when session closes
session$onSessionEnded(stopApp)
})




## _^x"(_^x"(_^x"(_^x"(_^x"(_^x"(_^x"(_^x"(_^x"(_^x"(
## END
