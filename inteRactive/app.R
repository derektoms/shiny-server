# Load libraries
library(plotly)
library(heatmaply)

# Global variables can go here


# Define the UI
ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            sliderInput('plotValue','Distribution mean',min=30,max=70,value=40)
            ),
        mainPanel(
            plotly::plotlyOutput('plot')
            )
        )
)


# Define the server code
server <- function(input, output) {
    
  output$plot <- plotly::renderPlotly({
    heatmaply::heatmaply(mtcars)
  })
}

# Return a Shiny app object
shinyApp(ui = ui, server = server)