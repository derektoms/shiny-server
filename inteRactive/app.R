# Load libraries
library(plotly)

# Global variables can go here
m <- 200


# Define the UI
ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            sliderInput('plotValue','Distribution mean',min=30,max=70,value=40)
            ),
        mainPanel(
            plotlyOutput('plot')
            )
        )
)


# Define the server code
server <- function(input, output) {
    meantoplot<- reactive({
        input$plotValue
    })
    
  output$plot <- renderPlotly({
    plot_ly(y = ~rnorm(meantoplot()),type="box")
  })
}

# Return a Shiny app object
shinyApp(ui = ui, server = server)