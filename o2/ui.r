library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

  # Application title
  titlePanel("Oxygen Diffusion"),

  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("depth",
                  "Media depth (cm):",
                  min = 0.1,
                  max = 1,
                  value = .3),
      sliderInput("O2_pressure",
                  "Oxygen concentration (%):",
                  min = 1,
                  max = 100,
                  value = 21),
      sliderInput("Qmet",
                  label=div(HTML("Cell metabolism (amol O<sub>2</sub>/s•cell):")),
                  min = 0,
                  max = 1000,
                  value = 86,
                  step=4),
      selectInput("elev",
                  "Location:",
                  list("Vancouver"=1,"Calgary"=0.87))
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      textOutput("caption"),
      plotOutput("wellPlot"),
      plotOutput("piePlot")
    )
  )
))

