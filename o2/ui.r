library(shiny)

# Define UI for application
shinyUI(fluidPage(

  # Application title
  titlePanel("Oxygen Diffusion"),

  sidebarLayout(
    # Culture parameters
    sidebarPanel(
      selectInput("cultureware",
                  "Culture vessel:",
                  list("T75"=75,"T25"=25,"24wp"=1.9,"6wp"=34.8,"96wp"=0.32,"35 mm"=9,"60 mm"=21,"100 mm"=55)),
    # Culture parameters
    
      selectInput("celltype",
                  "Cell type:",
                  list("Primary MEF"=75,"CHO"=88,"Primary hepatocyte"=600)),
    # Culture conditions

      sliderInput("O2_pressure",
                  "Oxygen concentration (%):",
                  min = 0,
                  max = 100,
                  value = 20),
      selectInput("elev",
                  "Location:",
                  list("Vancouver"=1,"Calgary"=0.87),
    # Advanced control
    
      sliderInput("depth",
                  "Media depth (cm):",
                  min = 0.1,
                  max = 1,
                  value = .3),
      sliderInput("Qmet",
                  label=div(HTML("Cell metabolism (amol O<sub>2</sub>/s•cell):")),
                  min = 0,
                  max = 1000,
                  value = 86,
                  step=4)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      textOutput("caption"),
      plotOutput("wellPlot"),
      plotOutput("piePlot")
    )
  )
))

