require(shiny)

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
      textInput("volume",
                "Media volume (ml):",
                value=15),
      selectInput("celltype",
                  "Cell type:",
                  list("MEF"=33,
                  "CHO"=88,
                  "HeLa"=12.5,
                  "Hepatocytes (rat)"=325,
                  "Erythrocytes (rabbit)"=0.02,
                  "Jurkat"=12,
                  "ESC (mouse)"=40
                  )),
      hr(),
      # Culture conditions
      sliderInput("O2_pressure",
                  "Oxygen concentration (%):",
                  min = 0,
                  max = 100,
                  value = 21),
      selectInput("elev",
                  "Location:",
                  list("Vancouver"=1,"Calgary"=0.87)),
      hr(),
      # Advanced control
      checkboxInput("advanced", "Advanced",TRUE),
      conditionalPanel(
        condition = "input.advanced == true",
        uiOutput("advanced")
    )),
    
    # Show a plot of the generated distribution
    mainPanel(
      column(6,
             verbatimTextOutput("oxy"),
             verbatimTextOutput("depth"),
             verbatimTextOutput("atmo"),
             verbatimTextOutput("Qcell"),
             plotOutput("piePlot")
      ),
      column(6,
             plotOutput("wellPlot")
      )
      
    )
  )
))
