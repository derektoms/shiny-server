########################################
#$#$#$#$#$#    HEADER      #$$#$#$#$#$#$
########################################
## last update: 29-09-2022
## Initial start for a bioreactor data visualization app

## Libraries


require(shiny)

########################################
#$#$#$#$#$#  USER INTERFACE   #$#$#$#$#$
########################################

## Define UI for application
shinyUI(fluidPage(

## Application title
  titlePanel("Myo Palate Bioreactor Data Visualization"),

## _^x"(_^x"(_^x"(_^x"(_^x"(_^x"(_^x"(_^x"(_^x"(_^x"(
  sidebarLayout(
## Data parameters
#### Next version, these need to be generated server-side based on the data table (e.g., factor levels of "cell.culture" and column names)
##
    sidebarPanel(
    selectInput("xaxis",
                label="Which parameter to display on the X axis?",
                choices = datacols,
                 selected = "duration"),
    selectInput("yaxis",
                  label="Which parameter to display on the Y axis?",
                  choices = datacols,
                  selected = "cellcount"),
    selectInput("cc3d",
                  label="Which 3D cell culture experiment to display?",
                  choices = list("CC3D015"='cc3d015', "CC3D016"='cc3d016', "CC3D017"='cc3d017'),
                  multiple = TRUE,
                  selected = c("cc3d015", "cc3d016", "cc3d017"))
      ),

## _^x"(_^x"(_^x"(_^x"(_^x"(_^x"(_^x"(_^x"(_^x"(_^x"(    
## Plotting
    mainPanel(
        plotOutput("scatterplot"),
        dataTableOutput("data2")
    )
  )
))