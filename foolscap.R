> filter(gpl,gpl=='GPL1261') %>% select(gpl,title)
# Source:   lazy query [?? x 2]
# Database: sqlite 3.19.3 [/Volumes/ULTRA/across_array/GEOmetadb.sqlite]
  gpl     title                                             
  <chr>   <chr>                                             
1 GPL1261 [Mouse430_2] Affymetrix Mouse Genome 430 2.0 Array
> filter(gpl,gpl=='GPL570') %>% select(gpl,title)
# Source:   lazy query [?? x 2]
# Database: sqlite 3.19.3 [/Volumes/ULTRA/across_array/GEOmetadb.sqlite]
  gpl    title                                                       
  <chr>  <chr>                                                       
1 GPL570 [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array


gse %>% 
  left_join(gse_gpl, copy = TRUE) %>% 
  filter(gpl == "GPL570") %>% 
  select(title,gse) %>%
  collect() %>% ## this was needed here to be able to filter a second time
  filter(str_detect(title, search_string))
 
 
  
  gse %>% 
    left_join(gse_gpl, copy = TRUE) %>% 
    filter(gpl == "GPL570" & str_detect(title, search_string)) %>% 
    select(title,gse) %>%
    collect() %>% ## this was needed here to actually be able to filter a second time
    filter(str_detect(title, search_string))
    
gse %>% 
    left_join(gse_gpl, copy = TRUE) %>% 
  select(title,gse) %>%
  filter(gpl == "GPL570" & str_detect(title, search_string)) %>% 
  
  collect() %>% ## this was needed here to actually be able to filter a second time
  filter(str_detect(title, search_string))
  
  gse %>% 
    left_join(gse_gpl, copy = TRUE) %>% 
    filter(str_detect(title, search_string))
    filter(gpl == "GPL570") %>% 
    select(title,gse) %>%
    collect() %>% ## this was needed here to be able to filter a second time
    
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

### Example


library(shiny)
library(DBI)
library(pool)

## set these up outside of the server function
## this could be where the database is first filtered based on GPL
## both can be stored (the other data really isn't necessary) so we
## limit calls to the server and use a local version of the table as
## a dataframe for each session

pool <- dbPool(
  drv = RMySQL::MySQL(),
  dbname = "shinydemo",
  host = "shiny-demo.csa7qlmguqrf.us-east-1.rds.amazonaws.com",
  username = "guest",
  password = "guest"
)

ui <- fluidPage(
  textInput("ID", "Enter your ID:", "5"),
  tableOutput("tbl"),
  numericInput("nrows", "How many cities to show?", 10),
  plotOutput("popPlot")
)

server <- function(input, output, session) {
  output$tbl <- renderTable({
    pool %>% tbl("City") %>%
      filter(ID == input$ID)
  })
  output$popPlot <- renderPlot({
    df <- pool %>% tbl("City") %>%
      head(as.integer(input$nrows)[1]) %>% collect()
    pop <- df$Population
    names(pop) <- df$Name
    barplot(pop)
  })
}

shinyApp(ui, server)

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#