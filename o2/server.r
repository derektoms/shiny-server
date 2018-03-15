require(shiny)
require(ggplot2)
require(ReacTran)
require(scales)
require(RColorBrewer)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

dataPlots<-reactive({
    
    ## All units are SI:
    # length -- m
    # time -- s
    # concentration -- mol/L
    # pressure -- Pa
    
    ## Constants
    k.H <- 107.100525e6 # L•Pa/mol, Henry's Constant
    R <- 8.314472e3 # Pa•m^3/K•mol, Universal Gas Constant
    D.l <- 2.861e-9 # m^2/s, Diffusion coefficient of oxygen in water at 37 °C \cite{Randers-Eichhorn1996}
    #D.l <- 2.18e-9 # m^2/s, Diffusion coefficient of oxygen in blood at 37 °C \cite{Goldstick1976}
    k.Q <- 2.64e-7 # mol/L, Respiratory constant \cite{FROESE1962}
    
    ## Variables
    Temp <- as.numeric(input$temp) + 273.15 # K (input from °C), Temperature
    a <- (as.numeric(input$volume)/as.numeric(input$surfA))/100 # m (from vol/SA), Depth of fluid
    #     a <- (input$depth)/100 # m, Depth of fluid (3 mm)

    p.O2.inc <- (input$O2_pressure)/100 # Partial pressure of O2 supplied to the incubator
    p.CO2.inc <- 0.05 # Partial pressure of CO2 supplied to the incubator
    p.H2O.inc <- 0.06 # Partial pressure of water vapour in the incubator (at 37 °C)

    p.atm <- as.numeric(input$elev)*101.325e3 # Pa, Atmospheric pressure; In Calgary: 88 kPa http://climate.weather.gc.ca

    OCR <- input$Qmet*1e-18 # mol O2/cell•s, Oxygen consumption for CHO \cite{Wagner2011,Jorjani1999}
    d.cell <- 2e9 # cells/m^2, Cell density \cite{Camire2007}
    
    ## Calculations
    # Oxygen partial pressure in the incubator
    p.O2 <- p.O2.inc*(1-p.CO2.inc-p.H2O.inc)*p.atm

    # Dissolved O2 at surface, using Henry's Law
    b <- 1000*p.O2/k.H # mol/m^3, conversion factor 1 m^3 = 1000 L

    # Cell monolayer oxygen consumption
    Q.cell <- OCR*d.cell # mol O2/s•m^2, Oxygen consumption

    ## Model space 
    # set up a grid with N equally sized boxes
    N <- 100
    zgrid <- setup.grid.1D(x.up=0,x.down=a,N=N) # model space is the media 'a' m deep, upstream is the cell layer at x=0

    z <- zgrid$x.mid # used for plotting concentrations (in the middle of each grid cell)

    ## Differential equation using the 1D transport function from the ReacTran package, where C is the concentration, C.down is the downstream boundary concentration (surface of the culture media), flux.up is the upstream flux (cellular O2 consumption), D is the diffusion coefficient of oxygen, and dx is the model space defined above

    diff <- function(t,Y,parms){
        tran = tran.1D(C=Y, C.down=b,flux.up=-Q.cell,D=D.l,dx=zgrid)
        list(dY = tran$dC,flux.up=tran$flux.up,flux.down=tran$flux.down)
    }

    ## Steady state solution
    std.out0 <- steady.1D(y=runif(zgrid$N),parms=NULL,func=diff,dimens=N)
    oxy <- data.frame(conc = as.numeric(std.out0$y),depth=as.numeric(z),area=as.numeric(1))

    ## Oxygen restriction: if the concentration of oxygen at the bottom of the vessel is below the respiratory constant, k.Q, the minimum at which cells can respire, a new steady state is obtained with C.up representing this minimum concentration. As a result, the cells experience an oxygen deficit, O2.def

    if (oxy[1,1]<k.Q){
        diff <- function(t,Y,parms){
            tran = tran.1D(C=Y, C.down=b,C.up=k.Q,D=D.l,dx=zgrid)
            list(dY = tran$dC,flux.up=tran$flux.up,flux.down=tran$flux.down)
        }
        std.out0 <- steady.1D(y=runif(zgrid$N),parms=NULL,func=diff,dimens=N)
        oxy <- data.frame(conc = as.numeric(std.out0$y),depth=as.numeric(z),area=as.numeric(1))
        d.C.max <- (D.l*(b-k.Q))/a # maximum oxygen concentration gradient
        O2.def <- (Q.cell-d.C.max)/Q.cell # oxygen deficit (%)
    } else {
        O2.def <- 0 # no oxygen deficit
    }

    O2.df <- data.frame(group=c("Available","Restricted"),value=c(100*(1-O2.def),100*O2.def)) # collect oxygen restriction data into a dataframe
    ## Pass calculated values out
    media.depth <<- a*1000
    cell.oxygen <<- oxy[1,1]
    
    ## Culture well plot
    well <<- ggplot(oxy,aes(x=area,y=depth,fill=conc)) + geom_tile() + scale_y_continuous(labels=function(x)x*100,limits=c(0,0.01)) + labs(fill=expression(O[2]~(mM)),x=NULL,y="depth (cm)",title="Culture Media Oxygen Profile") + scale_fill_gradientn(colours=brewer.pal(11,"RdYlGn"),limits=c(0,0.5),oob=squish, values=c(0,.25,1)) + scale_x_continuous(breaks=NULL,labels=NULL) + theme(legend.key.height=unit(3,'line'),legend.key.width=unit(2,'line'))
    
    # Oxygen restriction
    O2 <<- ggplot(O2.df,aes(x="",y=value,fill=group)) + geom_bar(width=1,stat="identity") + labs(x=NULL,y=NULL) + coord_polar("y",start=0) + scale_fill_manual(name=NULL,labels=c(expression(Normal~ O[2]~consumption),expression(O[2]~Restriction)),values=c("green","black")) + 
        ggtitle("Oxygen Restriction") +
    theme(legend.key.size=unit(2,'lines'),legend.text.align=0)
    
})

output$advanced <- renderUI({
    tagList(
    textInput("temp",
              "Incubation temperature °C):",
              value=37),
    sliderInput("surfA",
              label=div(HTML("Surface area (cm<sup>2</sup>):")),
              min = 0.1,
              max = 100,
              value = input$cultureware),
    sliderInput("Qmet",
              label=div(HTML("Cell metabolism (amol O<sub>2</sub>/s•cell):")),
              min = 0,
              max = 1000,
              value = input$celltype,
              step=4)
    )
})
            
  output$wellPlot <- renderPlot({
        dataPlots()
        well
  })

  output$piePlot <- renderPlot({
        dataPlots()
        O2
    })
    
  output$Qcell <- renderText({
      dataPlots()
      paste(input$Qmet/1e7)
  })  
  output$atmo <- renderText({
      dataPlots()
      paste("Atmosphere (atm) = ",input$elev)
  })
  output$depth <- renderText({
      dataPlots()
      paste("Depth (mm) = ",media.depth)
  })
  output$oxy <- renderText({
      dataPlots()
      paste("Oxygen at cell surface (mM) = ",format(cell.oxygen,digits=6))
  })
    
  ## Kill shinyApp when session closes
  session$onSessionEnded(stopApp)
})