## shiny app for user interface to generate CNV plots and browsable IGV files from EPIC idat files ##

library(shiny)
library(gridExtra)
source("./scripts/EPIC_functions v3.R")

inputChoicesEPIC <- "First select input directory"
refData <- NULL
levelN <- c(0,1,2,4,6,8,12)
x <- NULL

## user interface ##
ui <- navbarPage("UMCU Array R tools", 
                 tabPanel(tags$b("EPIC array"),              
                          # row layout with input and output definitions ----
                          fluidRow(
                            column(2,titlePanel(h4("Select folder")),actionButton("dir",tags$b("Input directory"),icon("folder-open"),style="color: #fff; background-color: #0088c7; border-color: #ffffff")),
                            column(4,titlePanel(h4("Select sample")),selectInput("samplesEPIC", NULL,choices = inputChoicesEPIC)),
                            column(2,titlePanel(h4("Plot CNVs")),actionButton("CNVplot", tags$b("Plot CNVs"), icon("paper-plane"), style="color: #fff; background-color: #0088c7; border-color: #ffffff")),
                            column(4,titlePanel(h4("Tumor cell percentage")),textOutput("tcp"))
                          ),
                          fluidRow(
                            column(5,sliderInput("level2N",label = "Set 2N level",min = -1.2,max = 1.2,value = 0,step = 0.01)),
                            column(5,sliderInput("level1N",label = "Set 1N level",min = -1.2,max = 1.2,value = -0.4,step = 0.01)),
                            column(2,checkboxInput("showCNlines",tags$b("Show CN lines"),value = TRUE)),
                            column(10, titlePanel(h4("predicted CDKN2A status")), textOutput("CDKN2A")),
                            column(2, titlePanel(h4("Export PDF")), downloadButton("createPDF", tags$b("generate PDF") , icon("paper-plane"), style="color: #fff; background-color: #0088c7; border-color: #ffffff"))
                          ),
                          fluidRow(
                            br(),
                            br(),
                            mainPanel(width = 12,plotOutput("plotEPIC")),
                            mainPanel(width = 12,plotOutput("plotPDF"))
                          )
                 ),
                 tags$style(HTML(".navbar-default .navbar-brand {color: #ffffff;}
                                 .navbar { background-color: #0088c7;}
                                 .navbar-default .navbar-nav > li > a {color:#ffffff;}
                                 "))
)
## server function ##
server <- function(input, output, session) {
  
  ## list the available samples in a chosen directory ##
  observeEvent(input$dir,ignoreInit = T,{
    chosenDir <<- choose.dir(parameters$dataFolder)
    if(!is.na(chosenDir)){
      updateSelectInput(session,"samplesEPIC",choices=c("All",loadEPICsamples(chosenDir)) )
    }
    
  })
  
  cnLevels <- reactiveValues()
  observe({
    #cnLevels$levelN <- c(0,1,2,4,6,8,12)
    cnLevels$tumorperc <- 2*(2^(input$level2N) -2^(input$level1N))
    cnLevels$cnLevels <- log2(2^(input$level2N) + ((levelN-2)*(2^(input$level2N) -2^(input$level1N))))
    if(is.null(x)){
      cnLevels$CDKN2A_status <- "First load data"
    }else if((which.min(abs(x@detail$ratio["CDKN2A/B"]- (cnLevels$cnLevels[1:3])))) == 1){
      cnLevels$CDKN2A_status <- "HOM DEL"
    } else if((which.min(abs(x@detail$ratio["CDKN2A/B"]- (cnLevels$cnLevels[1:3])))) == 2){
      cnLevels$CDKN2A_status <- "HET DEL"
    }else{
      cnLevels$CDKN2A_status <- which.min(abs(x@detail$ratio["CDKN2A/B"]- (cnLevels$cnLevels[1:3])))
    }
    
  })
  
  output$tcp <- renderText({
    #paste(input$level1N)
    paste(round(cnLevels$tumorperc*100,0),"%")
  })
  
  output$CDKN2A <- renderText({
    paste(cnLevels$CDKN2A_status)
  })
  
  ## If plot CNVs is clicked controls are loaded and CNV plots and IGV files are generated ##
  observeEvent(input$CNVplot,ignoreInit = T,{
    if(is.null(refData)){
      showModal(modalDialog("Loading Controls", footer=NULL))
      #refData <<- loadAndProcessControls(addCustomDetails = T)
      refData <<- readRDS(file = "T:/pathologie/PRL/Groep-Brosens/2. Anna Vera/15. ovaryNET/Ovary/data/ControlMetSet rds/controlMetSet.rds" )
      removeModal()
    }
    ## samples can either be processed 1 by 1 or all at the same time, depending on what is chosen in the app ##
    if(input$samplesEPIC == "All"){
      samples <- loadEPICsamples(chosenDir)
      for ( i in 1:length(samples)){ 
        showModal(modalDialog(paste0("Loading data and making plot for ",samples[i]), footer=NULL))
        patient <- loadAndProcessSamples(dataFolder = chosenDir,sampleName = samples[i])
        x <- segmentData(patient = patient,refData = refData,dataFolder = chosenDir)
        pdf(paste0(chosenDir,"/",sampleNames(patient),".pdf"),width=18,height = 7)
        CNV.genomeplot(x)
        dev.off()
        
        removeModal()
      }
    }
    
    if(input$samplesEPIC != "All"){
      showModal(modalDialog(paste0("Loading data and making plot for ",input$samplesEPIC), footer=NULL))
      patient <- loadAndProcessSamples(dataFolder = chosenDir,sampleName = input$samplesEPIC)
      x <<- segmentData(patient = patient,refData = refData,dataFolder = chosenDir)
      pdf(paste0(chosenDir,"/",sampleNames(patient),".pdf"),width=12,height = 7)
      CNV.genomeplot(x, main= c(sampleNames(patient), "\n", "This CNV tool has not been validated for diagnostic purposes"))
      dev.off()
      
      removeModal()
      output$plotEPIC <- renderPlot({
        CNV.genomeplot(x)
        min_x <- par("usr")[1]
        max_x <- par("usr")[2]
        xlabelPos <- 0.985*(max_x-min_x)
        
        if(input$showCNlines){
          for (l in 1:length(cnLevels$cnLevels)){
            lines(x=c(min_x,(0.97*(max_x-min_x))),y=rep(cnLevels$cnLevels[l],2),lwd=2,col="red",lty=2)
          }
          #abline(h=cnLevels$cnLevels,lwd=2,col="red",lty=2)
          text(rep(xlabelPos,length(cnLevels$cnLevels)),c(cnLevels$cnLevels),labels = c(paste0(levelN,"N")))
        }
        
      })
      
      output$createPDF <- downloadHandler(
        
        filename = function() {
          "plots.pdf"
        },
        
        content = function(file) {
          
          pdf(file)
          CNV.genomeplot(x)
          min_x <- par("usr")[1]
          max_x <- par("usr")[2]
          xlabelPos <- 0.955*(max_x-min_x)
          
          if(input$showCNlines){
            for (l in 1:length(cnLevels$cnLevels)){
              lines(x=c(min_x,(0.97*(max_x-min_x))),y=rep(cnLevels$cnLevels[l],2),lwd=2,col="red",lty=2)
            }
            text(rep(xlabelPos,length(cnLevels$cnLevels)),c(cnLevels$cnLevels),labels = c(paste0(levelN,"N")))
            text(round(cnLevels$tumorperc*100,0),"%")
            dev.off()
            
            # pdf(paste0(chosenDir,"/",sampleNames(patient),".pdf"),width=12,height = 7)
            #CNV.genomeplot(x, chr = "all", main= c(sampleNames(patient), "\n", "This CNV tool has not been validated for diagnostic purposes"))
          }
        })
      
    }
    
    
    
    
  })
  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)



