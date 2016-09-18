library(shiny)
library(Biostrings)
library(msa)
library(dplyr)
#######################
# Analyzing if several sequences contain the target sequence.
# This one can upload files, analyze it and download the result.

ui <- fluidPage(
    sidebarPanel(
      helpText("Input target sequence below."),
      textInput("text","DNA Sequences"),
      # Chick if you want reverse complement for the target sequence.
      checkboxInput("rev","Reverse complement", value = FALSE),
      # Input file, call it file1. Use "multiple = TRUE" to allow taking multiple files. Use "accept='.seq'" means we only take .seq file, which is the standard Sanger sequencing file containing only DNA strings.
      fileInput("file1",
                "Choose seq files from directory",
                multiple = TRUE,
                accept='.seq'),
      # We want to remove some low quality sequences.
      # helpText("Since the beginning and end of Sanger sequencing is messy, we want to chop these nucleotides out."),
      helpText("Chop first 49th nt and after 450nt to remove low quality sequences. Usually your target sequences should be located within this region"),
      numericInput("start","Choose where to start",value=50,step=1),
      numericInput("stop","Choose where to start",value=450,step=1),
      # Give the file a name.
      textInput("seqName","File name",value = "SeqInfo"),
      
      # Download CSV files.
      downloadButton('downloadData', 'DownloadCSV'),
      # Download the pdf with multiple sequence alignment.
      downloadButton('downloadPDF',"DownloadPDF")),
    
    mainPanel(
      h3('Contain the target sequence or not?'),
      h4('Your result:'),
      tableOutput("table")
    )
  )


server <- function(input, output) {
  
  datasetInput <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) {
      return(NULL)
    } else {
      inFile %>%
        rowwise() %>%
        do({
          DNASeq <- substring(readChar(.$datapath,nchar=input$stop),input$start)
          newtext <- if(input$rev) {
            reverseComplement(DNAString(input$text))
          }
          else {
            newtext <- input$text 
          }
          ContainT <- ifelse(grepl(newtext, DNASeq), TRUE, FALSE)
          DNASeqDF <- data.frame(.$name,nchar(DNASeq),ContainT,DNASeq,stringsAsFactors =F)
          colnames(DNASeqDF)=c("Name","nchar","Contains target string","DNA strings")
          DNASeqDF
        })
    }
  })
  
  output$table <- renderTable({
    datasetInput()
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {paste0(input$seqName,'.csv')},
    content = function(file) {
      write.csv(datasetInput(), file)
    })

  # Here we build the pdf file for multiple sequence alignmnet.
  output$downloadPDF = downloadHandler(
    # downloadHandler needs two arguments: filename and content.
     filename = paste0(input$seqName,'.pdf'),
     content = function(file) {
       newtext <- if(input$rev) {
         reverseComplement(DNAString(input$text))
       }
       else {
         newtext <- input$text 
       }
       
       df<-as.data.frame(datasetInput())
       DNA <- c(df[,5],newtext)
       names(DNA) <- c(df[,1],"Target Sequence")
       DNA <- unlist(DNA)
       
       msaPrettyPrint(
         msa(DNAStringSet(DNA),order="input")
         , file = 'report.pdf'
         , output="pdf"
         , showNames="left"
         , showLogo="top"
         , consensusColor="BlueRed"
         , logoColors="accessible area"
         , askForOverwrite=FALSE)
       file.rename("report.pdf",file)
     },
     contentType = 'application/pdf'
   )
  
  #you can look up the temdir with the following code
  #print(tempdir())
}

shinyApp(ui = ui, server = server)

#For the sample DNA sequence use:
# CGCGTCTTGTCGAACGAAGCCT