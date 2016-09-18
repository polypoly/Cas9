library(shiny)
library(Biostrings)
library(msa)
library(dplyr)
#######################
# Analyzing if several sequences contain the target sequence.
# This one can upload files, analyze it and download the result.

# Starting the ui function
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
      helpText("You can chop the nucleotides at two ends to remove low quality reads."),
      numericInput("start","Choose where to start",value=50,step=10),
      numericInput("stop","Choose where to stop",value=450,step=10),
      # Give the file a name.
      textInput("seqName","File name",value = "SeqInfo"),
      #
      selectInput("order", "Order of alignment in PDF",
                  choices = c("aligned", "input")),
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

# Starting the server function
server <- function(input, output) {
  
  # Read in the target sequence as targetInput().
  targetInput <- reactive({
    # First check if we need to reverse complement of the target sequence.
    newtext <- if(input$rev) {
      reverseComplement(DNAString(toupper(gsub("[[:space:]]", "", input$text))))
    }
    else {
      newtext <- toupper(gsub("[[:space:]]", "", input$text))
    }
  })
  
  # Read the data files (many input DNA sequences) in and make a data frame as datasetInput().
  datasetInput <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) {
      return(NULL)
    } else {
      inFile %>%
        rowwise() %>%
        do({
          DNASeq <- readChar(.$datapath,nchar=800)
          DNASeq <- gsub("[[:space:]]", "", DNASeq)
          DNASeq <- substring(DNASeq,first=input$start,last=input$stop)
          ContainT <- ifelse(grepl(targetInput(), DNASeq), TRUE, FALSE)
          DNASeqDF <- data.frame(.$name,nchar(DNASeq),ContainT,DNASeq,stringsAsFactors =F)
          colnames(DNASeqDF)=c("Name","Length","Contains target string","DNA strings")
          DNASeqDF
        })
    }
  })
  
  # Make the table that you see on the output page.
  output$table <- renderTable({
    datasetInput()
  })
  
  # This allows us to download csv files. 
  output$downloadData <- downloadHandler(
    filename = function() {paste0(input$seqName,'.csv')},
    content = function(file) {
      write.csv(datasetInput(), file)
    })

  # Here we build the pdf file for multiple sequence alignmnet.
  output$downloadPDF = downloadHandler(
    # DownloadHandler needs two arguments: filename and content.
     filename = function() {paste0(input$seqName,'.pdf')},
     content = function(file) {
       # Make the datasetInput a dataframe and add the target sequence to make a DNAStringSet. 
       df<-as.data.frame(datasetInput())
       DNA <- c(df[,4],targetInput())
       names(DNA) <- c(df[,1],"Target Sequence")
       
       # Generate multiple sequence alignment.
       msaPrettyPrint(
         msa(DNAStringSet(DNA),order=input$order)
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
  # You can look up the temdir with the following code
  #print(tempdir())
}

shinyApp(ui = ui, server = server)

# For the example DNA sequences use:
# GTTATCGCGTCTTGTCGAAC