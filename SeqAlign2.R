library(shiny)
library(Biostrings)
library(msa)
library(dplyr)
#######################
#Analyzing if several sequences contain the target sequence
#This one can upload files, analyze it and download the result

ui <- fluidPage(
  fluidPage(
    sidebarPanel(
      helpText("Input target sequence below"),
      textInput("text","DNA Sequences"),
      checkboxInput("rev","Reverse complement", value = FALSE),
      fileInput("file1",
                "Choose seq files from directory",
                multiple = TRUE,
                accept='.seq'),
      downloadButton('downloadData', 'DownloadCSV')),
    mainPanel(
      h3('Contain the target sequence or not?'),
      h4('Your result:'),
      tableOutput("table")#,
      #textOutput("filepath")
    )
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
          a <- substring(readChar(.$datapath,nchar=400),51)
          newtext <- if(input$rev) {
            reverseComplement(DNAString(input$text))
          }
          else {
            newtext <- input$text 
          }
          b <- ifelse(grepl(newtext, a), TRUE, FALSE)
          a <- data.frame(.$name,nchar(a),b,a,stringsAsFactors =F)
          colnames(a)=c("name","nchar","contains target string","DNA strings")
          a
        })
    }
    
  })
  output$table <- renderTable({
    datasetInput()
  })
  output$filepath <- renderText({tempdir()})
  print(tempdir())
  output$downloadData <- downloadHandler(
    filename = function() {paste0('test', '.csv')},
    content = function(file) {
      write.csv(datasetInput(), file)
    })
}

shinyApp(ui = ui, server = server)