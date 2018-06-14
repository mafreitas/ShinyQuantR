
#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#



# Define UI for application that draws a histogram
shinyUI(
  
  # Application title
  navbarPage("DataScienceTemplate",
    tabPanel("Load Data",
      sidebarLayout(
        # Sidebar
        sidebarPanel(
          fileInput('file1', 'Choose CSV File',
            accept=c('text/csv', 
              'text/comma-separated-values,text/plain', 
              '.csv')),
          checkboxInput('header', 'Header', TRUE),
          checkboxInput('mq', 'MaxQuant', FALSE),
          radioButtons('sep', 'Separator',
            c(Comma=',',Semicolon=';',Tab='\t'),
            '\t'),
          radioButtons('quote', 'Quote',
            c(None='','Double Quote'='"','Single Quote'="'"),
            '"')
        ),
        # MainPanel
        mainPanel(
          DT::dataTableOutput('rawdata')
        )
      )
    ),
    
    tabPanel("Transform Data",
      sidebarLayout(
        # Sidebar with a slider input
        sidebarPanel(
          
          uiOutput("selectdata"),
          
          
          textInput("remove","Min # of Values:",
            "3"),
          
          selectInput("transform","Log2 Transform:",
            c("None","log2")),
          
          selectInput("impute","Missing Value Imputation:",
            c("None","missforest")),
          
          selectInput("normalization","Normalization:",
            c("None","Trimmed Mean","Trimmed Median","Quantile"), selected = "None", multiple = FALSE),
          
          actionButton("process", "Perfrom Data Processing")
          
          
          
          
        ),
        # MainPanel
        mainPanel(
          fluidRow(
            column(8,DT::dataTableOutput('untransformed')),
            column(4,plotOutput('untransformed_plot'))
            
          ),
          fluidRow(
            column(8,DT::dataTableOutput('transformed')),
            column(4,plotOutput('transformed_plot'))
          )
          
        )
      )
    ),
    
    tabPanel("Select Data",
      sidebarLayout(
        # Sidebar
        sidebarPanel(
          uiOutput("selecta"),
          uiOutput("selectb")
        ),
        # MainPanel
        mainPanel(
          fluidRow(
            column(8,DT::dataTableOutput('groupa_data')),
            column(4,plotOutput('groupa_plot'))
            
          ),
          fluidRow(
            column(8,DT::dataTableOutput('groupb_data')),
            column(4,plotOutput('groupb_plot'))
          )
        )
      )
    ),
    
    tabPanel("Significance Test",
      sidebarLayout(
        # Sidebar
        sidebarPanel(
          selectInput("sigtest","Significance Test:",
            c("t-test","limmaFit"), multiple = FALSE),
          #numericInput("threshold", "PValue Threshold:", 0.05, min = 0, max = 1,step = 0.01)
          downloadButton('downloadData', label ="Download Result Set")
          
        ),
        # MainPanel
        mainPanel(
          fluidRow(
            DT::dataTableOutput('sig_data')
          ),
          fluidRow(
            plotOutput("plot1", # Equivalent to: click = clickOpts(id = "plot_click")
              brush = brushOpts(id = "plot1_brush")
            )
          ),
          fluidRow(
            DT::dataTableOutput('sig_data_selected')
          ),
          fluidRow(
            plotOutput("plot2")
          )
        )
      )
    )
  )
)





