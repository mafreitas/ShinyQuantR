#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

source("https://bioconductor.org/biocLite.R")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(shiny, DT, limma, 
               psych, ggplot2,
               genefilter,reshape2,
               stringr,missForest,shinyjs)

options(shiny.maxRequestSize=1000*1024^2) 

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  observeEvent(input$close, {
    js$closeWindow()
    stopApp()
  })
  
  
  datafile <- reactive({
    req(input$file1)
    
    data.return <- read.csv(
      input$file1$datapath,
      header = input$header,
      sep = input$sep,
      quote = input$quote
    )
    colnames(data.return) = make.names(colnames(data.return))
    print(colnames)
    if (!('Protein.IDs' %in% colnames(data.return)) ){
      colnames(data.return)[input$id_col] <- 'Protein.IDs' 
    }
        
    data.return
  })
  
  data <- reactive({
    req(datafile())
    
    data.return <- datafile()
    
    if (input$mq == T){
      data.return <- cbind(
        Protein.name = data.return$Protein.IDs,
        Protein.IDs = data.return$Protein.IDs,
        data.return[, grepl("LFQ" , names(data.return))]
      )
    }
    else {
      data.return = cbind(
        Protein.name = data.return$Protein.IDs,
        Protein.IDs = data.return$Protein.IDs,
        data.return
      )
    }
    
    newcolnames <- gsub("LFQ.intensity.", "", colnames(data.return))
    colnames(data.return) <- newcolnames
    data.return$Protein.name <-
      sub(".*\\|.*\\|", "", data.return$Protein.name)
    data.return$Protein.name <-
      sub("_.*", "", data.return$Protein.name)
    
    data.return$Protein.IDs <-
      str_match(data.return$Protein.IDs, "\\|([A-Z_0-9]*)\\|")[, 2]
    data.return$Protein.IDs <-
      paste0("<a href='http://www.uniprot.org/uniprot/",
        data.return$Protein.IDs,"' target='_blank'>",
        data.return$Protein.IDs,"</a>")
    data.return[data.return == 0] <- NA
    
    return(data.return)
  })
  
  data.selected <- reactive({
    data()[, input$include]
  })
  
  data.filtered <- reactive({
    data.selected()[rowSums(!is.na(data.selected())) >= input$remove, ]
  })
  
  protein.filtered <- reactive({
    data()[rowSums(!is.na(data.selected())) >= input$remove,2 ]
  })
  
  data.transformed = reactive({
    data.return <- data.filtered()
    if (input$transform == "None") {
      #data.return <- data.return
    } else if (input$transform == "log2") {
      data.return <- log2(data.return)
    }
    
    data.return
    
  })
  
  data.imputed = reactive({
    data.return <- data.transformed()
    if (input$impute == "None") {
      #data.return <- data.return
    } else if (input$impute == "missforest") {
      showModal(
        modalDialog(title = "Missing Value Imputation",
          "Please wait... This coule take a few minutes!")
      )
      data.return <- missForest(data.return)$ximp
      removeModal()
    }
    data.return
    
  })
  
  data.normalized = reactive({
    data.return <- data.imputed()
    if (input$normalization == "None") {
      #data.return <- data.return
    } else if (input$normalization == "Trimmed Mean") {
      # Calculate the trimmed mean of all expression values
      # A trimmed mean  is somewhere between the mean and median of the set of values.
      # We remove the top and bottom 2% of values (for example), and find the mean of the remaining values,
      trmean = apply(data.return, 2, mean, trim = 0.02)
      mean.of.trmeans = mean(trmean)
      data.return <- data.return / trmean * mean.of.trmeans
    } else if (input$normalization == "Trimmed Median") {
      # Divide all expression values by the mean for that batch, and multiply by a scaling factor,
      # such as the mean of the trimmed means.
      # Calculate the trimmed mean of all expression values
      trmed = apply(data.return, 2, median)
      med.of.trmedians = median(trmed)
      data.return <-  data.return / trmed * med.of.trmedians
    } else if (input$normalization == "Quantile") {
      data.return <- normalizeQuantiles(as.data.frame(data.return))
    }
    
    data.return
    
  })
  
  
  
  data.a <- reactive({
    data.normalized()[, input$groupa]
  })
  
  data.b <- reactive({
    data.normalized()[, input$groupb]
  })
  
  data.a.mean <- reactive({
    data.frame(mean = apply(data.normalized()[, input$groupa], 1, mean))
  })
  
  data.b.mean <- reactive({
    data.frame(mean = apply(data.normalized()[, input$groupb], 1, mean))
    
  })
  
  plotdata <- reactive({
    if (input$sigtest == "t-test") {
      m = as.matrix(cbind(data.a(), data.b()))
      f = factor(c(rep("groupa", ncol(data.a(
      ))), rep("groupb", ncol(data.b(
      )))))
      data.tmp = data.frame(
        Protein.IDs = protein.filtered(),
        logFC = data.b.mean()$mean - data.a.mean()$mean,
        pValue = rowttests(m, f)$p.value
      )
      data.tmp$neglogpValue = -log10(data.tmp$pValue)
      
    }
    else if (input$sigtest == "limmaFit") {
      m = as.matrix(cbind(data.a(), data.b()))
      fac <-
        factor(c(rep("groupA", ncol(data.a(
        ))), rep("groupB", ncol(data.b(
        )))))
      fit <- lmFit(m, design = model.matrix( ~ fac))
      
      fit <- eBayes(fit)
      data.tmp <- data.frame(
        Protein.IDs = protein.filtered(),
        topTable(fit, sort.by = 'none', n = Inf))
      data.tmp$neglogpValue = -log10(data.tmp$P.Value)
    }    
    data.tmp <- cbind(data.tmp,data.a(),data.b())
    return(data.tmp)
    
    
  })
  
  
  output$selectdata <- renderUI ({
    selectInput("include","Include Columns:",
      names(data()[,c(-1,-2)]), multiple = TRUE)
  })
  
  output$rawdata <- DT::renderDataTable({
    req(data)
    DT::datatable(
      format(data(),scientific = 1,digits=4),
      escape = FALSE,
      options = list(
        scrollY = '600px', 
        paging = FALSE,
        scrollX = TRUE
      ),
      rownames= FALSE
    )
  })
  
  output$untransformed <- DT::renderDataTable({
    input$process
    isolate( 
      DT::datatable(
        format(data.frame(Protein.IDs=data()$Protein.IDs,data.selected()), digits = 3, scientific = 5),  
        escape = FALSE,
        options = list(
          scrollY = '300px', 
          paging = FALSE,
          scrollX = TRUE
        ),
        rownames= FALSE
      ))
  })
  
  output$untransformed_plot <- renderPlot({
    input$process
    isolate( 
      boxplot(data.selected())
    )
  })
  
  output$transformed <- DT::renderDataTable({
    input$process
    isolate( 
      DT::datatable(
        format(data.frame(Protein.IDs=protein.filtered(),data.normalized()), digits = 3, scientific = 5),  
        escape = FALSE,
        options = list(
          scrollY = '300px', 
          paging = FALSE,
          scrollX = TRUE
        ),
        rownames= FALSE
      )
    )
  })
  
  output$transformed_plot <- renderPlot({
    input$process
    isolate( 
      boxplot(data.normalized())
    )
  })
  
  output$selecta <- renderUI ({
    
    selectInput("groupa","Group A:",
      names(data.normalized()), multiple = TRUE)
  })
  
  output$selectb <- renderUI ({
    selectInput("groupb","Group B:",
      names(data.normalized()), multiple = TRUE)
  })
  
  output$groupa_data <- DT::renderDataTable({
    
    DT::datatable(
      format(data.frame(Protein.IDs=protein.filtered(),data.a(),data.a.mean()), digits = 3, scientific = 5), 
      escape = FALSE,
      options = list(
        lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
        pageLength = 15,
        scrollY = '300px', 
        paging = FALSE,
        scrollX = TRUE
      ),
      rownames= FALSE
    )
  })
  
  output$groupb_data <- DT::renderDataTable({
    
    DT::datatable(
      format(data.frame(Protein.IDs=protein.filtered(),data.b(),data.b.mean()), digits = 3, scientific = 5), 
      escape = FALSE,
      options = list(
        lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
        pageLength = 15,
        scrollY = '300px', 
        paging = FALSE,
        scrollX = TRUE
      ),
      rownames= FALSE
    )
  })
  
  output$groupa_plot <- renderPlot({
    pairs.panels(data.a()) 
  })
  
  output$groupb_plot <- renderPlot({
    pairs.panels(data.b()) 
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('resultset-', Sys.Date(), '.csv', sep='')
    },
    content = function(file) {
      write.csv(plotdata(), file, row.names = FALSE)
    }
  )
  
  
  
  output$sig_data <- DT::renderDataTable({
    
    DT::datatable(
      format(plotdata(), digits = 3, scientific = 5), 
      escape = FALSE,
      options = list(
        scrollY = '300px', 
        paging = FALSE,
        scrollX = TRUE
      ),
      rownames= FALSE
    )
  })
  
  
  data.sig.selected <- reactive({
    data.return <- brushedPoints(plotdata(), input$plot1_brush)
    return(data.return)
  })
  
  data.sig.melted <- reactive({
    data.return <- data.sig.selected()

    data.return <- melt(data.return,id.vars='Protein.IDs',measure.vars=cbind(input$groupa,input$groupb))
    print (data.return)
    return(data.return)
  })
  
  output$sig_data_selected <- DT::renderDataTable({
    
    DT::datatable(
      format(data.sig.selected(), digits = 3, scientific = 5), 
      escape = FALSE,
      options = list(
        lengthMenu = list(c(5, 10, -1), c('5', '10', 'All')),
        pageLength = 5,
        buttons = list('copy', 'csv', 'excel', 'pdf', 'print'),
        scrollY = '300px', 
        paging = FALSE,
        scrollX = TRUE
      ),
      rownames= FALSE
    )
  })
  
  output$plot1 <- renderPlot({
    g = ggplot(data=plotdata(), aes(x=logFC, y=neglogpValue)) +
      geom_point(alpha=0.4, size=1.75) +
      xlab("log2 fold change") + ylab("-log10 p-value")
    
    g
  })
  
  output$plot2 <- renderPlot({
    g = ggplot(data=data.sig.melted(), aes(x = variable, y = value)) + 
      geom_line(aes(color = Protein.IDs, group = Protein.IDs))
  
    g
  })
  
  
  
})
