
source("theheartoftheall.R")

#Let the fun begin

shinyServer(function(input, output, session) {
  

  getAllSimulationResults <- eventReactive({input$startSim 
    session$clientData},{
    
    ParametersAndDosing     <- getAllParametersAndDosing(input)
    SimulationResults       <- getSimulationResults(ParametersAndDosing)
   # SimulationResults       <- convertUnits(SimulationResults,ParametersAndDosing$p)
  
    return(SimulationResults)
  })
  
  
  getAllPlotParameters <- reactive({
    plot_p                  <- getPlotParameters(input)
    
    return(plot_p)
  })
  
  getAllDatasets <- reactive({
    Datasets                <- loadData()
    return(Datasets)
  })

  output$pk_plot_plasma <- renderPlot({ 
   
    SimulationResults <- getAllSimulationResults()
    
    plot_p            <- getAllPlotParameters()
    pk_plot           <- makePKplot_plasma(SimulationResults,plot_p)
    Datasets          <- getAllDatasets()
    pk_plot_withdata  <- addPKdataToplot_plasma(pk_plot,Datasets,input)
    
    print(pk_plot) 
    
    })
  
  output$effect_plot <- renderPlot({ 
    
    SimulationResults <- getAllSimulationResults()
    
    plot_p            <- getAllPlotParameters()
    effect_plot       <- makeEffektplot(SimulationResults,plot_p)
    Datasets          <- getAllDatasets()
    effect_plot_withdata  <- addEffectdataToplot_effect(effect_plot,Datasets,input)
    
    print(effect_plot_withdata) 
  
  })
  
  # -----------DATASET OUTPUT-------------      
  
  observeEvent (input$show_rawdata, {
    
    SimulationResults <- getAllSimulationResults()
    plot_p            <- getAllPlotParameters()
    rawData           <- returnRawData(SimulationResults,plot_p)
    
    output$dataset <- renderTable(rawData)
  })
  
  output$downloadData <- downloadHandler(
    
    filename = function() { paste('rawData', '.csv', sep='') },
    content = function(file) {
      SimulationResults <- getAllSimulationResults()
      plot_p            <- getAllPlotParameters()
      write.csv(returnRawData(SimulationResults,plot_p), file)
    }
  )
  
})

