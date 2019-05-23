library(shiny)
#library(shinythemes)

shinyUI(fluidPage(
  #theme = "bootstrap_readable.css",
  br(),
  fluidRow(
    HTML("<h2>Target mediated drug disposition model </h2>"),
    #HTML("Text for introduction"),
    br(),
    #img(src='MiceAloneforShiny.png'),
    br(),
    br()
    ),
  fluidRow(
    sidebarLayout(  
      sidebarPanel(h4("Dose in mg:"),
                   sliderInput("dose","", 
                               min = 0, 
                               max = 600, 
                               value = 150),
                   sliderInput("interval", 
                               "Dosing interval in [days]:", 
                               min = 1, 
                               max = 100, 
                               value = 28),
                   sliderInput("n_doses", 
                               "Number of doses:", 
                               min = 1, 
                               max = 10, 
                               value = 10),
                   numericInput("ligandss", "Ligand concentration (fM)", 11, min = NA, max = NA, step = NA,
                                width = NULL),
                   sliderInput("kd", 
                               "Kd in [nM]:", 
                               min = 0.001, 
                               max = 1, 
                               value = 0.78, step= 0.001),
                   sliderInput("cl", 
                               "Clearance in [ml/day]:", #0.0068 L/hour
                               min = 1000*0.0022*24, 
                               max = 1000*0.0068*24, 
                               value = 1000*0.0068*24, step= 1000*0.0001*24),
                   sliderInput("targetTt12", 
                               "target turnover t1/2 in hours:", 
                               min = 0.1, 
                               max = 10, 
                               value = 0.5),
                   numericInput("Comp_factor", "complex elimination compared to mAb [factor]", 1, min = NA, max = NA, step = NA,
                                width = NULL),
                   actionButton("startSim", "Start Simulation", icon("paper-plane"), 
                                style="color: #fff; background-color: #FF0000; border-color: #2e6da4")
                   
      ),  
      mainPanel(  
        plotOutput("pk_plot"),
        plotOutput("pd_plot")#,height = "750px"
      )
    )
  )
))

