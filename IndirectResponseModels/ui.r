library(shiny)
library(markdown)
library(shinydashboard)

dashboardPage(skin = "black",
  dashboardHeader(title="Delayed Effects"),
  dashboardSidebar(
    radioButtons("effectModel", "Which effect model to use?", choices = c("Turnover model 3 - (Stimulation of build-up)","Turnover model 2 - (Inhibition of loss)"), selected = "Turnover model 3 - (Stimulation of build-up)"),
    sliderInput("dose", 
                "Dose in mg/kg",min = 0, 
                max = 50, 
                value = 1, step= 1),
    sliderInput("n_doses", 
                "how many doses",min = 1, 
                max = 10, 
                value = 1, step= 1),
    sliderInput("ii", 
                "dosing frequency (h)",min = 1, 
                max = 120, 
                value = 24, step= 1),
    actionButton("startSim", "Start Simulation", icon("paper-plane"), 
                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
    numericInput("untilTimeToPlot","Simulate until (hours)",min = 1,max = 1500,value = 700, step= 50),
    radioButtons("plasma_log", "Plot on a linear or log scale?", choices = c("Log","Linear"), selected = "Linear"),
    bookmarkButton(),
    downloadButton('downloadData', 'Download Simulation data'),
    radioButtons("Dataset", "Exercise", choices = c("Step 1 - find best model with one dose", "Step 2 - find best model with three doses","Step 3 - multiple dosing and efficacy"), selected = "Step 1 - find best model with one dose")
    
  ),
  dashboardBody(
          fluidRow(
            box(width = 8,title="PK", plotOutput("pk_plot_plasma",height = "200px"))
            
          ),
          fluidRow(
          box(width = 8,title="PD",plotOutput("effect_plot"))
          ),
          fluidRow(
            # box(title="PK Parameters important for all models",width = 3,
            #   sliderInput("Central", 
            #               "Central volume (ml/kg)",min = 50, 
            #               max = 500, 
            #               value =100, step= 1)
              # ,
              # sliderInput("Peripheral", 
              #             "Peripheral volume (ml/kg)",min = 1, 
              #             max = 100, 
              #             value =2, step= 1),
              # sliderInput("keC", 
              #             "elimination rate (1/h)",min = 0, 
              #             max = 0.5, 
              #             value =0.12, step= 0.01)
            #),
            box(width = 8,
                title = "Parameters for turnover models",
                sliderInput("Rss",
                            " Baseline effect",min = 1,
                            max = 100,
                            value =15, step= 1),
                sliderInput("kout",
                            "kout (1/h)",min = 0.01,
                            max = 1,
                            value =0.8, step= 0.01),
                sliderInput("Emax_turnover", 
                            "Emax",min = 1, 
                            max = 100, 
                            value =2, step= 0.1),
                sliderInput("EC50", 
                            "EC50 (ug/ml)",min = 0.01, 
                            max = 10, 
                            value =0.1, step= 0.1)
                
                # sliderInput("kin",
                #             " kin",min = 0,
                #             max = 10,
                #             value =0.5, step= 0.01),


                
            )
                               
          )
         

              )
     
    
)
