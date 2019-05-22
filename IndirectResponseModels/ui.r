library(shiny)
library(markdown)
library(shinydashboard)

dashboardPage(skin = "black",
  dashboardHeader(title="Delayed Effects"),
  dashboardSidebar(
    radioButtons("Dataset", "Exercise", choices = c("Step 1 - find best model with one dose", "Step 2 - find best model with three doses","Step 3 - multiple dosing and efficacy"), selected = "Step 1 - find best model with one dose"),
    radioButtons("effectModel", "Which effect model to use?", choices = c("Distributional delay", "Turnover model 2 - (Inhibition of loss)","Turnover model 3 - (Stimulation of build-up)"), selected = "Distributional delay"),
    sliderInput("dose", 
                "Dose in mg/kg",min = 0, 
                max = 50, 
                value = 1, step= 0.1),
    sliderInput("n_doses", 
                "how many doses",min = 1, 
                max = 10, 
                value = 1, step= 1),
    sliderInput("ii", 
                "dosing frequency (h)",min = 1, 
                max = 120, 
                value = 24, step= 1),
    numericInput("untilTimeToPlot","Simulate until (hours)",min = 1,max = 1500,value = 700, step= 50),
    radioButtons("plasma_log", "Plot on a linear or log scale?", choices = c("Log","Linear"), selected = "Linear"),
    bookmarkButton(),
    downloadButton('downloadData', 'Download Simulation data')
  ),
  dashboardBody(
          fluidRow(
            box(width = 4,title="PK", plotOutput("pk_plot_plasma")),
            box(width = 8,title="Effect",plotOutput("effect_plot"))
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
            box(title="Effect Parameters for all effect models",width = 4,
                sliderInput("EC50", 
                            "EC50 (ug/ml)",min = 0.01, 
                            max = 0.5, 
                            value =0.1, step= 0.01)
                
            ),
            box(width = 4,
                title = "Parameters for distributional model",
                sliderInput("kcp", 
                            "compartment exchange rate (1/h)",min = 0, 
                            max = 0.1, 
                            value =0.05, step= 0.01),
                sliderInput("Emax", 
                            "Emax",min = 20, 
                            max = 100, 
                            value =30, step= 1),
                sliderInput("E0", 
                            "E0",min = 0, 
                            max = 20, 
                            value =0, step= 1)
              
            ),
            box(width = 4,
                title = "Parameters for turnover model",
                sliderInput("Rss",
                           " Baseline effect",min = 1,
                           max = 100,
                           value =15, step= 1),
                # sliderInput("kin",
                #             " kin",min = 0,
                #             max = 10,
                #             value =0.5, step= 0.01),
                sliderInput("kout",
                            "kout (1/h)",min = 0.01,
                            max = 1,
                            value =0.8, step= 0.01),
                sliderInput("Emax_turnover", 
                            "Emax",min = 0.01, 
                            max = 1, 
                            value =0.5, step= 0.01)
                
            )
                               
          )
         

              )
     
    
)
