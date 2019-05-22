#Dose is in mg/kg
#it gets transfered into ug/kg
#Volume is in ml/kg therefore it is ug/ml


library (deSolve)
library (MASS)
library (lattice)
library(ggplot2)
library(scales)

#help functions
units <- function( number,from,to,molar_mass){
  #units converts the units
  
  # 10 3     E  3     kilo     k
  # 10 2     E  2     hecto    h
  # 10-0     E -0     base     b
  # 10-1     E -1     deci     d
  # 10-2     E -2     centi    c
  # 10-3     E -3     milli    m
  # 10-6     E -6     micro    u
  # 10-9     E -9     nano     n
  # 10-12    E-12     pico     p
  # 10-15    E-15     femto    f
  # other possible units
  # dalton
  # mol
  # grams
  
  #examples 
  #1. wieviel ng sind 50 g?
  #units(50,'b','n');  
  #1. wieviel mol sind 1 g stoff mit einem molekulargewicht von 1?
  #units(50,'b','n');  
  
  avo = 6.022*10^23; 
  
  #%%%%%%%%------------->%%%%%%%%%%%%%%%%%%%%%%%
  #convert to standart unit base
  switch (from,   
          k={
            standart = number*10^3;
          },
          h={
            standart = number*10^2;
          },
          b={
            standart = number*10^0;
          },
          d={
            standart = number*10^-1;
          },
          c={
            standart = number*10^-2;
          },
          m={
            standart = number*10^-3;
          },
          u={
            standart = number*10^-6;
          },
          n={
            standart = number*10^-9;
          },
          p={
            standart = number*10^-12;
          },
          f={
            standart = number*10^-15;
          },
          mol={#standart is grams       
            standart = number*molar_mass;
          },
          gram={   
            standart = number;
          },
          dalton={   
            standart = number/avo;   
          }
  )
  
  #%%%%%%%%<-------------%%%%%%%%%%%%%%%%%%%%%%%
  #convert to needed unit 
  
  switch (to,   
          k={
            standart = standart*10^-3;
          },
          h={
            standart = standart*10^-2;
          },
          b={
            standart = standart;
          },
          d={
            standart = standart*10^1;
          },
          c={
            standart = standart*10^2;
          },
          m={
            standart = standart*10^3;
          },
          u={
            standart = standart*10^6;
          },
          n={
            standart = standart*10^9;
          },
          p={
            standart = standart*10^12;
          },
          f={
            standart = standart*10^15;
          },
          mol={#standart is grams       
            standart = standart*(1/molar_mass);
          },
          gram={   
            standart = standart;
          },
          dalton={   
            standart = standart*avo; 
          }
  )
  
  return(standart)
  
}

#Data
DataEffect_1dose<-function(Datasets){
  # Format: time (hours), 0.1 mg/kg	0.3 mg/kg	1 mg/kg	3 mg/kg	5 mg/kg	15 mg/kg	50 mg/kg
  
  effect_1 = c(  #data is in pmol/L
    1.2, 12.15051208,
    6, 21.74665972,
    12, 29.18131362,
    24, 37.03490936,
    36, 40.94087811,
    48, 39.68344127,
    72, 39.47319387,
    96, 32.39906443,
    168, 21.30101371,
    336, 10.4949575
    )
  
  effect_1_Matrix = matrix(effect_1 ,nrow=20,ncol=2,byrow = TRUE)        # fill matrix by rows 
  
  #convert units
  
  
  #make groups
  
  sizes = dim(effect_1_Matrix)
  effect_1_Matrixwithgroups = c()
  for (i in 2:sizes[2]){
    effect_1_Matrixwithgroups <- rbind(effect_1_Matrixwithgroups,cbind(c(effect_1_Matrix[,1]),c(effect_1_Matrix[,i]),c(rep(i-1, length(effect_1_Matrix[,i])))))    
  }
  effect_1_Matrixwithgroups=data.frame(effect_1_Matrixwithgroups)
  
  #Human_Data77108_RO7062931_Plasma_0.1mg <- subset(Human_RO7062931PlasmaMatrixwithgroups,X3==1)
  effect_1_Matrixwithgroups1mg <- effect_1_Matrixwithgroups
  effect_1_Matrixwithgroups1mg$X3 = 1
  Datasets$DataEffect_1dose <- effect_1_Matrixwithgroups1mg
  
  return(Datasets)
}
DataEffect_15dose<-function(Datasets){
  # Format: time (hours), 0.1 mg/kg	0.3 mg/kg	1 mg/kg	3 mg/kg	5 mg/kg	15 mg/kg	50 mg/kg
  
  effect_15 = c(  #data is in pmol/L
    1.2,	11.59308016,
    6,	22.54121306,
    12,	30.58694867,
    24,	43.44980274,
    36,	42.83132646,
    48,	42.7092062,
    72,	44.93360365,
    96,	47.52414454,
    168,	44.31675663,
    336,	25.03251502
    )
  
  effect_15_Matrix = matrix(effect_15 ,nrow=20,ncol=2,byrow = TRUE)        # fill matrix by rows 
  
  #convert units
  
  
  #make groups
  
  sizes = dim(effect_15_Matrix)
  effect_15_Matrixwithgroups = c()
  for (i in 2:sizes[2]){
    effect_15_Matrixwithgroups <- rbind(effect_15_Matrixwithgroups,cbind(c(effect_15_Matrix[,1]),c(effect_15_Matrix[,i]),c(rep(i-1, length(effect_15_Matrix[,i])))))    
  }
  effect_15_Matrixwithgroups=data.frame(effect_15_Matrixwithgroups)
  
  #Human_Data77108_RO7062931_Plasma_0.1mg <- subset(Human_RO7062931PlasmaMatrixwithgroups,X3==1)
  effect_15_Matrixwithgroups1mg <- effect_15_Matrixwithgroups
  effect_15_Matrixwithgroups1mg$X3 = 2
  Datasets$DataEffect_15dose <- effect_15_Matrixwithgroups1mg
  
  return(Datasets)
}
DataEffect_45dose<-function(Datasets){
  # Format: time (hours), 0.1 mg/kg	0.3 mg/kg	1 mg/kg	3 mg/kg	5 mg/kg	15 mg/kg	50 mg/kg
  
  effect_45 = c(  #data is in pmol/L
    1.2,	12.76077941,
    6,	21.46423009,
    12,	30.67806824,
    24,	42.4966943,
    36,	45.49629581,
    48,	52.18673665,
    72,	48.34650472,
    96,	52.74243499,
    168,	50.76821404,
    336,	31.74671827
  )
  
  effect_45_Matrix = matrix(effect_45 ,nrow=20,ncol=2,byrow = TRUE)        # fill matrix by rows 
  
  #convert units
  
  
  #make groups
  
  sizes = dim(effect_45_Matrix)
  effect_45_Matrixwithgroups = c()
  for (i in 2:sizes[2]){
    effect_45_Matrixwithgroups <- rbind(effect_45_Matrixwithgroups,cbind(c(effect_45_Matrix[,1]),c(effect_45_Matrix[,i]),c(rep(i-1, length(effect_45_Matrix[,i])))))    
  }
  effect_45_Matrixwithgroups=data.frame(effect_45_Matrixwithgroups)
  
  #Human_Data77108_RO7062931_Plasma_0.1mg <- subset(Human_RO7062931PlasmaMatrixwithgroups,X3==1)
  effect_45_Matrixwithgroups1mg <- effect_45_Matrixwithgroups
  effect_45_Matrixwithgroups1mg$X3 = 3
  Datasets$DataEffect_45dose <- effect_45_Matrixwithgroups1mg
  
  return(Datasets)
}

loadData<-function(){
  Datasets <- list()
  Datasets <- DataEffect_1dose(Datasets)
  Datasets <- DataEffect_15dose(Datasets)
  Datasets <- DataEffect_45dose(Datasets)
  return(Datasets)
}

#Models
odes_DistributionalDelay <- function (t, A, p) {  # ODE system 
  

  ##assign the current values
  Dose_Central = A[1]
  Drug_Central = A[2]
  Drug_Peri = A[3]
  
  ##assign the parameters
  Central        <- p$Central
  Peripheral     <- p$Peripheral

  ka_Central <- p$ka
  ke_Central <- p$keC
  
  kcp <-p$kcp
  #kpc <-p$kpc
  
  ##reactions
  
  ReactionFlux1 = ka_Central*Dose_Central
  ReactionFlux2 = (ke_Central*Drug_Central)*Central
  ReactionFlux4 = (kcp*Drug_Central)*Central-(kcp*Drug_Peri)*Peripheral
  
  
  ## right hand side
  
  ddtDose_Central = -ReactionFlux1
  ddtDrug_Central =  1/Central*(ReactionFlux1 - ReactionFlux2 - ReactionFlux4)
  ddtAROB = 1/Peripheral*(ReactionFlux4)
  
  deriva <- list ( c (  ddtDose_Central, ddtDrug_Central, ddtAROB) )
  
  #print(deriva)
  
  return ( deriva )
}
odes_TurnoverModel_3_stimulation_of_buildup <- function (t, A, p) {  # ODE system 
  
  
  ##assign the current values
  Dose_Central = A[1]
  Drug_Central = A[2]
  Drug_Peri = A[3]
  Response = A[4]
  
  ##assign the parameters
  Central        <- p$Central
  Peripheral     <- p$Peripheral
  
  ka_Central <- p$ka
  ke_Central <- p$keC
  
  kcp <-p$kcp
  
  kin <- p$Rss*p$kout
  #kin <- p$kin
  kout <- p$kout
  #kpc <-p$kpc
  
  ##reactions
  
  ReactionFlux1 = ka_Central*Dose_Central
  ReactionFlux2 = (ke_Central*Drug_Central)*Central
  ReactionFlux4 = (kcp*Drug_Central)*Central-(kcp*Drug_Peri)*Peripheral
  
  
  
  ## right hand side
  
  ddtDose_Central = -ReactionFlux1
  ddtDrug_Central =  1/Central*(ReactionFlux1 - ReactionFlux2 - ReactionFlux4)
  ddtDrug_Peripheral = 1/Peripheral*(ReactionFlux4)
  ddtResponse = kin*(1+(p$Emax_turnover*Drug_Central/(p$EC50 +Drug_Central))) - kout*Response
  
  deriva <- list ( c (  ddtDose_Central, ddtDrug_Central, ddtDrug_Peripheral, ddtResponse) )
  
  #print(deriva)
  
  return ( deriva )
}
odes_TurnoverModel_2_inhibition_of_loss <- function (t, A, p) {  # ODE system 
  
  
  ##assign the current values
  Dose_Central = A[1]
  Drug_Central = A[2]
  Drug_Peri = A[3]
  Response = A[4]
  
  ##assign the parameters
  Central        <- p$Central
  Peripheral     <- p$Peripheral
  
  ka_Central <- p$ka
  ke_Central <- p$keC
  
  kcp <-p$kcp
  
  kin <- p$Rss*p$kout
  #kin <- p$kin
  kout <- p$kout
  #kpc <-p$kpc
  
  ##reactions
  
  ReactionFlux1 = ka_Central*Dose_Central
  ReactionFlux2 = (ke_Central*Drug_Central)*Central
  ReactionFlux4 = (kcp*Drug_Central)*Central-(kcp*Drug_Peri)*Peripheral
  
  
  
  ## right hand side
  
  ddtDose_Central = -ReactionFlux1
  ddtDrug_Central =  1/Central*(ReactionFlux1 - ReactionFlux2 - ReactionFlux4)
  ddtDrug_Peripheral = 1/Peripheral*(ReactionFlux4)
  ddtResponse = kin - kout*(1-(p$Emax_turnover*Drug_Central/(p$EC50 +Drug_Central)))*Response
  
  deriva <- list ( c (  ddtDose_Central, ddtDrug_Central, ddtDrug_Peripheral, ddtResponse) )
  
  #print(deriva)
  
  return ( deriva )
}

#inital conditions
getInitialConditions <- function (p){
  
  #A_init = c(0,0,0,0,0,p$mRNAss)
  
  A_init = c(0,0,0)
  
  return(A_init)
  
}
getInitialConditionsTurnoverModel <- function (p){
  
  #A_init = c(0,0,0,0)
  
  A_init = c(0,0,0,p$Rss)
  #A_init = c(0,0,0,p$kin/p$kout)
  
  return(A_init)
  
}

#parameters
getAllParametersAndDosing <-function(fromslider){
  
  
  p<-list()
  ## Fixed Parameters
  
   p$ka <- 1.5 #fromslider$ka  #in 1/hour
  
  ## Get all Parameters from the sliders
  
  #which model to use?
  p$effectModel <- fromslider$effectModel
  
  #p$keC <- fromslider$keC #1/h
  p$keC <- 0.04
  
  p$Emax = fromslider$Emax
  p$Emax_turnover = fromslider$Emax_turnover
  
  # p$Central  <- fromslider$Central #in ml/kg
  p$Central <- 100
  #p$Peripheral <- fromslider$Peripheral
  p$Peripheral <- 100
  
  p$kcp <- fromslider$kcp
  #p$kpc <- fromslider$kpc
  
  p$Rss <- fromslider$Rss
  #p$kin <- fromslider$kin
  
  p$kout <- fromslider$kout
  
  p$E0 = fromslider$E0
  #p$Emax = fromslider$Emax
  p$EC50 = fromslider$EC50
  
  ### Dosing and timepoints
  n_doses    <- fromslider$n_doses
  
  dose_cmt   <- 1
  ii         <- fromslider$ii

  dose_times <- seq (from = 0, by=ii, to=(n_doses*ii))
  doseAmount <- as.numeric(fromslider$dose)*1000 #in mg/kg to ug/kg
  dose_amts  <- c(rep (doseAmount, n_doses), 0)
  washouttime <- 4 #time of dosing time
  
  #build sequences with logarithmically after dosing
  times      <- seq(from=0, to=1500, by=1500/(500))  # Integration window and stepsize
 # times      <- c(seq(from=0, to=24, by=0.25), seq(from=24, to=365*24, by=24)) 
  times      <- sort(c(times,dose_times,dose_times+0.1,dose_times+0.2,dose_times+0.3))
  obs_c      <- c(1:2)  # Observation compartments
  #n_ind      <- 1
  
  dosing <- c()
  dosing$ii = ii
  dosing$dose_cmt = dose_cmt
  dosing$dose_times = dose_times
  dosing$dose_amts = dose_amts
  
  simulation <- c()
  simulation$times = times
  simulation$obs_c = obs_c
  
  ## Prepare return structure
  ParametersAndDosing <- c()
  ParametersAndDosing$simulation    = simulation
  ParametersAndDosing$dosing        = dosing
  ParametersAndDosing$p             = p
  
  return(ParametersAndDosing)
  
}
getPlotParameters<- function(fromslider){
  
  #Parameters that influence the visualization are seperated 
  #to model parameters to ensure that the model is not 
  #simulated again if they are changed
  plot_p <- list()
  plot_p$untilTimeToPlot <- fromslider$untilTimeToPlot
   if(fromslider$plasma_log== "Linear"){  
    plot_p$plasma_log = FALSE
  }else{
    plot_p$plasma_log = TRUE
  }

  return(plot_p)
}

#call the model and get back simulation results
getSimulationResults <- function(ParametersAndDosing){
  
  Simulation_setting          <- ParametersAndDosing$simulation
  Dosing                      <- ParametersAndDosing$dosing  
  p                           <- ParametersAndDosing$p   
 
  #Define dosing events
  dosing <- data.frame(var = c("Dose_Central"), time = Dosing$dose_times,value = Dosing$dose_amts, method = c("add"))
  
  if(p$effectModel=="Distributional delay"){
    A_init   <- getInitialConditions(p)  # Initial state of ODE system
    
    #Make species names vector
    yini <- c(Dose_Central = A_init[1], Drug_Central = A_init[2], AROB = A_init[3])
    #Make the integration
    #'method' should be one of “lsoda”, “lsode”, “lsodes”, “lsodar”, “vode”, “daspk”, “euler”, “rk4”, “ode23”, “ode45”, “radau”, “bdf”, “bdf_d”, “adams”, “impAdams”, “impAdams_d”, “iteration”
    des_out <- ode(time = Simulation_setting$times, func=odes_DistributionalDelay, parms = p,y = yini, events= list(data = dosing))
    
    #write output in a dataframe
    dat_ind <- c()
    for (k in 1:length(A_init)) {
      concentrations <- des_out[,(k+1)]
      dat_ind <- rbind (dat_ind, cbind(t=des_out[,1], comp=k, ipred=concentrations))
    }
    comb_dat = data.frame(dat_ind)
    
    toplot <-list()
    toplot$PK_plasma <- subset(comb_dat,comp==2) 
    toplot$PK_peri  <- subset(comb_dat,comp==3) 
    
    toplot$Effect <- toplot$PK_peri
    toplot$Effect$ipred <- p$E0 + (p$Emax-p$E0)*(toplot$PK_peri$ipred)/(p$EC50 +(toplot$PK_peri$ipred))
    
  }else if (p$effectModel=="Turnover model 3 - (Stimulation of build-up)"){
    
    A_init   <- getInitialConditionsTurnoverModel(p)  # Initial state of ODE system
    
    #Make species names vector
    yini <- c(Dose_Central = A_init[1], Drug_Central = A_init[2], AROB = A_init[3], Response = A_init[4])
    #Make the integration
    #'method' should be one of “lsoda”, “lsode”, “lsodes”, “lsodar”, “vode”, “daspk”, “euler”, “rk4”, “ode23”, “ode45”, “radau”, “bdf”, “bdf_d”, “adams”, “impAdams”, “impAdams_d”, “iteration”
    des_out <- ode(time = Simulation_setting$times, func=odes_TurnoverModel_3_stimulation_of_buildup, parms = p,y = yini, events= list(data = dosing))#, method="vode"
    
    #write output in a dataframe
    dat_ind <- c()
    for (k in 1:length(A_init)) {
      concentrations <- des_out[,(k+1)]
      dat_ind <- rbind (dat_ind, cbind(t=des_out[,1], comp=k, ipred=concentrations))
    }
    comb_dat = data.frame(dat_ind)
    
    toplot <-list()
    toplot$PK_plasma <- subset(comb_dat,comp==2) 
    toplot$PK_peri  <- subset(comb_dat,comp==3) 
    toplot$Effect  <- subset(comb_dat,comp==4) 
    
  }else if (p$effectModel=="Turnover model 2 - (Inhibition of loss)"){
    
    A_init   <- getInitialConditionsTurnoverModel(p)  # Initial state of ODE system
    
    #Make species names vector
    yini <- c(Dose_Central = A_init[1], Drug_Central = A_init[2], AROB = A_init[3], Response = A_init[4])
    #Make the integration
    #'method' should be one of “lsoda”, “lsode”, “lsodes”, “lsodar”, “vode”, “daspk”, “euler”, “rk4”, “ode23”, “ode45”, “radau”, “bdf”, “bdf_d”, “adams”, “impAdams”, “impAdams_d”, “iteration”
    des_out <- ode(time = Simulation_setting$times, func=odes_TurnoverModel_2_inhibition_of_loss, parms = p,y = yini, events= list(data = dosing))#, method="vode"
    
    #write output in a dataframe
    dat_ind <- c()
    for (k in 1:length(A_init)) {
      concentrations <- des_out[,(k+1)]
      dat_ind <- rbind (dat_ind, cbind(t=des_out[,1], comp=k, ipred=concentrations))
    }
    comb_dat = data.frame(dat_ind)
    
    toplot <-list()
    toplot$PK_plasma <- subset(comb_dat,comp==2) 
    toplot$PK_peri  <- subset(comb_dat,comp==3) 
    toplot$Effect  <- subset(comb_dat,comp==4) 
    
  }
         
    
  return(toplot)
  
}

#convert the units from the simulation to the desired output units
convertUnits <- function(SimulationResults,p){
  
  #convert Plasma from nmol/ml to ug/ml for plotting
  SimulationResults$PK_plasma$ipred = units(units(SimulationResults$PK_plasma$ipred,'mol','gram',p$MW),'n','u') 
  SimulationResults$PK_plasma_UB$ipred = units(units(SimulationResults$PK_plasma_UB$ipred,'mol','gram',p$MW),'n','u') 
  SimulationResults$PK_plasma_LB$ipred = units(units(SimulationResults$PK_plasma_LB$ipred,'mol','gram',p$MW),'n','u') 
  
  #convert Liver cons from nmol/g to ug/g for plotting
  SimulationResults$PK_liver$ipred = units(units(SimulationResults$PK_liver$ipred,'mol','gram',p$MW),'n','u') 
  SimulationResults$PK_liver_UB$ipred = units(units(SimulationResults$PK_liver_UB$ipred,'mol','gram',p$MW),'n','u') 
  SimulationResults$PK_liver_LB$ipred = units(units(SimulationResults$PK_liver_LB$ipred,'mol','gram',p$MW),'n','u') 
  
  #convert Liver cons from hours to days for plotting
  SimulationResults$PK_liver$t = SimulationResults$PK_liver$t
  SimulationResults$PK_liver_UB$t = SimulationResults$PK_liver_UB$t
  SimulationResults$PK_liver_LB$t = SimulationResults$PK_liver_LB$t
  
  #convert Kidney cons from nmol/g to ug/g for plotting
  SimulationResults$CKidney$ipred = units(units(SimulationResults$CKidney$ipred,'mol','gram',p$MW),'n','u') 
  SimulationResults$CKidney_UB$ipred = units(units(SimulationResults$CKidney_UB$ipred,'mol','gram',p$MW),'n','u') 
  SimulationResults$CKidney_LB$ipred = units(units(SimulationResults$CKidney_LB$ipred,'mol','gram',p$MW),'n','u') 
  
  #convert Kidney cons from hours to days for plotting
  SimulationResults$CKidney$t = SimulationResults$CKidney$t
  SimulationResults$CKidney_UB$t = SimulationResults$CKidney_UB$t
  SimulationResults$CKidney_LB$t = SimulationResults$CKidney_LB$t
  
  return(SimulationResults)
}

#Plotting functions
base_breaks <- function(n = 100){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
  
  #....“I don't even know what that means
  #No one knows what it means, but it's provocative....
  #....No, it's not, it's gross 
  #It's gets the people going...
  
}
base_breaks_log <- function(n = 20){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}
base_breaks_lin <- function(n = 20){
  function(x) {
    allticks<-axisTicks(range(x, na.rm = TRUE), log = FALSE, n = n)
    allPositiveticks = allticks[allticks>=0]
    return(allPositiveticks)
  }
}

makePKplot_plasma <- function(SimulationResults,plot_p){
  
  
  #make a prediction dataframe to show the confidence bands
 # predframe <- data.frame(SimulationResults$PK_plasma$t,lwr=SimulationResults$PK_plasma_LB$ipred,upr=SimulationResults$PK_plasma_UB$ipred)
  predframe_central <- data.frame(x=SimulationResults$PK_plasma$t, y=SimulationResults$PK_plasma$ipred, comp=1)
  predframe_peripheral <- data.frame(x=SimulationResults$PK_peri$t, y=SimulationResults$PK_peri$ipred,comp=2)
  
  #build the plot
  
  xlabel <- "Time in hours"
  ylabel <- "concentration (ug/ml)"
  
  
  # myplot_step1 <- ggplot(predframe_central,na.rm=TRUE)+
  # geom_line(size=2, na.rm=TRUE,aes(y=SimulationResults$PK_plasma$ipred, x=SimulationResults$PK_plasma$t, colour = "#ed4354"))+
  # #geom_ribbon(aes(ymin=lower, ymax=upper, x=SimulationResults$PK_plasma$t, fill = "Range"), alpha = 0.3, colour = "#ed4354")+
  # scale_colour_manual("",values="#ed4354",guide=FALSE)+
  # scale_fill_manual("",values="#ed4354",guide=FALSE)+
  # #scale_x_continuous(limits=c(0,plot_p$untilTimeToPlot))+
  #   theme_bw()+
  #   xlab(xlabel) +
  #   ylab(ylabel)+
  #   theme(axis.title.y = element_text(size = rel(1.5), angle = 90, vjust=2),axis.title.x = element_text(size = rel(1.5), vjust=-0.5),legend.position = c(0.9, 0.95),legend.justification=c(1, 1))
  # 
  
  myplot_step1 <- ggplot(data=predframe_central,aes(x=x, y=y,group=as.factor(comp),color=as.factor(comp),linetype=as.factor(comp)))+
    geom_line(size=2)+
    xlab(xlabel) +
    ylab(ylabel)+
    scale_color_manual(values=c("#ED4354", "#72BF43", "#ED4354", "#ED4354", "#72BF43"),name="Simulated\nconcentrations",breaks=c("1", "2"),labels=c("Central compartment", "Peripheral compartment"))+
    scale_linetype_manual(values=c(1, 2),name="Simulated\nconcentrations",breaks=c("1", "2"),labels=c("Central compartment", "Peripheral compartment"))+
    theme_bw()+
    theme(axis.title.y = element_text(size = rel(1.5), angle = 90, vjust=2),axis.title.x = element_text(size = rel(1.5), vjust=-0.5),legend.position = c(0.85, 0.98),legend.justification=c(1, 1))
  
    myplot_step1 <- myplot_step1 + geom_line(data=predframe_peripheral,aes(x=x, y=y,group=as.factor(comp),color=as.factor(comp),linetype=as.factor(comp)),size=2)
  
  
  
  # myplot_step1 <- ggplot(data=data.frame(SimulationResults$PK_plasma),aes(x=t, y=ipred,group=as.factor(comp),colour=as.factor(comp)), na.rm=TRUE)+
  #   geom_line(size=2, na.rm=TRUE)+ 
  #   geom_ribbon(data=predframe,aes(ymin=lwr,ymax=upr),alpha=0.3)+
  #   #geom_line(data=SimulationResults$Ab.mean,aes(x=t,y=x),color="#ed4354",size=2, na.rm=TRUE)+       
  #   theme_bw()+
  #   xlab(xlabel) +
  #   ylab(ylabel)+
  #   scale_color_manual(values=c("#ed4354","#000000"),name="Simulated\nconcentrations",labels=c("Concentration in blood"),guide=FALSE)+
  #   #scale_y_continuous(trans = log_trans(), breaks = base_breaks(),
  #   #                   labels = prettyNum, limits=c(0.0001,NA))+
  #   scale_x_continuous(limits=c(0,plot_p$untilTimeToPlot))+
  #   theme(axis.title.y = element_text(size = rel(1.5), angle = 90, vjust=2),axis.title.x = element_text(size = rel(1.5), vjust=-0.5),legend.position = c(0.9, 0.95),legend.justification=c(1, 1))
  # 
  if (plot_p$plasma_log) {
    myplot_step1 <- myplot_step1+
      scale_y_continuous(trans = log_trans(), breaks = base_breaks_log(),
                         labels = prettyNum)+
      coord_cartesian(ylim = c(0.01,1000), xlim =c(0,plot_p$untilTimeToPlot)) 
  }else{
    myplot_step1 <- myplot_step1+
      scale_y_continuous(breaks = base_breaks_lin(),
                         labels = prettyNum)+
      coord_cartesian(xlim =c(0,plot_p$untilTimeToPlot)) 
  }
  
  
  
  
}
makeEffektplot <- function(SimulationResults,plot_p){
  
  
  #make a prediction dataframe to show the confidence bands
  # predframe <- data.frame(SimulationResults$PK_plasma$t,lwr=SimulationResults$PK_plasma_LB$ipred,upr=SimulationResults$PK_plasma_UB$ipred)
  predframe_effect <- data.frame(x=SimulationResults$Effect$t, y=SimulationResults$Effect$ipred, comp=1)
  
  #build the plot
  
  xlabel <- "Time in hours"
  ylabel <- "effect"
  
  
  # myplot_step1 <- ggplot(predframe_central,na.rm=TRUE)+
  # geom_line(size=2, na.rm=TRUE,aes(y=SimulationResults$PK_plasma$ipred, x=SimulationResults$PK_plasma$t, colour = "#ed4354"))+
  # #geom_ribbon(aes(ymin=lower, ymax=upper, x=SimulationResults$PK_plasma$t, fill = "Range"), alpha = 0.3, colour = "#ed4354")+
  # scale_colour_manual("",values="#ed4354",guide=FALSE)+
  # scale_fill_manual("",values="#ed4354",guide=FALSE)+
  # #scale_x_continuous(limits=c(0,plot_p$untilTimeToPlot))+
  #   theme_bw()+
  #   xlab(xlabel) +
  #   ylab(ylabel)+
  #   theme(axis.title.y = element_text(size = rel(1.5), angle = 90, vjust=2),axis.title.x = element_text(size = rel(1.5), vjust=-0.5),legend.position = c(0.9, 0.95),legend.justification=c(1, 1))
  # 
  
  myplot_step1 <- ggplot(data=predframe_effect,aes(x=x, y=y))+
    geom_line(size=2)+
    xlab(xlabel) +
    ylab(ylabel)+
  #  scale_color_manual(values=c("#ED4354", "#72BF43", "#ED4354", "#ED4354", "#72BF43"),name="Simulated\nconcentrations",breaks=c("1", "2"),labels=c("Central compartment", "Peripheral compartment"))+
  #  scale_linetype_manual(values=c(1, 1),name="Simulated\neffect",breaks=c("1", "2"),labels=c("Central compartment", "Peripheral compartment"))+
    theme_bw()+
    theme(axis.title.y = element_text(size = rel(1.5), angle = 90, vjust=2),axis.title.x = element_text(size = rel(1.5), vjust=-0.5),legend.position = c(0.85, 0.98),legend.justification=c(1, 1))
  
 # myplot_step1 <- myplot_step1 + geom_line(data=predframe_peripheral,aes(x=x, y=y,group=as.factor(comp),color=as.factor(comp),linetype=as.factor(comp)),size=2)
  
  
  
  # myplot_step1 <- ggplot(data=data.frame(SimulationResults$PK_plasma),aes(x=t, y=ipred,group=as.factor(comp),colour=as.factor(comp)), na.rm=TRUE)+
  #   geom_line(size=2, na.rm=TRUE)+ 
  #   geom_ribbon(data=predframe,aes(ymin=lwr,ymax=upr),alpha=0.3)+
  #   #geom_line(data=SimulationResults$Ab.mean,aes(x=t,y=x),color="#ed4354",size=2, na.rm=TRUE)+       
  #   theme_bw()+
  #   xlab(xlabel) +
  #   ylab(ylabel)+
  #   scale_color_manual(values=c("#ed4354","#000000"),name="Simulated\nconcentrations",labels=c("Concentration in blood"),guide=FALSE)+
  #   #scale_y_continuous(trans = log_trans(), breaks = base_breaks(),
  #   #                   labels = prettyNum, limits=c(0.0001,NA))+
  #   scale_x_continuous(limits=c(0,plot_p$untilTimeToPlot))+
  #   theme(axis.title.y = element_text(size = rel(1.5), angle = 90, vjust=2),axis.title.x = element_text(size = rel(1.5), vjust=-0.5),legend.position = c(0.9, 0.95),legend.justification=c(1, 1))
  # 
  if (plot_p$plasma_log) {
    myplot_step1 <- myplot_step1+
      scale_y_continuous(trans = log_trans(), breaks = base_breaks_log(),
                         labels = prettyNum)+
      coord_cartesian(ylim = c(1,120), xlim =c(0,plot_p$untilTimeToPlot)) 
  }else{
    myplot_step1 <- myplot_step1+
      scale_y_continuous(breaks = base_breaks_lin(),
                         labels = prettyNum)+
      coord_cartesian(xlim =c(0,plot_p$untilTimeToPlot)) 
  }
  
  
  
  
}
addPKdataToplot_plasma <- function(myplot_step1,Datasets,input){
  
  #1mg/kg dose plotting
  if (input$Dataset=="Step 1"){   
    myplot_step2 <- myplot_step1+ 
      geom_point(data=data.frame(Datasets$DataPK_1dose), aes(x=X1,y=X2,group=as.factor(X3),shape=as.factor(X3),linetype=NULL,color=NULL),size=5, na.rm=TRUE)               
  }else{
    myplot_step2 <- myplot_step1
  }
  
  
  if (input$Dataset=="Step 1"){
    myplot_step_end <- myplot_step2+  
      scale_shape_manual(values=c(1),name="Measured\nconcentrations",breaks=c("1"),
                         labels=c("Dose 1"))
      # +scale_fill _manual(values=c(1, 2, 3, 4, 5, 6, 7 ,8 , 9, 10, 11, 12, 13 , 14, 15, 16),name="Measured\nconcentrations",breaks=c("1", "2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"),
      #                   labels=c("1 mg/kg Shionogi LNA (mouse)", "5 mg/kg Shionogi LNA (mouse)", "10 mg/kg Shionogi LNA (mouse)",
      #                            "0.3 mg/kg Roche naked LNA (rat)", "1 mg/kg Roche naked LNA (rat)", "5 mg/kg Roche naked LNA (rat)",
      #                            "0.45 mg/kg Roche GalNac-LNA (rat)","1.5 mg/kg Roche GalNac-LNA (rat)", "7.5 mg/kg Roche GalNac-LNA (rat)",
      #                            "0.1 mg/kg Roche GalNac-LNA (cyno)","0.3 mg/kg Roche GalNac-LNA (cyno)","1 mg/kg Roche GalNac-LNA (cyno)",
      #                            "3 mg/kg Roche GalNac-LNA (cyno)","5 mg/kg Roche GalNac-LNA (cyno)","15 mg/kg Roche GalNac-LNA (cyno)",
      #                            "50 mg/kg Roche GalNac-LNA (cyno)"))       
  }
  else{
    myplot_step_end <- myplot_step2
  }
  
  
  #print the final plot
  return(myplot_step_end)
  
}
addEffectdataToplot_effect <- function(myplot_step1,Datasets,input){
  
  
  #1mg/kg dose plotting
  if (input$Dataset=="Step 1 - find best model with one dose" | input$Dataset=="Step 2 - find best model with three doses"){   
    myplot_step2 <- myplot_step1+ 
      geom_point(data=data.frame(Datasets$DataEffect_1dose), aes(x=X1,y=X2,group=as.factor(X3),shape=as.factor(X3),linetype=NULL,color=NULL),size=5, na.rm=TRUE)               
  }else{
    myplot_step2 <- myplot_step1
  }
  
  if (input$Dataset=="Step 2 - find best model with three doses"){   
    myplot_step3 <- myplot_step2+ 
      geom_point(data=data.frame(Datasets$DataEffect_15dose), aes(x=X1,y=X2,group=as.factor(X3),shape=as.factor(X3),linetype=NULL,color=NULL),size=5, na.rm=TRUE)+
      geom_point(data=data.frame(Datasets$DataEffect_45dose), aes(x=X1,y=X2,group=as.factor(X3),shape=as.factor(X3),linetype=NULL,color=NULL),size=5, na.rm=TRUE)
  }else{
    myplot_step3 <- myplot_step2
  }
  
  if (input$Dataset=="Step 3 - multiple dosing and efficacy"){   
    myplot_step4 <- myplot_step3+
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=47, ymax=52, alpha=0.2, fill="green")+
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=52, ymax=Inf, alpha=0.2, fill="red")
    coord_cartesian(ylim = c(0,60))
  }else{
    myplot_step4 <- myplot_step3
  }
  
  if (input$Dataset=="Step 1 - find best model with one dose"| input$Dataset=="Step 2 - find best model with three doses" ){
    myplot_step_end <- myplot_step4+  
      scale_shape_manual(values=c(1,2,3),name="Measured\nconcentrations",breaks=c("1","2","3"),
                         labels=c("1 mg/kg", "15 mg/kg", "45 mg/kg"))
    # +scale_fill _manual(values=c(1, 2, 3, 4, 5, 6, 7 ,8 , 9, 10, 11, 12, 13 , 14, 15, 16),name="Measured\nconcentrations",breaks=c("1", "2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"),
    #                   labels=c("1 mg/kg Shionogi LNA (mouse)", "5 mg/kg Shionogi LNA (mouse)", "10 mg/kg Shionogi LNA (mouse)",
    #                            "0.3 mg/kg Roche naked LNA (rat)", "1 mg/kg Roche naked LNA (rat)", "5 mg/kg Roche naked LNA (rat)",
    #                            "0.45 mg/kg Roche GalNac-LNA (rat)","1.5 mg/kg Roche GalNac-LNA (rat)", "7.5 mg/kg Roche GalNac-LNA (rat)",
    #                            "0.1 mg/kg Roche GalNac-LNA (cyno)","0.3 mg/kg Roche GalNac-LNA (cyno)","1 mg/kg Roche GalNac-LNA (cyno)",
    #                            "3 mg/kg Roche GalNac-LNA (cyno)","5 mg/kg Roche GalNac-LNA (cyno)","15 mg/kg Roche GalNac-LNA (cyno)",
    #                            "50 mg/kg Roche GalNac-LNA (cyno)"))       
  }
  else{
    myplot_step_end <- myplot_step4
  }
  
  
  #print the final plot
  return(myplot_step_end)
  
}
returnRawData <- function(SimulationResults,plot_p){
  
  #
  #howlong <- max(c(plot_p$untilTimeToPlot,plot_p$untilTimeToPlot_liver*24,plot_p$untilTimeToPlot_kidney*24))
  
  newdata_PK_plasma <- subset(SimulationResults$PK_plasma)
  newdata_PK_effect <- subset(SimulationResults$Effect)
  
  datatable_temp <- t(rbind(newdata_PK_plasma$t, newdata_PK_plasma$ipred, newdata_PK_effect$ipred))
  dimnames(datatable_temp) <- list(rownames(datatable_temp, do.NULL = FALSE, prefix = "col"),c("Time (h)","Plasma cons (ug/ml)", "Effect"))
  
  #predframe <- data.frame(x=SimulationResults$PK_plasma$t, y=SimulationResults$PK_plasma$ipred, lower=SimulationResults$PK_plasma_LB$ipred, upper=SimulationResults$PK_plasma_UB$ipred)
  #datatable_temp <- datatable
  #select the lines that are interesting
  #datatable_temp <- table(t(rbind(head(SimulationResults$PK_plasma$t), head(SimulationResults$PK_plasma$ipred))))
  #dimnames(datatable_temp) = list( 
  #  +   c("row1", "row2"),         # row names 
  #  +   colnames(datatable_temp, do.NULL = FALSE, prefix = "col")) # column names 
  #print(datatable_temp)
  # datatable_temp[13]<-NULL #remove the too long columns 4,5,13 and 14
  # datatable_temp[13]<-NULL
  # datatable_temp[13]<-NULL
  # datatable_temp[4]<-NULL #remove the original time (keep only the time in hours)
  # datatable_temp[4]<-NULL
  #make a prediction dataframe to show the confidence bands
  # predframe <- data.frame(SimulationResults$PK_plasma$t,lwr=SimulationResults$PK_plasma_LB$ipred,upr=SimulationResults$PK_plasma_UB$ipred)
  
  
  return(datatable_temp)
  
}

##########To be commented out when used in a shiny app

# input<-list()
# 
# input$Dataset <-1
# 
# #dosing
# input$n_doses    <- 5
# input$ii  <- 24
# input$interval   <- 7
# input$dose       <- 0.1
# 
# #For PK
# input$ka         <- 8.5  #in 1/hour
# input$Central     <- 35 #in ml
# input$Peripheral <- 2
# input$keC <- 0.05
# input$kcp <- 0.1
# #input$kpc <- 0.1
# 
# #forEffect
# input$E0 = 0
# # input$Emax = 100
# input$EC50 = 20
# 
# input$Rss = 0.5
# input$kout = 2
# 
# input$plasma_log <- 'Linear'
# input$untilTimeToPlot <- 300
# 
# input$effectModel<- "Turnover model 2 - (Inhibition of loss)"
# ParametersAndDosing     <- getAllParametersAndDosing(input)
# SimulationResults       <- getSimulationResults(ParametersAndDosing)
# #SimulationResults       <- convertUnits(SimulationResults,ParametersAndDosing$p)
# plot_p                  <- getPlotParameters(input)
# Datasets                <- loadData()
# pk_plot                 <- makePKplot_plasma(SimulationResults,plot_p)
# effect_plot             <- makeEffektplot(SimulationResults,plot_p)
# pk_plot_withdata        <- addPKdataToplot_plasma(pk_plot,Datasets,input)
# #pk_plot_liver           <- makePKplot_liver(SimulationResults,plot_p)
# #pk_plot_liver_withdata  <- addPKdataToplot_liver(pk_plot_liver,Datasets,input)
# #pk_plot_kidney           <- makePDplot_kidney(SimulationResults,plot_p)
# #pk_plot_kidney_withdata  <- addPKdataToplot_kidney(pk_plot_kidney,Datasets,input)

##############end of commenting out





