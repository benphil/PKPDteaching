library(shiny)
library (deSolve)
library (MASS)
library (lattice)
library(ggplot2)
library(plyr)    ## for the ldply and ddply aggregation functions
library(scales)
library(dplyr) 

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

#pdf(file=paste(Dose,"mg.",TVKD,"nM Kd.","k_Tar=",k_Tar,".k_mAb=",k_mAb,".Kcomplex=",k_comp,".",SIM_TITLE,
#               ".","PK, target and complex plots with PI using TMDD model",".n=",n.ind,".",SIM_DAT,".pdf", sep=""))


# ## Target and complex profile (Plasma)
# pl <- ddply(out,"time",function(d){
#   data.frame(low = quantile(d$ligandP,probs=0.025),
#              mid = quantile(d$ligandP,probs=0.5),
#              hi = quantile(d$ligandP,probs=0.975))
# })
# 
# cmxpl <- ddply(out,"time",function(d){
#   data.frame(lowcmx = quantile(d$complexP,probs=0.025), #No need to divide by VC as complex is already as conc
#              midcmx = quantile(d$complexP,probs=0.5),
#              hicmx = quantile(d$complexP,probs=0.975))
# })
# Tarpl <- cbind(pl, cmxpl)

shinyServer(function(input, output, session) {
  
  simulate <- eventReactive({input$startSim 
    session$clientData},{
    
    #print("Simulate")
    ### Dose and time Settings
    
    n_doses    <- input$n_doses-1
    
    ii         <- input$interval*24
    dose_times <- seq (from = 0, by=ii, to=n_doses*ii)
    
    SIM_TITLE <- "mAB multiple dose" # give each simualtion a unique title
    Version <- "1" #Simulation version number
    SIM_DAT <- Sys.Date()
    
    ####################### Parameters to adjust for simulations
    #### mAb binding kinetics
    TVKD <- input$kd #0.78/10 # nM 0.78nM as default
    ####Target dynamics
    BASE_CONC <- (input$ligandss/1000)/1000 #0.011/1000  ## nM. 11fM as default value
    HL_Tar <- input$targetTt12 #0.5 # hours Target turnover half life 0.5 based on Biesma et al
    Comp_factor <- input$Comp_factor # fold higher complex elimination compared to mAb.
    #### Dose
    Dose <- input$dose # mg. total dose of mAb
    YTE_factor <- 3   # set to 1 for non YTE sims. set to 3 for 3x decrease in CL
    TAU <- ii #24*7*4*3   #Dosing interval (days)
    
    ##############  Do not alter anyting below this line!
    #### 
    Kon  <- 82   # nM-1 s-1 - fixed at estimate for IMA-026 and Koff derived. Keep fixed and just vary KD as this is OK for this level
    # of investigation given that the model is essntially not be sensitive to Koff
    set.seed(1234567)
    
    washouttime <- 2 #time of dosing time
    #SIM_DURATION <- ii*24*n_doses+ii*24*washouttime
    
    SIM_DURATION <- 24*7*4*6*3 #Duration of simulation in hours
    DT <- 10  # time step for simulations (hours)
    ## If you wish to give multiple doses
    
    FIRST_DOSE <- 0 #time of first dose (hours)
    LAST_DOSE <- FIRST_DOSE +24*7*4*8   #time of last dose (hours)
    Time_labels <- seq(FIRST_DOSE,SIM_DURATION,by=TAU)/24
    #
    n.ind <- 5 # number of individuals in simulation ~ >500 for final simulations but wil take time so use ~30 to check model is working first
    
    
    ## SC
    ## Define dosing events longitudinally (like in nonmem)
    Mw <- 146 #kDa - molecular weight
    
    
    dose_times <- seq (from = 0, by=ii, to=(n_doses*ii))
    times      <- seq(from=0, to=SIM_DURATION, by=(SIM_DURATION)/(300))  # Integration window and stepsize
    times      <- sort(c(times,dose_times,dose_times+0.1,dose_times+0.2,dose_times+0.3))
    
    dose.events <- data.frame(var="depot",
                              time=dose_times, #If you wish to give multiple doses
                              #time=FIRST_DOSE,                        #single dose - just gives dose at time specified
                              value=Dose*1000/Mw, method="add") # dose express in nMol
    
    ## ## parameter values
    ## PK assumed to be equivalent to IMA-026
    TVCL <- input$cl/(1000*24)  # from ml/day to L/hr 
    TVVC <- 3.22        # L - central volume
    TVVP <- 4.06        # L - peripheral volume
    TVQ  <- 0.019      # L/hr
    TVKA <- 0.2/24       # hr-1 - absorption rate constant
    BIOAV <- 0.75       #F
    
    #### Target parameters
    Kdeg <- 0.693/HL_Tar #hr-1 
    BASE_P <- BASE_CONC*TVVC # nMoles
    Rsyn_P <- BASE_P * Kdeg 
    Kel_mAB <- TVCL/TVVC
    Kcomplex <- Kel_mAB * Comp_factor
    ## Calculate maximal theoretical accumulation of complex
    Max_acc_ratio <- Kdeg / Kcomplex #maximum complex accumulation ratio
    Max_acc_pM <- signif(BASE_CONC*(Kdeg / Kcomplex)*1000,3) #maximum complex accumulation concentration
    
    ## PKPD with TMDD in plasma
    pkeq <- function(t, y, p) {    
      
      #### PK param  
      CL <- p[1]
      VC <- p[2]
      VP <- p[3]
      Q <- p[4]
      KA <- p[5]
      KD <- p[6]
      Koff <- Kon*KD
      
      dy1 <- - KA* y[1]         #depot drug amount
      
      dy2 <- BIOAV*KA*y[1] - (CL/VC)*y[2] -  Kon*(y[4]/VC)*(y[2]/VC)*VC + Koff*y[5]*VC - (Q/VC)*y[2] + (Q/VP)*y[3]  # drug amount in plasma amount
      
      dy3 <- (Q/VC)*y[2] - (Q/VP)*y[3] # drug amount in peripheral cmt
      
      dy4 <- Rsyn_P - Kdeg*y[4] - Kon*(y[4]/VC)*(y[2]/VC)*VC + Koff*y[5]*VC #free target amount in plasma
      
      dy5 <- Kon*(y[4]/VC)*(y[2]/VC) - Koff*y[5] - Kcomplex*y[5]# complex conc in plasma
      
      list(c(dy1, dy2, dy3, dy4, dy5))
    }
    
    ## define matrix of parameters (see desolve guide)
    
    # params <- mvrnorm(n.ind, mu = log(c(TVCL,TVVC,TVVP,TVQ,TVKA,TVKD)),
    #                   Sigma = matrix(c(0.04, 0,    0,    0,    0,    0,
    #                                    0,    0.01, 0,    0,    0,    0,
    #                                    0,    0,    0.01, 0,    0,    0,
    #                                    0,    0,    0,    0.01, 0,    0,
    #                                    0,    0,    0,    0,    0.01, 0,
    #                                    0,    0,    0,    0,    0,    0.25), nrow=6))
    # 
    
    #params <- exp(params) #for multiple simulations
    params <- t(c(TVCL,TVVC,TVVP,TVQ,TVKA,TVKD)) #for single simulation
    
    ## Initial values of state variables
    yini <- c(depot = 0, plasma = 0, peripheral = 0,ligandP = BASE_P, complexP = 0)
    ## observation times
    #times <- seq(0,SIM_DURATION,DT) # sequence starting at 0 ending at defined duration going up by x hours
    
    start.time <- Sys.time()
    
    ## loop through matrix row by row (i.e. by ID) output data.frame as you go
    out <- ldply(1:nrow(params),function(i){
      
      outi <- ode(func = pkeq, times = times,
                  y = yini, parms = params[i,],events=list(data=dose.events), method="lsode")
      
      outi <- as.data.frame(outi)
      outi$plasma <- outi$plasma/params[i, 2]  # divides amount by individual V1 to get conc
      outi$ligandP <- outi$ligandP/params[i, 2]  # divides amount by individual V1 to get conc
      outi$SIM <- i
      outi
    })
    
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(c("Simulation time:",time.taken))
    
    allstuffforplotting<-list()
    
    allstuffforplotting$out <- out
    allstuffforplotting$Time_labels <- Time_labels
    allstuffforplotting$Kdeg <- Kdeg
    allstuffforplotting$TVCL <- TVCL
    allstuffforplotting$TVVC <- TVVC
    allstuffforplotting$Kcomplex <- Kcomplex
    #allstuffforplotting$mAbpl <- mAbpl
    allstuffforplotting$SIM_DURATION <- SIM_DURATION
    allstuffforplotting$BASE_CONC <- BASE_CONC
    allstuffforplotting$Max_acc_pM <- Max_acc_pM
    
    return(allstuffforplotting) 
  })
  
  output$pk_plot <- renderPlot({
    
    allstuffforplotting<-simulate()
    
     out <- allstuffforplotting$out
     Time_labels <- allstuffforplotting$Time_labels
     Kdeg <- allstuffforplotting$Kdeg 
     TVCL <-   allstuffforplotting$TVCL
     TVVC <- allstuffforplotting$TVVC
     Kcomplex <- allstuffforplotting$Kcomplex
     mAbpl <- allstuffforplotting$mAbpl
     SIM_DURATION <- allstuffforplotting$SIM_DURATION
     BASE_CONC <- allstuffforplotting$BASE_CONC 
     Max_acc_pM <- allstuffforplotting$Max_acc_pM
    
    
    ## plot the results
    ## mAb profile (Plasma)
    # mAbpl <- ddply(out,"time",function(d){
    #   data.frame(lowpl = quantile(d$plasma,probs=0.025),
    #              midpl = quantile(d$plasma,probs=0.5),
    #              hipl = quantile(d$plasma,probs=0.975))
    # })
    
    ## Plot mAb plasma levels
    ##get key infor into pdf title and ensure that numbers are no too long
    k_Tar <- signif(Kdeg,3)
    k_mAb <- signif(TVCL/TVVC,3)
    k_comp <- signif(Kcomplex,3)
    
    mAbpl<-mutate(out,midpl=plasma,lowpl=plasma,hipl=plasma) #for single simulation plotting
    
    myplot_step1 <-ggplot(mAbpl,aes(x=time/24,y=midpl)) +
     # geom_ribbon(aes(ymin=lowpl,ymax=hipl),alpha=0.4,fill = "red") +
      #scale_x_continuous(breaks=c(0,7,14,21,28,56,28*3,28*4,28*5,28*6), limits=c(0,SIM_DURATION/24)) + #Option to put in specific points
      scale_x_continuous(breaks=Time_labels, limits=c(0,SIM_DURATION/24)) +
      geom_line(size=2,color="red") + 
      scale_y_log10(breaks=c(10,30,100,300,600), limits=c(6,600)) +
      labs(x="Time (Days)", y="mAb Concentration (nM)", 
           title ="Predicted mAb profile in Plasma") + theme_bw()
    
    # ## Target and complex profile (Plasma)
    # pl <- ddply(out,"time",function(d){
    #   data.frame(low = quantile(d$ligandP,probs=0.025),
    #              mid = quantile(d$ligandP,probs=0.5),
    #              hi = quantile(d$ligandP,probs=0.975))
    # })
    # 
    # cmxpl <- ddply(out,"time",function(d){
    #   data.frame(lowcmx = quantile(d$complexP,probs=0.025), #No need to divide by VC as complex is already as conc
    #              midcmx = quantile(d$complexP,probs=0.5),
    #              hicmx = quantile(d$complexP,probs=0.975))
    # })
    # Tarpl <- cbind(pl, cmxpl)
    # 
    # myplot_step1 <- ggplot(Tarpl,aes(x=time/24,y=mid*1000)) +
    #   geom_ribbon(aes(ymin=low*1000,ymax=hi*1000),alpha=0.4, fill = "red") +
    #   geom_ribbon(aes(ymin=lowcmx*1000,ymax=hicmx*1000),alpha=0.4,fill = "blue") +
    #   geom_line(aes(x=time/24,y=midcmx*1000)) +
    #   #scale_x_continuous(breaks=c(0,7,14,21,28,56,28*3,28*4,28*5,28*6), limits=c(0,SIM_DURATION/24)) + #Option to put in specific points
    #   scale_x_continuous(breaks=Time_labels, limits=c(0,SIM_DURATION/24)) +
    #   geom_line() + 
    #  # scale_y_log10(breaks=sort(c(0.0001,0.001,0.01,0.1,1,3,10,30,100,300,1000,3000,10000,30000,100000,Max_acc_pM)), limits=c(0.001,100000)) +
    #   geom_line(aes(y=Max_acc_pM),lty=2,size=1, col ="red") + 
    #   annotate("text", x=(SIM_DURATION*0.9/24), y=Max_acc_pM*1.25, label="max acc (pM)", col= "red") +
    #   geom_line(aes(y=BASE_CONC*1000),lty=2,size=1, col ="red") + 
    #   annotate("text", x=(SIM_DURATION*0.9/24), y=BASE_CONC*1300, label="Baseline (pM)", col= "red") +
    #   geom_line(aes(y=BASE_CONC*100),lty=2,size=1, col ="dark green") +
    #   annotate("text", x=(SIM_DURATION*0.9/24), y=BASE_CONC*80, label="90% supression", col= "dark green") +
    #   labs(x="Time (Days)", y="Free ligand or complex Concentration in Plasma (pM)", 
    #        title ="Predicted Free Ligand or complex profile in Plasma (pM)") + theme_bw()+
    #   scale_y_continuous(trans = log_trans(), breaks = base_breaks_log(),limits=c(BASE_CONC*1300/100,NA),
    #                      labels = prettyNum)
    
  
    #print the final plot
    print(myplot_step1)
   
    
  })
  
  output$pd_plot <- renderPlot({
    
    allstuffforplotting<-simulate()
    
    out <- allstuffforplotting$out
    Time_labels <- allstuffforplotting$Time_labels
    Kdeg <- allstuffforplotting$Kdeg 
    TVCL <-   allstuffforplotting$TVCL
    TVVC <- allstuffforplotting$TVVC
    Kcomplex <- allstuffforplotting$Kcomplex
    #mAbpl <- allstuffforplotting$mAbpl
    SIM_DURATION <- allstuffforplotting$SIM_DURATION
    BASE_CONC <- allstuffforplotting$BASE_CONC 
    Max_acc_pM <- allstuffforplotting$Max_acc_pM
    
    # 
    # ## Target and complex profile (Plasma)
    # pl <- ddply(out,"time",function(d){
    #   data.frame(low = quantile(d$ligandP,probs=0.025),
    #              mid = quantile(d$ligandP,probs=0.5),
    #              hi = quantile(d$ligandP,probs=0.975))
    # })
    # 
    # cmxpl <- ddply(out,"time",function(d){
    #   data.frame(lowcmx = quantile(d$complexP,probs=0.025), #No need to divide by VC as complex is already as conc
    #              midcmx = quantile(d$complexP,probs=0.5),
    #              hicmx = quantile(d$complexP,probs=0.975))
    # })
    # Tarpl <- data.frame(cbind(pl, cmxpl))
    
    Tarpl<-mutate(out,mid=ligandP,low=ligandP,hi=ligandP,midcmx=complexP,lowcmx=complexP,hicmx=complexP) #for single simulation plotting
    #print(Tarpl)
    
    myplot_step1 <- ggplot(data=Tarpl,aes(x=time/24,y=mid*1000)) +
      #geom_ribbon(aes(ymin=low*1000,ymax=hi*1000),alpha=0.4, fill = "red") +
      #geom_ribbon(aes(ymin=lowcmx*1000,ymax=hicmx*1000),alpha=0.4,fill = "blue") +
      geom_line(aes(x=time/24,y=midcmx*1000),size=2,color="green") +
      #scale_x_continuous(breaks=c(0,7,14,21,28,56,28*3,28*4,28*5,28*6), limits=c(0,SIM_DURATION/24)) + #Option to put in specific points
      scale_x_continuous(breaks=Time_labels, limits=c(0,SIM_DURATION/24)) +
      geom_line(size=2,color="blue") +
      # scale_y_log10(breaks=sort(c(0.0001,0.001,0.01,0.1,1,3,10,30,100,300,1000,3000,10000,30000,100000,Max_acc_pM)), limits=c(0.001,100000)) +
      geom_line(aes(y=Max_acc_pM),lty=2,size=1, col ="red") +
      annotate("text", x=(SIM_DURATION*0.9/24), y=Max_acc_pM*1.25, label="max acc (pM)", col= "red") +
      geom_line(aes(y=BASE_CONC*1000),lty=2,size=1, col ="red") +
      annotate("text", x=(SIM_DURATION*0.9/24), y=BASE_CONC*1300, label="Baseline (pM)", col= "red") +
      geom_line(aes(y=BASE_CONC*100),lty=2,size=1, col ="dark green") +
      annotate("text", x=(SIM_DURATION*0.9/24), y=BASE_CONC*80, label="90% supression", col= "dark green") +
      labs(x="Time (Days)", y="Free ligand or complex Concentration in Plasma (pM)",
           title ="Predicted Free Ligand or complex profile in Plasma (pM)") + theme_bw()+
      scale_y_continuous(trans = log_trans(), breaks = base_breaks_log(),limits=c(BASE_CONC*1300/100,NA),
                         labels = prettyNum)
    
    
  
    #print the final plot
    print(myplot_step1)
    
    
  })
  
})

