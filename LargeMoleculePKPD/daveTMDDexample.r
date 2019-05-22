## Full TMDD model in central compartment
# Author: Dave Fairman
# Date: 07-March-2017
# Version 1
#rm(list=ls(all=TRUE))
#FOLDER <- "DMDG course"                     #Alter to specify your own location
#base.dir <- "."                #Alter to specify your own location
#setwd(file.path(base.dir, FOLDER))
SIM_TITLE <- "mAB multiple dose" # give each simualtion a unique title
Version <- "1" #Simulation version number
SIM_DAT <- Sys.Date()

####################### Parameters to adjust for simulations
#### mAb binding kinetics
TVKD <- 0.78/10 # nM 0.78nM as default
####Target dynamics
BASE_CONC <- 0.011/1000  ## nM. 11pM as default value
HL_Tar <- 0.5 # hours Target turnover half life 0.5 based on Biesma et al
Comp_factor <- 30 # fold higher complex elimination compared to mAb.
#### Dose
Dose <- 150 # mg. total dose of mAb
YTE_factor <- 3   # set to 1 for non YTE sims. set to 3 for 3x decrease in CL
TAU <- 24*7*4*3   #Dosing interval (hours)

##############  Do not alter anyting below this line!
#### 
Kon  <- 82   # nM-1 s-1 - fixed at estimate for IMA-026 and Koff derived. Keep fixed and just vary KD as this is OK for this level
# of investigation given that the model is essntially not be sensitive to Koff
set.seed(1234567)
SIM_DURATION <- 24*7*4*6 #Duration of simulation in hours
DT <- 1  # time step for simulations (hours)
## If you wish to give multiple doses

FIRST_DOSE <- 0 #time of first dose (hours)
LAST_DOSE <- FIRST_DOSE +24*7*4*8   #time of last dose (hours)
Time_labels <- seq(FIRST_DOSE,SIM_DURATION,by=TAU)/24
#
n.ind <- 30 # number of individuals in simulation ~ >500 for final simulations but wil take time so use ~30 to check model is working first
# library
library(deSolve)   ## DE solver
library(ggplot2)   ##  for plotting
library(MASS)     ## for mvrnorm function
library(plyr)    ## for the ldply and ddply aggregation functions
## SC
## Define dosing events longitudinally (like in nonmem)
Mw <- 146 #kDa - molecular weight

dose.events <- data.frame(var="depot",
                          time=seq(FIRST_DOSE,LAST_DOSE,by=TAU), #If you wish to give multiple doses
                          #time=FIRST_DOSE,                        #single dose - just gives dose at time specified
                          value=Dose*1000/Mw, method="add") # dose express in nMol

## ## parameter values
## PK assumed to be equivalent to IMA-026
TVCL <- 0.0068/YTE_factor       # L/hr 
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

params <- mvrnorm(n.ind, mu = log(c(TVCL,TVVC,TVVP,TVQ,TVKA,TVKD)),
                  Sigma = matrix(c(0.04, 0,    0,    0,    0,    0,
                                   0,    0.01, 0,    0,    0,    0,
                                   0,    0,    0.01, 0,    0,    0,
                                   0,    0,    0,    0.01, 0,    0,
                                   0,    0,    0,    0,    0.01, 0,
                                   0,    0,    0,    0,    0,    0.25), nrow=6))

params <- exp(params)

## Initial values of state variables
yini <- c(depot = 0, plasma = 0, peripheral = 0,ligandP = BASE_P, complexP = 0)
## observation times
times <- seq(0,SIM_DURATION,DT) # sequence starting at 0 ending at defined duration going up by x hours

## loop through matrix row by row (i.e. by ID) output data.frame as you go
out <- ldply(1:nrow(params),function(i){
  
  outi <- ode(func = pkeq, times = times,
              y = yini, parms = params[i,],events=list(data=dose.events))
  
  outi <- as.data.frame(outi)
  outi$plasma <- outi$plasma/params[i, 2]  # divides amount by individual V1 to get conc
  outi$ligandP <- outi$ligandP/params[i, 2]  # divides amount by individual V1 to get conc
  outi$SIM <- i
  outi
})

## plot the results
## mAb profile (Plasma)
mAbpl <- ddply(out,"time",function(d){
  data.frame(lowpl = quantile(d$plasma,probs=0.025),
             midpl = quantile(d$plasma,probs=0.5),
             hipl = quantile(d$plasma,probs=0.975))
})

## Plot mAb plasma levels
##get key infor into pdf title and ensure that numbers are no too long
k_Tar <- signif(Kdeg,3)
k_mAb <- signif(TVCL/TVVC,3)
k_comp <- signif(Kcomplex,3)

pdf(file=paste(Dose,"mg.",TVKD,"nM Kd.","k_Tar=",k_Tar,".k_mAb=",k_mAb,".Kcomplex=",k_comp,".",SIM_TITLE,
               ".","PK, target and complex plots with PI using TMDD model",".n=",n.ind,".",SIM_DAT,".pdf", sep=""))

ggplot(mAbpl,aes(x=time/24,y=midpl)) +
  geom_ribbon(aes(ymin=lowpl,ymax=hipl),alpha=0.4,fill = "red") +
  #scale_x_continuous(breaks=c(0,7,14,21,28,56,28*3,28*4,28*5,28*6), limits=c(0,SIM_DURATION/24)) + #Option to put in specific points
  scale_x_continuous(breaks=Time_labels, limits=c(0,SIM_DURATION/24)) +
  geom_line() + 
  scale_y_log10(breaks=c(10,30,100,300,600), limits=c(6,600)) +
  labs(x="Time (Days)", y="mAb Concentration (nM)", 
       title ="Predicted mAb profile in Plasma") + theme_bw()

## Target and complex profile (Plasma)
pl <- ddply(out,"time",function(d){
  data.frame(low = quantile(d$ligandP,probs=0.025),
             mid = quantile(d$ligandP,probs=0.5),
             hi = quantile(d$ligandP,probs=0.975))
})

cmxpl <- ddply(out,"time",function(d){
  data.frame(lowcmx = quantile(d$complexP,probs=0.025), #No need to divide by VC as complex is already as conc
             midcmx = quantile(d$complexP,probs=0.5),
             hicmx = quantile(d$complexP,probs=0.975))
})
Tarpl <- cbind(pl, cmxpl)

myplot <- ggplot(Tarpl,aes(x=time/24,y=mid*1000)) +
  geom_ribbon(aes(ymin=low*1000,ymax=hi*1000),alpha=0.4, fill = "red") +
  geom_ribbon(aes(ymin=lowcmx*1000,ymax=hicmx*1000),alpha=0.4,fill = "blue") +
  geom_line(aes(x=time/24,y=midcmx*1000)) +
  #scale_x_continuous(breaks=c(0,7,14,21,28,56,28*3,28*4,28*5,28*6), limits=c(0,SIM_DURATION/24)) + #Option to put in specific points
  scale_x_continuous(breaks=Time_labels, limits=c(0,SIM_DURATION/24)) +
  geom_line() + 
  scale_y_log10(breaks=c(0.0001,0.001,0.01,0.1,1,3,Max_acc_pM,10), limits=c(BASE_CONC,13)) +
  geom_line(aes(y=Max_acc_pM),lty=2,size=1, col ="red") + 
  annotate("text", x=(SIM_DURATION*0.9/24), y=Max_acc_pM*1.25, label="max acc (pM)", col= "red") +
  geom_line(aes(y=BASE_CONC*1000),lty=2,size=1, col ="red") + 
  annotate("text", x=(SIM_DURATION*0.9/24), y=BASE_CONC*1300, label="Baseline (pM)", col= "red") +
  geom_line(aes(y=BASE_CONC*100),lty=2,size=1, col ="dark green") +
  annotate("text", x=(SIM_DURATION*0.9/24), y=BASE_CONC*80, label="90% supression", col= "dark green") +
  labs(x="Time (Days)", y="Free ligand or complex Concentration in Plasma (pM)", 
       title ="Predicted Free Ligand or complex profile in Plasma") + theme_bw()

dev.off()