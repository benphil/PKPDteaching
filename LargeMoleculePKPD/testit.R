draw_params <- function (theta) {
  p <- list()  # Parameter list 
  p$CL <- theta[1]
  p$V  <- theta[2]
  p$KA <- theta[3]
  p$S1 <- p$V 
  return(p)
}

### ODE system
des <- function (t, A, p) {  # ODE system 
  p$K20 <- p$CL/p$V
  dAdt_1 <- -p$KA*A[1] 
  dAdt_2 <- -p$K20*A[2]  +p$KA*A[1]
  return ( list ( c (  dAdt_1, dAdt_2 ) ) )
}

### Perform numerical integration, collect data
num_int_wrapper<- function (i, times, A_init, des, p_ind) {
  des_out <- lsoda(A_init, times, des, p_ind)
  dat_ind <- c()
  for (j in 1:length(A_init)) {
    dat_ind <- rbind (dat_ind, cbind(i, t=des_out[,1], comp=j, ipred=des_out[,(j+1)]))
  }
  return(data.frame(dat_ind))
}

vanBuerenData<-function(){
  
  #SF_ugml2nM = 1e6/par.MW_mAbs; % ug/ml -> nM
  
  # Format: time (days), 2mg/kg, 2mg/kg, 20mg/kg, 20kg,mg, 40kg/mg, 40kg/mg
  
  LammertsVanBuerrenData = c( 0.004,  63.7,  30.6,	587.2,	449.7,	798.1,	929.5,
     0.02,	36.1,	43.1,	492.2,	631.0,	160.0,	376.2,
     0.04,	55.1,	50.4,	534.0,	446.6,	985.4,	930.6,
     0.08,	46.6,	42.1,	336.3,	313.2,	636.1,	668.2,
     0.1,     39.8,	39.2,	415.8,	393.8,	1161.1,	983.8,
     0.3,     40.6,	30.9,	414.4,	387.3,	683.7,	713.3,
     0.5,     29.3,	31.7,	549.4,	422.1,	1077.8,	628.2,
     1,       18.3,	18.1,	274.8,	247.4,	692.6,	1029.5,
     2,       10.6,	11.5,	232.2,	230.7,	451.7,	571.8,
     3,       5.1,     5.7,     234.5,	211.5,	626.0,	655.8,
     6,       1.0,     2.1,     200.6,	151.8,	409.6,	333.1,
     8,       0.4,     0.8,     180.5,	145.5,	319.3,	357.3,
     11,      0.2,     0.4,     111.528,	86.4,	261.6,	249.8,
     15,      0.1,     0.2,     79.640,	68.2,	192.4,	154.0,
     18,      0.1,     0.1,     74.940,	64.9,	165.2,	84.3,
     22,      NaN,     NaN,     38.036,	28.9,	125.1,	68.4,
     26,      NaN,     NaN,     13.485,	15.9,	84.1,	28.0,
     30,      NaN,     NaN,     5.103,	11.3,	38.4,	13.4,
     36,      NaN,     NaN,     1.266,	1.9,     13.3,	1.5)
  
  LammertsVanBuerrenMatrix = matrix(LammertsVanBuerrenData,nrow=19,ncol=7,byrow = TRUE)        # fill matrix by rows 
  
  #convert in mg/ml from ug/ml
  LammertsVanBuerrenMatrix <- LammertsVanBuerrenMatrix/1000
  LammertsVanBuerrenMatrix[,1] <-LammertsVanBuerrenMatrix[,1]*1000 #but not time
  
  
  
  return(data.frame(LammertsVanBuerrenMatrix))
}



  ### Dose and time Settings
  A_init     <- rep(0, 2)  # Initial state of ODE system
  n_doses    <- 3
  dose_cmt   <- 1
  ii         <- 24
  dose_times <- seq (from = 0, by=ii, to=n_doses*ii)
  dose_amts  <- c(rep (100, n_doses), 0)
  times      <- seq(from=0, to=ii*n_doses, by=.5)  # Integration window and stepsize
  obs_c      <- c(1:2)  # Observation compartments
  n_ind      <- 1
  
  
  ### Combine Parameters
  theta <- c(1, 1, 0.1) 
  p_ind   <- draw_params (theta=theta)
  
  
  comb_dat <- c()
  for (k in 1:(length(dose_times)-1)) {
    
    if (k > 1) {
      A_upd <- dat_ind[dat_ind$t==tail(time_window,1),]$ipred
    } else {
      A_upd <- A_init
    }
    A_upd[dose_cmt] <- A_upd[dose_cmt] + dose_amts[k]
    time_window <- times[(times > dose_times[k]) & (times <= dose_times[k+1])]
    
    dat_ind <- num_int_wrapper (1, time_window, A_upd, des, p_ind)
    
    #des_out <- lsoda(A_init, times, des, p_ind)
    #dat_ind <- data.frame(des_out)
    
    comb_dat <- rbind (comb_dat, dat_ind)
  }

xlabel <- "time in hours"
ylabel <- "concentration in mg/ml"

dataset <- subset(comb_dat,comp!=1)

expDaten = vanBuerenData()

# p<-ggplot(data=data.frame(dataset), aes(x=t, y=ipred,color=as.factor(comp),group=as.factor(comp))) + 
#         geom_line(size=2)+
#         xlab(xlabel) + 
#         ylab(ylabel)
# 
# p+geom_point(data=expDaten, aes(x=X1,y=X7,group=NULL,color=NULL))

p=ggplot(data=data.frame(dataset),aes(x=t/24, y=ipred))+
  geom_line(size=2,color="blue")+
  geom_hline(yintercept = 0.1,linetype=2,color="red") +
  #geom_line(data=c.mean,aes(x=t/24,y=ipred,size=2),,color="blue")+       
  xlab(xlabel) +
  ylab(ylabel)+
  #scale_color_manual(values=c("#f5949e","#ed4354","#55c508","#bcb7b8","#000000","#000000", "#ED4354", "#ED4354", "#72BF43"),name="Simulated\nconcentrations",breaks=c("1", "2", "3", "4", "5"),labels=c("Central compartment", "Peripheral compartment","Mean central","Mean peripheral","LLoQ"))+
  #scale_linetype_manual(values=c(1, 1, 1, 1, 2),name="Simulated\nconcentrations",breaks=c("1", "2","3","4","5"),labels=c("Central compartment", "Peripheral compartment","Mean central","Mean peripheral","LLoQ"))+
  #theme_bw()+
  theme(axis.title.y = element_text(size = rel(1.5), angle = 90, vjust=2),axis.title.x = element_text(size = rel(1.5), vjust=-0.5),legend.position = c(0.85, 0.98),legend.justification=c(1, 1))

print(p)


    
