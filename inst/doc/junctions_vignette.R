## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#    sim_inf_chrom(pop_size,
#                  initial_heterozygosity,
#                  total_runtime,
#                  morgan,
#                  markers,
#                  seed)

## ---- single_sim_inf, fig.width = 7, fig.height = 7----------------------
  N = 100 #population size
  H_0 = 0.5 #initial heterozygosity
  maxT = 1000 #run time
  C = 1 #number of recombinations per meiosis
  numMarkers = -1 #we first ignore the effect of imposing randomly distributed markers

  v <- sim_inf_chrom(N, H_0, maxT, C, numMarkers, 42)

  plot(v$avgJunctions, type = "l", xlab = "Generations",
        ylab = "Number of Junctions", main = "Example Infinite Chromosome")

## ---- repl_sim_inf-------------------------------------------------------
  number_replicates <- 100
  v <- c();
  for(r in 1:number_replicates) {
    v2 <- sim_inf_chrom(N, H_0, maxT, C, numMarkers, sample(1e9,1))  
    v <- rbind(v, as.numeric(v2$avgJunctions))
  }
  v <- colMeans(v) #mean across replicates

## ---- plot_sim_inf, fig.width = 7, fig.height = 7------------------------
  clarity <- seq(1, maxT, length.out = 50) #we plot not all points, for clarity
  plot(v[clarity]~clarity, lwd = 2,
       xlab = "Generations",
       ylab = "Number of Junctions",
       main = "Average behaviour Infinite chromosome", pch = 16) 

  t <- 0:maxT
  predicted <- number_of_junctions(N = N, H_0 = H_0, C = C, t = t)
  lines(predicted~t,col="blue") 

## ---- single_sim_random, fig.width = 7, fig.height = 7-------------------
N = 100 #population size
H_0 = 0.5 #initial heterozygosity
maxT = 1000 #run time
C = 1 #number of recombinations per meiosis
numMarkers = 1000 #1000 markers

#single example run
v <- sim_inf_chrom(N, H_0, maxT, C, numMarkers, 42)
plot(v$avgJunctions, type = "l", xlab = "Generations",
     ylab = "Number of Junctions", main = "Example Infinite Chromosome")

lines(v$detectedJunctions, col = "blue")
legend("bottomright", c("Real number","Number detected"),
       lty = 1, col = c("black", "blue"))

## ---- repl_sim_random----------------------------------------------------
meanJ <- c()
detectJ <- c()
for(r in 1:number_replicates) {
  v2 <- sim_inf_chrom(N, H_0, maxT, C, numMarkers, sample(1e9,1))
  meanJ <- rbind(meanJ, as.numeric(v2$avgJunctions))
  detectJ <- rbind(detectJ, as.numeric(v2$detectedJunctions))
}
meanJ <- colMeans(meanJ)
detectJ <- colMeans(detectJ)

## ---- plot_sim_random, fig.width = 7, fig.height = 7---------------------
clarity <- seq(1,maxT,length.out=50) #we plot not all points, for clarity
plot(meanJ[clarity]~clarity, xlab = "Generations",
     ylab = "Number of Junctions",
     main = "Average behaviour Infinite chromosome",
     pch=16) 

points(detectJ[clarity]~clarity, pch = 17, col = "blue")

t <- 0:maxT
predicted <- number_of_junctions(N = N, H_0 = H_0, C = C, t = t)
lines(predicted~t,lwd=2) 

K <- max(detectJ) #now substitute K with the observed maximum detected.
pred <- K - K *(1-H_0*C/K)^t
lines(pred~t,col="blue",lwd=2) 

legend("bottomright",
       c("Real number","Detected",
         "Predicted Real","Predicted Detected"),
       pch=c(16,17,NA,NA),
       lty=c(NA,NA,1,1),
       col=c("black","blue",
             "black","blue"),
       lwd=2)

## ---- eval = FALSE-------------------------------------------------------
#  sim_fin_chrom(pop_size,
#                initial_heterozygosity,
#                total_runtime,
#                morgan,
#                seed,
#                R)

## ---- single_sim_fin, fig.width = 7, fig.height = 7----------------------
  R = 100 #chromosome size
  N = 100 #population size 
  H_0 = 0.5 #initial heterozygosity
  C <- 1 #number of recombinations per meiosis
  maxT <- 1000

  #single example run
  v <- sim_fin_chrom(pop_size = N, 
                     initial_heterozygosity = H_0, 
                     total_runtime = maxT, 
                     morgan = C, 
                     seed = 42,
                     R = R)
  plot(v$avgJunctions, type="l", xlab="Generations",
     ylab="Number of Junctions", main="Example Finite Chromosome")

## ---- repl_sim_fin-------------------------------------------------------
v <- c();
for(r in 1:number_replicates) {
  v2 <- sim_fin_chrom(N,H_0, maxT, C, sample(1e9,1), R)
  v <- rbind(v, as.numeric( v2$avgJunctions))
}
v <- colMeans(v)

## ---- plot_sim_fin, fig.width = 7, fig.height = 7------------------------
clarity <- seq(1, 1000, length.out=50) #we plot not all points, for clarity
plot(v[clarity]~clarity, lwd = 2, 
     xlab = "Generations", 
     ylab = "Number of Junctions", 
     main = "Average behaviour Finite Chromosome", 
     pch  = 16) 

t <- 0:maxT
predicted <- number_of_junctions(N = N, R = R, 
                                 H_0 = H_0, C = C, 
                                 t)
lines(predicted~t,col="blue") 
legend("bottomright", c("Simulated","Predicted"),
       pch = c(16, NA),
       lty = c(NA, 1),
       col = c("black", "blue"))

## ---- equations----------------------------------------------------------
J <- number_of_junctions(N = 100, R = 1000, H_0 = 0.5, C = 1, t = 200)
J <- tail(J,1)
J

## ------------------------------------------------------------------------
  time_estim <- estimate_time(J = J, N = 100, R = 1000, H_0 = 0.5, C = 1)
  time_estim 

## ------------------------------------------------------------------------
time_error(J = J, N = 100, R = 1000, H_0 = 0.5, C = 1, time_estim, relative = TRUE)

## ------------------------------------------------------------------------
time_error(J = J, N = 100, R = 1000, H_0 = 0.5, C = 1, time_estim, relative = FALSE)

## ------------------------------------------------------------------------
calculate_mat(N = 100, R = 1000, H_0 = 0.5, C = 1)

## ---- fig.width = 7, fig.height = 7--------------------------------------
maxT <- 1000
t <- 0:maxT
J <- number_of_junctions(N = 100, 
                         R = 1000,
                         H_0 = 0.5,
                         C = 1,
                         t = 0:maxT)
par(mar=c(4,5,2,3))
plot(J~t, 
     type = "l",
     xlab = "Generations",
     ylab = "Number of Junctions",
     xlim = c(0, maxT))
#vertical line that indicates the upper limit
abline(v = calculate_mat(N = 100, R = 1000, H_0 = 0.5, C = 1), lty = 2)
par(new = T)

v <- time_error(J = J, 
                     N = 100, 
                     R = 1000, 
                     H_0 = 0.5, 
                     C = 1, 
                     0:(maxT - 1), #to avoid an error at t = maxT
                     relative = TRUE)
plot(v, 
     col = "red", type = "l",
     xlim = c(0, maxT), 
     ylim = c(0, 1), 
     xaxt = "n", yaxt = "n",
     xlab = "", ylab = "")
axis(4)
mtext("Relative error",side=4, line = 2)
legend("bottomright", 
       c("Number of junctions", "Relative Error", "t_MAX"),
       lty = c(1, 1, 2),
       col = c("black", "red", "black"))

