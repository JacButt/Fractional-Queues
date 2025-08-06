######################################################
# This code is used to generate plots of our functional limits for the GI/Gi/1 queue model, using the formulation in terms
# of inverse subordinators. We assume that alpha_1 = alpha_2 = alpha for our two fractional Poisson Processes in this case, with arrival "rates" 
# lambda and departures "rates" mu

set.seed(1) # Used for debugging purposes

t_step = 0.01 #time partition size
x_step = 0.1 # space partition size
n_steps = 10^5 #number of steps in the algorithm

#####################################

#Initialising vectors for the Subordinators and inverse-subordinators of the arrival and departure processes
#To be filled later in the code

Sub_Arrival <- rep(0, n_steps) 
Sub_Departure <- rep(0, n_steps)
Inv_Arrival <- rep(0, n_steps)
Inv_Departure <- rep(0, n_steps)

# Initialising vectors for spatial and time coordinates

T <- seq(from = 0, to = n_steps * t_step, by = t_step) 
X <- seq(from = 0, to = n_steps * x_step, by = x_step) 


alpha = 0.7 #shape parameter
lambda = 4
mu = 1

#######################################
# algorithm to generate subordinators XA and XD.
########################################

k = 0 # dummy variable for the loop
xa = 0 # initial values for the subordinators at time 0
xd = 0

# While loop to fill in values of the subordinators. Done by following the CMS formula


while (k < n_steps){
  k = k + 1
  uA <- runif(1)
  vA <- runif(1)
  uD <- runif(1)
  vD <- runif(1)
  phiA = pi * (vA - 1/2)
  phiD = pi * (vD - 1/2)
  L =  t_step^(1/alpha) * (1 + (tan((pi * alpha)/2))^2)^(1/(2 * alpha))
  xStepA = L * (sin(alpha * pi * vA)/(cos(phiA))^(1/alpha)) * (cos(phiA - alpha * pi * vA)/(-log(uA)))^((1-alpha)/alpha)
  xa = xa + xStepA
  xStepD = L * (sin(alpha * pi * vD)/(cos(phiD))^(1/alpha)) * (cos(phiD - alpha * pi * vD)/(-log(uD)))^((1-alpha)/alpha)
  xd = xd + xStepD
  Sub_Arrival[k] = xa
  Sub_Departure[k] = xd
}

################################################################
# two algorithms to generate the inverse subordinators YA and YD from XA, XD respectively. These are generated directly using the definition of inverse functions
#############################################################

###
# Initialising a while loop to fill in the values of Inv_Arrivals, by inverting space-time and using Sub_Arrivals.
# For each spatial coordinate j in X we calculate the minimum point in T for which Sub_Arrivals > j and 
# set that value as Inv_Arrivals[j]

j = 1
k = 1
lim_Arrival = n_steps
while (j < n_steps + 1){
  if (k >= n_steps){
    lim_Arrival = j
    break
  }
  else if (Sub_Arrival[k] < X[j]){
    k = k + 1
  }
  else {
    Inv_Arrival[j] = T[k]
    j = j + 1
  }
}

# Same procedure as before but using Sub_Departures to get Inv_Departures

j = 1
k = 1
lim_Departure = n_steps
while (j < n_steps + 1){
  if (k >= n_steps){
    lim_Departure = j
    break
  }
  else if (Sub_Departure[k] < X[j]){
    k = k + 1
  }
  else {
    Inv_Departure[j] = T[k]
    j = j + 1
  }
}

lim_final = min(lim_Arrival, lim_Departure) - 1 # Cut-off point for the final loop


#######################################################################
# Final section that fills in values of the functional limit against time using the values of lambda, mu and the inverse subordinators Inv_Arrival, Inv_Departure
#######################################################################

Q <- rep(0, lim_final) # vector initialisation for the functional limit given in the final line of (3.20)

k = 1 # dummy variable for while loop
R = 0 # Value of the functional limit before reflection
inf_R = 0 # Variable tracking the infimum of R so we can reflect it across the axis
R_diff = 0 # The reflection of R across the axis, giving the actual value for the functional limit

####
# This while loop fills the vector Q. Throughout each pass of the loop, the value of the infimum of R
# up to that point is stored, so that the functional limit can be calculated using R - R_Diff

while (k < lim_final + 1){
  R = (lambda^(alpha)) * Inv_Arrival[k] - (mu^(alpha)) * Inv_Departure[k]
  if(R < inf_R){
    inf_R = R
  }
  R_diff = R - inf_R
  Q[k] = R_diff
  k = k + 1
}

# Finally, we plot the values of the functional limit against time. Note that due to the inversion of the subordinator
# before, we have to use the spatial coordinates from before as the time coordinates here 

plot(X[0:lim_final], unlist(Q), type="l", xlab = "time", ylab = "X") 

