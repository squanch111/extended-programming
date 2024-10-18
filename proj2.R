# Contribution:
# Team Members:
# - Member1 (S2702040)
# - Member2 (S2663848)
# - Member3 (S2678971)

# - [Member1 and Member2]: Developed the infection modeling code and bootstrapping functions, 
#   as well as performed data processing (70%).
# - [Member3]: Worked on the plotting functions and added uncertainty analysis with bootstrapping,as well as comments(30%).

# Overview:
# This code performs simulation-based modeling to estimate infection trajectories from observed Covid-19 death data.
# Using a lognormal distribution to model the delay between infection and death, the code optimizes the estimated 
# infection times for each case, and simulates the number of infections per day.
# Additionally, bootstrapping is used to quantify uncertainty in the infection estimates 
# once the method has converged. The method is using Poisson-simulated data to replace the real death data multiple times.
# Finally this code will show three plots: 
## 1. Real Death Trajectory VS Simulated Death Trajectory
## 2. [After Bootstrapping] Poisson-simulated Death Trajectory VS Simulated Death Trajectory  # This plot can show the uncertainty
## 3. The Simulated Incidence Trajectory(Before and After Bootstrapping), Actual Death Trajectory 
##    and Lockdown Start Day are plotted for visual comparison
# The code also tracks the execution time of the simulation.


# 1. Load and process the data
# Input: "engcov.txt" - data file containing time and number of deaths per day
# Output: 'data' - The first 150 rows of processed data
# Purpose: Load raw data into memory and prepare it for analysis
# engcov <- read.csv("~/Desktop/engcov.txt", sep="")
data <- engcov
data <- data[1:150, ] 

# 2. Set the model parameters
# Input: manually set parameters (mean, standard deviation, lockdown date, etc.)
# Output: 'duration_probs' - Normalized lognormal distribution used to simulate the time from infection to death
# Purpose: Set the parameters of the lognormal distribution and generate the probability distribution of the time from infection to death
meanlog <- 3.152
sdlog <- 0.451
lockdown_day <- 84  # The date of UK lockdown
duration_probs <- dlnorm(1:80, meanlog, sdlog)
duration_probs <- duration_probs / sum(duration_probs)  # Standardization

# 3. Initialize time from infection to death and time of death
# Input: 'data$julian' - date of death per day, 'data$nhs' - number of deaths per day
# Output: 'death' - expanded death date vector; 't0' -initial time of infection for each death case
# Purpose: Calculate the initial infection time for each case and prepare the input for simulation optimization
death <- rep(data$julian, data$nhs)  # Expand the number of deaths into individuals
n <- 29442 # Total number of deaths
infection_to_death <- sample(1:80, n, replace=TRUE, prob=duration_probs)
t0 <- death - infection_to_death

# 4. Simulate and optimize process core functions
# Function: deconv
# Inputs:
# - t: A vector containing the dates of death.
# - deaths: A vector containing the number of deaths on each day.
# - n.rep: The number of iterations to perform for optimization or bootstrapping (default = 100).
# - bs: A logical flag indicating whether to perform bootstrapping (default = FALSE).
# - t0: (Optional) A vector of initial guesses for infection times. If NULL, the function will initialize these internally.
# Outputs:
# - A list containing:
#   - P: A vector of the error (P-value) after each iteration.
#   - inft: A matrix of daily infection counts at each iteration.
#   - t0: The optimized infection times for each individual case.
# Purpose:
# This function optimizes the infection times by minimizing the error between the simulated and actual daily death counts.
# The error is calculated based on a modified Pearson statistic. Bootstrapping can be used to estimate uncertainty in 
# the infection counts.
deconv <- function(t, deaths, n.rep=100, bs=FALSE, t0=NULL) {
  colors <- rainbow(n.rep) # Initialize colors for drawing
  P_history <- numeric(n.rep)  # Save the history for each iteration
  inft <- matrix(0, nrow=310, ncol=n.rep)   # Save a history of the number of daily infections
  #If do not use Bootstrap, draw the actual death curve
  if (bs == FALSE){
    plot(data$julian, data$nhs, type='h', col="blue", lwd=2, xlab="Day", ylab="Deaths",
         main="Real and Simulated Deaths Over Iterations") # plot the real deaths
    legend("topright", legend=c("Actual Deaths"), col=c("blue"), cex=1.2, lty=1, lwd=2)
  }
  # If no initial infection time t0 is passed in, t0 is initialized with the given t and deaths
  if (is.null(t0)) {
    t0 <- death - infection_to_death # Initial guess time of infection
  }
  # Calculate actual death count
  di <- tabulate(death, nbins=310) 
  
  # Perform 100 iterations
  for (rep in 1:n.rep) {
    # Randomly shuffle the order of t0
    indices <- sample(seq_along(t0))  # Randomize the index order to generate a random out-of-order sequence number vector of length t0
    
    # Regenerate a new infection-to-death time vector and calculate simulated time of death in simulated_deaths
    infection_to_death <- sample(1:80, n, replace=TRUE, prob=duration_probs)
    # Keep t0 within reasonable limits
    t0 <- pmin(pmax(t0, 1), 310) 
    simulated_deaths <- t0 + infection_to_death

    # Calculate the daily simulated death count
    dsi <- tabulate(simulated_deaths, nbins=310)
    
    # Add bootstrap
    # use Poisson distribution
    if (bs){
      # Generate Poisson-distributed random numbers for 310 days
      di_bs <- rpois(length(di), lambda=di) 
      current_sum <- sum(di_bs)
      # Adjust the Poisson data to match actual death counts
      if (sum(di_bs) < n) {
        # Randomly choose days to increase the death counts
        adjustment_di <- sample(1:310, (n - sum(di_bs)), replace=TRUE)  
        adjustment_table <- table(adjustment_di)
        # Add the counts from adjustment_table to the corresponding days in di
        di_bs[as.numeric(names(adjustment_table))] <- di_bs[as.numeric(names(adjustment_table))] + as.numeric(adjustment_table)
      } 
      else if (sum(di_bs) > n) {
        # Randomly choose days to decrease the death counts
        # only from days where deaths occurred
        adjustment_di <- sample(which(di_bs > 0), (sum(di_bs) - n), replace=TRUE)  
        adjustment_table <- table(adjustment_di)
        # Decrease in di, ensuring no day is decreased by more than 1
        di_bs[as.numeric(names(adjustment_table))] <- di_bs[as.numeric(names(adjustment_table))] - pmin(adjustment_table, 1)
        }

      if(rep == 1){
        plot(1:310, di, type='h', col=rgb(0.4,0.1,0,0.1), lwd=2, xlab="Day", ylab="Deaths",
             main="Real and Simulated Deaths Over Iterations(Bootstrap)") # plot the death generated by Poisson distribution
        legend("topright", legend=c("Actual Deaths(Poisson)"), col=rgb(0.4,0.1,0,0.1), lty=1, lwd=2)
      }
      p <- sum((di_bs - dsi)^2 / pmax(1, dsi),na.rm=TRUE) # Avoid multiple computing di
    }
    else{   
      # Calculate the current error p-value
      p <- sum((di - dsi)^2 / pmax(1, dsi),na.rm=TRUE) 
    }

    # Adjust the step size according to the number of iterations
    if (rep > 75) {
      step_size <- sample(c(-2, -1, 1, 2), length(t0), replace=TRUE)
    } else if (rep > 50) {
      step_size <- sample(c(-4, -2, -1, 1, 2, 4), length(t0), replace=TRUE)
    } else {
      step_size <- sample(c(-8, -4, -2, -1, 1, 2, 4, 8), length(t0), replace=TRUE)
    }
    
    # Iterate over each element in t0 and apply the step size
    for (i in indices) {
      step <- step_size[i]
      former_death_day <- simulated_deaths[i]
      current_death_day <-  pmin(pmax(t0[i] + step, 1), 310) + infection_to_death[i]
      # Updated simulated death
      dsi[former_death_day] <- dsi[former_death_day] - 1
      dsi[current_death_day] <- dsi[current_death_day] + 1
      # Calculate the new p-value
      new_p <- sum((di - dsi)^2 / pmax(1, dsi),na.rm=TRUE)
      # Compare the new P value and decide whether to accept the new t0
      if (new_p < p) {
        t0[i] <- pmin(pmax(t0[i] + step, 1), 310)
        p <- new_p
      } else {
        dsi[former_death_day] <- dsi[former_death_day] + 1
        dsi[current_death_day] <- dsi[current_death_day] - 1
      }
    }
    # Save the P-value and the number of daily infections
    P_history[rep] <- p
    inft[, rep] <- tabulate(t0, nbins=310)
    # Output iteration information
    cat("Iteration:", rep, "p value:", p, "\n")
    # The simulated death data is overlaid with a different color after each iteration
    if (bs){
      lines(1:310, dsi[1:310], col=rgb(0.8, 0.2, 0.2), lwd=1)
    }
    else{
      lines(1:211, dsi[1:211], col=colors[rep], lwd=1)  # Use different colors to draw simulated death data
    }
    
  }
  
  # Return error history and daily number of infections
  return(list(P=P_history, inft=inft, t0=t0))
  
}

# 5. Perform the simulation and record the time
# Input: Pass in data for simulation and optimization
# Output: runtime, optimized results (number of infections)
# Objective: To evaluate the operational efficiency of the model and obtain the final distribution of the number of infections
start_time = Sys.time() # Record start time
system.time(result_optimized <- deconv(data$julian, data$nhs, bs=FALSE))
end_time = Sys.time()  # End of record time
total_time = end_time-start_time # Calculate total run time
cat("Execution time:", total_time, "mins\n")  # Output execution time

# Extract actual death data and final simulated death data
final_simulated_deaths <- result_optimized$inft[, ncol(result_optimized$inft)]  # Simulated final infection rate

#Poisson distribution simulation after convergence of p and t0 (substituting the convergent t0)
system.time(bootstrap_result <- deconv(data$julian, data$nhs, bs=TRUE, t0=result_optimized$t0))
final_simulated_deaths_bs <- bootstrap_result$inft[, ncol(bootstrap_result$inft)]

# plot
plot(1:310, final_simulated_deaths[1:310],type='l' ,col="red", lwd=2, xlab="Day", ylab="Counts",
     main="Actual Deaths vs Simulated Incidence")  #Simulated infection curve
lines(data$julian, data$nhs, col="blue", lwd=2 ) # True death curve
lines(1:310, final_simulated_deaths_bs[1:310],type='l' ,col="green", lwd=2)# simulation after strapping
# Mark the first day of lockdown (Day 84)
abline(v=84, col="black", lty=2, lwd=2)

# add legend
legend("topright", legend=c("Simulated Incidence", "Actual Deaths", "Lockdown Start", "Simulated Incidence After Bootstrapping"),
       col=c("red", "blue", "black","green"), lty=c(1,1,2,1),
       lwd=c(2,2,2,2),
       border=NA,
       y.intersp=1.2,
       box.lwd=0,
       cex = 0.8)



