## EXPERIMENTS (new tidier source!)


## Initialise variables ------------------

log_n_particles_vals <- 4:8
n_particles_vals <- 2 ^ (log_n_particles_vals)

n_leaves <- n_particles_vals[1]

n_reps <- 100

## Generate synthetic data ---------------

delta <- 0.1
sigma <- 0.1
n_obs <- n_particles_vals[length(n_particles_vals)]

ou_data <- ou_sim(n_obs, delta, sigma)


## Generate immortal trajectory ----------

std_n_particles <- 100

std_smc <- standard_SMC(std_n_particles, ou_data, ou_emission_density, ou_transition_sam, ou_initial_sam, sigma = sigma, delta = delta)
imtl_index <- sample(1:std_n_particles, 1, prob = std_smc$weights[nrow(stdSMC$weights), ])

std_anc <- std_smc$ancestry
class(std_anc) <- 'genealogy'
std_anc@N <- as.integer(std_n_particles)
std_anc@Ngen <- as.integer(n_obs - 1)

imtl_anc <- traceAncestry(std_anc, sampl = imtl_index)

imtl_states <- numeric(n_obs)
for (t in 1:(n_obs - 1)) {
  imtl_states[t] <- std_smc$positions[t, imtl_anc[t]]
}
imtl_states[n_obs] <- std_smc$positions[n_obs, imtl_index]


## Run simulations -----------------------

treeht_iters <- function(data, j, n_particles, imtl_states){
  smc_sam <- conditional_SMC(n_particles, data, imtl_states, ou_emission_density, ou_transition_sam, ou_initial_sam, sigma = sigma, delta = delta)

  anc <- smc_sam$ancestry
  class(anc) <- 'genealogy'
  anc@N <- as.integer(n_particles)
  anc@Ngen <- as.integer(n_obs - 1)

  ancestrySize(anc, sampl=sample(1:n_particles, n_leaves, replace=F))$treeHeight
}

no_cores <- future::availableCores()
registerDoParallel(makeCluster(no_cores, type='FORK', outfile="debug_file.txt"))

for (i in 1:length(n_particles_vals)){
  tree_height <- foreach(j = 1:n_reps, .combine = c)  %dorng% SMC_treeht_reps(ou_data, j, n_particles_vals[i], imtl_states)
  write.table(t(tree_height), file="treeht_out.csv", sep=",", append=TRUE, row.names=FALSE, col.names=FALSE)
}

stopImplicitCluster()


## Process results -----------------------

tree_height <- read.csv("treeht_out.csv", header=FALSE)

mean_tree_ht <- apply(tree_height, 1, mean)
var_tree_ht <- apply(tree_height, 1, var)
se_tree_ht <- (var_tree_ht / n_reps) ^ 0.5

plot(n_particles_vals, mean_tree_ht / n_particles_vals, type = 'b', pch = 16, col = 2, lwd = 2, xlab = "number of particles", ylab = "mean tree height / no. particles", main = paste("Tree height profile: conditional SMC, no. leaves=", n_leaves))
lines(n_particles_vals, (mean_tree_ht - se_tree_ht) / n_particles_vals, col = 2, lty = 2)
lines(n_particles_vals, (mean_tree_ht + se_tree_ht) / n_particles_vals, col = 2, lty = 2)
legend("topright", c("mean", "+/- 1 std error"), lty = c(1, 2), lwd = c(2,1), col = 2, pch = c(16, NA))





