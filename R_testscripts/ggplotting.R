## ggplot2 plots :)

library(ggplot2)
par(mfrow = c(1, 1))

## WITHOUT GGPLOT2 -------------------
# par(mfrow = c(1, 1))
#
# tree_height0 <- read.csv("treeht_out_0sd.csv", header=FALSE)
# tree_height1 <- read.csv("treeht_out_1sd.csv", header=FALSE)
# tree_height2 <- read.csv("treeht_out_2sd.csv", header=FALSE)
# tree_height3 <- read.csv("treeht_out_3sd.csv", header=FALSE)
#
# mean_tree_ht0 <- apply(tree_height0, 1, mean)
# var_tree_ht0 <- apply(tree_height0, 1, var)
# se_tree_ht0 <- (var_tree_ht0 / n_reps) ^ 0.5
#
# mean_tree_ht1 <- apply(tree_height1, 1, mean)
# var_tree_ht1 <- apply(tree_height1, 1, var)
# se_tree_ht1 <- (var_tree_ht1 / n_reps) ^ 0.5
#
# mean_tree_ht2 <- apply(tree_height2, 1, mean)
# var_tree_ht2 <- apply(tree_height2, 1, var)
# se_tree_ht2 <- (var_tree_ht2 / n_reps) ^ 0.5
#
# mean_tree_ht3 <- apply(tree_height3, 1, mean)
# var_tree_ht3 <- apply(tree_height3, 1, var)
# se_tree_ht3 <- (var_tree_ht3 / n_reps) ^ 0.5
#
# plot(n_particles_vals, mean_tree_ht0 / n_particles_vals, ylim=c(0.1, 2.3), type = 'b', pch = 16, col = 'red', lwd = 2, xlab = "number of particles", ylab = "mean tree height / no. particles", main = paste("Conditional SMC, no. leaves=", n_leaves))
# lines(n_particles_vals, (mean_tree_ht0 - se_tree_ht0) / n_particles_vals, col = 'red', lty = 2)
# lines(n_particles_vals, (mean_tree_ht0 + se_tree_ht0) / n_particles_vals, col = 'red', lty = 2)
#
# lines(n_particles_vals, mean_tree_ht1 / n_particles_vals, type = 'b', pch = 16, col = 'purple', lwd = 2)
# lines(n_particles_vals, (mean_tree_ht1 - se_tree_ht1) / n_particles_vals, col = 'purple', lty = 2)
# lines(n_particles_vals, (mean_tree_ht1 + se_tree_ht1) / n_particles_vals, col = 'purple', lty = 2)
#
# lines(n_particles_vals, mean_tree_ht2 / n_particles_vals, type = 'b', pch = 16, col = 'blue', lwd = 2)
# lines(n_particles_vals, (mean_tree_ht2 - se_tree_ht2) / n_particles_vals, col = 'blue', lty = 2)
# lines(n_particles_vals, (mean_tree_ht2 + se_tree_ht2) / n_particles_vals, col = 'blue', lty = 2)
#
# lines(n_particles_vals, mean_tree_ht3 / n_particles_vals, type = 'b', pch = 16, col = 'green4', lwd = 2)
# lines(n_particles_vals, (mean_tree_ht3 - se_tree_ht3) / n_particles_vals, col = 'green4', lty = 2)
# lines(n_particles_vals, (mean_tree_ht3 + se_tree_ht3) / n_particles_vals, col = 'green4', lty = 2)
#
# legend("topright", c("MAP", "MAP+1SD", "MAP+2SD", "MAP+3SD"), lty = 1, lwd = 2, col = c('red', 'purple', 'blue', 'green4'), pch = 16)


## USE THIS ONE --------------------

n_reps <- 100

tree_height0 <- read.csv("Sim_data/n_16/treeht_out_0sd.csv", header=FALSE)
tree_height1 <- read.csv("Sim_data/n_16/treeht_out_1sd.csv", header=FALSE)
tree_height2 <- read.csv("Sim_data/n_16/treeht_out_2sd.csv", header=FALSE)
tree_height3 <- read.csv("Sim_data/n_16/treeht_out_3sd.csv", header=FALSE)

mean_tree_ht0 <- apply(tree_height0, 1, mean, na.rm=T)
var_tree_ht0 <- apply(tree_height0, 1, var, na.rm=T)
se_tree_ht0 <- (var_tree_ht0 / n_reps) ^ 0.5

mean_tree_ht1 <- apply(tree_height1, 1, mean, na.rm=T)
var_tree_ht1 <- apply(tree_height1, 1, var, na.rm=T)
se_tree_ht1 <- (var_tree_ht1 / n_reps) ^ 0.5

mean_tree_ht2 <- apply(tree_height2, 1, mean, na.rm=T)
var_tree_ht2 <- apply(tree_height2, 1, var, na.rm=T)
se_tree_ht2 <- (var_tree_ht2 / n_reps) ^ 0.5

mean_tree_ht3 <- apply(tree_height3, 1, mean, na.rm=T)
var_tree_ht3 <- apply(tree_height3, 1, var, na.rm=T)
se_tree_ht3 <- (var_tree_ht3 / n_reps) ^ 0.5

df <- data.frame(npart = rep(seq(from = 256, to = 4096, by = 256), 4),
                 meanht = c(mean_tree_ht0, mean_tree_ht1, mean_tree_ht2, mean_tree_ht3),
                 seht = c(se_tree_ht0, se_tree_ht1, se_tree_ht2, se_tree_ht3),
                 immortal = rep(c("MAP", "MAP+1SD", "MAP+2SD", "MAP+3SD"), each=length(mean_tree_ht0)))

ggplot(df, aes(npart, meanht / npart)) +
  geom_point(aes(col = immortal)) +
  geom_ribbon(aes(ymin = (meanht - seht) / npart,
                  ymax = (meanht + seht) / npart, fill = immortal),
                  alpha=0.3) +
  xlab("number of particles") +
  ylab("mean tree height / n. particles") +
  ggtitle("Tree height profile: conditional SMC, n. leaves=16")

# # ------------------
# # using median and quatiles instead of mean +/- SE
#
# tree_height0 <- read.csv("Sim_data/n_16/treeht_out_0sd.csv", header=FALSE)
# tree_height1 <- read.csv("Sim_data/n_16/treeht_out_1sd.csv", header=FALSE)
# tree_height2 <- read.csv("Sim_data/n_16/treeht_out_2sd.csv", header=FALSE)
# tree_height3 <- read.csv("Sim_data/n_16/treeht_out_3sd.csv", header=FALSE)
#
# med_tree_ht0 <- apply(tree_height0, 1, median)
# upper_tree_ht0 <- apply(tree_height0, 1, quantile, probs=0.9)
# lower_tree_ht0 <- apply(tree_height0, 1, quantile, probs=0.1)
#
# med_tree_ht1 <- apply(tree_height1, 1, median)
# upper_tree_ht1 <- apply(tree_height1, 1, quantile, probs=0.9)
# lower_tree_ht1 <- apply(tree_height1, 1, quantile, probs=0.1)
#
# med_tree_ht2 <- apply(tree_height2, 1, median)
# upper_tree_ht2 <- apply(tree_height2, 1, quantile, probs=0.9)
# lower_tree_ht2 <- apply(tree_height2, 1, quantile, probs=0.1)
#
# med_tree_ht3 <- apply(tree_height3, 1, median)
# upper_tree_ht3 <- apply(tree_height3, 1, quantile, probs=0.9)
# lower_tree_ht3 <- apply(tree_height3, 1, quantile, probs=0.1)
#
# df <- data.frame(npart = rep(seq(from = 256, to = 4096, by = 256), 4),
#                  median = c(med_tree_ht0, med_tree_ht1, med_tree_ht2, med_tree_ht3),
#                  upper = c(upper_tree_ht0, upper_tree_ht1, upper_tree_ht2, upper_tree_ht3),
#                  lower = c(lower_tree_ht0, lower_tree_ht1, lower_tree_ht2, lower_tree_ht3),
#                  immortal = rep(c("MAP", "MAP+1SD", "MAP+2SD", "MAP+3SD"), each=length(med_tree_ht0)))
#
# ggplot(df, aes(npart, median / npart)) +
#   geom_point(aes(col = immortal)) +
#   geom_ribbon(aes(ymin = lower / npart,
#                   ymax = upper / npart, fill = immortal),
#               alpha=0.3) +
#   xlab("number of particles N") +
#   ylab("median tree height / N") +
#   ggtitle("Tree height profile: conditional SMC, n=16")
