## Make ggplot2 plots seen in report

library(ggplot2)
par(mfrow = c(1, 1))


## Conditional SMC (for each plot) --------------------

n_reps <- 1000
n_obs <- 6 * 4096

tree_height0 <- read.csv("treeht_cond_n2_0sd.csv", header=FALSE)
tree_height1 <- read.csv("treeht_cond_n2_1sd.csv", header=FALSE)
tree_height2 <- read.csv("treeht_cond_n2_2sd.csv", header=FALSE)
tree_height3 <- read.csv("treeht_cond_n2_3sd.csv", header=FALSE)

tree_height0[is.na(tree_height0)] <- n_obs
tree_height1[is.na(tree_height1)] <- n_obs
tree_height2[is.na(tree_height2)] <- n_obs
tree_height3[is.na(tree_height3)] <- n_obs

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
  ggtitle("Tree height profile: conditional SMC, n. leaves=2")


## Standard SMC ----------------------------------

treeht_std2 <- read.csv("Sim_data/180604/treeht_std_n2.csv", header=FALSE)
treeht_std16 <- read.csv("Sim_data/180604/treeht_std_n16.csv", header=FALSE)

mean_std2 <- apply(treeht_std2, 1, mean, na.rm=T)
var_std2 <- apply(treeht_std2, 1, var, na.rm=T)
se_std2 <- (var_std2 / n_reps) ^ 0.5

mean_std16 <- apply(treeht_std16, 1, mean, na.rm=T)
var_std16 <- apply(treeht_std16, 1, var, na.rm=T)
se_std16 <- (var_std16 / n_reps) ^ 0.5

df <- data.frame(npart = rep(seq(from = 256, to = 4096, by = 256), 2),
                 meanht = c(mean_std2, mean_std16),
                 seht = c(se_std2, se_std16),
                 n.leaves = as.factor(rep(c(2, 16), each=length(mean_std2))))

# df <- data.frame(npart = seq(from = 256, to = 4096, by = 256),
#                  meanht = mean_std2,
#                  seht = se_std2,
#                  n.leaves = rep(2, length(mean_std2)))

ggplot(df, aes(npart, meanht / npart, group=n.leaves)) +
  geom_line(aes(col = n.leaves)) +
  geom_ribbon(aes(ymin = (meanht - seht) / npart,
                  ymax = (meanht + seht) / npart, fill = n.leaves),
              alpha=0.3) +
  xlab("number of particles") +
  ylab("mean tree height / n. particles") +
  ggtitle("Tree height profile: standard SMC")
