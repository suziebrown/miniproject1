## make standard dataset for use in all runs

ou_data <- ou_sim(4096 * 6, 0.1, 0.1)
write.table(ou_data, sep = ',', row.names = FALSE, col.names = FALSE, file = "oudata.csv")
