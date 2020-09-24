library("pROC")
library("ggplot2")

dataDir <- "~/Projects/afni_python/"
data_raw <- read.delim(paste0(dataDir, "vCAT_006_task.csv"), sep = ",")

num_runs <- 4

# clean gaps, fill stim_ID
for (i in 1:dim(data_raw)[1]) {
  if (is.na(data_raw[i, 16])) {
    data_raw <- data_raw[-i, ]
  }
  if (is.na(data_raw[i, 17])) {
    data_raw[i, 18] <- "SceneFace"
  }
}

# num trials per run
num_runTrial <- dim(data_raw)[1] / num_runs

# separate into runs
c <- 1
cc <- num_runTrial
for (i in 1:num_runs) {
  assign(paste0("data_run", i), data_raw[c:cc, 16:22])
  c <- cc + 1
  cc <- cc + num_runTrial
}


# make TF for scene-face Hit, Miss
out_hit <- paste0(dataDir, "tf_hit.txt")
out_miss <- paste0(dataDir, "tf_miss.txt")

for (i in 1:num_runs) {

  # determine append
  if (i == 1) {
    h_ap <- F
  } else {
    h_ap <- T
  }

  # run data
  h_df <- get(paste0("data_run", i))

  # get start, duration of hits
  ind_sf_hit <- which(h_df$stim_ID == "SceneFace" & h_df$trial_acc == 1)
  start_hit <- round(h_df$StimOnset[ind_sf_hit], 1)
  dur_hit <- round(h_df$resp_trial_1.rt[ind_sf_hit], 2)

  # marry, account for NR
  if (length(start_hit) == 0) {
    mar_hit <- "*"
  } else {
    mar_hit <- paste0(start_hit, ":", dur_hit)
  }

  # write
  cat(mar_hit, "\n", file = out_hit, append = h_ap, sep = "\t")

  # same for misses
  ind_sf_miss <- which(h_df$stim_ID == "SceneFace" & h_df$trial_acc == 0 & h_df$resp_trial_1.keys != "None")
  start_miss <- round(h_df$StimOnset[ind_sf_miss], 1)
  dur_miss <- round(h_df$resp_trial_1.rt[ind_sf_miss], 2)
  if (length(start_miss) == 0) {
    mar_miss <- "*"
  } else {
    mar_miss <- paste0(start_miss, ":", dur_miss)
  }
  cat(mar_miss, "\n", file = out_miss, append = h_ap, sep = "\t")
}


# look at set behavior
data_set1 <- rbind(data_run1, data_run2)
data_set2 <- rbind(data_run3, data_run4)

ind_sf_s1 <- grep("SceneFace", data_set1[, 3])
test_acc1 <- data_set1[ind_sf_s1, 5]
plot(test_acc1, main = "Set 1 Scene-Face", ylab = "Accuracy", xlab = "Trial")

ind_sf_s2 <- grep("SceneFace", data_set2[, 3])
test_acc2 <- data_set2[ind_sf_s2, 5]
plot(test_acc2, main = "Set 2 Scene-Face", ylab = "Accuracy", xlab = "Trial")

roc_set1 <- roc(test_acc1, 1:length(test_acc1))
roc_set2 <- roc(test_acc2, 1:length(test_acc2))

a <- ggroc(list(set1 = roc_set1, set2 = roc_set2)) + ggtitle("Face-Scene Accuracy")
aa <- a + xlab("False Positive Rate") + ylab("True Positive Rate")
aa
