


args <- commandArgs()
dataDir <- args[6]
outDir <- args[7]
subjStr <- args[8]
numRuns <- args[9]

# dataDir <- "/Users/nmuncy/Projects/learn_mvpa/vCAT_data/vCAT_003"
# outDir <- "/Users/nmuncy/Projects/afni_python/sub-003/ses-S1"
# subjStr <- "vCAT_003"
# numRuns <- 2

for(run in 1:numRuns){
  
  # get data
  data_raw <- read.delim(paste0(dataDir,"/", subjStr, "_", run, "_localizer.csv"), sep = ",")
  
  # determine append
  if (run == 1) {
    h_ap <- F
  } else {
    h_ap <- T
  }
  
  # set out files
  out_scene <- paste0(outDir, "/tf_loc_scene.txt")
  out_face <- paste0(outDir, "/tf_loc_face.txt")
  out_num <- paste0(outDir, "/tf_loc_num.txt")

  # determine index, onset for scene, face, integer
  ind_scene <- grep("scene_img", data_raw$stim)
  ons_scene <- round(data_raw$StimOnset[ind_scene], 1)
  
  ind_face <- grep("face_img", data_raw$stim)
  ons_face <- round(data_raw$StimOnset[ind_face], 1)
  
  ind_num <- grep("^[0-9]", data_raw$stim)
  ons_num <- round(data_raw$StimOnset[ind_num], 1)
  
  # write out
  cat(ons_scene, "\n", file = out_scene, append = h_ap, sep = "\t")
  cat(ons_face, "\n", file = out_face, append = h_ap, sep = "\t")
  cat(ons_num, "\n", file = out_num, append = h_ap, sep = "\t")
}
