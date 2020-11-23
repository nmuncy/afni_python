

args <- commandArgs()
dataDir <- args[6]
outDir <- args[7]
subjStr <- args[8]
numRuns <- args[9]
phase <- args[10]


# # For testing
# subjStr <- "vCAT_005"
# dataDir <- paste0("/Users/nmuncy/Projects/learn_mvpa/vCAT_data/", subjStr)
# outDir <- paste0("/Users/nmuncy/Projects/afni_python/",subjStr,"/ses-S1")
# numRuns <- 4
# phase <- "task"
# run <- 1


for(run in 1:numRuns){

  # determine append
  if (run == 1) {
    h_ap <- F
  } else {
    h_ap <- T
  }

  if(phase == "loc"){

    # get data
    data_raw <- read.delim(paste0(dataDir,"/", subjStr, "_", run, "_localizer.csv"), sep = ",")

    # set out files
    #   3 timing files - for scene, face, and numbers
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

  if(phase == "task"){

    # get data
    data_raw <- read.delim(paste0(dataDir, "/", subjStr, "_simp_", phase, run, ".csv"), sep=",")

    # replace NA with NR for BL types
    data_raw$k_img[is.na(data_raw$k_img)] <- "NR"
    data_raw$stim[is.na(data_raw$stim)] <- "NR"

    # determine onsets for type (BL, cond)
    #   face vs scene
    #   responded, correct vs incorrect
    #   event and preceding event
    for(type in c("cond", "BL")){

      # get indices
      ind_hold <- grep(type, data_raw$trialtype)

      # start vectors
      for(i in c("face", "scene")){
        for(j in c("prec", "event")){
          for(k in c("cor", "icor")){
            assign(paste0("ons_",type,"_",i,"_",j,"_",k), vector())
          }
        }
      }

      # fill vectors
      #   get previous vector
      #   assign (append) new, rounded value
      for(i in ind_hold){
        if(data_raw$stim[i-1] == "face1" || data_raw$stim[i-1] == "face2"){
          if(data_raw$resp[i-1] != "None" && data_raw$acc[i-1] == 1){

            hold_event <- get(paste0("ons_",type,"_face_event_cor"))
            assign(paste0("ons_",type,"_face_event_cor"), c(hold_event, round(data_raw$onset[i],1)))

            hold_prec <- get(paste0("ons_",type,"_face_prec_cor"))
            assign(paste0("ons_",type,"_face_prec_cor"), c(hold_prec, round(data_raw$onset[i-1],1)))

          }else if(data_raw$resp[i-1] != "None" && data_raw$acc[i-1] == 0){

            hold_event <- get(paste0("ons_",type,"_face_event_icor"))
            assign(paste0("ons_",type,"_face_event_icor"), c(hold_event, round(data_raw$onset[i],1)))

            hold_prec <- get(paste0("ons_",type,"_face_prec_icor"))
            assign(paste0("ons_",type,"_face_prec_icor"), c(hold_prec, round(data_raw$onset[i-1],1)))
          }
        }else if(data_raw$stim[i-1] == "scene1" || data_raw$stim[i-1] == "scene2"){
          if(data_raw$resp[i-1] != "None" && data_raw$acc[i-1] == 1){

            hold_event <- get(paste0("ons_",type,"_scene_event_cor"))
            assign(paste0("ons_",type,"_scene_event_cor"), c(hold_event, round(data_raw$onset[i],1)))

            hold_prec <- get(paste0("ons_",type,"_scene_prec_cor"))
            assign(paste0("ons_",type,"_scene_prec_cor"), c(hold_prec, round(data_raw$onset[i-1],1)))

          }else if(data_raw$resp[i-1] != "None" && data_raw$acc[i-1] == 0){

            hold_event <- get(paste0("ons_",type,"_scene_event_icor"))
            assign(paste0("ons_",type,"_scene_event_icor"), c(hold_event, round(data_raw$onset[i],1)))

            hold_prec <- get(paste0("ons_",type,"_scene_prec_icor"))
            assign(paste0("ons_",type,"_scene_prec_icor"), c(hold_prec, round(data_raw$onset[i-1],1)))
          }
        }
      }

      # write out
      for(i in c("face", "scene")){
        for(j in c("prec", "event")){
          for(k in c("cor", "icor")){
            hold_out <- get(paste0("ons_",type,"_",i,"_",j,"_",k))
            if(length(hold_out) == 0){
              hold_out <- "*"
            }
            out_file <- paste0(outDir, "/tf_Study_", substr(type,1,1), substr(i,1,1), substr(j,1,1), substr(k,1,1),".txt")
            cat(hold_out, "\n", file = out_file, append = h_ap, sep = "\t")
          }
        }
      }
    }

  }
}



### Code for building each part individually
#     Highly repetitive

# # get indices
# ind_cond <- grep("COND", data_raw$stim)
# ind_bl <- grep("BL", data_raw$trialtype)
#
# # start vectors for cond
# #   face/scene * cond/precede * i/cor
# for(i in c("face", "scene")){
#   for(j in c("cond", "prec")){
#     for(k in c("cor", "icor")){
#       assign(paste0("ons_",i,"_",j,"_",k), vector())
#     }
#   }
# }
#
# # determine onsets for cond
# #   face versus scene
# #   responded, correct vs incorrect
# #   condition and fixed preceding condition
# for(i in ind_cond){
#   if(data_raw$stim[i-1] == "face1" || data_raw$stim[i-1] == "face2"){
#     if(data_raw$resp[i-1] != "None" && data_raw$acc[i-1] == 1){
#       ons_face_cond_cor <- c(ons_face_cond_cor, round(data_raw$onset[i],1))
#       ons_face_prec_cor <- c(ons_face_prec_cor, round(data_raw$onset[i-1],1))
#     }else if(data_raw$resp[i-1] != "None" && data_raw$acc[i-1] == 0){
#       ons_face_cond_icor <- c(ons_face_cond_icor, round(data_raw$onset[i],1))
#       ons_face_prec_icor <- c(ons_face_prec_icor, round(data_raw$onset[i-1],1))
#     }
#   }else if(data_raw$stim[i-1] == "scene1" || data_raw$stim[i-1] == "scene2"){
#     if(data_raw$resp[i-1] != "None" && data_raw$acc[i-1] == 1){
#       ons_scene_cond_cor <- c(ons_scene_cond_cor, round(data_raw$onset[i],1))
#       ons_scene_prec_cor <- c(ons_scene_prec_cor, round(data_raw$onset[i-1],1))
#     }else if(data_raw$resp[i-1] != "None" && data_raw$acc[i-1] == 0){
#       ons_scene_cond_icor <- c(ons_scene_cond_icor, round(data_raw$onset[i],1))
#       ons_scene_prec_icor <- c(ons_scene_prec_icor, round(data_raw$onset[i-1],1))
#     }
#   }
# }
#
# # start vectors for bl - same as above
# #   face/scene * i/corr * fixed/bl
# for(i in c("face", "scene")){
#   for(j in c("cor", "icor")){
#     for(k in c("bl", "blprec")){
#       assign(paste0("ons_",i,"_",k,"_",j), vector())
#     }
#   }
# }
#
# # determine onsets for bl - same as above
# #   face versus scene
# #   responded, corr vs incorrect
# for(i in ind_bl){
#   if(data_raw$stim[i-1] == "face1" || data_raw$stim[i-1] == "face2"){
#     if(data_raw$resp[i-1] != "None" && data_raw$acc[i-1] == 1){
#       ons_face_blprec_cor <- c(ons_face_blprec_cor, round(data_raw$onset[i-1],1))
#       ons_face_bl_cor <- c(ons_face_bl_cor, round(data_raw$onset[i],1))
#     }else if(data_raw$resp[i-1] != "None" && data_raw$acc[i-1] == 0){
#       ons_face_blprec_icor <- c(ons_face_blprec_icor, round(data_raw$onset[i-1],1))
#       ons_face_bl_icor <- c(ons_face_bl_icor, round(data_raw$onset[i],1))
#     }
#   }else if(data_raw$stim[i-1] == "scene1" || data_raw$stim[i-1] == "scene2"){
#     if(data_raw$resp[i-1] != "None" && data_raw$acc[i-1] == 1){
#       ons_scene_blprec_cor <- c(ons_scene_blprec_cor, round(data_raw$onset[i-1],1))
#       ons_scene_bl_cor <- c(ons_scene_bl_cor, round(data_raw$onset[i],1))
#     }else if(data_raw$resp[i-1] != "None" && data_raw$acc[i-1] == 0){
#       ons_scene_blprec_icor <- c(ons_scene_blprec_icor, round(data_raw$onset[i-1],1))
#       ons_scene_bl_icor <- c(ons_scene_bl_icor, round(data_raw$onset[i],1))
#     }
#   }
# }

