#example for HPV6
#we present one example for E1 gene
#we also add the code for E2, since it is composed of two regions

setwd("D:/HPV6/")
#we import the alignment sequences for HPV6
my_sequences <- readLines("hpv6_aligned_clean.fasta")
type="HPV6"
risk="low"

#we import different libraries
library(stringr)
library(data.table)
library(dplyr)

#we import the annotation file for HPV6, and focus on E1 gene
gff_file = read.table(file = "coordinates.txt", header=T)
gff_sub = gff_file[gff_file$type==type,]
gff_sub = gff_sub[gff_sub$gene=="E1",]

#summary_table initialization
summary_apobec = data.frame(id_seq=character(), type=character(),risk=character(),
                            nb_mut_apo=integer(), nb_tcw=integer(), nb_mut_ct_no_tcw=integer(), nb_c_no_tcw=integer(), 
                            ratio_apo=double(), ratio_non_apo=double(),fold_change=double(),pval=double())

nb_tcw = 0
nb_no_tcw = 0
#firstline = reference sequence
#nb_tcw and nb_no_tcw define before and not include in the loop because done just one time (based on the ref)
for(my_seq in 1:length(my_sequences)){
  print(my_seq)
  nb_mut_tcw = 0
  nb_mut_no_tcw = 0
  #even numbers correspond to sequences
  if(my_seq%%2 == 0){
    #we read the reference genome
    if(my_seq == 2){
      ref_split <- strsplit(my_sequences, "")[[2]]
      for(i in 3:length(ref_split)){
        if(i>=gff_sub$start && i<=gff_sub$end){
          pos1 = ref_split[i-2]
          pos2 = ref_split[i-1]
          pos3 = ref_split[i]
          motif = paste(ref_split[i-2], ref_split[i-1], ref_split[i], "")
          motif = gsub(" ", "", motif, fixed = TRUE)
          #we check if we have a TCW motif
          if(motif=="tca" || motif=="tct"){
            nb_tcw=nb_tcw+1
          }else{
            if(pos2=="c" && motif!="tca"){
              nb_no_tcw=nb_no_tcw+1
            }else if(pos2=="c" && motif!="tct"){
              nb_no_tcw=nb_no_tcw+1
            }
          }
        }
      }
    #we look for other sequences than the reference
    }else{
      seq_split <- strsplit(my_sequences, "")[[my_seq]]
      for(i in 3:length(ref_split)){
        if(i>=gff_sub$start && i<=gff_sub$end){
          pos1_ref = ref_split[i-2]
          pos2_ref = ref_split[i-1]
          pos3_ref = ref_split[i]
          pos1_seq = seq_split[i-2]
          pos2_seq = seq_split[i-1]
          pos3_seq = seq_split[i]
          motif_ref = paste(pos1_ref, pos2_ref, pos3_ref, "")
          motif_ref = gsub(" ", "", motif_ref, fixed = TRUE)
          motif_seq = paste(pos1_seq, pos2_seq, pos3_seq, "")
          motif_seq = gsub(" ", "", motif_seq, fixed = TRUE)
          #we check if we have a TCW motif in reference
          if(motif_ref=="tca" || motif_ref=="tct"){
            #if yes, we check if we have a c>t mutation
            if(motif_seq=="tta" || motif_seq=="ttt"){
              nb_mut_tcw = nb_mut_tcw+1
            }
          }else{
            #this is not a c>t mutation in TCW
            if(pos2_ref=="c" && motif_ref!="tca"){
              if(pos2_seq=="t"){
                nb_mut_no_tcw=nb_mut_no_tcw+1
              }
            }else if(pos2_ref=="c" && motif_ref!="tct"){
              if(pos2_seq=="t"){
                nb_mut_no_tcw=nb_mut_no_tcw+1
              }
            }
          }
        }
      }
    }
    
    #we perform a fisher test to check for the significance
    my_fisher = fisher.test(matrix(c((nb_mut_tcw),((nb_tcw)-(nb_mut_tcw)),(nb_mut_no_tcw),((nb_no_tcw)-(nb_mut_no_tcw))),2,2, byrow=TRUE)) 
    my_p_val = my_fisher$p.value
    
    #we calculate the ratios and the ratioC>T
    ratio_apo = (nb_mut_tcw/nb_tcw*100)
    ratio_non_apo = (nb_mut_no_tcw/nb_no_tcw*100)
    fold_change = ratio_apo/ratio_non_apo
    
    #we count the observations for each sequence
    summary_apobec[nrow(summary_apobec) + 1,] = list(my_sequences[my_seq-1], risk, type,
                                                     nb_mut_tcw, nb_tcw,
                                                     nb_mut_no_tcw, nb_no_tcw, ratio_apo,
                                                     ratio_non_apo, fold_change,my_p_val)
  }
}

#we consider a ratioC>T >=2 as a gene sequence enriched in APOBEC mutations
summary_apobec=summary_apobec[summary_apobec$fold_change>=2,]
summary_apobec=summary_apobec[complete.cases(summary_apobec),]



########
#SPECIFICALLY FOR E2 (because two rows in gff, not one)
#we import the annotation file for HPV6, and focus on E2 gene
gff_file = read.table(file = "coordinates.txt", header=T)
gff_sub = gff_file[gff_file$type==type,]
gff_sub = gff_sub[gff_sub$gene=="E2",]

#summary_table initialization
summary_apobec = data.frame(id_seq=character(), type=character(),risk=character(),
                            nb_mut_apo=integer(), nb_tcw=integer(), nb_mut_ct_no_tcw=integer(), nb_c_no_tcw=integer(), 
                            ratio_apo=double(), ratio_non_apo=double(),fold_change=double(),pval=double())


nb_tcw = 0
nb_no_tcw = 0
#firstline = reference sequence
#nb_tcw and nb_no_tcw define before and not include in the loop because done just one time (based on the ref)
for(my_seq in 1:length(my_sequences)){
  print(my_seq)
  nb_mut_tcw = 0
  nb_mut_no_tcw = 0
  #even numbers correspond to sequences
  if(my_seq%%2 == 0){
    #we read the reference genome
    if(my_seq == 2){
      ref_split <- strsplit(my_sequences, "")[[2]]
      for(i in 3:length(ref_split)){
        #we check the two regions of E2
        if((i>=gff_sub$start[1] && i<=gff_sub$end[1]) || (i>=gff_sub$start[2] && i<=gff_sub$end[2])){
          pos1 = ref_split[i-2]
          pos2 = ref_split[i-1]
          pos3 = ref_split[i]
          motif = paste(ref_split[i-2], ref_split[i-1], ref_split[i], "")
          motif = gsub(" ", "", motif, fixed = TRUE)
          #we check if we have a TCW motif
          if(motif=="tca" || motif=="tct"){
            tcw_motif[nrow(tcw_motif) + 1,] = list(my_sequences[my_seq-1], risk, type, i) ####2
            nb_tcw=nb_tcw+1
          }else{
            if(pos2=="c" && motif!="tca"){
              c_no_tcw[nrow(c_no_tcw) + 1,] = list(my_sequences[my_seq-1], risk, type, i) ####4
              nb_no_tcw=nb_no_tcw+1
            }else if(pos2=="c" && motif!="tct"){
              c_no_tcw[nrow(c_no_tcw) + 1,] = list(my_sequences[my_seq-1], risk, type, i) ####4
              nb_no_tcw=nb_no_tcw+1
            }
          }
        }
      }
    #we look for other sequences than the reference
    }else{
      seq_split <- strsplit(my_sequences, "")[[my_seq]]
      for(i in 3:length(ref_split)){
        if((i>=gff_sub$start[1] && i<=gff_sub$end[1]) || (i>=gff_sub$start[2] && i<=gff_sub$end[2])){
          pos1_ref = ref_split[i-2]
          pos2_ref = ref_split[i-1]
          pos3_ref = ref_split[i]
          pos1_seq = seq_split[i-2]
          pos2_seq = seq_split[i-1]
          pos3_seq = seq_split[i]
          motif_ref = paste(pos1_ref, pos2_ref, pos3_ref, "")
          motif_ref = gsub(" ", "", motif_ref, fixed = TRUE)
          motif_seq = paste(pos1_seq, pos2_seq, pos3_seq, "")
          motif_seq = gsub(" ", "", motif_seq, fixed = TRUE)
          #we check if we have a TCW motif in reference
          if(motif_ref=="tca" || motif_ref=="tct"){
            #if yes, we check if we have a c>t mutation
            if(motif_seq=="tta" || motif_seq=="ttt"){
              nb_mut_tcw = nb_mut_tcw+1
            }
          }else{
            #this is not a c>t mutation in TCW
            if(pos2_ref=="c" && motif_ref!="tca"){
              if(pos2_seq=="t"){
                nb_mut_no_tcw=nb_mut_no_tcw+1
              }
            }else if(pos2_ref=="c" && motif_ref!="tct"){
              if(pos2_seq=="t"){
                nb_mut_no_tcw=nb_mut_no_tcw+1
              }
            }
          }
        }
      }
    }
    
    #we perform a fisher test to check for the significance
    my_fisher = fisher.test(matrix(c((nb_mut_tcw),((nb_tcw)-(nb_mut_tcw)),(nb_mut_no_tcw),((nb_no_tcw)-(nb_mut_no_tcw))),2,2, byrow=TRUE)) 
    my_p_val = my_fisher$p.value
    
    #we calculate the ratios and the ratioC>T
    ratio_apo = (nb_mut_tcw/nb_tcw*100)
    ratio_non_apo = (nb_mut_no_tcw/nb_no_tcw*100)
    fold_change = ratio_apo/ratio_non_apo
    
    #we count the observations for each sequence
    summary_apobec[nrow(summary_apobec) + 1,] = list(my_sequences[my_seq-1], risk, type,
                                                     nb_mut_tcw, nb_tcw,
                                                     nb_mut_no_tcw, nb_no_tcw, ratio_apo,
                                                     ratio_non_apo, fold_change,my_p_val)
  }
}

#we consider a ratioC>T >=2 as a gene sequence enriched in APOBEC mutations
summary_apobec=summary_apobec[summary_apobec$fold_change>=2,]
summary_apobec=summary_apobec[complete.cases(summary_apobec),]


