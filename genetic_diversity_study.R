library(stringr)
library(dplyr)
setwd("F:/HPV/genetic_diversity")

#we get all pileup files within the directory
filenames = dir(pattern="*.pileup")

#we generate an empty dataframe
data_summary=data.frame(sample=character(),genome=character(),nb_pos=integer(),
                        nb_lack=integer(),position_nucl=integer(),type_mut=character())

#for each pileup file
for(my_file in filenames){
  #we open the file as a tsv format
  file_name = sub("\\..*", "", my_file)
  my_data = readr::read_tsv(my_file, col_names = F)
  my_data = my_data[,c(1:5)]
  
  #we import a file that indicates the different HPV observed for each sample 
  my_references = read.table("reference_observed.txt", header = T, sep = "\t")
  my_references = my_references[my_references$sample==file_name,]
  my_references = c(my_references$references)
  
  final_df = data.frame(1)
  
  #for each type in the sample
  for(one_ref in my_references){
    print(one_ref)
    my_data_focus = my_data[my_data$X1==one_ref,]
    ref_genome = my_data_focus$X1[1]
    
    #we set all variables to 0
    lacking_data=0
    my_data_focus$prop_lack=0
    my_data_focus$match=0
    my_data_focus$mut_A=0
    my_data_focus$mut_T=0
    my_data_focus$mut_G=0
    my_data_focus$mut_C=0
    my_data_focus$insert=0
    my_data_focus$del=0
  
    #we count the number and type of observation at each position of the genome
    for(i in 1:nrow(my_data_focus)){
      #we consider a lacking position when the deapth of coverage is less than 10X
      if(my_data_focus$X4[i]<10){
        lacking_data=lacking_data+1
      } else{
        my_data_focus$prop_lack[i] = 0
        my_data_focus$match[i] = str_count(my_data_focus$X5[i], "\\.")
        my_data_focus$match[i] = my_data_focus$match[i] + str_count(my_data_focus$X5[i], "\\,")
        my_data_focus$mut_A[i] = str_count(my_data_focus$X5[i], "A|a")
        my_data_focus$mut_T[i] = str_count(my_data_focus$X5[i], "T|t")
        my_data_focus$mut_G[i] = str_count(my_data_focus$X5[i], "G|g")
        my_data_focus$mut_C[i] = str_count(my_data_focus$X5[i], "C|c")
        my_data_focus$insert[i] = str_count(my_data_focus$X5[i], "\\+")
        my_data_focus$del[i] = str_count(my_data_focus$X5[i], "\\*")
        #my_data_focus$remo_del[i] = str_count(my_data_focus$X5[i], "\\-")
      }
    }
    
    #we calculate the percentage of each observation
    my_data_focus$prop_match= my_data_focus$match/my_data_focus$X4*100
    my_data_focus$prop_mut_A = my_data_focus$mut_A/my_data_focus$X4*100
    my_data_focus$prop_mut_T = my_data_focus$mut_T/my_data_focus$X4*100
    my_data_focus$prop_mut_G = my_data_focus$mut_G/my_data_focus$X4*100
    my_data_focus$prop_mut_C = my_data_focus$mut_C/my_data_focus$X4*100
    my_data_focus$prop_insert = my_data_focus$insert/my_data_focus$X4*100
    my_data_focus$prop_del = my_data_focus$del/my_data_focus$X4*100

    #for each position, we define the consensus base
    for(i in 1:nrow(my_data_focus)){
      #if the depth is less than 10X, it is a lacking position
      if(my_data_focus$X4[i]<10){
        my_data_focus$prop_lack[i] = 100
        my_data_focus$prop_match[i] = 0
        my_data_focus$prop_mut_A[i] = 0
        my_data_focus$prop_mut_T[i] = 0
        my_data_focus$prop_mut_G[i] = 0
        my_data_focus$prop_mut_C[i] = 0
        my_data_focus$prop_insert[i] = 0
        my_data_focus$prop_del[i] = 0
      }

      my_data_focus$obs[i] = sum(my_data_focus[i,c(14:20)]>50)
      
      #if the proportion of insertion is > 50%, we add the additional bases
      if(my_data_focus$prop_insert[i]>50){
        test = my_data_focus$X5[i]
        test2 = str_split(test, "\\+")
        test3=as.data.frame(test2, col.names = "observ")
        for(j in 1:nrow(test3)){
          test3$observ[j]=gsub("\\,","",test3$observ[j])
          test3$observ[j]=gsub("\\.","",test3$observ[j])
          test3$observ[j]=gsub("\\$","",test3$observ[j])
          test3$observ[j]=toupper(test3$observ[j])
        }
        test3$number=1
        #we check the different insertion detected
        test4 = test3 %>%
          group_by(observ) %>%
          summarise(total=sum(number)) %>%
          top_n(1, total)
        #we consider the insertion the most present
        insert_seq = test4$observ[1]
        insert_seq = gsub('[[:digit:]]+', '', insert_seq)
        my_data_focus$prop_match[i] = 0
        my_data_focus$prop_mut_A[i] = 0
        my_data_focus$prop_mut_T[i] = 0
        my_data_focus$prop_mut_G[i] = 0
        my_data_focus$prop_mut_C[i] = 0
        my_data_focus$prop_insert[i] = test4$total[1]/my_data_focus$X4[i]*100
        my_data_focus$prop_del[i] = 0        
      }
      
      #if the proportion of bases is >50% of match, we consider a match
      #else, it is a mutation
      if(my_data_focus$prop_match[i]>50){
        my_data_focus$sequence[i]=my_data_focus$X3[i]
      }else if(my_data_focus$prop_mut_A[i]>50){
        my_data_focus$sequence[i]="A"
      }else if(my_data_focus$prop_mut_T[i]>50){
        my_data_focus$sequence[i]="T"
      }else if(my_data_focus$prop_mut_G[i]>50){
        my_data_focus$sequence[i]="G"
      }else if(my_data_focus$prop_mut_C[i]>50){
        my_data_focus$sequence[i]="C"
      }else if(my_data_focus$prop_insert[i]>50){
        my_data_focus$sequence[i]=paste(my_data_focus$X3[i], insert_seq, sep="")
      }else if(my_data_focus$prop_del[i]>50){
        my_data_focus$sequence[i]="-"
      }else{
        my_data_focus$sequence[i]="N"
      }
    }
    
    #we generate a sub_data table
    my_sub_data = my_data_focus[,c(1,2,3,22)]

    #for each position, we check if the consensus base is different to the reference
    #if it is the case, we add the position in a dataframe
    for(my_position in 1:nrow(my_sub_data)){
      if(my_sub_data$X3[my_position] != my_sub_data$sequence[my_position]){
        if(my_sub_data$sequence[my_position] == "N"){

        }else{
          data_summary[nrow(data_summary) + 1,] = list(file_name, one_ref, dim(my_sub_data)[1],lacking_data,
                                                       my_position, "mutation")
        }
      }
    }
  }
}
    
#write the table into a file
write.table(x = data_summary, file = "LIST_MUTATIONS.txt",
            sep="\t", col.names = F, row.names = F, quote = F)


#import the table if you disconnected the session
#data_summary=read.table("LIST_MUTATIONS.txt", sep="\t", header=T)

#count the number of mutations per sample and type
mutation_count = data_summary %>%
  group_by(sample, genome) %>%
  summarise(COUNT=n())

#import a file that indicates the HPV types for each genbank accession number
HPV_ref = read.table(file = "HPV_ref_type.txt", header=T, sep="\t")

#we list the different genomes
list_genomes=unique(mutation_count$genome)

#we perform the correspondance to add the HPV type
for(i in 1:nrow(mutation_count)){
  for(j in 1:nrow(HPV_ref)){
    if(mutation_count$genome[i]==HPV_ref$genome[j]){
      mutation_count$type_HPV[i]=as.character(HPV_ref$HPV_type[j])
    }
  }
}

#we indicate the different genbank identifiers per category of risk
#HIGH RISK
high_risk=c("K02718.1", "X05015.1", "J04353.1", "M12732.1",
            "X74477.1","M62849.1", "X74479.1", "M62877.1",
            "X74481.1","X74483.1","D90400.1", "X77858.1",
            "U31794.1", "DQ080079.1")

#LOW RISK
low_risk=c("X00203.1", "M14119.1", "M73236.1", "AJ620205.1",
           "U31788.1", "U37488.1", "U31793.1", "X94164.1",
           "AJ620209.1", "AF436128.1", "AY057438.1")

#INTERMEDIATE RISK
intermediate_risk=c("X74482.1", "D21208.1", "U21941.1", "X94165.1",
                    "AB027021.1")

#UNKNOWN RISK
unknown_risk=c("X74475.1", "X74476.1", "AY395706.1", "AB027020.1",
               "AF436130.1", "AF151983.1", "AF293960.1", "AF349909.1",
               "AJ400628.2", "DQ080081.1", "GQ244463.1")

#we add the risk category in the table
for(i in 1:nrow(mutation_count)){
  if(mutation_count$genome[i] %in% high_risk){
    mutation_count$risk[i]="high"
  }else if(mutation_count$genome[i] %in% low_risk){
    mutation_count$risk[i]="low"
  }else if(mutation_count$genome[i] %in% intermediate_risk){
    mutation_count$risk[i]="intermediate"
  }else{
    mutation_count$risk[i]="unknown"
  }
}

library(ggplot2)

#we produce a plot that shows the number of mutations per HPV type
P1 = ggplot(data = mutation_count, aes(x = type_HPV, y=COUNT, fill=risk))+
  geom_boxplot(outlier.shape = NA)+
  scale_x_discrete(limits=c("HPV16","HPV18","HPV31","HPV33","HPV35","HPV39","HPV45","HPV51","HPV52","HPV56","HPV58", "HPV59", "HPV66", "HPV68",
                            "HPV53","HPV67","HPV70","HPV73","HPV82",
                            "HPV6","HPV11","HPV42","HPV43","HPV44","HPV54","HPV61","HPV72","HPV81","HPV89","HPV90",
                            "HPV32","HPV34","HPV62","HPV69","HPV74","HPV83","HPV84","HPV86","HPV87","HPV101","HPV114"))+
  scale_fill_manual(values=c("#c00000", "#ffd966", "#fff2cc", "#bacee2"))+
  theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        #       text=element_text(size = 16, family="Arial"),
        axis.text.x=element_text(colour="black", size = 11, angle=45),
        axis.text.y=element_text(colour="black", size = 16),
        legend.position = "right")

#show the plot
P1

#we produce a plot that shows the number of mutations per risk category
P2 = ggplot(data=mutation_count, aes(x=risk, y=COUNT, fill=risk))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_manual(values=c("#c00000", "#ffd966", "#fff2cc", "#bacee2"))+
  ylim(0,500)+
  theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        #       text=element_text(size = 16, family="Arial"),
        axis.text.x=element_text(colour="black", size = 11, angle=45),
        axis.text.y=element_text(colour="black", size = 16),
        legend.position = "right")

#show the plot
P2

#generate sub tables for low and high risks
#then perform a mann whitney u test to check if one category has more mutations
t1=mutation_count[mutation_count$risk=="low",]
t2=mutation_count[mutation_count$risk=="high",]
wilcox.test(t1$COUNT,t2$COUNT)

#write the table
write.table(x = mutation_count, file = "compare.txt",
            sep="\t", col.names = F, row.names = F, quote = F)


#import the gff file to perform analysis at the gene level
gff_file = readr::read_tsv(file = "HPV_annotation.gff")

#we update the data, especially when two genes overlapped
gff_file = gff_file[gff_file$V3=="gene",]
gff_file$V4 = as.integer(gff_file$V4)
gff_file$V5 = as.integer(gff_file$V5)
gff_file = gff_file[((gff_file$V5)-(gff_file$V4))>10,]
gff_file=gff_file[,c(1,4,5,10,11)]

gff_file$clean=""
for(i in 2:nrow(gff_file)){
  if(gff_file$V1[i]==gff_file$V1[i-1]){
    if(gff_file$V4[i]<=gff_file$V5[i-1]){
      if(gff_file$V10[i]=="E4"){
        temp_val=gff_file$V5[i-1]
        gff_file$V5[i-1]=((gff_file$V4[i])-1)
        gff_file[nrow(gff_file) + 1,] = list(gff_file$V1[i], ((gff_file$V5[i])+1), temp_val, "E2", gff_file$V11[i],"")
      }else if(gff_file$V10[i]=="E2" && gff_file$V10[i-1]=="E4"){
        temp_val=gff_file$V5[i-1]
        gff_file$V4[i]=temp_val+1
        gff_file$clean[i]="remove"
      } else{
      temp_val=gff_file$V5[i-1]
      gff_file$V5[i-1]=((gff_file$V4[i])-1)
      gff_file[nrow(gff_file) + 1,] = list(gff_file$V1[i], ((gff_file$V5[i-1])+1), temp_val, paste(gff_file$V10[i-1],gff_file$V10[i], sep="-"), gff_file$V11[i],"")
      gff_file$V4[i]=temp_val+1
      }
      
    }
  }
}

gff_file=gff_file[gff_file$clean=="",]
gff_file$size=((gff_file$V5)-(gff_file$V4))

#we save the updated annotation file
write.table(gff_file, "coordinates3.txt", sep="\t")

#we import the novel gff file if needed
#gff_file = read.table("coordinates3.txt", header = T, sep="\t")

#we indicate the location (i.e. the gene) for each mutation observed
data_summary$location=""
data_summary$size=NA
for(i in 1:nrow(data_summary)){
  for(j in 1:nrow(gff_file)){
    if(data_summary$genome[i]==gff_file$V1[j]){
      if(data_summary$position_nucl[i] >= gff_file$V4[j] && data_summary$position_nucl[i] <= gff_file$V5[j]){
        data_summary$location[i]=gff_file$V10[j]
        data_summary$size[i]=gff_file$size[j]
      }
    }
  }
}

data_summary$number=1

#we indicates the risk category where a mutation was observed
for(i in 1:nrow(data_summary)){
  if(data_summary$genome[i] %in% high_risk){
    data_summary$risk[i]="high"
  }else if(data_summary$genome[i] %in% low_risk){
    data_summary$risk[i]="low"
  }else if(data_summary$genome[i] %in% intermediate_risk){
    data_summary$risk[i]="intermediate"
  }else{
    data_summary$risk[i]="unknown"
  }
}

#we count the number of mutations per sample, reference and gene
summarize_mut = data_summary %>%
  group_by(sample,genome, location, size) %>%
  summarise(total=sum(number)) 

#we normalize the number of mutations with the size of the gene
summarize_mut$prop = summarize_mut$total/summarize_mut$size

#we indicates the risk category
for(i in 1:nrow(summarize_mut)){
  if(summarize_mut$genome[i] %in% high_risk){
    summarize_mut$risk[i]="high"
  }else if(summarize_mut$genome[i] %in% low_risk){
    summarize_mut$risk[i]="low"
  }else if(summarize_mut$genome[i] %in% intermediate_risk){
    summarize_mut$risk[i]="intermediate"
  }else{
    summarize_mut$risk[i]="unknown"
  }
}

#we exclude mutations not located or a gene, or mutations that overlapped two genes to avoid any bias
summarize_mut = summarize_mut[summarize_mut$location!="L2-L1",]
summarize_mut = summarize_mut[summarize_mut$location!="E1-E2",]
summarize_mut = summarize_mut[summarize_mut$location!="E6-E7",]
summarize_mut = summarize_mut[summarize_mut$location!="",]
summarize_mut = summarize_mut[summarize_mut$location!="E4-16",]

#we plot the ratio of mutations per gene and per risk category
plot3 = ggplot(data = summarize_mut, aes(x=location, y=prop, fill=risk))+
  geom_boxplot(outlier.shape = NA)+
  ylim(0,0.06)+
  scale_fill_manual(values=c("#c00000", "#ffd966", "#fff2cc", "#bacee2"))+
  theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        #       text=element_text(size = 16, family="Arial"),
        axis.text.x=element_text(colour="black", size = 11, angle=45),
        axis.text.y=element_text(colour="black", size = 16),
        legend.position = "right")
  
#show the plot
plot3

#focus on the E6 gene for low and high risks
summarize_mut_low = summarize_mut[summarize_mut$risk=="low",]
summarize_mut_high = summarize_mut[summarize_mut$risk=="high",]
summarize_mut_high_E6 = summarize_mut_high[summarize_mut_high$location=="E6",]
summarize_mut_low_E6 = summarize_mut_low[summarize_mut_low$location=="E6",]
#we check for the difference and calculate the medians
wilcox.test(summarize_mut_high_E6$prop, summarize_mut_low_E6$prop)
median(summarize_mut_high_E6$prop)
median(summarize_mut_low_E6$prop)
#we combine the two sub datasets
combine = rbind(summarize_mut_high_E6, summarize_mut_low_E6)

##we plot the ratio of mutations betwee low and risk risks for E6
plot4 = ggplot(data = combine, aes(x=risk, y=prop, fill=risk))+
  geom_boxplot(outlier.shape = NA)+
  ylim(0,0.03)+
  scale_fill_manual(values=c("#c00000",  "#fff2cc"))+
  theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 11, angle=45),
        axis.text.y=element_text(colour="black", size = 16),
        legend.position = "right")

#show the plot
plot4
