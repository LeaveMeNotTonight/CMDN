library(tidyverse)
library(CluMSID)
library(CluMSIDdata)
library(grid)
library(OrgMassSpecR)
library(pheatmap)
library(reshape2)
library(MSMSsim)
library(msentropy)
library(readxl)
#rm(list=ls())

#1.  Differential expression metabolic feature space
S9Metabolic_space<-read.csv("Area_PG1.csv")
DEM_space<-subset(S9Metabolic_space,S9Metabolic_space$FC>2&S9Metabolic_space$p<0.05&S9Metabolic_space$S.N.average>=10)
rownames(DEM_space$Average.Mz)<-NULL
write.csv(DEM_space,"DEM_space.csv")
# Analyze parent compound
ParentCompoundList<-S9Metabolic_space %>%
  filter(Metabolite.name != "Unknown")
ParentCompoundChangedList<-ParentCompoundList %>%
  filter(p<=0.05&FC<1)
ParentCompoundUnchangedList<-ParentCompoundList %>%
  filter(!(p<=0.05&FC<1))
write.csv(ParentCompoundChangedList,"ParentCompoundChangedList.csv")
#2.  Rule-based MS1 biotransformation analysis
#2.1.1 Import metabolism rule & compound list; generate rule-based metabolite MS1 information;
#Meta_P1_rule<-read.csv("Phase 1 reaction rule all.csv")
Meta_P1_rule<-read_excel("C:/ZHD/UNC/UNC_Research/Pesticide_Project/Pesticide-S9-New/Phase 1 reaction rule Final (1+2).xlsx", sheet = 1)
CompoundList<-read.csv("PG1 parent.csv")
n_P1_rule <- nrow(Meta_P1_rule)
n_CompoundList <- nrow(CompoundList)
CMList <- CompoundList[rep(1:n_CompoundList, each = n_P1_rule), ]
rownames(CMList) <- NULL
Meta_P1List<-replicate(n_CompoundList,Meta_P1_rule,simplify = FALSE)
Meta_P1List_df <- do.call(rbind, Meta_P1List)
RB_MS1_List <- cbind(CMList, Meta_P1List_df)
RB_MS1_List$RB_Metabolite_Name <- paste(RB_MS1_List$Compound.Name, RB_MS1_List$Reaction_Type, sep = ":")
RB_MS1_List$RB_Metabolite_mz <-RB_MS1_List$m.z + RB_MS1_List$Mass_Difference
RB_MS1_List_filtered <- RB_MS1_List %>%
  filter(ifelse(!grepl("S", Molecular_Formula), !grepl("S", Reaction_Type), TRUE))
RB_MS1_List_filtered <- RB_MS1_List_filtered %>%
  filter(ifelse(!grepl("Cl", Molecular_Formula), !grepl("Cl", Reaction_Type), TRUE))
RB_MS1_List_filtered <- RB_MS1_List_filtered %>%
  filter(ifelse(!grepl("F", Molecular_Formula), !grepl("F", Reaction_Type), TRUE))
RB_MS1_List_filtered <- RB_MS1_List_filtered %>%
  filter(ifelse(!grepl("Br", Molecular_Formula), !grepl("Br", Reaction_Type), TRUE))
write.csv(RB_MS1_List_filtered,"RuleBased_MS1_FeatureList.csv")

 #2.2 Biotransformer prediction list; Add the biotransformation result to the MS1 pseudo list
RB_MS1_List<-read_excel("PG1_RuleBased+combo+biotransformer_MS1_FeatureList.xlsx")

#find the RB-list metabolites in DEM space
RB_MS1_List.expanded.total <- tibble()

for(i in 1:nrow(RB_MS1_List))
{
  #i <- 43
  RB_MS1_List.i <- RB_MS1_List[i , ]
  target.i <- RB_MS1_List$RB_Metabolite_mz[i]
  DEM_space.i <- DEM_space %>% 
    filter(.data = . , 
           abs((Average.Mz - target.i) / target.i) <= 5e-6)
  match.i <- nrow(DEM_space.i)
  
  if(match.i == 0)
  {
    RB_MS1_List.expanded.i <- RB_MS1_List.i %>% 
      mutate(.data = . , 
             Average.Mz = NA , 
             Average.Rt = NA , 
             Average.AID = NA)
  }
  if(match.i == 1)
  {
    RB_MS1_List.expanded.i <- RB_MS1_List.i %>% 
      mutate(.data = . , 
             Average.Mz = DEM_space.i$Average.Mz , 
             Average.Rt = DEM_space.i$Average.Rt.min. , 
             Average.AID = DEM_space.i$Alignment.ID)
  }
  if(match.i > 1)
  {
    RB_MS1_List.expanded.i <- rbind(RB_MS1_List.i , RB_MS1_List.i[rep(1, match.i-1), ]) %>% 
      mutate(.data = . , 
             Average.Mz = DEM_space.i$Average.Mz , 
             Average.Rt = DEM_space.i$Average.Rt.min. , 
             Average.AID = DEM_space.i$Alignment.I)
  }
  
  if(i == 1)
  {
    RB_MS1_List.expanded.total <- RB_MS1_List.expanded.i
  }
  else
  {
    RB_MS1_List.expanded.total <- rbind(RB_MS1_List.expanded.total , RB_MS1_List.expanded.i)
  }
}

Integrated_RB_MS1_List<-subset(RB_MS1_List.expanded.total,Average.Mz>=0)
Integrated_RB_MS1_List$mz_defect<-(abs(Integrated_RB_MS1_List$Average.Mz - Integrated_RB_MS1_List$RB_Metabolite_mz) / Integrated_RB_MS1_List$RB_Metabolite_mz)*1e6
Integrated_RB_MS1_List$RT_shift<-(Integrated_RB_MS1_List$Average.Rt-Integrated_RB_MS1_List$rt)
dup_names <- Integrated_RB_MS1_List$RB_Metabolite_Name[duplicated(Integrated_RB_MS1_List$RB_Metabolite_Name)]
for (RB_Metabolite_Name in dup_names) {
  indices <- which(Integrated_RB_MS1_List$RB_Metabolite_Name == RB_Metabolite_Name)
  num_duplicates <- length(indices)
  if (num_duplicates > 1) {
    Integrated_RB_MS1_List$RB_Metabolite_Name[indices] <- paste(RB_Metabolite_Name, seq_len(num_duplicates), sep = "-")
  }
}
Integrated_RB_MS1_List_unique <- Integrated_RB_MS1_List %>%
  distinct(Average.AID, .keep_all = TRUE)
write.csv(Integrated_RB_MS1_List,"Integrated_RB_MS1_List_plusBiotransformer.csv")
#write.csv(Integrated_RB_MS1_List_unique,"Integrated_RB_MS1_List_plusBiotransformer_unique.csv")
#Integrated_RB_MS1_List_duplicates <- Integrated_RB_MS1_List[duplicated(Integrated_RB_MS1_List$Average.AID) | duplicated(Integrated_RB_MS1_List$Average.AID, fromLast = TRUE), ]
#Integrated_RB_MS1_List_duplicates_unique <- Integrated_RB_MS1_List_duplicates %>%
#  distinct(Average.AID, .keep_all = TRUE)


#2.4 Build Inclusion list for PRM MS2 collection.

il.total <- Integrated_RB_MS1_List %>% 
  arrange(. , Average.Rt) %>% 
  mutate(. , `inclusion list batch` = rep(1:15 , 2000)[1:nrow(Integrated_RB_MS1_List)])

colnames.template <- read_csv(file = 'C:/ZHD/UNC/UNC_Research/Pesticide_Project/Pesticide-S9-Metabolism/021424PG8-S9-Analysis/MS2InclusionListTemplate.csv') %>% 
  colnames()
inclusionlist.folderpath <- 'C:/ZHD/UNC/UNC_Research/Pesticide_Project/Pesticide-S9-New/PG1-New/PG1-IL/'

for(i in levels(factor(il.total$`inclusion list batch`)))
{
  #i <- 3
  il.i <- il.total %>%
    filter(. , 
           `inclusion list batch` == i)
  
  overlaptotal.i <- c()
  for(k in 1:nrow(il.i))
  {
    #k <- 3
    rtmed.i <- il.i$Average.Rt[k]
    
    overlap.i <- il.i[-k , ] %>% 
      filter(. , 
             abs(Average.Rt - rtmed.i) <= 0.1) %>% 
      nrow()
    overlaptotal.i <- c(overlaptotal.i , overlap.i)
  }
  
  export.i <- data.frame(V1 = il.i$Average.Mz , 
                         V2 = NA , 
                         V3 = 'Chemical formula' , 
                         V4 = '+ H' , 
                         V5 = 1 , 
                         V6 = 'Positive' , 
                         V7 = il.i$Average.Rt - 0.25 , 
                         V8 = il.i$Average.Rt + 0.25 , 
                         V9 = NA , 
                         V10 = 'NCE' , 
                         V11 = NA , 
                         V12 = paste0(il.i$`xcms feature id` , 
                                      '_mz' , 
                                      round(il.i$Average.Mz , 3) , 
                                      '_rt' , 
                                      round(il.i$Average.Rt , 2)))
  colnames(export.i) <- colnames.template
  
  write_csv(x = export.i , 
            file = paste0(inclusionlist.folderpath , 
                          'PG1-' , 
                          max(overlaptotal.i) , 
                          '_InclusionList_' , 
                          i , 
                          '.csv') , 
            na = '')
  
}
#2.5 Build Inclusion list for STD
parentCompound<-read.csv("PG1 parent.csv")
il.total <- parentCompound %>% 
  arrange(. , RT) %>% 
  mutate(. , `inclusion list batch` = rep(1:5 , 2000)[1:nrow(parentCompound)])

colnames.template <- read_csv(file = 'C:/ZHD/UNC/UNC_Research/Pesticide_Project/Pesticide-S9-Metabolism/021424PG8-S9-Analysis/MS2InclusionListTemplate.csv') %>% 
  colnames()
inclusionlist.folderpath <- 'C:/ZHD/UNC/UNC_Research/Pesticide_Project/Pesticide-S9-New/PG1-New/PG1STD-IL/'

for(i in levels(factor(il.total$`inclusion list batch`)))
{
  #i <- 3
  il.i <- il.total %>%
    filter(. , 
           `inclusion list batch` == i)
  
  overlaptotal.i <- c()
  for(k in 1:nrow(il.i))
  {
    #k <- 3
    rtmed.i <- il.i$RT[k]
    
    overlap.i <- il.i[-k , ] %>% 
      filter(. , 
             abs(RT - rtmed.i) <= 0.5) %>% 
      nrow()
    overlaptotal.i <- c(overlaptotal.i , overlap.i)
  }
  
  export.i <- data.frame(V1 = il.i$m.z , 
                         V2 = NA , 
                         V3 = 'Chemical formula' , 
                         V4 = '+ H' , 
                         V5 = 1 , 
                         V6 = 'Positive' , 
                         V7 = il.i$RT - 0.25 , 
                         V8 = il.i$RT + 0.25 , 
                         V9 = NA , 
                         V10 = 'NCE' , 
                         V11 = NA , 
                         V12 = paste0(il.i$`xcms feature id` , 
                                      '_mz' , 
                                      round(il.i$m.z , 3) , 
                                      '_rt' , 
                                      round(il.i$RT , 2)))
  colnames(export.i) <- colnames.template
  
  write_csv(x = export.i , 
            file = paste0(inclusionlist.folderpath , 
                          'PG8STD-' , 
                          max(overlaptotal.i) , 
                          '_InclusionList_' , 
                          i , 
                          '.csv') , 
            na = '')
  
}


#3. MSMS networking
# Clean the annotated parent compounds first; make sure no duplicate, maximize the matched compound coverage and their MS/MS availability
MSDialFile<-read.csv("Area_PG1_PRM.csv")
MSDialFileWithPRM<-MSDialFile %>% filter(MS.MS.spectrum != "")
annotatedMSDialFile<-subset(MSDialFileWithPRM,MSDialFileWithPRM$Metabolite.name!="Unknown")
annotatedMSDialFile <- annotatedMSDialFile[!duplicated(annotatedMSDialFile$Metabolite.name), ]# 删除重复值
write.csv(annotatedMSDialFile,"annotatedMSDialFile.csv")
ParentCompoundFeatureList<- data.frame(matrix(nrow = nrow(annotatedMSDialFile),  ncol = 0))
ParentCompoundFeatureList$Alignment.ID <-annotatedMSDialFile$Alignment.ID
ParentCompoundFeatureList$Average.Mz <-annotatedMSDialFile$Average.Mz
ParentCompoundFeatureList$Average.Rt.min. <-annotatedMSDialFile$Average.Rt.min.
ParentCompoundFeatureList$Metabolite.name <-annotatedMSDialFile$Metabolite.name  
write.csv(ParentCompoundFeatureList,"ParentCompoundFeatureList.csv")
# To assign new alignment ID after PRM run
MSDialFileWithPRM <- MSDialFileWithPRM %>% rename(
  Average.Mz.New = Average.Mz,
  Average.Rt.min.New = Average.Rt.min.,
  Alignment.ID.New = Alignment.ID
)

Integrated_RB_MS1_List_merged <- Integrated_RB_MS1_List %>%
  rowwise() %>%
  mutate(
    Average.Mz.New = map_dbl(Average.Mz, ~ MSDialFileWithPRM$Average.Mz.New[which.min(abs(MSDialFileWithPRM$Average.Mz.New - .x))]),
    Average.Rt.min.New = map_dbl(Average.Rt, ~ MSDialFileWithPRM$Average.Rt.min.New[which.min(abs(MSDialFileWithPRM$Average.Rt.min.New - .x))]),
    
    Alignment.ID.New = map2_dbl(Average.Mz, Average.Rt, ~ {
      mz_diff = abs(MSDialFileWithPRM$Average.Mz.New - .x) / .x
      rt_diff = abs(MSDialFileWithPRM$Average.Rt.min.New - .y)
      combined_idx = which.min(mz_diff + rt_diff)
      MSDialFileWithPRM$Alignment.ID.New[combined_idx]
    })
  ) %>%
#    Average.Mz.New = MSDialFileWithPRM$Average.Mz.New[which.min(abs(MSDialFileWithPRM$Average.Mz.New - Average.Mz))],
#    Average.Rt.min.New = MSDialFileWithPRM$Average.Rt.min.New[which.min(abs(MSDialFileWithPRM$Average.Rt.min.New - Average.Rt))],
#    Alignment.ID.New = MSDialFileWithPRM$Alignment.ID.New[combined_idx]
#  ) %>%
#  ungroup() %>%
filter(abs((Average.Mz.New - Average.Mz)/Average.Mz) < 5e-6 & abs(Average.Rt.min.New - Average.Rt) < 0.1)
#Integrated_RB_MS1_List_merged$RB_Metabolite_Name <- paste(Integrated_RB_MS1_List_merged$Compound.Name, Integrated_RB_MS1_List_merged$Reaction_Type, sep = ":")

write.csv(Integrated_RB_MS1_List_merged,"Integrated_RB_MS1_List_merged.csv")
Integrated_RB_MS1_List_merged_duplicates <- Integrated_RB_MS1_List_merged[duplicated(Integrated_RB_MS1_List_merged$Alignment.ID.New) | duplicated(Integrated_RB_MS1_List_merged$Alignment.ID.New, fromLast = TRUE), ]
Integrated_RB_MS1_List_merged_duplicates_unique <- Integrated_RB_MS1_List_merged %>%
  distinct(Alignment.ID.New, .keep_all = TRUE)

MS2FeatureList <- data.frame(matrix(nrow = nrow(Integrated_RB_MS1_List_merged),  ncol = 0)) # This is the final metabolite list!!!
MS2FeatureList$Alignment.ID <-Integrated_RB_MS1_List_merged$Alignment.ID.New
MS2FeatureList$Average.Mz <-Integrated_RB_MS1_List_merged$Average.Mz.New
MS2FeatureList$Average.Rt.min. <-Integrated_RB_MS1_List_merged$Average.Rt.min.New
MS2FeatureList$Metabolite.name <-Integrated_RB_MS1_List_merged$RB_Metabolite_Name

dup_names <- MS2FeatureList$Metabolite.name[duplicated(MS2FeatureList$Metabolite.name)]
for (Metabolite.name in dup_names) {
  indices <- which(MS2FeatureList$Metabolite.name == Metabolite.name)
  num_duplicates <- length(indices)
  if (num_duplicates > 1) {
    MS2FeatureList$Metabolite.name[indices] <- paste(Metabolite.name, seq_len(num_duplicates), sep = "-")
  }
}
MS2FeatureList_Merged<-rbind(ParentCompoundFeatureList,MS2FeatureList)
write.csv(MS2FeatureList_Merged,"MS2FeatureList_Merged.csv")
MS2FeatureList_Merged_unique <- MS2FeatureList_Merged %>%
  distinct(Alignment.ID, .keep_all = TRUE)

#add MS.MS.spectrum to the merged PRM list
names(MSDialFileWithPRM)[names(MSDialFileWithPRM) == "Alignment.ID.New"] <- "Alignment.ID"
MS2FeatureList_PRMmerged <- MS2FeatureList_Merged %>%
  left_join(MSDialFileWithPRM %>% select(Alignment.ID, MS.MS.spectrum), by = "Alignment.ID")
MS2FeatureList_Final<-subset(MS2FeatureList_PRMmerged,MS2FeatureList_PRMmerged$MS.MS.spectrum!=0)

#MSMS list generation
df <- data.frame(name=MS2FeatureList_Final$Metabolite.name,data_column =MS2FeatureList_Final$MS.MS.spectrum,AID=MS2FeatureList_Final$Alignment.ID)
duplicated(df$name)
duplicated(df$AID)
split_data <- function(data, name) {
  split_data <- strsplit(data, " ")[[1]]
  result <- do.call(rbind, lapply(split_data, function(x) {
    values <- strsplit(x, ":")[[1]]
    matrix(c(as.numeric(values[1]), as.numeric(values[2])), nrow = 1, ncol = 2)
  }))
  attr(result, "name") <- name
  return(result)
}

list_of_MSMSSpectrum <- lapply(1:nrow(df), function(i) {
  split_data(df$data_column[i], df$name[i])
})
DistMatrix<- as.data.frame(matrix(nrow = length(df$name), ncol = length(df$name)))
rownames(DistMatrix) <- df$name
colnames(DistMatrix) <- df$name
Maxmz<-max(MS2FeatureList_Final$Average.Mz)

#msentropy_similarity(list_of_MSMSSpectrum[[1]], list_of_MSMSSpectrum[[1]], ms2_tolerance_in_da = 0.02)

for (i in 1:length(df$name)){
  for (j in 1:length(df$name)){
    #i=1
    #j=1
    DistMatrix[i,j]<-msentropy_similarity(list_of_MSMSSpectrum[[i]], list_of_MSMSSpectrum[[j]], 
                                          ms2_tolerance_in_da = 0.02,
                                          clean_spectra = TRUE,
                                          min_mz = 0,
                                          max_mz = 1000,
                                          noise_threshold = 0.01,
                                          max_peak_num = 100,
                                          weighted = TRUE,)
#      SpectrumSimilarity(list_of_MSMSSpectrum[[i]], list_of_MSMSSpectrum[[j]], t = 0.25, b = 10, 
#                   top.label = NULL, bottom.label = NULL, 
#                   xlim = c(100, 1000))
  }
}

#DistMatrix[is.na(DistMatrix)] <- 0
write.csv(DistMatrix,"DistMatrix.csv")
#Correlation Analysis
pheatmap(DistMatrix,show_colnames = FALSE, show_rownames = FALSE)
heatmap(as.matrix(DistMatrix))

DistMatrix_clean<-as.matrix(DistMatrix)
HCplot(DistMatrix_clean, type = "heatmap", 
       cexRow = 0.1, cexCol = 0.1,
       margins = c(6,6))


#generate similarity network file (edge)
points <- rownames(DistMatrix)
DistMatrix_df <- data.frame(Source = character(),
                            Target = character(),
                            Distance = numeric(),
                            stringsAsFactors = FALSE)
n <- length(df$AID)
for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    DistMatrix_df <- rbind(DistMatrix_df, data.frame(Source = points[i],
                                                     Target = points[j],
                                                     Distance = DistMatrix[i, j]))
  }
}

DistMatrix_edge<-subset(DistMatrix_df,DistMatrix_df$Distance!=0)
write.csv(DistMatrix_edge,"DistMatrix_edge.csv")
DistMatrix_edge0<-DistMatrix_edge
DistMatrix_edge0$Source_DR <- gsub(":.*", "", DistMatrix_edge0$Source)
DistMatrix_edge0$Target_DR <- gsub(":.*", "", DistMatrix_edge0$Target)
DistMatrix_edge_DR_0 <- subset(DistMatrix_edge0,DistMatrix_edge0$Source_DR==DistMatrix_edge0$Target_DR)
DistMatrix_edge_DR_0 <- DistMatrix_edge_DR_0 %>%
  mutate(`Distance_Type` = ifelse(!grepl(":", Source) | !grepl(":", Target), "P-M", "M-M"))

#3.1 De-redundancy; to remove the MS/MS from one feature
MS2FeatureList <- MS2FeatureList %>%
  mutate(Metabolite_origin = sub(":.*", "", Metabolite.name))
# Step 1: Merge based on condition 1 (Metabolite.name == Source and Metabolite_origin == Target)
result_1 <- MS2FeatureList %>%
  left_join(DistMatrix_edge, by = c("Metabolite.name" = "Source", "Metabolite_origin" = "Target")) %>%
  select(Metabolite.name, Metabolite_origin, Distance_1 = Distance)

# Step 2: Merge based on condition 2 (Metabolite.name == Target and Metabolite_origin == Source)
result_2 <- MS2FeatureList %>%
  left_join(DistMatrix_edge, by = c("Metabolite.name" = "Target", "Metabolite_origin" = "Source")) %>%
  select(Metabolite.name, Metabolite_origin, Distance_2 = Distance)

# Step 3: Merge the two results and select the non-empty Distance value
result <- result_1 %>%
  left_join(result_2, by = c("Metabolite.name", "Metabolite_origin")) %>%
  mutate(Distance = coalesce(Distance_1, Distance_2)) %>%
  select(-Distance_1, -Distance_2)
result <- result %>%
  mutate(Distance = replace(Distance, is.na(Distance), 0))

# Step 4: Merge table with AID
MS2FeatureList <- MS2FeatureList %>%
  left_join(result %>% select(Metabolite.name, Distance), by = "Metabolite.name")

# Step 5: Delete all duplicate metadata_annotation and keep the one with the highest distance
MS2FeatureList_unique <- MS2FeatureList %>%
  group_by(Alignment.ID) %>%     
  slice_max(order_by = Distance) %>%  
  ungroup()  
MS2FeatureList_unique <- MS2FeatureList_unique %>%
  distinct(Alignment.ID, .keep_all = TRUE) %>% 
  arrange(Alignment.ID)                       
ParentCompoundFeatureList_rename<-ParentCompoundFeatureList %>%
  rename(Parent_MZ = Average.Mz, Parent_RT = Average.Rt.min.)
MS2FeatureList_unique_plusMZRT<-  MS2FeatureList_unique %>%
  left_join(ParentCompoundFeatureList_rename, by = c("Metabolite_origin" = "Metabolite.name"))

MS2FeatureList_unique_plusMZRT$RT_shift<-(MS2FeatureList_unique_plusMZRT$Average.Rt.min.-MS2FeatureList_unique_plusMZRT$Parent_RT)
write.csv(MS2FeatureList_unique_plusMZRT,"MS2FeatureList_unique_plusMZRT.csv")

MS2FeatureList_unique <- MS2FeatureList_unique %>%
  select(1:4)
duplicated(MS2FeatureList_unique$Alignment.ID)
MS2FeatureList_unique_Merged<-rbind(ParentCompoundFeatureList,MS2FeatureList_unique)
MS2FeatureList_unique_Merged <- MS2FeatureList_unique_Merged %>%
  group_by(Alignment.ID) %>%
  filter(!(n() > 1 & grepl(":", Metabolite.name))) %>%
  slice(1) %>% 
  ungroup()
duplicated(MS2FeatureList_unique_Merged$Alignment.ID)
write.csv(MS2FeatureList_unique_Merged,"MS2FeatureList_unique_Merged.csv")

# Rerun matrix
#add MS.MS.spectrum to the merged PRM list
names(MSDialFileWithPRM)[names(MSDialFileWithPRM) == "Alignment.ID.New"] <- "Alignment.ID"
MS2FeatureList_unique_Merged_New <- MS2FeatureList_unique_Merged %>%
  left_join(MSDialFileWithPRM %>% select(Alignment.ID, MS.MS.spectrum), by = "Alignment.ID")
MS2FeatureList_Final_Filtered<-subset(MS2FeatureList_unique_Merged_New,MS2FeatureList_unique_Merged_New$MS.MS.spectrum!=0)
duplicated(MS2FeatureList_Final_Filtered$Alignment.ID)

#MSMS list generation
df1 <- data.frame(name=MS2FeatureList_Final_Filtered$Metabolite.name,data_column =MS2FeatureList_Final_Filtered$MS.MS.spectrum,AID=MS2FeatureList_Final_Filtered$Alignment.ID)
duplicated(df1$name)
duplicated(df1$AID)
split_data <- function(data, name) {
  split_data <- strsplit(data, " ")[[1]]
  result <- do.call(rbind, lapply(split_data, function(x) {
    values <- strsplit(x, ":")[[1]]
    matrix(c(as.numeric(values[1]), as.numeric(values[2])), nrow = 1, ncol = 2)
  }))
  attr(result, "name") <- name
  return(result)
}

list_of_MSMSSpectrum <- lapply(1:nrow(df1), function(i) {
  split_data(df1$data_column[i], df1$name[i])
})
#Cos-similarity Distance Matrix
DistMatrix_New<- as.data.frame(matrix(nrow = length(df1$name), ncol = length(df1$name)))
rownames(DistMatrix_New) <- df1$name
colnames(DistMatrix_New) <- df1$name
Maxmz<-max(MS2FeatureList_Final_Filtered$Average.Mz)

#msentropy_similarity(list_of_MSMSSpectrum[[1]], list_of_MSMSSpectrum[[1]], ms2_tolerance_in_da = 0.02)

for (i in 1:length(df1$name)){
  for (j in 1:length(df1$name)){
    #i=1
    #j=1
    DistMatrix_New[i,j]<-msentropy_similarity(list_of_MSMSSpectrum[[i]], list_of_MSMSSpectrum[[j]], 
                                          ms2_tolerance_in_da = 0.02,
                                          clean_spectra = TRUE,
                                          min_mz = 0,
                                          max_mz = 1000,
                                          noise_threshold = 0.01,
                                          max_peak_num = 100,
                                          weighted = TRUE,)
    #      SpectrumSimilarity(list_of_MSMSSpectrum[[i]], list_of_MSMSSpectrum[[j]], t = 0.25, b = 10, 
    #                   top.label = NULL, bottom.label = NULL, 
    #                   xlim = c(100, 1000))
  }
}

#DistMatrix[is.na(DistMatrix)] <- 0
write.csv(DistMatrix_New,"DistMatrix_New.csv")
#Correlation Analysis
pheatmap(DistMatrix,show_colnames = FALSE, show_rownames = FALSE)
heatmap(as.matrix(DistMatrix_New))

#generate similarity network file (edge)
points <- rownames(DistMatrix_New)
DistMatrix_df_New <- data.frame(Source = character(),
                            Target = character(),
                            Distance = numeric(),
                            stringsAsFactors = FALSE)
n <- length(df1$AID)
for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    DistMatrix_df_New <- rbind(DistMatrix_df_New, data.frame(Source = points[i],
                                                     Target = points[j],
                                                     Distance = DistMatrix_New[i, j]))
  }
}

DistMatrix_edge_New<-subset(DistMatrix_df_New,DistMatrix_df_New$Distance!=0)
write.csv(DistMatrix_edge_New,"DistMatrix_edge_New.csv")

#3.2 Remove Non-protonated adducts/ ISF based on similar RT and spectra
DistMatrix_edge_Final <- DistMatrix_edge_New %>%
  left_join(MS2FeatureList_Final %>% select(Metabolite.name, Average.Rt.min.), 
            by = c("Source" = "Metabolite.name"))
colnames(DistMatrix_edge_Final)[colnames(DistMatrix_edge_Final) == "Average.Rt.min."] <- "Source_RT"
DistMatrix_edge_Final <- DistMatrix_edge_Final %>%
  left_join(MS2FeatureList_Final %>% select(Metabolite.name, Average.Rt.min.), 
            by = c("Target" = "Metabolite.name"))
colnames(DistMatrix_edge_Final)[colnames(DistMatrix_edge_Final) == "Average.Rt.min."] <- "Target_RT"
DistMatrix_edge_Final$rt_defect<-abs(DistMatrix_edge_Final$Source_RT-DistMatrix_edge_Final$Target_RT)
DistMatrix_edge_Final_ISF<-subset(DistMatrix_edge_Final,DistMatrix_edge_Final$Distance>=0.5&DistMatrix_edge_Final$rt_defect<=0.05)
#DistMatrix_edge_Final_nonISF<- DistMatrix_edge_Final %>%
#  filter(!(Distance >= 0.5 & rt_defect <= 0.05))
write.csv(DistMatrix_edge_Final_ISF,"DistMatrix_edge_Final_ISF.csv")

#Dimensionality reduction; assign origin compound for both source and target
DistMatrix_edge1<-DistMatrix_edge_Final
DistMatrix_edge1$Source_DR <- gsub(":.*", "", DistMatrix_edge1$Source)
DistMatrix_edge1$Target_DR <- gsub(":.*", "", DistMatrix_edge1$Target)
DistMatrix_edge_DR <- subset(DistMatrix_edge1,DistMatrix_edge1$Source_DR==DistMatrix_edge1$Target_DR)
DistMatrix_edge_DR <- DistMatrix_edge_DR %>%
  mutate(`Distance_Type` = ifelse(!grepl(":", Source) | !grepl(":", Target), "P-M", "M-M"))
write.csv(DistMatrix_edge_DR,"DistMatrix_edge_DR.csv")
DistMatrix_edge_DR_cytoscape <- DistMatrix_edge_DR
DistMatrix_edge_DR_cytoscape$Distance <-round(DistMatrix_edge_DR_cytoscape$Distance, 3)
DistMatrix_edge_DR_cytoscape <- DistMatrix_edge_DR_cytoscape %>%
  mutate(edge_ID = row_number())
DistMatrix_edge_DR_cytoscape <- DistMatrix_edge_DR_cytoscape %>%
  mutate(HCM = ifelse(Distance >= 0.5 & Distance_Type == "P-M", "HC", "LC"))
write.csv(DistMatrix_edge_DR_cytoscape,"DistMatrix_edge_DR_cytoscape.csv")
#HC Metabolite
highconfidenceMetabolite<-filter(DistMatrix_edge_DR_cytoscape,Distance>=0.5&Distance_Type=="P-M")
highconfidenceMetabolite <- highconfidenceMetabolite %>%
  mutate(
    temp = ifelse(grepl(":", Source), Source, Target), 
    Source = ifelse(grepl(":", Source), Target, Source), 
    Target = temp 
  ) %>%
  select(-temp)

HC_Metabolite_list<- MS2FeatureList %>%
  inner_join(highconfidenceMetabolite, by = c("Metabolite.name" = "Target"))
HC_Metabolite_list <- HC_Metabolite_list %>%
  mutate(Ion_Form = "[M+H]+")
matched_rows <- DistMatrix_edge_Final_ISF %>%
  filter(Source %in% HC_Metabolite_list$Metabolite.name & 
           Target %in% HC_Metabolite_list$Metabolite.name)

HC_Metabolite_list <- HC_Metabolite_list %>%
  mutate(Ion_Form = ifelse(Metabolite.name %in% matched_rows$Source | 
                             Metabolite.name %in% matched_rows$Target, 
                           "Non-protonated adduct/ISF", Ion_Form))

write.csv(highconfidenceMetabolite,"highconfidenceMetabolite.csv")
write.csv(HC_Metabolite_list,"HC_Metabolite_list.csv")

#Correlation Analysis
pheatmap(DistMatrix,show_colnames = FALSE, show_rownames = FALSE)
heatmap(as.matrix(DistMatrix))
