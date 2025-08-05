library(tidyverse)
library(keras)
library(Retip)
library(rcdk)
library(dplyr)
#>Starts parallel computing
prep.wizard()

# import excel file for training and testing data
RP2 <- readxl::read_excel("PG_DB.xlsx", sheet = "Training", col_types = c("text", 
                                                                                 "text", "text", "numeric"))
# import excel file for external validation set
RP_ext <- readxl::read_excel("PG_DB.xlsx", sheet = "Training", col_types = c("text", 
                                                                                  "text", "text", "numeric"))
RP2 <- RP2 %>%
  mutate(SMILES = str_trim(SMILES, side = "right"))
#> or use HILIC database included in Retip
#HILIC <- HILIC

#> Calculate Chemical Descriptors from CDK
descs <- getCD(RP2)

#> Clean dataset from NA and low variance value
db_rt <- proc.data(descs)


#> Plot chem space, the first value is your library with Chemical descriptors calculated, 
#> the second one is your target that can be a included database 
#> or your favourite one that you have uploaded inside Retip

chem.space(db_rt,t="CHEBI")


#> Split in training and testing using caret::createDataPartition
set.seed(101)
inTraining <- caret::createDataPartition(db_rt$XLogP, p = .8, list = FALSE)
training <- db_rt[ inTraining,]
testing  <- db_rt[-inTraining,]

#> Train Model

xgb <- fit.xgboost(training)

rf  <- fit.rf(training)

brnn <- fit.brnn(training)

#keras <- fit.keras(training,testing)

#lightgbm <- fit.lightgbm(training,testing)

#> first you have to put the testing daframe and then the name of the models you have computed
stat <- get.score(testing,xgb,rf,brnn)

#> first you have to put the testing daframe and then the name of the models you want to visualize. 
#> Last value is the title of the graphics.
p.model(testing, m=xgb,title = "XGBoost-Pesticide Database")


#> import dataset

pathogen_box <- readxl::read_excel("PG_DB AITest.xlsx", col_types = c("text", 
                                                                               "text", "text"))
pathogen_box <- pathogen_box %>%
  mutate(SMILES = str_trim(SMILES, side = "right"))
#> compute Chemical descriptors
pathogen_box_desc <- getCD(pathogen_box)

#> perform the RT spell
pathogen_box_pred <- RT.spell(training,pathogen_box_desc,model=xgb)
write.csv(pathogen_box_pred,"pathogen_box_pred1_PDB.csv")
