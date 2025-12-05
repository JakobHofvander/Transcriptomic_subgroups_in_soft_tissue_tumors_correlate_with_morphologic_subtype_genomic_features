##
install.packages("viridis")
install.packages("caret")
library(viridis)
library(pheatmap)
library(janitor)
library(tidyverse)
library(data.table)
library(car)
library(randomForest)
library(caret)
library(survminer)
library(survival)

sessionInfo()

## Get the data first .... # you can optimize this after the presentation and not now 
### look at the bottom to evaluate ntrees to use
### use tuneRF to estimate mtry = refers to the number of features to consider when making a split at each node of the decision tree.
### clean this up in the next steps!


##### read in the transcriptomic data
TmV_25 = fread("/Users/jhofvand/Dropbox/BMC/RNAseq_Fresh/clustering_all_samples/Annottaions_Figge_2024/tSNEA_704/TmV_25_Ano2.csv")

## remove diagnosis with less than 3 samples from training set
training = TmV_25 %>% 
  group_by(Cluster2) %>% 
  mutate(n = n()) %>% 
  filter(n > 3) %>% 
  ungroup() %>% 
  dplyr::select(!n) %>% 
  #dplyr::select(!Cluster2) %>% 
  dplyr::select(!Genome) %>% 
  select(-contains("-")) # %>% column_to_rownames("UID")

### create the validation set (I mean prediction set)
### same here - remove dig < 4
validation = fread("/Users/jhofvand/Dropbox/BMC/RNAseq_Fresh/clustering_all_samples/Annottaions_Figge_2024/tSNEA_704/Silhouette_per_diagnosis.tsv") %>% 
  group_by(Cluster) %>%
  mutate(n = n()) %>% 
  filter(n > 3) %>% 
  filter(`Silhouette Score` <= quantile(`Silhouette Score`, 0.25)) %>%
  ungroup() %>% dplyr::select(UID)

## spindle cell is not a real diagnosis, put them all in validation so that a diagnosis can be predicted
spindle = fread("/Users/jhofvand/Dropbox/BMC/RNAseq_Fresh/clustering_all_samples/Annottaions_Figge_2024/tSNEA_704/Silhouette_per_diagnosis.tsv") %>% 
  filter(Cluster == "Spindle cell") %>% dplyr::select(UID)
validation = rbind(validation, spindle) %>% distinct()

## combine to get expression per gene for both sets
validation = inner_join(training, validation) %>% column_to_rownames("UID")
training = anti_join(training, validation) %>% column_to_rownames("UID")

## transform to factor
training$Cluster2 = as.factor(training$Cluster2)

##### Define training control for cross-validation
train_control <- trainControl(method = "cv", 
                              number = 10,
                              repeats = 5)
sessionInfo()
# Train the Random Forest model and tune mtry and ntree using 10-fold cross-validation
rf_model <- train(Cluster2 ~ ., data = training,
                  method = "rf",
                  trControl = train_control,
                  tuneLength = 5)

setwd("~/Dropbox/BMC/RNAseq_Fresh/clustering_all_samples/Annottaions_Figge_2024/tSNEA_704/random_forrest/")
save(rf_model,file = "rf_model.RData")
print(rf_model)
# load("rf_model.RData")

#### use the model top predict other diagnosis 
row = rownames(validation)
p = predict(rf_model, newdata=validation) %>% as.data.frame() %>% rownames_to_column()
p$rowname = row
colnames(p) = c("UID", "prediction_rf_model")

pred1=predict(rf_model, newdata=validation, type = "prob") %>% as.data.frame() %>% rownames_to_column(var = "UID")
pred2 = inner_join(p, pred1) 

diff = TmV_25 %>% dplyr::select(UID, Cluster2) %>% 
  inner_join(pred2)

diff= diff %>% 
  dplyr::select(!Cluster2) %>% 
  #column_to_rownames("UID") %>% 
  pivot_longer(!UID:prediction_rf_model) %>% 
  group_by(UID, prediction_rf_model) %>% 
  slice_max(value, n = 1) %>% 
  inner_join(sil2 , by = "UID") %>% arrange(desc(diff)) %>% 
  dplyr::select(UID, value) %>% distinct()

new = inner_join(diff, pred2) 
new %>% clipr::write_clip()


### intersect med sil score
sil = read_xlsx("~/Dropbox/Jakob+Figge/Transcriptomics/Sil_score.xlsx") 
full_join(sil, new, by = "UID") %>% clipr::write_clip()

### only diff diagnosis
read_xlsx("~/Dropbox/Jakob+Figge/Transcriptomics/Sil_score_v2.xlsx") %>% 
  dplyr::select(UID:diff_sil_vs_median) %>% 
  inner_join(new) %>% filter(Cluster != prediction_rf_model) %>% 
  arrange(desc(diff_sil_vs_median)) %>% 
  clipr::write_clip()


##### some qc comparison stuff ..
nd = read_xlsx("~/Dropbox/Jakob+Figge/Transcriptomics/Sil_score_v2.xlsx") %>% 
  dplyr::select(UID:diff_sil_vs_median) %>% 
  inner_join(new) %>% filter(Cluster != prediction_rf_model) %>% 
  dplyr::select(UID) %>% mutate(group = "New_diag")

sd = read_xlsx("~/Dropbox/Jakob+Figge/Transcriptomics/Sil_score_v2.xlsx") %>% 
  dplyr::select(UID:diff_sil_vs_median) %>% 
  inner_join(new) %>% filter(Cluster == prediction_rf_model) %>% 
  dplyr::select(UID) %>% mutate(group = "Same_diag")

a = rbind(nd, sd)

#### plot some qc metrics for the model
## OOB error rate
setwd("~/Dropbox/BMC/RNAseq_Fresh/clustering_all_samples/Annottaions_Figge_2024/tSNEA_704/random_forrest/")
load("rf_model.RData")
print(rf_model)

#### Plot feature importance
importance <- varImp(rf_model, scale = FALSE)
plot(importance)
importance %>% head()
####
#confusion matrix
# ROC cureve?

library(pROC)
# Select a parameter setting
selectedIndices <- rf_model$pred$mtry == 96
# Plot:
plot.roc(rf_model$pred$obs[selectedIndices],
         rf_model$pred$M[selectedIndices])

### cant be shown because when you train the model you should use, savePredictions = T # redo if requested

