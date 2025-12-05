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

# read in expression data
TmV_25 = fread("source_data/TmV_25_Ano.csv")

## remove diagnosis (Cluster) with less than 3 samples from training set
training = TmV_25 %>% 
  group_by(Cluster) %>% 
  mutate(n = n()) %>% 
  filter(n > 3) %>% 
  ungroup() %>% 
  dplyr::select(!n)

## create the prediction set based on sil score
validation = fread("source_data/Silhouette_per_diagnosis.tsv") %>% 
  inner_join(training) %>% 
  dplyr::select(lab_no, Cluster, Silhouette_Score) %>% 
  group_by(Cluster) %>%
  mutate(n = n()) %>% 
  filter(n > 3) %>% 
  filter(Silhouette_Score <= quantile(Silhouette_Score, 0.25)) %>%
  ungroup() %>% dplyr::select(lab_no)

## Spindle tum NOS is not a real diagnosis, put them all in validation so that a correct diagnosis can be predicted
spindle = TmV_25 %>% 
  filter(Cluster == "Spindle tum NOS") %>% dplyr::select(lab_no)
validation = rbind(validation, spindle) %>% distinct()

## combine to get expression per gene for both sets
validation = inner_join(training, validation) %>% column_to_rownames("lab_no")
training = anti_join(training, validation) %>% column_to_rownames("lab_no")

## transform to factor
training$Cluster = as.factor(training$Cluster)

##### Define training control for cross-validation
train_control <- trainControl(method = "cv", 
                              number = 10,
                              repeats = 5)
# Train the Random Forest model and tune mtry and ntree using 10-fold cross-validation
rf_model <- train(Cluster ~ ., data = training,
                  method = "rf",
                  trControl = train_control,
                  tuneLength = 5)

save(rf_model,file = "rf_model.RData")
print(rf_model)
# load("rf_model.RData")

# use the model to predict other diagnosis 
row = rownames(validation)
p = predict(rf_model, newdata=validation) %>% as.data.frame() %>% rownames_to_column()
p$rowname = row
colnames(p) = c("lab_no", "prediction_rf_model")

pred1=predict(rf_model, newdata=validation, type = "prob") %>% as.data.frame() %>% rownames_to_column(var = "lab_no")
pred2 = inner_join(p, pred1) 

