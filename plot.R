#亚组分析
#ps1+2  seed3 3+4 seed1 k1
#GTR  seed3 STR seed25 k3
#radiotherapy_n  seed6 radiotherapy_y seed3
#Tumor location_c  seed5 k3 Tumor location_nc seed1 k3
#Number of tumors mul样本过少
#Size of tumor <med seed5 k5 >=med seed4 k5
#Ki-67 <med seed2 k5 >=med seed5 k5
#age <65 seed2 k5 >=65 seed5 k5
#gender female seed6 k5 male seed1 k5
p=pipline(6,K=5)
pos_train=p[[9]][[2]][p[[9]][[2]]['Model']=="Radiomics-Clinical",][1:2]
pos_train$subgroup='Without postoperative radiotherapy'
pos_test=p[[10]][[2]][p[[10]][[2]]['Model']=="Radiomics-Clinical",][1:2]
pos_test$subgroup='Without postoperative radiotherapy'
pos_train_auc=p[[9]][[1]][1]
pos_test_auc=p[[10]][[1]][1]

p=pipline(3,K=5)
neg_train=p[[9]][[2]][p[[9]][[2]]['Model']=="Radiomics-Clinical",][1:2]
neg_train$subgroup='Postoperative radiotherapy'
neg_test=p[[10]][[2]][p[[10]][[2]]['Model']=="Radiomics-Clinical",][1:2]
neg_test$subgroup='Postoperative radiotherapy'
neg_train_auc=p[[9]][[1]][1]
neg_test_auc=p[[10]][[1]][1]

train_data=rbind(pos_train,neg_train)
test_data=rbind(pos_test,neg_test)
train_auc=c(pos_train_auc,neg_train_auc)
test_auc=c(pos_test_auc,neg_test_auc)

curve_colors <- c( "Without postoperative radiotherapy" = "#B54764","Postoperative radiotherapy" = "#7895C1")
# 使用 ggplot2 绘制 ROC 曲线
ggplot(train_data, aes(x = x, y = y, color = subgroup)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  labs(title = "Postoperative radiotherapy subgroup in training cohort", 
       x = "False Positive Rate", 
       y = "True Positive Rate") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),# 将标题居中
        legend.position = c(0.65, 0.25),           # 图例放置在内部（x=0.8，y=0.8）位置  
        legend.title = element_blank()) + 
  scale_color_manual(name = "Subgroup",   
                     values = curve_colors,  
                     labels = c(paste0("Without postoperative radiotherapy (AUC=", train_auc[1],")"), paste0("Postoperative radiotherapy (AUC=", train_auc[2],")")))

# 使用 ggplot2 绘制 ROC 曲线
ggplot(test_data, aes(x = x, y = y, color = subgroup)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  labs(title = "Postoperative radiotherapy subgroup in testing cohort", 
       x = "False Positive Rate", 
       y = "True Positive Rate") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),# 将标题居中
        legend.position = c(0.65, 0.25),           # 图例放置在内部（x=0.8，y=0.8）位置  
        legend.title = element_blank()) + 
  scale_color_manual(name = "Subgroup",   
                     values = curve_colors,  
                     labels = c(paste0("Without postoperative radiotherapy (AUC=", test_auc[1],")"), paste0("Postoperative radiotherapy (AUC=", test_auc[2],")")))




#Radiomics-Clinical模型，Radiomics模型，Clinical baseline模型,Inspection模型的roc曲线对比图
curve_colors <- c( "Radiomics-Clinical" = "#B54764",
                   "Radiomics" = "#EDC3A5", "Clinical baseline" = "#ABD1BC",
                   "Inspection" = "#A8CBDF")#"Clinical-Radiomics" = "#992224","DIABLO_t1" = "#E3625D","DIABLO_t1c" = "#EF8B67", "DIABLO_t2" = "#F0C284", "DIABLO_t2f" = "#FCB6A5",
p=pipline(651)
train_p=p[[9]][[2]]
test_p=p[[10]][[2]]
auc_train=p[[9]][[1]]
auc_test=p[[10]][[1]]

# 使用 ggplot2 绘制 ROC 曲线
ggplot(train_p, aes(x = x, y = y, color = Model)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  labs(title = "ROC Curves for Train Data by Models", 
       x = "False Positive Rate", 
       y = "True Positive Rate") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +  # 将标题居中
  scale_color_manual(values = curve_colors) +  # 手动设置曲线的颜色
  annotate("text", x = 0.5, y = 0.3, label = paste("Radiomics-Clinical:", auc_train[1]), color = curve_colors["Radiomics-Clinical"], hjust = 0)  +
  annotate("text", x = 0.5, y = 0.25, label = paste("Radiomics:", auc_train[2]), color = curve_colors["Radiomics"], hjust = 0) +
  annotate("text", x = 0.5, y = 0.2, label = paste("Inspection:", auc_train[3]), color = curve_colors["Inspection"], hjust = 0) +
  annotate("text", x = 0.5, y = 0.15, label = paste("Clinical baseline:", auc_train[4]), color = curve_colors["Clinical baseline"], hjust = 0) 

ggplot(test_p, aes(x = x, y = y, color = Model)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  labs(title = "ROC Curves for Test Data by Models", 
       x = "False Positive Rate", 
       y = "True Positive Rate") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +  # 将标题居中
  scale_color_manual(values = curve_colors) +  # 手动设置曲线的颜色
  annotate("text", x = 0.5, y = 0.3, label = paste("Radiomics-Clinical:", auc_test[1]), color = curve_colors["Radiomics-Clinical"], hjust = 0)  +
  annotate("text", x = 0.5, y = 0.25, label = paste("Radiomics:", auc_test[2]), color = curve_colors["Radiomics"], hjust = 0) +
  annotate("text", x = 0.5, y = 0.2, label = paste("Inspection:", auc_test[3]), color = curve_colors["Inspection"], hjust = 0) +
  annotate("text", x = 0.5, y = 0.15, label = paste("Clinical baseline:", auc_test[4]), color = curve_colors["Clinical baseline"], hjust = 0) 




#生存分析曲线
sur_train=as.data.frame(p[[11]])
sur_test=as.data.frame(p[[12]])

fit <- survfit(Surv(as.numeric(sur_train$V3),V1==1) ~ risk, data = sur_train)
library(survminer)

ggsurvplot(fit, 
           data = sur_train, 
           pval = TRUE, 
           conf.int = TRUE, 
           #risk.table = TRUE, 
           legend.labs = c("High Risk", "Low Risk"), 
           title = "Survival Analysis by Risk Group of Train Data")

fit <- survfit(Surv(as.numeric(sur_test$V3),V1==1) ~ risk, data = sur_test)
ggsurvplot(fit, 
           data = sur_test, 
           pval = TRUE, 
           conf.int = TRUE, 
           #risk.table = TRUE, 
           legend.labs = c("High Risk", "Low Risk"), 
           title = "Survival Analysis by Risk Group of Test Data")
