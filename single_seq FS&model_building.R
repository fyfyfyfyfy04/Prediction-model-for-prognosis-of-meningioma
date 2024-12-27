library(survival)
library(dplyr)
library(glmnet)
library(ggpubr)
library(ROCR)
library(caret)
library(smotefamily)
library(mlr3verse)
library(mlr3learners)
library(mlr3extralearners)

set.seed(11)
#划分数据集
sample=createDataPartition(data$os, p = 0.7, list = FALSE)
train <- data[sample, ]
test <- data[-sample, ]
#标准化
mean_train = colMeans(train[,2:1052])
sd_train <- sapply(train[,2:1052], sd, na.rm = TRUE)
train[,2:1052] = scale(train[,2:1052])
test[,2:1052]= scale(test[,2:1052],center = mean_train,scale = sd_train)

newx=t(train)
colnames(newx)=newx[1,]
newx = newx[-1,]
x=t(newx)
x=data.frame(x)
x <- sapply(x, as.numeric)
x=data.frame(x)
x <- x[sapply(x, is.numeric) | sapply(x, is.factor)]
smote<-SMOTE(x[, -c(1052,1053,1054)],x[,1052])$data #训练集过采样

y <- smote[, 1052]
x <- smote[,-1052]

x_test <- test[, -c(1053:1055)]
y_test <- test[, 1053:1054]

#单变量cox
#设置p值的阈值
pfilter <- 0.05   
#新建空白数据框
uniresult <- data.frame()
#单因素COX回归分析中p值＜0.05的基因，其分析结果输入到之前新建的空白数据框uniresult中
for(i in colnames(train[,2:1052])){   
  unicox <- coxph(Surv(time = os.date, event = os) ~ train[,i], data = train)  
  unisum<- summary(unicox)   
  pvalue <- round(unisum$coefficients[,5],3) 
  if(pvalue<pfilter){ 
    uniresult <- rbind(uniresult,
                       cbind(feature=i,
                             HR=unisum$coefficients[,2],
                             L95CI=unisum$conf.int[,3],
                             H95CI=unisum$conf.int[,4],
                             pvalue=unisum$coefficients[,5]
                       ))
  }
}   
x = subset(x, select = uniresult$feature)
x_test = subset(x_test, select = uniresult$feature)
x_test = cbind(test$X,x_test)

#构建lasso回归模型
model_lasso <- glmnet(x, unlist(y),family="binomial", nlambda=1000, alpha=1)
print(model_lasso)
plot(model_lasso, xvar="lambda", label=TRUE)

y=as.numeric(y)
cv_fit <- cv.glmnet(x=data.matrix(x), y=data.matrix(y),family="binomial",alpha =1,nfolds = 10)#nlambda=1000,,max_iter=1000
plot(cv_fit)

lasso_min = glmnet(x, unlist(y), lambda=cv_fit$lambda.min)
lasso_1se = glmnet(x, unlist(y), lambda=cv_fit$lambda.1se)

choose_gene_min=rownames(lasso_min$beta)[as.numeric(lasso_min$beta)!=0]
choose_gene_1se=rownames(lasso_1se$beta)[as.numeric(lasso_1se$beta)!=0]
length(choose_gene_min)
length(choose_gene_1se)

newx=t(x_test)
colnames(newx)=newx[1,]
newx = newx[-1,]
x_test=t(newx)
#输出lasso回归模型预测值
lasso.prob.train <- predict(cv_fit, newx = as.matrix(x) , s=c(cv_fit$lambda.min,cv_fit$lambda.1se),type = "response")
lasso.prob.test <- predict(cv_fit, newx = x_test , s=c(cv_fit$lambda.min,cv_fit$lambda.1se),type = "response")
re_train=cbind(y ,lasso.prob.train)
re_test=cbind(y_test ,lasso.prob.test)
head(re_train)

#auc
pred_min_train <- prediction(re_train[,2], re_train[,1])
auc_min_train = performance(pred_min_train,"auc")@y.values[[1]]
auc_min_train
pred_1se_train <- prediction(re_train[,3], re_train[,1])
auc_1se_train = performance(pred_1se_train,"auc")@y.values[[1]]
auc_1se_train


pred_min_test <- prediction(re_test[,3], re_test[,1])
auc_min_test = performance(pred_min_test,"auc")@y.values[[1]]
auc_min_test
pred_1se_test <- prediction(re_test[,4], re_test[,1])
auc_1se_test = performance(pred_1se_test,"auc")@y.values[[1]]
auc_1se_test

#roc曲线绘图
train_pred <- performance(pred_min_train,"tpr","fpr")
test_pred <- performance(pred_min_test,"tpr","fpr")
plot(train_pred,colorize=F, col="#0072B5FF",lwd = 2,main = "ROC Curve for T2f") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
lines(test_pred@x.values[[1]],test_pred@y.values[[1]],col = "#BC3C29FF",lwd = 2 )
text(0.8,0.25, labels = paste0("train auc = ",round(performance(pred_1se_train,"auc")@y.values[[1]],3)))
text(0.79,0.15, labels = paste0("test auc = ",round(performance(pred_1se_test,"auc")@y.values[[1]],3)))
legend("topleft", legend = c("train", "test"), col = c("#0072B5FF", "#BC3C29FF"),   
       lwd = 2)  
#保存筛选特征数据集
sel=data[choose_gene_1se]
sel=cbind(data['X'],sel,data['os'])
write.csv(sel,'D:/file/radiomics/result/update/os/t2f_sel_feature.csv')
