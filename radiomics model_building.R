library(mixOmics)
library(BiocParallel)
library(caret)
library(base)
library(dplyr)
library(tidyr) 
library(pROC)


data_scale=function(train,test,end_col){
    mean_train = colMeans(train[,2:end_col])
    sd_train <- sapply(train[,2:end_col], sd, na.rm = TRUE)
    train[,2:end_col] = scale(train[,2:end_col])
    test[,2:end_col]= scale(test[,2:end_col],center = mean_train,scale = sd_train)
    return(list(train,test))
}

radiomics_model <- function(seed,K) {
    '''
    利用mixOmics函数包中的DIABLO方法融合四个序列的影像数据构建影像组学模型
    输入参数:
    seed:随机种子
    K:SMOTE方法中K值
    输出参数:
    r_result_train:训练集预测结果
    r_result_test:测试集预测结果
    auc_train:训练集AUC
    auc_test:测试集AUC
    roc_train:训练集ROC绘图数据
    roc_test:测试集ROC绘图数据
    '''
    set.seed(seed)
    sample=createDataPartition(t1_data$os, p = 0.7, list = FALSE)

    #分别对四个序列的影像数据进行数据集划分以及标准化
    t1_train <- t1_data[sample, ]
    t1_test <- t1_data[-sample, ]
    mean_train = colMeans(t1_train[,2:16])
    sd_train <- sapply(t1_train[,2:16], sd, na.rm = TRUE)
    t1_train[,2:16] = scale(t1_train[,2:16])
    t1_test[,2:16]= scale(t1_test[,2:16],center = mean_train,scale = sd_train)

    t1c_train <- t1c_data[sample, ]
    t1c_test <- t1c_data[-sample, ]
    mean_train = colMeans(t1c_train[,2:20])
    sd_train <- sapply(t1c_train[,2:20], sd, na.rm = TRUE)
    t1c_train[,2:20] = scale(t1c_train[,2:20])
    t1c_test[,2:20]= scale(t1c_test[,2:20],center = mean_train,scale = sd_train)

    t2_train <- t2_data[sample, ]
    t2_test <- t2_data[-sample, ]
    mean_train = colMeans(t2_train[,2:44])
    sd_train <- sapply(t2_train[,2:44], sd, na.rm = TRUE)
    t2_train[,2:44] = scale(t2_train[,2:44])
    t2_test[,2:44]= scale(t2_test[,2:44],center = mean_train,scale = sd_train)

    t2f_train <- t2f_data[sample, ]
    t2f_test <- t2f_data[-sample, ]
    mean_train = colMeans(t2f_train[,2:26])
    sd_train <- sapply(t2f_train[,2:26], sd, na.rm = TRUE)
    t2f_train[,2:26] = scale(t2f_train[,2:26])
    t2f_test[,2:26]= scale(t2f_test[,2:26],center = mean_train,scale = sd_train)

    Y=t1_train$os
    Y_test=t1_test$os
    if(length(unique(Y_test))<2){
    return(T)
    }

    t1_train=a(t1_train,15)
    t1c_train=a(t1c_train,19)
    t2_train=a(t2_train,43)
    t2f_train=a(t2f_train,25)
    train=cbind(t1_train,t1c_train,t2_train,t2f_train)
    smote<-SMOTE_(train,Y,K=K)$data #smote过采样，进行部分代码修改，便于后续取出原数据集样本
    t1_train=smote[,1:15]
    t1c_train=smote[,16:34]
    t2_train=smote[,35:77]
    t2f_train=smote[,78:102]
    Y=smote[,103]

    t1_test=a(t1_test,15)
    t1c_test=a(t1c_test,19)
    t2_test=a(t2_test,43)
    t2f_test=a(t2f_test,25)

    data=list(t1=t1_train,t1c=t1c_train,t2=t2_train,t2f=t2f_train)
    test=list(t1=t1_test,t1c=t1c_test,t2=t2_test,t2f=t2f_test)

    # for square matrix filled with 0.1s
    design = matrix(0.3, ncol = length(data), nrow = length(data), 
                    dimnames = list(names(data), names(data)))
    diag(design) = 0 # set diagonal to 0s
    list.keepX = list (t1 = c(15,15,5,10), 
                        t1c = c(15,15,5,5),
                        t2 = c(25,5,40,5),
                        t2f = c(5,15,25,5))
    final.diablo.model = block.splsda(X = data, Y = Y, ncomp = 5, # set the optimised DIABLO model
                                    keepX = list.keepX, design = design)
    predict.diablo = predict(final.diablo.model, newdata = test)
    pred_test=predict.diablo[["AveragedPredict"]][,,5][,2]#WeightedPredict
    class_test=as.data.frame(predict.diablo[["AveragedPredict.class"]][["max.dist"]][,5])[,1]
    roc_test=roc(Y_test,pred_test)
    auc_test=roc_test[['auc']]


    predict.diablo = predict(final.diablo.model, newdata = data)
    pred_train=predict.diablo[["AveragedPredict"]][,,5][,2]
    class_train=as.data.frame(predict.diablo[["AveragedPredict.class"]][["max.dist"]][,5])[,1]

    ori_sam_index=which(nchar(rownames(smote)) > 3)
    r_result_train=as.data.frame(cbind(pred_train[ori_sam_index],class_train[ori_sam_index],Y[ori_sam_index]))
    r_result_test=as.data.frame(cbind(pred_test,class_test,Y_test))
    r_result_train$ID=rownames(r_result_train)
    r_result_test$ID=rownames(r_result_test)
    colnames(r_result_train)=c('r_pred','r_class','r_real','ID')
    colnames(r_result_test)=c('r_pred','r_class','r_real','ID')
    roc_train=roc(r_result_train$r_real,as.numeric(r_result_train$r_pred))
    auc_train=roc_train[['auc']]

    return(list(r_result_train,r_result_test,auc_train,auc_test,roc_train,roc_test))
    }