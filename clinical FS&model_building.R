library(survival)
library(rms)
library(regplot)
library(riskRegression)
library(pec)
library(caret)
library(timeROC)
library(tidyr)
library(pROC)
library(ROCR)

clinical_model=function(seed,K) {
    '''
    对临床基线数据进行特征筛选并建模，使用方法单变量及多变量逻辑回归
    输入参数：
    seed:随机种子
    K:SMOTE采样参数k
    输出数据：
    c_result_train:训练集预测结果
    c_result_test:测试集预测结果
    train_auc:训练集AUC
    test_auc:测试集AUC
    roc_train:训练集ROC曲线绘图数据
    roc_test:测试集ROC曲线绘图数据
    '''
    set.seed(seed)
    sample=createDataPartition(clinic_data$os, p = 0.7, list = FALSE)
    #划分数据集
    train <- clinic_data[sample, ]
    test <- clinic_data[-sample, ]

    train_data = train[1:20]
    Y = as.factor(train$os)
    smote<-SMOTE_(train_data,Y,K=K)$data #smote过采样，进行部分代码修改，便于后续取出原数据集样本
    train = smote
    train$class = as.numeric(train$class)

    # 定义一个处理函数  
    round_if_not_integer <- function(x) {  
    if (!is.integer(x)) {  
        return(round(x))  # 四舍五入  
    }  
    return(x)  # 如果是整数，就返回原值  
    }  

    # 应用处理函数到数据框的每一个元素  
    train[1:20] <- as.data.frame(lapply(train[1:20], function(column) {  
    sapply(column, round_if_not_integer)  
    }))  

    data=train
    #设置p值的阈值
    pfilter <-0.05   
    #新建空白数据框
    uniresult <- data.frame()
    #单因素逻辑回归分析中p值＜0.05的基因，其分析结果输入到之前新建的空白数据框uniresult中
    for(i in colnames(train[,1:20])){  
    logic <- glm(class ~ unlist(train[i]), data= train, family = binomial())
    unisum<- summary(logic) 
    if(dim(unisum$coefficients)[1]<2){
        next
    }
    pvalue <- round(unisum$coefficients[2,4],2) 
    if(pvalue<pfilter){ 
        uniresult <- rbind(uniresult,
                            cbind(feature=i,
                                pvalue=unisum$coefficients[2,4]
                            ))
    }
    }   

    #如果单因素逻辑回归分析中p值＜0.05的基因大于0，则进行多因素逻辑回归分析
    if(dim(uniresult)[1]!=0){
    col=c(uniresult[,1],'class')
    data = subset(data, select = col)

    #进行多因素逻辑回归分析
    fea=as.formula(paste("class ~", paste(colnames(data)[1:(ncol(data)-1)], collapse = "+")))
    logic_ <- glm(fea, data= train, family = binomial())

    multisum <- summary(logic_)
    #提取所有基因的多因素COX回归分析结果至multiresult对象中
    feature <- rownames(multisum$coefficients)[-1]
    pvalue <- multisum$coefficients[-1,4]
    multiresult <- data.frame(feature=feature,
                                pvalue=pvalue)
    multiresult <- multiresult[multiresult$pvalue<pfilter,]

    if(dim(multiresult)[1]!=0){
        col=c(multiresult[,1],'class')
        data = subset(data, select = col)
    }

    }
    #使用筛选后的特征构建逻辑回归模型
    fea=as.formula(paste("class ~", paste(colnames(data)[1:(ncol(data)-1)], collapse = "+")))
    logic_final <- glm(fea, data= train, family = binomial())

    #预测值
    train$pre <- predict(logic_final, newdata = train,type ="response")
    test$pre <- predict(logic_final, newdata = test,type ="response")
    if(length(unique(test$os))<2){
    return(T)
    }
    roc_test=roc(test$os,test$pre)
    test_auc=roc_test[["auc"]]

    #分类结果
    class_train=ifelse(train$pre > 0.5, 1, 0)
    class_test=ifelse(test$pre > 0.5, 1, 0)

    ori_sam_index=which(nchar(rownames(smote)) > 3)
    c_result_train=as.data.frame(cbind(train$pre[ori_sam_index],class_train[ori_sam_index]))
    c_result_test=as.data.frame(cbind(test$pre,class_test))
    c_result_train$ID=rownames(train)[ori_sam_index]
    c_result_test$ID=rownames(test)
    colnames(c_result_train)=c('c_pred','c_class','ID')
    colnames(c_result_test)=c('c_pred','c_class','ID')

    roc_train=roc(train$class[ori_sam_index],as.numeric(train$pre[ori_sam_index]))
    train_auc=roc_train[["auc"]]

    return(list(c_result_train,c_result_test,train_auc,test_auc,roc_train,roc_test))
}

test_model=function(seed,K) {
    '''
    对血液及尿液检验数据进行特征筛选并建模，使用方法单变量及多变量逻辑回归
    输入参数：
    seed:随机种子
    K:SMOTE采样参数k
    输出数据：
    c_result_train:训练集预测结果
    c_result_test:测试集预测结果
    train_auc:训练集AUC
    test_auc:测试集AUC
    roc_train:训练集ROC曲线绘图数据
    roc_test:测试集ROC曲线绘图数据
    '''
    set.seed(seed)
    sample=createDataPartition(test_data$os, p = 0.7, list = FALSE)

    train <- test_data[sample, ]
    test <- test_data[-sample, ]


    mean_train = colMeans(train[1:39])
    sd_train <- sapply(train[1:39], sd, na.rm = TRUE)
    train[1:39] = scale(train[1:39])
    test[1:39]= scale(test[1:39],center = mean_train,scale = sd_train)

    train_data = train[1:50]
    Y = as.factor(train$os)
    smote<-SMOTE_(train_data,Y,K=K)$data
    train = smote
    train$class = as.numeric(train$class)

    # 定义一个处理函数  
    round_if_not_integer <- function(x) {  
    if (!is.integer(x)) {  
        return(round(x))  # 四舍五入  
    }  
    return(x)  # 如果是整数，就返回原值  
    }  

    # 应用处理函数到数据框的每一个元素  
    train[40:50] <- as.data.frame(lapply(train[40:50], function(column) {  
    sapply(column, round_if_not_integer)  
    }))  
    data=train

    #设置p值的阈值
    pfilter <-0.05   
    #新建空白数据框
    uniresult <- data.frame()
    #单因素逻辑回归分析中p值＜0.05的基因，其分析结果输入到之前新建的空白数据框uniresult中
    for(i in colnames(train[,1:50])){  
    logic <- glm(class ~ unlist(train[i]), data= train, family = binomial())
    unisum<- summary(logic)
    if(dim(unisum$coefficients)[1]<2){
        next
    }
    pvalue <- round(unisum$coefficients[2,4],2) 
    if(pvalue<pfilter){ 
        uniresult <- rbind(uniresult,
                            cbind(feature=i,
                                #HR=unisum$coefficients[,2],
                                #L95CI=unisum$conf.int[,3],
                                #H95CI=unisum$conf.int[,4],
                                pvalue=unisum$coefficients[2,4]
                            ))
    }
    }   

    if(dim(uniresult)[1]!=0){
    col=c(uniresult[,1],'class')
    data = subset(data, select = col)

    fea=as.formula(paste("class ~", paste(colnames(data)[1:(ncol(data)-1)], collapse = "+")))
    logic_ <- glm(fea, data= train, family = binomial())

    multisum <- summary(logic_)
    #提取所有基因的多因素COX回归分析结果至multiresult对象中
    feature <- rownames(multisum$coefficients)[-1]
    pvalue <- multisum$coefficients[-1,4]
    multiresult <- data.frame(feature=feature,
                                pvalue=pvalue)
    multiresult <- multiresult[multiresult$pvalue<pfilter,]

    if(dim(multiresult)[1]!=0){
        col=c(multiresult[,1],'class')
        data = subset(data, select = col)}
    }



    fea=as.formula(paste("class ~", paste(colnames(data)[1:(ncol(data)-1)], collapse = "+")))
    logic_final <- glm(fea, data= train, family = binomial())

    train$pre <- predict(logic_final, newdata = train,type ="response")
    test$pre <- predict(logic_final, newdata = test,type ="response")
    if(length(unique(test$os))<2){
    return(T)
    }
    roc_test=roc(test$os,test$pre)
    test_auc=roc_test[["auc"]]

    class_train=ifelse(train$pre > 0.5, 1, 0)
    class_test=ifelse(test$pre > 0.5, 1, 0)

    ori_sam_index=which(nchar(rownames(train)) > 3)
    t_result_train=as.data.frame(cbind(train$pre[ori_sam_index],class_train[ori_sam_index]))
    t_result_test=as.data.frame(cbind(test$pre,class_test))
    t_result_train$ID=rownames(train)[ori_sam_index]
    t_result_test$ID=rownames(test)
    tr=test_data[sample, ]
    os_date_train=tr[match(t_result_train$ID,rownames(train_data)),]$os_date
    os_date_test=tr[match(t_result_test$ID,rownames(test_data)),]$os_date
    t_result_train$os_date=os_date_train
    t_result_test$os_date=os_date_test
    colnames(t_result_train)=c('t_pred','t_class','ID','os_date')
    colnames(t_result_test)=c('t_pred','t_class','ID','os_date')

    roc_train=roc(train$class[ori_sam_index],as.numeric(train$pre[ori_sam_index]))
    train_auc=roc_train[["auc"]]

    return(list(t_result_train,t_result_test,train_auc,test_auc,roc_train,roc_test))
}