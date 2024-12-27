pipline <- function(seed,K) {
    '''
    将影像组学建模，临床数据建模以及最终投票模型构建打包为一个函数
    输入参数：
    seed:随机种子
    K:smote参数k值
    返回参数：
    result:包含各模型的auc值，各模型roc绘图数据，生存分析绘图数据
    '''
    set.seed(seed)
    radiomics_pred=radiomics_model(seed,K)
    if(isTRUE(radiomics_pred)){
    return(T)
    }
    r_train=radiomics_pred[[1]]
    r_test=radiomics_pred[[2]]
    r_train_auc=radiomics_pred[[3]]
    r_test_auc=radiomics_pred[[4]]
    r_roc_train=radiomics_pred[[5]]
    r_roc_test=radiomics_pred[[6]]
    test_pred=test_model(seed,K)
    if(isTRUE(test_pred)){
    return(T)
    }
    t_train=test_pred[[1]]
    t_test=test_pred[[2]]
    t_train_auc=test_pred[[3]]
    t_test_auc=test_pred[[4]]
    t_roc_train=test_pred[[5]]
    t_roc_test=test_pred[[6]]
    clinical_pred=clinical_model(seed,K)
    if(isTRUE(clinical_pred)){
    return(T)
    }
    c_train=clinical_pred[[1]]
    c_test=clinical_pred[[2]]
    c_train_auc=clinical_pred[[3]]
    c_test_auc=clinical_pred[[4]]
    c_roc_train=clinical_pred[[5]]
    c_roc_test=clinical_pred[[6]]

    both_sample_train <- Reduce(intersect, list(r_train$ID,t_train$ID, c_train$ID))
    both_sample_test <- Reduce(intersect, list(r_test$ID,t_test$ID, c_test$ID))
    if(length(both_sample_train)==0 | length(both_sample_test)==0){
    return(T)
    }

    train_vote=vote_model(r_train,t_train,c_train,both_sample_train)
    test_vote=vote_model(r_test,t_test,c_test,both_sample_test)
    if(isTRUE(train_vote) | isTRUE(test_vote)){
    return(T)
    }
    train_auc=train_vote[[2]]
    test_auc=test_vote[[2]]
    train_roc=train_vote[[1]]
    test_roc=test_vote[[1]]
    sur_train=train_vote[[3]]
    sur_test=test_vote[[3]]

    train_roc_data=roc_plot_data(train_roc,r_roc_train,t_roc_train,c_roc_train)
    test_roc_data=roc_plot_data(test_roc,r_roc_test,t_roc_test,c_roc_test)

    result=list(r_train_auc,r_test_auc,t_train_auc,t_test_auc,c_train_auc,c_test_auc,train_auc,test_auc,train_roc_data,test_roc_data,sur_train,sur_test)

    return(result)
}
