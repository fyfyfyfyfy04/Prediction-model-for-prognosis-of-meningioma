# 定义一个函数来计算多数投票结果  
majority_vote <- function(row) {  
  # 计算每个标签的数量  
  vote <- table(row)  
  # 返回投票数量最多的标签  
  return(as.numeric(names(vote)[which.max(vote)]))  
}  

#输出最后预测值
pred_vote <- function(row) {  
  if(row[4]==1){
    return(max(row[1:3]))
  } else{
    return(min(row[1:3]))
  }
}  

# 投票模型
vote_model <- function(r_test,t_test,c_test,both_sample) {
    '''
    输入参数：
    r_test:影像组学预测结果
    t_test:检验模型预测结果
    c_test:临床基线模型预测结果
    both_sample:三个数据集的共有样本
    输出：
    final_roc:最终预测roc曲线绘图数据
    final_auc:最终预测auc值
    risk:预测结果（将样本分为高风险及低风险人群）
    sur_data:生存分析绘图数据
    '''
    r_test=r_test %>%  
    filter(ID %in% both_sample) 
    t_test=t_test %>%  
    filter(ID %in% both_sample)
    c_test=c_test %>%  
    filter(ID %in% both_sample)
    if(length(unique(r_test$r_real))<2){
    return(T)
    }

    data=merge(r_test, t_test, by = "ID", all = TRUE)
    data=merge(data, c_test, by = "ID", all = TRUE)

    pred_class = as.data.frame(cbind(data$r_class,data$t_class,data$c_class))

    # 对每一行应用函数  
    pred_class$final_class <- apply(pred_class, 1, majority_vote) 

    pred=as.data.frame(cbind(data$r_pred,data$t_pred,data$c_pred,pred_class$final_class))
    pred$final_pred <- apply(pred, 1, pred_vote) 
    final_roc=roc(as.factor(data$r_real),as.numeric(pred$final_pred))
    final_auc=final_roc[['auc']]
    risk=ifelse(pred_class$final_class==1, "High Risk", "Low Risk")
    sur_data=cbind(data$r_real,risk,as.numeric(data$os_date))

    return(list(final_roc,final_auc,sur_data))
}