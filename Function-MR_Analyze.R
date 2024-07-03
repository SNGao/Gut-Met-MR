
# MR分析及其可视化（老函数，在新一轮分析中已经被淘汰）
## 针对已经匹配过的Data开展基本的MR分析
## 可自定义输入想要的method类型
MR_ana_visual <- function(mr_dat = NULL,
                          save_path = NULL,
                          output_name = NULL,
                          methods = NULL){
  file = paste(c(paste(c(save_path,output_name), collapse = '/'),'.pdf'), collapse = '')
  pdf(file)
  generate_odds_ratios(mr_res = mr(mr_dat, method_list = methods))
  print(mr_scatter_plot(mr_results = mr(mr_dat,method_list = methods), mr_dat))
  print(mr_heterogeneity(mr_dat))
  print(mr_pleiotropy_test(mr_dat))
  print(mr_funnel_plot(singlesnp_results = mr_singlesnp(mr_dat)))
  print(mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(mr_dat)))
  dev.off()
}


# 检查是否存在异常值（MR-PRESSO 逐渐被淘汰）
outlier_check <- function(data = NULL){
  result <- mr_presso(
    BetaOutcome ="beta.outcome", 
    BetaExposure = "beta.exposure", 
    SdOutcome ="se.outcome", 
    SdExposure = "se.exposure", 
    OUTLIERtest = TRUE,DISTORTIONtest = TRUE, NbDistribution = 500,
    data = data,
    SignifThreshold = 0.05)
  print(result)
  data_2 <- data[-c(which(result$`MR-PRESSO results`$`Outlier Test`$Pvalue<0.05)),]
  return(data_2)
}

# 循环MR-PRESSO分析，剔除异常SNP
## 需要输入初步的MR-PRESSO分析结果（在不逐一剔除SNP时，不需要本分析）
## 当异常值存在时，则依次开展MR-PRESSO循环检验，直到结果符合
## 如果single P值小于<0.05则全部剔除，反之则
PRESSO_check <- function(data = NULL,
                         result = NULL){
  while (result$`MR-PRESSO results`$`Global Test`$Pvalue<0.05){
    # 提取P_value最小的序数
    Outlier_list_p = result$`MR-PRESSO results`$`Outlier Test`$Pvalue
    Outlier_list_p = as.numeric(gsub('<', '', Outlier_list_p))
    min_count <- sum(Outlier_list_p == min(Outlier_list_p)) #计算最小值的个数
    
    data = data[-c(which.min(Outlier_list_p)), ]
    
    result <- mr_presso(
      BetaOutcome ="beta.outcome", 
      BetaExposure = "beta.exposure", 
      SdOutcome ="se.outcome", 
      SdExposure = "se.exposure", 
      OUTLIERtest = TRUE,DISTORTIONtest = FALSE, NbDistribution = 500,
      data = data,
      SignifThreshold = 0.05)
  }
  return(data)
}


# 生成Harmonized文件
## 期间对多个SNP开展PRESSO剔除检验
PRESSO_Check_loop <- function(
                              folder_path = NULL,
                              out_2 = NULL)
{
  result_total = data.frame()
  file_names <- list.files(folder_path)
  #print(file_names)
  
  for (file_name in file_names){ #
    print(file_name)
    file_path = paste(c(folder_path, file_name), collapse = '')
    dat = read.csv(file_path);  dat$X = NULL
    
    mr_dat = Two_MR(dat, out_2)
    if (dim(mr_dat)[1] == 0){ next }
    
    if(dim(mr_dat)[1]>=4){
      result <- mr_presso(
        BetaOutcome ="beta.outcome", 
        BetaExposure = "beta.exposure", 
        SdOutcome ="se.outcome", 
        SdExposure = "se.exposure", 
        OUTLIERtest = TRUE,DISTORTIONtest = FALSE, NbDistribution = 500,
        data = mr_dat,
        SignifThreshold = 0.1)
      
      file_name_new = paste(c(folder_path, paste(c('PRESSO_',file_name), collapse = '')), collapse = '')
      print(result$`MR-PRESSO results`$`Global Test`$Pvalue)
      
      result_temp = data.frame(file_name, result$`MR-PRESSO results`$`Global Test`$Pvalue)
      colnames(result_temp) = c('exposure', 'PRESSO_p_val')
      result_total = rbind(result_total, result_temp)
      
      save_path = paste(c(folder_path, 'Total_PRESSO_results.csv'), collapse = '')
      write.csv(result_total, save_path)
      
      if(result$`MR-PRESSO results`$`Global Test`$Pvalue < 0.05){
        new_data = PRESSO_check(mr_dat, result)
        write.csv(new_data, file_name_new)
      }
      if(result$`MR-PRESSO results`$`Global Test`$Pvalue >= 0.05){
        write.csv(mr_dat, file_name_new)
      }
      
    } else {
      file_name_new = paste(c(folder_path, paste(c('PRESSO_NULL_',file_name), collapse = '')), collapse = '')
      write.csv(mr_dat, file_name_new)
    }
  }
}
