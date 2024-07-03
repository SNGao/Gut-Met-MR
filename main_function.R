### 运行所需的自定义函数
## 从本地读取Exposure数据
## 修改读入table的变量名字
## 对exp数据进行预处理
## 开展二样本MR初步分析
## MR分析及其可视化
## 检查是否存在异常值
## 在线提取缺失的eaf值


# 从本地读取Exposure数据(UKB)
dat_read_ukb <- function(filepath = NULL){
  exp <- read_delim(file=filepath, delim="\t")
  temp = strsplit(exp$variant,':')
  temp2 = data.frame(do.call(rbind, temp))
  colnames(temp2) = c('chr','SNP','other_allele','effect_allele')
  temp2$SNP = paste('rs', temp2$SNP,sep='')
  
  exp = cbind(temp2, exp)
  print(colnames(exp))
  print('我需要的colname label')
  print(c('SNP','effect_allele','other_allele','freq','p','beta','SE'))
  return(data.frame(exp))
}

# 从本地读取Exposure数据
dat_read <- function(filepath = NULL){
  exp <- fread(filepath, header = T)
  exp$V1 = NULL
  print(colnames(exp))
  print('我需要的colname label')
  print(c('SNP','effect_allele','other_allele','freq','p','beta','SE'))
  return(exp)
}

# 修改读入table的变量名字
rename_dat <- function(data){
  list = colnames(data)
  print(list)
  if ('OR' %in% list){
    data$beta = log(data$OR)
  }
  
  if ('Zscore' %in% list){
  }
  
  data$direction = NULL
  for (i in c(1:dim(data)[2])){
    cor_name = list[i]
    if(cor_name %in% c('other_allele','ALLELE0','ref','A2','REF',
                       'Allele2','ref.allele','OTHER_ALLELE')){ colnames(data)[i]='other_allele' }
    else if(cor_name %in% c('effect_allele','ALLELE1','alt','A1','ALT',
                            'Allele1','eff.allele','EFFECT_ALLELE')){ colnames(data)[i]='effect_allele' }
    else if(cor_name %in% c('A1FREQ','freq','af_alt','MAF','FRQ_A_40463','EAF',
                            'effect_allele_frequency','Freq1','Freq','EAF_HRC',
                            'eaf.exposure')){ colnames(data)[i]='freq' }
    else if(cor_name %in% c('P_LINREG','p','pval','P','p_value','P.Neff','PVAL',
                            'P.weightedSumZ','P-value','Pval','pval.exposure')){ colnames(data)[i]='p' }
    else if(cor_name %in% c('BETA','beta','Effect','LogOR','stdBeta','ES',
                            'Beta','est','beta.exposure')){ colnames(data)[i]='beta' }
    else if(cor_name %in% c('SE','sebeta','standard_error','StdErr','se', 
                            'StdErrLogOR','Markername','se.exposure')){ colnames(data)[i]='SE' }
    else if(cor_name %in% c('SNP','rsids','rsid','snp','MarkerName','ID',
                            'rsID','variant_id')){ colnames(data)[i]='SNP' }
  }
  print(colnames(data))
  data$effect_allele = toupper(data$effect_allele)
  data$other_allele = toupper(data$other_allele)
  return(data.frame(data))
}

# 对exp数据进行预处理(自动判断是否有满足p_value的SNP，并自动计算eaf值)
exp_P_LD <- function(data = NULL,
                     p_value = 5e-8,
                     clump_kb = 10000,
                     clump_r2 = 0.01,
                     file_p_eaf = NULL){
  
  data_subset <- as.data.frame(subset(data, p<p_value))
  if (dim(data_subset)[1] == 0){
    print('是时候换下一个数据集，或放宽SNP筛选标准了')
    return('NULL')
  }
  if (!('freq' %in% colnames(data_subset))){
    data_subset = snp_add_eaf(data_subset)
  }
  colnames(data_subset)[colnames(data_subset) == 'eaf'] = 'freq'
  write.csv(data_subset, file="temp.csv")
  
  ## 生成符合格式的数据
  exp_dat <- read_exposure_data(filename = "temp.csv",
                                sep = ",",snp_col = "SNP",
                                eaf_col = 'freq',
                                pval_col = "p",
                                beta_col = "beta",se_col = "SE",
                                effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele",clump = FALSE)

  ## 去除连锁不平衡
  exp_dat_clumped <- clump_data(exp_dat,
                                clump_kb = clump_kb,
                                clump_r2 = clump_r2)
  # write.csv(exp_dat_clumped, file_p_eaf)
}

# 开展二样本MR初步分析
Two_MR <- function(exp_dat = dat,
                   out_dat = out_2){
  exp_out_dat <- merge(exp_dat, out_dat,by.x="SNP", by.y='SNP')
  if (dim(exp_out_dat)[1] > 0){
  write.csv(exp_out_dat, 'temp2.csv')
  out_dat2<-read_outcome_data(snps = exp_dat$SNP,
                              filename = "temp2.csv",sep = ",",
                              snp_col = "SNP",beta_col = "beta",
                              eaf='freq',
                              se_col = "SE",effect_allele_col = "effect_allele",
                              other_allele_col = "other_allele",pval_col = "p")
  
  mr_dat<-harmonise_data(exposure_dat = exp_dat,
                         outcome_dat = out_dat2)
  return(mr_dat)
  } else {
    print('无匹配的SNP')
    return(exp_out_dat)
    }
}


# MR分析及其可视化
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

## 检查是否存在异常值
outlier_check <- function(data = NULL){
  result <- mr_presso(
    BetaOutcome ="beta.outcome", 
    BetaExposure = "beta.exposure", 
    SdOutcome ="se.outcome", 
    SdExposure = "se.exposure", 
    OUTLIERtest = TRUE,DISTORTIONtest = TRUE, NbDistribution = 1000,
    data = data,
    SignifThreshold = 0.05)
  print(result)
  data_2 <- data[-c(which(result$`MR-PRESSO results`$`Outlier Test`$Pvalue<0.05)),]
  return(data_2)
}

## 当异常值存在时，则依次开展MR-PRESSO循环检验，直到结果符合
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
      OUTLIERtest = TRUE,DISTORTIONtest = FALSE, NbDistribution = 1000,
      data = data,
      SignifThreshold = 0.05)
  }
  return(data)
}


## 在线提取缺失的eaf值(可能存在个别差异，但整体差异不大)
snp_add_eaf <- function(dat, build = "37", pop = "EUR")
{
  stopifnot(build %in% c("37","38"))
  stopifnot("SNP" %in% names(dat))
  
  # Create and get a url
  server <- ifelse(build == "37","http://grch37.rest.ensembl.org","http://rest.ensembl.org")
  pop <- paste0("1000GENOMES:phase_3:",pop)
  
  snp_reverse_base <- function(x)
  {
    x <- str_to_upper(x)
    #stopifnot(x %in% c("A","T","C","G"))
    switch(x,"A"="T","T"="A","C"="G","G"="C")
  }
  
  res_tab <- lapply(1:nrow(dat), function(i)
  {
    print(paste0("seaching for No.", i, " SNP"))
    dat_i <- dat[i,]
    
    ext <- paste0("/variation/Homo_sapiens/",dat_i$SNP, "?content-type=application/json;pops=1")
    url <- paste(server, ext, sep = "")
    res <- httr::GET(url)
    
    # Converts http errors to R errors or warnings
    httr::stop_for_status(res)
    
    # Convert R objects from JSON
    res <- httr::content(res)
    res_pop <- jsonlite::fromJSON(jsonlite::toJSON(res))$populations
    
    # Filter query results based on population set
    res_pop <- try(res_pop[res_pop$population == pop,])
    if("try-error" %in% class(res_pop))
    {
      print(paste0("There is not information for population ",pop))
      queried_effect_allele <- "NR"
      queried_other_allele <- "NR"
      queried_eaf <- -1
    }
    else
    {
      if(nrow(res_pop)==0)
      {
        print(paste0("There is not information for population ",pop))
        queried_effect_allele <- "NR"
        queried_other_allele <- "NR"
        queried_eaf <- -1
      }
      else
      {
        queried_effect_allele <- res_pop[1,"allele"][[1]]
        queried_other_allele <- res_pop[2,"allele"][[1]]
        queried_eaf <- res_pop[1,"frequency"][[1]]    
      }
    }
    
    effect_allele <- ifelse("effect_allele.exposure" %in% names(dat),
                            dat_i$effect_allele.exposure,
                            dat_i$effect_allele)
    
    other_allele <- ifelse("effect_allele.exposure" %in% names(dat),
                           dat_i$other_allele.exposure,
                           dat_i$other_allele)
    
    if("effect_allele.exposure" %in% names(dat))
    {
      name_output <- unique(c(names(dat), "eaf.exposure","reliability.exposure"))
    }
    else
    {
      name_output <- unique(c(names(dat), "eaf","reliability.exposure"))
    }
    
    len_effect_allele <- nchar(effect_allele)
    len_other_allele <- nchar(other_allele)
    
    if(len_effect_allele==1&len_other_allele==1)
    {
      if((queried_effect_allele==effect_allele & queried_other_allele==other_allele)|
         (queried_effect_allele==other_allele & queried_other_allele==effect_allele))
      {
        dat_i$eaf.exposure <- ifelse(effect_allele == queried_effect_allele,
                                     queried_eaf,
                                     1-queried_eaf)
        dat_i$eaf <- dat_i$eaf.exposure 
        dat_i$reliability.exposure <- "high"
      }
      else
      {
        r_queried_effect_allele <- snp_reverse_base(queried_effect_allele)
        r_queried_other_allele <- snp_reverse_base(queried_other_allele)
        if((r_queried_effect_allele==effect_allele & r_queried_other_allele==other_allele)|
           (r_queried_effect_allele==other_allele & r_queried_other_allele==effect_allele))
        {
          dat_i$eaf.exposure <- ifelse(effect_allele == r_queried_effect_allele,
                                       queried_eaf,
                                       1-queried_eaf)
          dat_i$eaf <- dat_i$eaf.exposure 
          dat_i$reliability.exposure <- "high"
        }
        else
        {
          dat_i$eaf.exposure <- ifelse(effect_allele == queried_effect_allele,
                                       queried_eaf,
                                       1-queried_eaf)
          dat_i$eaf <- dat_i$eaf.exposure 
          dat_i$reliability.exposure <- "low"
        }
      }
    }
    
    else
    {
      # To identify the potential DEL/ INS
      short_allele <- ifelse(len_effect_allele==1,
                             effect_allele,
                             other_allele)
      short_allele_eaf <- ifelse(short_allele == queried_effect_allele, 
                                 queried_eaf, 
                                 1-queried_eaf)
      dat_i$eaf.exposure <- ifelse(effect_allele == short_allele,
                                   short_allele_eaf,
                                   1-short_allele_eaf)
      dat_i$eaf <- dat_i$eaf.exposure 
      dat_i$reliability.exposure <- "low"
    }
    
    dat_i[name_output]
  })
  
  return(do.call(rbind, res_tab))
}

## 循环计算F+MAF值
cal_F_MAF_check <- function(folder_path = NULL){
  
  file_names = list.files(folder_path)
  Check_names = list()
  Sum = 0
  
  for (file_name in file_names){
    file_path = paste(c(folder_path, file_name), collapse = '')
    dat = read.csv(file_path);  dat$X = NULL
    
    ## 补充计算F-statistic
    dat['F_stat'] = (dat$beta.exposure^2)/(dat$se.exposure^2)
    dat['F_stat_check'] = ifelse(dat$F_stat>10, 0, 1) # 1表示属于弱工具量，需要被剔除
    
    ## 补充计算MAF值
    dat['MAF'] = ifelse(dat$eaf.exposure<=0.5, dat$eaf.exposure, 1-dat$eaf.exposure)
    dat['MAF_check'] = ifelse(dat$MAF>0.01, 0, 1) # 1表示低MAF，需要被剔除
    
    ## 生成备份数据
    file_name2 = paste(c(folder_path, paste(c('check_',file_name), collapse = '')), collapse = '')
    #write.csv(dat, file_name2)
    
    ## 生成最终数据
    check_list = dat$F_stat_check == 0 & dat$MAF_check == 0
    dat_final = dat[check_list, ]
    file_name3 = paste(c(folder_path, paste(c('final_',file_name), collapse = '')), collapse = '')
    write.csv(dat_final, file_name3)
    
    Sum = Sum +sum(dim(dat)[1] != dim(dat_final)[1])
    if (dim(dat)[1] != dim(dat_final)[1]){
      print(file_name)
      Check_names = append(Check_names, file_name)}
  }
  print('被矫正的样本数量')
  print(Sum)
  Check_file = t(data.frame(Check_names)); colnames(Check_file) = 'Metabolites'
  save_path = paste(c(folder_path, '被矫正的代谢物.csv'), collapse = '')
  write.csv(Check_file, save_path)
}

## Harmonized + PRESSO剔除多效性变量
PRESSO_Check_loop <- function(folder_path = NULL,
                         out_2 = NULL){
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
      OUTLIERtest = TRUE,DISTORTIONtest = FALSE, NbDistribution = 1000,
      data = mr_dat,
      SignifThreshold = 0.05)
    
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
