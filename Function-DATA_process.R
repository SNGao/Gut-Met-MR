
# 对exp数据进行预处理
## 自动判断是否有满足p_value的SNP，并自动计算eaf值
## 已经关闭文件导出功能
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

# 根据结局和暴露，对SNP开展harmonize，返回协调数据集
## 若无协调SNP，则返回空集，并输出Warning
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