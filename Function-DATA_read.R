
# 从本地读取Exposure数据(UKB)
## 该函数仅当UKB数据下载到本地后使用
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
## 输出所需要的Colname，及时修正
dat_read <- function(filepath = NULL){
  exp <- fread(filepath, header = T)
  exp$V1 = NULL
  print(colnames(exp))
  print('我需要的colname label')
  print(c('SNP','effect_allele','other_allele','freq','p','beta','SE'))
  return(exp)
}

# 修改读入table的变量名字，为后续分析服务
## 衔接dat_read函数
## 注意，考虑的名称修改存在局限，可以补充增加自定义的函数名字（该功能尚未实现，后续添加）
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
    if(cor_name %in% c('other_allele','ALLELE0','ref','A2',
                       'Allele2','ref.allele','OTHER_ALLELE')){ colnames(data)[i]='other_allele' }
    else if(cor_name %in% c('effect_allele','ALLELE1','alt','A1',
                            'Allele1','eff.allele','EFFECT_ALLELE')){ colnames(data)[i]='effect_allele' }
    else if(cor_name %in% c('A1FREQ','freq','af_alt','MAF','FRQ_A_40463','EAF',
                            'effect_allele_frequency','Freq1','Freq','EAF_HRC',
                            'eaf.exposure')){ colnames(data)[i]='freq' }
    else if(cor_name %in% c('P_LINREG','p','pval','P','p_value','P.Neff',
                            'P.weightedSumZ','P-value','Pval','pval.exposure')){ colnames(data)[i]='p' }
    else if(cor_name %in% c('BETA','beta','Effect','LogOR','stdBeta',
                            'Beta','est','beta.exposure')){ colnames(data)[i]='beta' }
    else if(cor_name %in% c('SE','sebeta','standard_error','StdErr',
                            'StdErrLogOR','Markername','se.exposure')){ colnames(data)[i]='SE' }
    else if(cor_name %in% c('SNP','rsids','rsid','snp','MarkerName',
                            'rsID','variant_id')){ colnames(data)[i]='SNP' }
  }
  print(colnames(data))
  data$effect_allele = toupper(data$effect_allele)
  data$other_allele = toupper(data$other_allele)
  return(data.frame(data))
}


