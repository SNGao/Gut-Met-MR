library(TwoSampleMR)
library(data.table)

## 输入outcome的ieu编号，自动产生以BioAge为Exposure的两样本MR结果
MR_ieu_BA_exp <- function(Exposure, ieu_label){
  exposure_dat = Exposure #非ieu数据库的结果，由自身计算产生
  outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = ieu_label)
  dat <- harmonise_data(
    exposure_dat = exposure_dat,
    outcome_dat = outcome_dat
  ); mr(dat)
}

## 输入Exposure的ieu编号，自动产生以BioAge为OutCome的两样本MR结果
MR_ieu_BA_out <- function(OutCome, ieu_label){
  exposure_dat <- extract_instruments(outcomes = ieu_label, p1=5e-8) # 参数可以自定义
  
  temp <- merge(exposure_dat,OutCome, by.x = "SNP", by.y = 'SNP')
  write.csv(temp, 'outcome_temp.csv')
  outcome_dat <- read_outcome_data(snps = exposure_dat$SNP,
                                   filename = "outcome_temp.csv",sep = ",",
                                   snp_col = "SNP",beta_col = "beta",eaf = 'freq',
                                   se_col = "SE",effect_allele_col = "effect_allele",
                                   other_allele_col = "other_allele",pval_col = "p")
  dat <- harmonise_data(
    exposure_dat = exposure_dat,
    outcome_dat = outcome_dat
  ); mr(dat)
  #mr_scatter_plot(mr_results = mr(dat,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),dat)
} 

## 输入输出均为ieu编号
Exp_Out_ieu <- function(exp_label, out_label){
  exp_dat = extract_instruments(outcomes = exp_label, p1=5e-7)
  out_dat = extract_outcome_data(snps = exp_dat$SNP, outcomes = out_label)
  dat <- harmonise_data(
    exposure_dat = exp_dat,
    outcome_dat = out_dat
  );  
  mr(dat)
  #mr_scatter_plot(mr_results = mr(dat,method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")),dat)
}
