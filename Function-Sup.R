### 替换奇怪的字符们
replace_.. = function(str){
  str_new = sub('  ', ' ', sub('\\.', ' ', str))
  return(str_new)
}

replace__ = function(str){
  str_new = sub('_', ' ', str)
  return(str_new)
}


### 本地开展CLUMP分析
exp_P_LD_local <- function(data = NULL,
                           p_value = 5e-6,
                           clump_kb = 10000,
                           clump_r2 = 0.01,
                           output_file_name = 'output'){
  
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
  SNP_clumped <- ld_clump_local(
    dat = dplyr::tibble(rsid=exp_dat$SNP, pval=exp_dat$pval.exposure, id=exp_dat$id.exposure),
    clump_kb = 10000,
    clump_r2 = 0.01,
    clump_p = p_value,
    bfile = '04-敏感性分析/data_maf0.01_rs_ref/EUR',
    plink_bin = genetics.binaRies::get_plink_binary()
  )
  exp_dat_clumped = exp_dat[exp_dat$SNP %in% SNP_clumped$rsid,]
  
  output_path = paste0(c(here(), '/04-敏感性分析/output/', 'Clumped_SNP_', output_file_name, '.csv'),
                       collapse = '')
  
  write.csv(exp_dat_clumped, 
            output_path)
  print('The file has been saved')
  
  return(exp_dat_clumped)
}


## Input the ID-name you want (like: 'midas_s_113782')
select_FASTA <- function(sintax_sequences,
                         seq.name = 'midas_s_113782',
                         new.name = NULL,
                         rename = FALSE,
                         save.folder = '04_Formal_Results/Phylogenetic tree/selected_FASTA/') {
  full.names = names(sintax_sequences)
  Num = which(full.names == full.names[grepl(seq.name, full.names)])
  if (rename == FALSE) {
    if (is.null(new.name)) {
      writeXStringSet(sintax_sequences[Num],
                      paste0(save.folder, seq.name, '.fasta'))
      cat('The file (orginal name) has been saved')
    } else {
      writeXStringSet(sintax_sequences[Num],
                      paste0(save.folder, seq.name,'(', new.name,').fasta')) # for backup
      writeXStringSet(sintax_sequences[Num],
                      paste0(save.folder, new.name, '.fasta'))
      cat('The two files (orginal and new name) have been saved')
    }
  } else {
    seq.tmp = sintax_sequences[Num]
    writeXStringSet(seq.tmp,
                    paste0(save.folder, seq.name,'(', new.name,').fasta')) # for backup
    names(seq.tmp) = paste0(seq.name, '-', new.name)
    writeXStringSet(seq.tmp,
                    paste0(save.folder, new.name, '.fasta'))
    cat('The name in the file has been modified \n')
    cat('The two files (orginal and new name) have been saved\n')
  }
}
