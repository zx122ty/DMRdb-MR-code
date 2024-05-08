MR_analysis = function(exposure_file, outcome_id, outcome_file, exposure_sample, exposure_file_address){
  exposure_id  <- gsub('.csv','',exposure_file)
  exposure_dat = fread(paste0(exposure_file_address, exposure_file))
  exposure_dat$id.exposure =  exposure_id
  
  if(length(exposure_sample)!=0 && !is.na(exposure_sample)[1] && exposure_sample[1]>0){
    exposure_dat$R2<- 2*(1-exposure_dat$eaf.exposure)*exposure_dat$eaf.exposure*(exposure_dat$beta.exposure)^2
    exposure_dat$F<- (exposure_dat$R2)/(1-exposure_dat$R2)*(exposure_sample-2)
    exposure_dat <- exposure_dat[exposure_dat$F > 10,]
  }
  
  if(nrow(exposure_dat) > 0){
    dat <- df_batch[df_batch$SNP %in% exposure_dat$SNP,]
    dat = as.data.frame(dat)
    
    if(nrow(dat)>0){
      outcome_dat <- TwoSampleMR::format_data(
        dat,
        snps = exposure_dat$SNP,
        type = "outcome",
        snp_col = "SNP",
        beta_col = "beta",
        se_col = "se",
        eaf_col = "eaf",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        pval_col = "pval",
        samplesize_col = "samplesize",
      )
      
      
      dat <- TwoSampleMR::harmonise_data(exposure_dat = exposure_dat,
                                         outcome_dat = outcome_dat,
                                         action = 3)
      
      dat = dat[dat$pval.outcome>=1e-05,]
      
      result1 <- TwoSampleMR::mr(dat, method_list = 'mr_ivw')
      
      if(nrow(result1)>0){
        result = TwoSampleMR::mr(dat)
        rslt <- '...'
        if(!dir.exists(paste0(rslt, exposure_id,'~',outcome_id))){dir.create(paste0(rslt, exposure_id,'~', outcome_id))}
        rslt <- paste0(rslt, exposure_id,'~', outcome_id)
        fwrite(result1, file = paste0(rslt, '/mr_result_ivw.csv'),row.names = F)
        fwrite(result, file = paste0(rslt, '/mr_result.csv'),row.names = F)
        fwrite(dat, file = paste0(rslt, '/dat.csv'),row.names = F)
        
        plotall(rslt, result, dat)
        
      }
    }
  }
}