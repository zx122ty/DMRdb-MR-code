plotall = function(file, mr_results, dat){

  # 0  change name------------
  dat$outcome = mr_results$outcome[1]
  dat$id.outcome = mr_results$id.outcome[1]
  dat$exposure = mr_results$exposure[1]
  dat$id.exposure = mr_results$id.exposure[1]
  write.csv(dat, row.names = F, file = paste0(file, '/dat.csv'))
  
  
  # 1 generate_odds_ratios-----------
  mr_result2<-generate_odds_ratios(mr_results)
  write.csv(mr_result2, row.names = F, file = paste0(file, '/mr_result2.csv'))
  
  
  # 2 mr_pleiotropy_test------------
  ple <- mr_pleiotropy_test(dat)
  write.csv(ple,paste0(file, '/ple.csv'), row.names = FALSE)
  
  
  # 3 mr_scatter_plot-------------
  A3<-function (mr_results, dat) 
  {
    requireNamespace("plyr", quietly = TRUE)
    mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), 
                         function(d) {
                           d <- plyr::mutate(d)
                           d <- subset(d, mr_keep)
                           index <- d$beta.exposure < 0
                           d$beta.exposure[index] <- d$beta.exposure[index] * 
                             -1
                           d$beta.outcome[index] <- d$beta.outcome[index] * 
                             -1
                           mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & 
                                             id.outcome == d$id.outcome[1])
                           
                           
                           mrres$a <- 0
                           if ("MR Egger" %in% mrres$method) {
                             temp <- mr_egger_regression(d$beta.exposure, 
                                                         d$beta.outcome, d$se.exposure, d$se.outcome, 
                                                         default_parameters())
                             mrres$a[mrres$method == "MR Egger"] <- temp$b_i
                           }
                           if ("MR Egger (bootstrap)" %in% mrres$method) {
                             temp <- mr_egger_regression_bootstrap(d$beta.exposure, 
                                                                   d$beta.outcome, d$se.exposure, d$se.outcome, 
                                                                   default_parameters())
                             mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
                           }
                           
                           
                           ggplot2::ggplot(data = d, ggplot2::aes(x = beta.exposure, 
                                                                  y = beta.outcome)) + ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome - 
                                                                                                                             se.outcome, ymax = beta.outcome + se.outcome), 
                                                                                                              colour = "grey", width = 0) + ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure - 
                                                                                                                                                                                   se.exposure, xmax = beta.exposure + se.exposure), 
                                                                                                                                                                    colour = "grey", height = 0) + ggplot2::geom_point() + 
                             ggplot2::geom_abline(data = mrres, ggplot2::aes(intercept = a, 
                                                                             slope = b, colour = method), show.legend = TRUE) + 
                             ggplot2::scale_colour_manual(values = c("#a6cee3", 
                                                                     "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", 
                                                                     "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", 
                                                                     "#6a3d9a", "#ffff99", "#b15928")) + ggplot2::labs(colour = "MR Test", 
                                                                                                                       x = paste("SNP effect on", d$exposure[1]), y = paste("SNP effect on", 
                                                                                                                                                                            d$outcome[1])) + ggplot2::theme(legend.position = "right", 
                                                                                                                                                                                                            legend.direction = "vertical") + ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1))+
                             theme_bw() 
                         })
    mrres
  }
  p <- A3(mr_results, dat)
  pdf(paste0(file, "/scatter.pdf"), width = 8, height = 6 )
  print(p[[1]])
  dev.off()
  
 
  # 4 mr_singlesnp-----------
  mr_outcome_single <- mr_singlesnp(dat)
  A4 = function (singlesnp_results, exponentiate = FALSE) 
  {
    res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), 
                       function(d) {
                         d <- plyr::mutate(d)
                         if (sum(!grepl("All", d$SNP)) < 2) {
                           return(blank_plot("Insufficient number of SNPs"))
                         }
                         levels(d$SNP)[levels(d$SNP) == "All - Inverse variance weighted"] <- "All - IVW"
                         levels(d$SNP)[levels(d$SNP) == "All - MR Egger"] <- "All - Egger"
                         am <- grep("All", d$SNP, value = TRUE)
                         d$up <- d$b + 1.96 * d$se
                         d$lo <- d$b - 1.96 * d$se
                         d$tot <- 0.01
                         d$tot[d$SNP %in% am] <- 1
                         d$SNP <- as.character(d$SNP)
                         nom <- d$SNP[!d$SNP %in% am]
                         nom <- nom[order(d$b)]
                         d <- rbind(d, d[nrow(d), ])
                         d$SNP[nrow(d) - 1] <- ""
                         d$b[nrow(d) - 1] <- NA
                         d$up[nrow(d) - 1] <- NA
                         d$lo[nrow(d) - 1] <- NA
                         d$SNP <- ordered(d$SNP, levels = c(am, "", nom))
                         xint <- 0
                         if (exponentiate) {
                           d$b <- exp(d$b)
                           d$up <- exp(d$up)
                           d$lo <- exp(d$lo)
                           xint <- 1
                         }
                         ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + 
                           ggplot2::geom_vline(xintercept = xint, linetype = "dotted") + 
                           ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, 
                                                                xmax = up, size = as.factor(tot), colour = as.factor(tot)), 
                                                   height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) + 
                           ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% 
                                                                                 "")), colour = "grey") + ggplot2::scale_colour_manual(values = c("black", 
                                                                                                                                                  "red")) + ggplot2::scale_size_manual(values = c(0.3, 
                                                                                                                                                                                                  1)) + ggplot2::theme(legend.position = "none", 
                                                                                                                                                                                                                       axis.text.y = ggplot2::element_text(size = 8), 
                                                                                                                                                                                                                       axis.ticks.y = ggplot2::element_line(size = 0), 
                                                                                                                                                                                                                       axis.title.x = ggplot2::element_text(size = 8)) + 
                           ggplot2::labs(y = "", x = paste0("MR effect size for\n'", 
                                                            d$exposure[1], "' on '", d$outcome[1], "'"))+
                           theme_bw()
                       })
    res
  }
  p4 <- A4(mr_outcome_single)
  pdf(paste0(file,"/single_forest.pdf"), width = 8, height = 6)
  print(p4[[1]])
  dev.off()
  
  
  # 5 funnel plot------------------
  A5 = function (singlesnp_results) 
  {
    res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), 
                       function(d) {
                         d <- plyr::mutate(d)
                         if (sum(!grepl("All", d$SNP)) < 2) {
                           return(blank_plot("Insufficient number of SNPs"))
                         }
                         am <- grep("All", d$SNP, value = TRUE)
                         d$SNP <- gsub("All - ", "", d$SNP)
                         am <- gsub("All - ", "", am)
                         ggplot2::ggplot(subset(d, !SNP %in% am), ggplot2::aes(y = 1/se, 
                                                                               x = b)) + ggplot2::geom_point() + ggplot2::geom_vline(data = subset(d, 
                                                                                                                                                   SNP %in% am), ggplot2::aes(xintercept = b, colour = SNP)) + 
                           ggplot2::scale_colour_manual(values = c("#a6cee3", 
                                                                   "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", 
                                                                   "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", 
                                                                   "#6a3d9a", "#ffff99", "#b15928")) + ggplot2::labs(y = expression(1/SE[IV]), 
                                                                                                                     x = expression(beta[IV]), colour = "MR Method") + 
                           ggplot2::theme(legend.position = "top", legend.direction = "vertical")+
                           
                           theme_bw()
                       })
    res
  }
  p5 <- A5(mr_outcome_single)
  pdf(paste0(file,"/funnel.pdf"), width = 8, height = 6 )
  print(p5[[1]])
  dev.off()
 
  
  # 6 Leave-one-out plot----------------------
  mr_outcome_loo <- mr_leaveoneout(dat)
  A6 = function (leaveoneout_results) 
  {
    res <- plyr::dlply(leaveoneout_results, c("id.exposure", 
                                              "id.outcome"), function(d) {
                                                d <- plyr::mutate(d)
                                                if (sum(!grepl("All", d$SNP)) < 3) {
                                                  return(blank_plot("Insufficient number of SNPs"))
                                                }
                                                d$up <- d$b + 1.96 * d$se
                                                d$lo <- d$b - 1.96 * d$se
                                                d$tot <- 1
                                                d$tot[d$SNP != "All"] <- 0.01
                                                d$SNP <- as.character(d$SNP)
                                                nom <- d$SNP[d$SNP != "All"]
                                                nom <- nom[order(d$b)]
                                                d <- rbind(d, d[nrow(d), ])
                                                d$SNP[nrow(d) - 1] <- ""
                                                d$b[nrow(d) - 1] <- NA
                                                d$up[nrow(d) - 1] <- NA
                                                d$lo[nrow(d) - 1] <- NA
                                                d$SNP <- ordered(d$SNP, levels = c("All", "", nom))
                                                print(
                                                  
                                                  ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + ggplot2::geom_vline(xintercept = 0, 
                                                                                                                         linetype = "dotted") + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, 
                                                                                                                                                                                     xmax = up, size = as.factor(tot), colour = as.factor(tot)), 
                                                                                                                                                                        height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) + 
                                                    ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% 
                                                                                                          "")), colour = "grey") + ggplot2::scale_colour_manual(values = c("black", 
                                                                                                                                                                           "red")) + ggplot2::scale_size_manual(values = c(0.3, 
                                                                                                                                                                                                                           1)) + ggplot2::theme(legend.position = "none", axis.text.y = ggplot2::element_text(size = 8), 
                                                                                                                                                                                                                                                axis.ticks.y = ggplot2::element_line(size = 0), 
                                                                                                                                                                                                                                                axis.title.x = ggplot2::element_text(size = 8)) + 
                                                    ggplot2::labs(y = "", x = paste0("MR leave-one-out sensitivity analysis for\n'", 
                                                                                     d$exposure[1], "' on '", d$outcome[1], "'"))+
                                                    
                                                    theme_bw()
                                                  
                                                  
                                                )
                                                
                                              })
    res
  }
  p6 <- A6(mr_outcome_loo)
  pdf(paste0(file,"/leave-one-out.pdf"), width = 8, height = 6 )
  print(p6[[1]])
  dev.off()
  
  
  # 7 mr_heterogeneity-------------------
  H7<-mr_heterogeneity(dat)
  write.csv(H7, paste0(file,'/H.csv'), row.names = FALSE)
  
  
  # 8 ivw_radial-------------------
  if(T){
    if (exists("RadialMR_r") && exists("RadialMR_r1")) {
      result <- try({
        pdf(paste0(file,"/radial.pdf"), width = 12, height = 8)
        on.exit(dev.off())
        plot_radial(c(RadialMR_r, RadialMR_r1))
      }, silent = TRUE)
      
      if (inherits(result, "try-error")) {
        cat("An error occurred in plotting RadialMR_r and RadialMR_r1: ", result, "\nSkipping plot.\n")
      }
      
    } else if (exists("RadialMR_r")) {
      
      result <- try({
        pdf(paste0(file,"/radial.pdf"), width = 12, height = 8)
        on.exit(dev.off())
        plot_radial(RadialMR_r)
      }, silent = TRUE)
      
      if (inherits(result, "try-error")) {
        cat("An error occurred in plotting RadialMR_r: ", result, "\nSkipping plot.\n")
      }
      
    } else if (exists("RadialMR_r1")) {
      
      result <- try({
        pdf(paste0(file,"/radial.pdf"), width = 12, height = 8)
        on.exit(dev.off())
        plot_radial(RadialMR_r1)
      }, silent = TRUE)
      
      if (inherits(result, "try-error")) {
        cat("An error occurred in plotting RadialMR_r1: ", result, "\nSkipping plot.\n")
      }
      
    } else {
      cat("Neither RadialMR_r nor RadialMR_r1 exists, skipping plot.\n")
    }
    
  }
 

  # 9 mr_presso-------------------
  if (T) {
    tryCatch({
      presso_result <- mr_presso(BetaOutcome ="beta.outcome", 
                                 BetaExposure = "beta.exposure", 
                                 SdOutcome ="se.outcome", 
                                 SdExposure = "se.exposure", 
                                 OUTLIERtest = TRUE, 
                                 DISTORTIONtest = TRUE, 
                                 data = dat, 
                                 NbDistribution = 1000, 
                                 SignifThreshold = 0.05)
      

      vector3 = c(43, 20, 12)

      index = presso_result[["MR-PRESSO results"]][["Distortion Test"]][["Outliers Indices"]]
      vector <- rep(0, length(dat$SNP))
      for(i in index){vector[i] <- 1}
      
      
      SNP <- dat$SNP
      RSSobs <- presso_result[["MR-PRESSO results"]][["Outlier Test"]]$RSSobs
      if(length(RSSobs) == 0) {
        RSSobs <- rep(NA, length(SNP))
      }
      
      Pval <- presso_result[["MR-PRESSO results"]][["Outlier Test"]]$Pvalue
      if(length(Pval) == 0) {
        Pval <- rep(NA, length(SNP))
      }

      summary  <- data.frame(SNP = SNP, RSSobs = RSSobs, Pval = Pval)
      summary$Outliers = vector

      globalRSSobs = presso_result[["MR-PRESSO results"]][["Global Test"]][["RSSobs"]]
      global = presso_result[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
      distortion_test_pval= presso_result[["MR-PRESSO results"]][["Distortion Test"]][["Pvalue"]]
      
      all <- dplyr::tibble(
        SNP = "All",
        RSSobs = globalRSSobs,
        Pval = global,
        Outliers = NA,
        Distortion_Test_Pval = if(length(distortion_test_pval) == 0) NA else distortion_test_pval
      )


      vector2 <- rep(NA, length(dat$SNP))
      Distortion_Test = data.frame("Distortion_Test_Pval"=vector2)

      summary = cbind(summary,Distortion_Test)
      all = rbind(all, summary)
      write.csv(all, paste0(file,'/mrpresso.csv'), row.names = FALSE,na = "")
      
      
    }, error = function(e) {
      cat("An error has occurred: ", e$message, "\n")
    })
    
  }
  
  
}