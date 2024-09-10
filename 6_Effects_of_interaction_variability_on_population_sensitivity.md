# Effects of intereaction variability on population sensitivity tested by multiple regression approaches

## Load required packages

``` r
sessionInfo() #save session information
```

    ## R version 3.6.3 (2020-02-29)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19045)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=Japanese_Japan.932  LC_CTYPE=Japanese_Japan.932   
    ## [3] LC_MONETARY=Japanese_Japan.932 LC_NUMERIC=C                  
    ## [5] LC_TIME=Japanese_Japan.932    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] afex_0.28-0        emmeans_1.5.2-1    car_3.0-10         carData_3.0-4     
    ##  [5] lme4_1.1-25        Matrix_1.2-18      ggeffects_1.2.0    visreg_2.7.0      
    ##  [9] glmmTMB_1.0.2.1    igraph_1.2.6       rmarkdown_2.5      patchwork_1.1.1   
    ## [13] ggsci_2.9          ggplot2_3.3.2      RColorBrewer_1.1-2 tidyr_1.1.2       
    ## [17] dplyr_1.0.2       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] splines_3.6.3       statmod_1.4.35      cellranger_1.1.0    yaml_2.2.1         
    ##  [5] numDeriv_2016.8-1.1 pillar_1.4.6        lattice_0.20-38     glue_1.4.2         
    ##  [9] digest_0.6.27       minqa_1.2.4         colorspace_2.0-3    sandwich_3.0-0     
    ## [13] plyr_1.8.6          htmltools_0.5.2     pkgconfig_2.0.3     haven_2.3.1        
    ## [17] purrr_0.3.4         xtable_1.8-4        mvtnorm_1.1-1       scales_1.1.1       
    ## [21] openxlsx_4.2.3      rio_0.5.16          tibble_3.0.4        generics_0.1.0     
    ## [25] farver_2.0.3        ellipsis_0.3.2      TH.data_1.0-10      withr_2.5.0        
    ## [29] TMB_1.7.18          cli_3.4.1           survival_3.1-8      magrittr_2.0.3     
    ## [33] crayon_1.3.4        readxl_1.3.1        estimability_1.3    evaluate_0.14      
    ## [37] fansi_1.0.3         nlme_3.1-150        MASS_7.3-53         forcats_0.5.0      
    ## [41] foreign_0.8-75      tools_3.6.3         data.table_1.13.2   hms_0.5.3          
    ## [45] lifecycle_1.0.3     multcomp_1.4-14     stringr_1.4.0       munsell_0.5.0      
    ## [49] zip_2.1.1           compiler_3.6.3      rlang_1.1.0         grid_3.6.3         
    ## [53] nloptr_1.2.2.2      rstudioapi_0.11     labeling_0.4.2      boot_1.3-24        
    ## [57] lmerTest_3.1-3      gtable_0.3.0        codetools_0.2-16    abind_1.4-5        
    ## [61] curl_4.3            reshape2_1.4.4      R6_2.5.0            zoo_1.8-8          
    ## [65] knitr_1.30          fastmap_1.1.0       utf8_1.1.4          stringi_1.4.6      
    ## [69] parallel_3.6.3      Rcpp_1.0.11         vctrs_0.6.1         tidyselect_1.1.0   
    ## [73] xfun_0.19           coda_0.19-4

``` r
library(ggplot2); packageVersion("ggplot2") #3.3.2
```

    ## [1] '3.3.2'

``` r
library(dplyr); packageVersion("dplyr") #1.0.2
```

    ## [1] '1.0.2'

``` r
library(tidyr); packageVersion("tidyr") #1.1.2
```

    ## [1] '1.1.2'

``` r
library(patchwork); packageVersion("patchwork") #1.1.1
```

    ## [1] '1.1.1'

``` r
library(RColorBrewer); packageVersion("RColorBrewer") #1.1.2
```

    ## [1] '1.1.2'

``` r
library(glmmTMB); packageVersion("glmmTMB") #1.0.2.1
```

    ## [1] '1.0.2.1'

``` r
library(visreg); packageVersion("visreg") #2.7.0
```

    ## [1] '2.7.0'

``` r
library(ggeffects); packageVersion("ggeffects") #1.2.0
```

    ## [1] '1.2.0'

``` r
library(lme4); packageVersion("lme4") #1.1.25
```

    ## [1] '1.1.25'

``` r
library(car); packageVersion("car") #3.0.10
```

    ## [1] '3.0.10'

``` r
library(emmeans); packageVersion("emmeans") #1.5.2.1
```

    ## [1] '1.5.2.1'

``` r
library(afex); packageVersion("afex") #0.28.0
```

    ## [1] '0.28.0'

## Data compilation

``` r
all_multsmap_coefs <- read.csv("./processed_data/all_regulsmap_coefs240123.csv")
all_Ts <- read.csv("./processed_data/all_Ts240123.csv")
Abundance_TrEff <- read.csv("./processed_data/Abundance_TrEff.csv", header=TRUE) #load population sensitivity data

all_multsmap_coefs$Recipient <- #define necessary columns
  factor(all_multsmap_coefs$Recipient, 
         levels=unique(all_multsmap_coefs$Recipient))
all_multsmap_coefs$Donor <- 
  factor(all_multsmap_coefs$Donor, 
         levels=c("Phytopla1", "Phytopla1_1", "Roti1", "Roti1_1", "Roti1_2", 
                  "Zoopla1", "Zoopla1_1", "Zoopla1_2", "Zoopla1_3", "Zoopla1_4", 
                  "Mp1", "Mp1_1", "Det1", "Herb1", "Herb1_1", "Herb1_2", 
                  "Pred.P1", "Pred.P1_1", "Pred.P1_2", "Pred.P1_3", "Pred.B1", 
                  "Pred.S1", "Pred.S1_1", "Pred.S1_2", "Pred.S1_3", "Moll1", "0"))
all_multsmap_coefs$Treatment <- factor(all_multsmap_coefs$Treatment, levels=c("Cont", "Fipro", "Pent", "Joint"))
all_multsmap_coefs <- 
  all_multsmap_coefs %>%
  unite(col="Recipient_Donor", Recipient, Donor, remove=FALSE)
```

## Extract interaction density-dependence

Density-dependence in per capita interaction effects were detemined by
regressing S-map coefs by recipient density at each time point.

``` r
all_multsmap_coefs_intersp <- 
  all_multsmap_coefs %>%
  separate(Donor, into=c("Donor_", "delay"), sep="_", remove=FALSE) %>% 
  filter(Recipient!=Donor_ & r_d!="c_0") %>% 
  select(-Donor_, - delay)%>% 
  unite("Year_Tank", Year, Tank, remove=FALSE) %>%
  ungroup(.)
```

    ## Warning: Expected 2 pieces. Missing pieces filled with `NA` in 9800 rows [1, 2, 3, 4, 5,
    ## 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].

``` r
all_multsmap_coefs_intersp_ab <- left_join(mutate(all_multsmap_coefs_intersp, Year=as.character(Year)), all_Ts)
```

    ## Joining, by = c("Year_Tank", "Week", "Recipient")

``` r
ddepend_lmres_list <-
  all_multsmap_coefs_intersp_ab %>%
  unite("Recipient_Donor", Recipient_Donor) %>%
  split(., f=.$Recipient_Donor) %>%
  lapply(function(data) as.data.frame(summary(lm(smap_coef ~ scAbundance, data))$coefficients)) #regressions in parallel

ddepend_lmres <-
  ddepend_lmres_list %>% do.call(rbind, .) %>%
  mutate(Recipient_Donor=rep(names(ddepend_lmres_list), each=2),
         Source=rep(c("Intercept", "scAbundance"), times=nrow(.)/2))

ddepend_lmres <- 
  ddepend_lmres %>%
  separate(Recipient_Donor, into=c("Recipient", "Donor"), sep="_") %>%
  filter(Source=="scAbundance") %>%
  mutate(SignDD=ifelse(Estimate>0, "PositiveDD", "NegativeDD"))

head(ddepend_lmres) #display data sheet
```

    ##        Estimate   Std. Error    t value     Pr(>|t|) Recipient     Donor      Source
    ## 1 -0.0212262846 0.0036509319 -5.8139361 2.025812e-08      Det1       Mp1 scAbundance
    ## 2 -0.0419480400 0.0049520150 -8.4709032 2.916973e-15      Det1 Phytopla1 scAbundance
    ## 3 -0.0057660151 0.0024092079 -2.3933240 1.749947e-02      Det1     Roti1 scAbundance
    ## 4 -0.0002354927 0.0003287232 -0.7163860 4.746326e-01     Herb1      Det1 scAbundance
    ## 5 -0.0009535404 0.0003550246 -2.6858436 7.875216e-03     Herb1       Mp1 scAbundance
    ## 6 -0.0001063029 0.0003216795 -0.3304621 7.414150e-01     Herb1 Phytopla1 scAbundance
    ##       SignDD
    ## 1 NegativeDD
    ## 2 NegativeDD
    ## 3 NegativeDD
    ## 4 NegativeDD
    ## 5 NegativeDD
    ## 6 NegativeDD

## Calculate mean interaction strength and interaction temporal variability in the controls

``` r
# calculate smap coef statistics

all_multsmap_coefs_smrz <- 
  all_multsmap_coefs %>%
  group_by(Treatment, Recipient, Donor, r_d, Recipient_Donor, Year, Tank) %>% #summarise per year per tank
  summarise(smap_coef_mean=mean(smap_coef, na.rm=TRUE), 
            smap_coef_SD=sd(smap_coef, na.rm=TRUE)) %>% 
  group_by(Treatment, Recipient, Donor, r_d, Recipient_Donor) %>% #summarise over years and tanks
  summarise(smap_coef_mean=mean(abs(smap_coef_mean), na.rm=TRUE), 
            smap_coef_SD=mean(smap_coef_SD, na.rm=TRUE))
```

    ## `summarise()` regrouping output by 'Treatment', 'Recipient', 'Donor', 'r_d',
    ## 'Recipient_Donor', 'Year' (override with `.groups` argument)`summarise()` regrouping
    ## output by 'Treatment', 'Recipient', 'Donor', 'r_d' (override with `.groups` argument)

``` r
# remove unnecessary rows

all_multsmap_coefs_smrz_intersp <- 
  all_multsmap_coefs_smrz %>%
  separate(Donor, into=c("Donor_", "delay"), sep="_", remove=FALSE) %>% 
  filter(Recipient!=Donor_ & r_d!="c_0") %>% #omit intraspecific effects and intercept
  select(-Donor_, - delay)%>% #omit unnecessary columns
  ungroup(.)
```

    ## Warning: Expected 2 pieces. Missing pieces filled with `NA` in 164 rows [1, 3, 4, 5, 6, 7,
    ## 8, 11, 12, 13, 14, 15, 20, 21, 22, 24, 25, 26, 27, 28, ...].

``` r
# calculate diffs in smap coefs (unnecessary for analysis but used to compile data)

smapcoef_TrEff <- #difference in smap coefs between the pesticide treatments and the controls 
  all_multsmap_coefs_smrz_intersp %>%
  ungroup(.) %>%
  select(-r_d, -smap_coef_SD) %>%
  spread(value=smap_coef_mean, key=Treatment) %>% 
  mutate(Fipro_Cont=Fipro-Cont, 
         Pent_Cont=Pent-Cont,
         Joint_Cont=Joint-Cont) %>%
  select(-Cont, -Fipro, -Pent, -Joint) %>% 
  gather(key=Treatment_Cont, value=IS_response, -(1:3)) %>% 
  mutate(Treatment_Cont=factor(Treatment_Cont, levels=c("Fipro_Cont", "Pent_Cont", "Joint_Cont"))) %>%
  left_join(., select(filter(all_multsmap_coefs_smrz_intersp, Treatment=="Cont"), Recipient_Donor, smap_coef_SD, smap_coef_mean)) 
```

    ## Joining, by = "Recipient_Donor"

``` r
# assemble (merge) data

smapcoef_TrEff <- 
  smapcoef_TrEff %>%
  left_join(., Abundance_TrEff) %>%
  rename(R_Abundance_response_nRR_0.1=Abundance_response_nRR_0.1) %>%
  left_join(., Abundance_TrEff, by=c("Donor"="Recipient", "Year"="Year", "Treatment_Cont"="Treatment_Cont")) %>%
  rename(D_Abundance_response_nRR_0.1=Abundance_response_nRR_0.1)
```

    ## Joining, by = c("Recipient", "Treatment_Cont")

``` r
smapcoef_TrEff_SignDD <- left_join(smapcoef_TrEff, ddepend_lmres)
```

    ## Joining, by = c("Recipient", "Donor")

``` r
smapcoef_TrEff_SignDD <- 
  smapcoef_TrEff_SignDD %>% 
  mutate(abs_R_Abundance_response_nRR_0.1=abs(R_Abundance_response_nRR_0.1), 
         abs_D_Abundance_response_nRR_0.1=abs(D_Abundance_response_nRR_0.1), 
         log_abs_DD=log(abs(Estimate)))

# display data sheet

head(as.data.frame(smapcoef_TrEff_SignDD))
```

    ##   Recipient Donor Recipient_Donor Treatment_Cont  IS_response smap_coef_SD smap_coef_mean
    ## 1 Phytopla1 Roti1 Phytopla1_Roti1     Fipro_Cont -0.006422403   0.06765366     0.03356245
    ## 2 Phytopla1 Roti1 Phytopla1_Roti1     Fipro_Cont -0.006422403   0.06765366     0.03356245
    ## 3 Phytopla1 Roti1 Phytopla1_Roti1     Fipro_Cont -0.006422403   0.06765366     0.03356245
    ## 4 Phytopla1  Det1  Phytopla1_Det1     Fipro_Cont -0.030756020   0.08478359     0.16869603
    ## 5 Phytopla1  Det1  Phytopla1_Det1     Fipro_Cont -0.030756020   0.08478359     0.16869603
    ## 6 Phytopla1  Det1  Phytopla1_Det1     Fipro_Cont -0.030756020   0.08478359     0.16869603
    ##   Year R_Abundance_response_nRR_0.1 abs_Abundance_response_nRR_0.1.x
    ## 1 2017                   0.16204186                       0.16204186
    ## 2 2018                   0.06312511                       0.06312511
    ## 3 2019                   2.72402931                       2.72402931
    ## 4 2017                   0.16204186                       0.16204186
    ## 5 2018                   0.06312511                       0.06312511
    ## 6 2019                   2.72402931                       2.72402931
    ##   D_Abundance_response_nRR_0.1 abs_Abundance_response_nRR_0.1.y   Estimate  Std. Error
    ## 1                    0.3195754                        0.3195754 0.06052547 0.005231589
    ## 2                   -0.5466622                        0.5466622 0.06052547 0.005231589
    ## 3                   -0.6312718                        0.6312718 0.06052547 0.005231589
    ## 4                   -0.7515908                        0.7515908 0.02543747 0.007142959
    ## 5                   -0.3493756                        0.3493756 0.02543747 0.007142959
    ## 6                   -0.3972443                        0.3972443 0.02543747 0.007142959
    ##     t value     Pr(>|t|)      Source     SignDD abs_R_Abundance_response_nRR_0.1
    ## 1 11.569233 2.323623e-24 scAbundance PositiveDD                       0.16204186
    ## 2 11.569233 2.323623e-24 scAbundance PositiveDD                       0.06312511
    ## 3 11.569233 2.323623e-24 scAbundance PositiveDD                       2.72402931
    ## 4  3.561195 4.549211e-04 scAbundance PositiveDD                       0.16204186
    ## 5  3.561195 4.549211e-04 scAbundance PositiveDD                       0.06312511
    ## 6  3.561195 4.549211e-04 scAbundance PositiveDD                       2.72402931
    ##   abs_D_Abundance_response_nRR_0.1 log_abs_DD
    ## 1                        0.3195754  -2.804691
    ## 2                        0.5466622  -2.804691
    ## 3                        0.6312718  -2.804691
    ## 4                        0.7515908  -3.671532
    ## 5                        0.3493756  -3.671532
    ## 6                        0.3972443  -3.671532

## Likelihood ratio tests and multiple regressions

``` r
theme_set(theme_test())

Group.list <- data.frame(Group=c("Phytopla1", "Roti1", "Zoopla1", "Mp1", "Det1", "Herb1", "Pred.P1", "Pred.B1", "Pred.S1", "Moll1"), 
                         Function=c("Producer", "Prey", "Prey", "Producer", "Prey", 
                                    "Prey", "Predator", "Predator", "Predator", "Prey"),
                         No.connections=c(4, 7, 2, 5, 4, 5, 3, 1, 2, 5), #unnecessary but maybe useful for readers interested in this property
                         No.causes=c(3, 5, 1, 3, 2, 3, 3, 1, 1, 0), #unnecessary but maybe useful for readers interested in this property
                         No.effects=c(3, 3, 1, 3, 3, 3, 1, 0, 1, 5)) #unnecessary but maybe useful for readers interested in this property
smapcoef_TrEff_SignDD_ <- left_join(smapcoef_TrEff_SignDD, Group.list, by=c("Recipient"="Group"))

# Do likelihood ratio tests repeatedly

Trt_Cont <- c("Fipro_Cont", "Pent_Cont", "Joint_Cont")
LRT_res_list <- NULL
for (i in Trt_Cont) {
  LRT_res <- mixed(log(abs_R_Abundance_response_nRR_0.1+1) ~ 
                     scale(log_abs_DD, scale=FALSE)*SignDD + scale(smap_coef_SD, scale=FALSE) + scale(smap_coef_mean, scale=FALSE) + Function + #centring is necessary for type-3-like tests
                     (1|Recipient_Donor) + (1|Year), 
                   subset(smapcoef_TrEff_SignDD_, Treatment_Cont==i), method="LRT", type=3, 
                   control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-6)), 
                   progress=FALSE)
  LRT_res_list <- c(LRT_res_list, list(LRT_res))
  ISDD_LRT_res <- as.data.frame(anova(LRT_res))
  ISDD_LRT_res[,2] <- round(ISDD_LRT_res[,2], digits=2)
  ISDD_LRT_res[,4] <- round(ISDD_LRT_res[,4], digits=2)
  write.csv(ISDD_LRT_res, paste("./tables/ISDD_LRT_", i, ".csv", sep=""))
}
```

    ## Contrasts set to contr.sum for the following variables: SignDD, Function, Recipient_Donor

    ## Numerical variables NOT centered on 0: log_abs_DD, smap_coef_SD, smap_coef_mean
    ## If in interactions, interpretation of lower order (e.g., main) effects difficult.

    ## REML argument to lmer() set to FALSE for method = 'PB' or 'LRT'

    ## Contrasts set to contr.sum for the following variables: SignDD, Function, Recipient_Donor

    ## Numerical variables NOT centered on 0: log_abs_DD, smap_coef_SD, smap_coef_mean
    ## If in interactions, interpretation of lower order (e.g., main) effects difficult.

    ## REML argument to lmer() set to FALSE for method = 'PB' or 'LRT'

    ## Contrasts set to contr.sum for the following variables: SignDD, Function, Recipient_Donor

    ## Numerical variables NOT centered on 0: log_abs_DD, smap_coef_SD, smap_coef_mean
    ## If in interactions, interpretation of lower order (e.g., main) effects difficult.

    ## REML argument to lmer() set to FALSE for method = 'PB' or 'LRT'

``` r
names(LRT_res_list) <- Trt_Cont
LRT_res_list
```

    ## $Fipro_Cont
    ## Mixed Model Anova Table (Type 3 tests, LRT-method)
    ## 
    ## Model: log(abs_R_Abundance_response_nRR_0.1 + 1) ~ scale(log_abs_DD, 
    ## Model:     scale = FALSE) * SignDD + scale(smap_coef_SD, scale = FALSE) + 
    ## Model:     scale(smap_coef_mean, scale = FALSE) + Function + (1 | Recipient_Donor) + 
    ## Model:     (1 | Year)
    ## Data: subset
    ## Data: smapcoef_TrEff_SignDD_
    ## Data: Treatment_Cont == i
    ## Df full model: 11
    ##                                    Effect df     Chisq p.value
    ## 1        scale(log_abs_DD, scale = FALSE)  1      0.35    .554
    ## 2                                  SignDD  1      0.03    .858
    ## 3      scale(smap_coef_SD, scale = FALSE)  1 26.78 ***   <.001
    ## 4    scale(smap_coef_mean, scale = FALSE)  1      0.28    .594
    ## 5                                Function  2      2.88    .237
    ## 6 scale(log_abs_DD, scale = FALSE):SignDD  1      0.41    .523
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '+' 0.1 ' ' 1
    ## 
    ## $Pent_Cont
    ## Mixed Model Anova Table (Type 3 tests, LRT-method)
    ## 
    ## Model: log(abs_R_Abundance_response_nRR_0.1 + 1) ~ scale(log_abs_DD, 
    ## Model:     scale = FALSE) * SignDD + scale(smap_coef_SD, scale = FALSE) + 
    ## Model:     scale(smap_coef_mean, scale = FALSE) + Function + (1 | Recipient_Donor) + 
    ## Model:     (1 | Year)
    ## Data: subset
    ## Data: smapcoef_TrEff_SignDD_
    ## Data: Treatment_Cont == i
    ## Df full model: 11
    ##                                    Effect df   Chisq p.value
    ## 1        scale(log_abs_DD, scale = FALSE)  1    0.68    .409
    ## 2                                  SignDD  1    0.07    .796
    ## 3      scale(smap_coef_SD, scale = FALSE)  1 7.62 **    .006
    ## 4    scale(smap_coef_mean, scale = FALSE)  1    0.25    .620
    ## 5                                Function  2    1.23    .541
    ## 6 scale(log_abs_DD, scale = FALSE):SignDD  1    1.06    .304
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '+' 0.1 ' ' 1
    ## 
    ## $Joint_Cont
    ## Mixed Model Anova Table (Type 3 tests, LRT-method)
    ## 
    ## Model: log(abs_R_Abundance_response_nRR_0.1 + 1) ~ scale(log_abs_DD, 
    ## Model:     scale = FALSE) * SignDD + scale(smap_coef_SD, scale = FALSE) + 
    ## Model:     scale(smap_coef_mean, scale = FALSE) + Function + (1 | Recipient_Donor) + 
    ## Model:     (1 | Year)
    ## Data: subset
    ## Data: smapcoef_TrEff_SignDD_
    ## Data: Treatment_Cont == i
    ## Df full model: 11
    ##                                    Effect df     Chisq p.value
    ## 1        scale(log_abs_DD, scale = FALSE)  1      1.61    .205
    ## 2                                  SignDD  1      0.12    .727
    ## 3      scale(smap_coef_SD, scale = FALSE)  1 21.90 ***   <.001
    ## 4    scale(smap_coef_mean, scale = FALSE)  1      0.01    .921
    ## 5                                Function  2   9.96 **    .007
    ## 6 scale(log_abs_DD, scale = FALSE):SignDD  1   7.50 **    .006
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '+' 0.1 ' ' 1

``` r
# Display multiple regression results

## Fipro

Vulnerability.Ta_FC.TMB <- glmmTMB(log(abs_R_Abundance_response_nRR_0.1+1) ~ 
                                     log_abs_DD*SignDD + smap_coef_SD + smap_coef_mean + Function + 
                                     (1|Recipient_Donor) + (1|Year), REML=TRUE, 
                                   subset(smapcoef_TrEff_SignDD_, Treatment_Cont=="Fipro_Cont"), 
                                   control=glmmTMBControl(optimizer=optim,
                                                          optArgs=list(method="BFGS")))
p1<- plot(visreg(Vulnerability.Ta_FC.TMB, "smap_coef_mean", re.form=NA, plot=FALSE), band=TRUE, partial=TRUE, gg=TRUE, overlay=TRUE, 
          xlab="Mean interaction strength", ylab="")
```

    ## Warning:   Note that you are attempting to plot a 'main effect' in a model that contains an
    ##   interaction.  This is potentially misleading; you may wish to consider using the 'by'
    ##   argument.

    ## Conditions used in construction of plot
    ## log_abs_DD: -4.206996
    ## SignDD: NegativeDD
    ## smap_coef_SD: 0.03164694
    ## Function: Prey
    ## Recipient_Donor: Det1_Mp1
    ## Year: 2018

``` r
p2<- plot(visreg(Vulnerability.Ta_FC.TMB, "smap_coef_SD", re.form=NA, plot=FALSE), band=TRUE, partial=TRUE, gg=TRUE, overlay=TRUE, 
          xlab="Interaction temporal variability", ylab="")
```

    ## Warning:   Note that you are attempting to plot a 'main effect' in a model that contains an
    ##   interaction.  This is potentially misleading; you may wish to consider using the 'by'
    ##   argument.

    ## Conditions used in construction of plot
    ## log_abs_DD: -4.206996
    ## SignDD: NegativeDD
    ## smap_coef_mean: 0.05984084
    ## Function: Prey
    ## Recipient_Donor: Det1_Mp1
    ## Year: 2018

``` r
p3<- plot(visreg(Vulnerability.Ta_FC.TMB,  "log_abs_DD", by="SignDD", re.form=NA, plot=FALSE), band=TRUE, partial=TRUE, gg=TRUE, overlay=TRUE, 
          xlab="ln(Interaction density dependence)", ylab="ln(Recipient sensitivity + 1)")
p3 <- p3 + scale_y_continuous(limits=c(-0.4, 2.1)) + 
  theme(legend.position=c(0.5, 0.93), legend.justification="center", legend.direction="horizontal") +
  scale_color_manual(name=NULL, labels=c("Negative IDD", "Positive IDD"), values=c("#FF4E37", "#008DFFFF")) +
  scale_fill_manual(name=NULL, labels=c("Negative IDD", "Positive IDD"), values=alpha(c("#FF4E37", "#008DFFFF"), 0.3)) + 
  ggtitle("Model for I vs. C")
```

    ## Scale for 'colour' is already present. Adding another scale for 'colour', which will
    ## replace the existing scale.

    ## Scale for 'fill' is already present. Adding another scale for 'fill', which will
    ## replace the existing scale.

``` r
p4<- plot(visreg(Vulnerability.Ta_FC.TMB, "Function", re.form=NA, plot=FALSE), band=TRUE, partial=TRUE, gg=TRUE, overlay=TRUE, 
          xlab="Recipient function", ylab="")
```

    ## Warning:   Note that you are attempting to plot a 'main effect' in a model that contains an
    ##   interaction.  This is potentially misleading; you may wish to consider using the 'by'
    ##   argument.

    ## Conditions used in construction of plot
    ## log_abs_DD: -4.206996
    ## SignDD: NegativeDD
    ## smap_coef_SD: 0.03164694
    ## smap_coef_mean: 0.05984084
    ## Recipient_Donor: Det1_Mp1
    ## Year: 2018

``` r
(fig.fipro <- p3+p2+p1+p4 + plot_layout(nrow=1, guides="keep"))
```

![](6_Effects_of_interaction_variability_on_population_sensitivity_files/figure-markdown_github/multiple%20regression-1.png)

``` r
## Pent

Vulnerability.Ta_PC.TMB <- glmmTMB(log(abs_R_Abundance_response_nRR_0.1+1) ~ 
                                     (log_abs_DD)*SignDD + smap_coef_SD + smap_coef_mean + Function +
                                     (1|Recipient_Donor) + (1|Year), REML=TRUE, 
                                   subset(smapcoef_TrEff_SignDD_, Treatment_Cont=="Pent_Cont"), 
                                   control=glmmTMBControl(optimizer=optim,
                                                          optArgs=list(method="BFGS")))
p1<- plot(visreg(Vulnerability.Ta_PC.TMB, "smap_coef_mean", re.form=NA, plot=FALSE), band=TRUE, partial=TRUE, gg=TRUE, overlay=TRUE, 
          xlab="Mean interaction strength", ylab="")
```

    ## Warning:   Note that you are attempting to plot a 'main effect' in a model that contains an
    ##   interaction.  This is potentially misleading; you may wish to consider using the 'by'
    ##   argument.

    ## Conditions used in construction of plot
    ## log_abs_DD: -4.206996
    ## SignDD: NegativeDD
    ## smap_coef_SD: 0.03164694
    ## Function: Prey
    ## Recipient_Donor: Det1_Mp1
    ## Year: 2018

``` r
p2<- plot(visreg(Vulnerability.Ta_PC.TMB, "smap_coef_SD", re.form=NA, plot=FALSE), band=TRUE, partial=TRUE, gg=TRUE, overlay=TRUE, 
          xlab="Interaction temporal variability", ylab="")
```

    ## Warning:   Note that you are attempting to plot a 'main effect' in a model that contains an
    ##   interaction.  This is potentially misleading; you may wish to consider using the 'by'
    ##   argument.

    ## Conditions used in construction of plot
    ## log_abs_DD: -4.206996
    ## SignDD: NegativeDD
    ## smap_coef_mean: 0.05984084
    ## Function: Prey
    ## Recipient_Donor: Det1_Mp1
    ## Year: 2018

``` r
p3<- plot(visreg(Vulnerability.Ta_PC.TMB,  "log_abs_DD", by="SignDD", re.form=NA, plot=FALSE), band=TRUE, partial=TRUE, gg=TRUE, overlay=TRUE, 
          xlab="ln(Interaction density dependence)", ylab="ln(Recipient sensitivity + 1)") + theme(legend.position="none") + ggtitle("Model for H vs. C")
p4<- plot(visreg(Vulnerability.Ta_PC.TMB, "Function", re.form=NA, plot=FALSE), band=TRUE, partial=TRUE, gg=TRUE, overlay=TRUE, 
          xlab="Recipient function", ylab="")
```

    ## Warning:   Note that you are attempting to plot a 'main effect' in a model that contains an
    ##   interaction.  This is potentially misleading; you may wish to consider using the 'by'
    ##   argument.

    ## Conditions used in construction of plot
    ## log_abs_DD: -4.206996
    ## SignDD: NegativeDD
    ## smap_coef_SD: 0.03164694
    ## smap_coef_mean: 0.05984084
    ## Recipient_Donor: Det1_Mp1
    ## Year: 2018

``` r
(fig.pent <- p3+p2+p1+p4 + plot_layout(nrow=1, guides="keep"))
```

![](6_Effects_of_interaction_variability_on_population_sensitivity_files/figure-markdown_github/multiple%20regression-2.png)

``` r
## Joint

Vulnerability.Ta_JC.TMB <- glmmTMB(log(abs_R_Abundance_response_nRR_0.1+1) ~ 
                                     (log_abs_DD)*SignDD + smap_coef_SD + smap_coef_mean + Function +
                                     (1|Recipient_Donor) + (1|Year), REML=TRUE, 
                                   subset(smapcoef_TrEff_SignDD_, Treatment_Cont=="Joint_Cont"), 
                                   control=glmmTMBControl(optimizer=optim,
                                                          optArgs=list(method="BFGS")))
p1<- plot(visreg(Vulnerability.Ta_JC.TMB, "smap_coef_mean", re.form=NA, plot=FALSE), band=TRUE, partial=TRUE, gg=TRUE, overlay=TRUE, 
          xlab="Mean interaction strength", ylab="")
```

    ## Warning:   Note that you are attempting to plot a 'main effect' in a model that contains an
    ##   interaction.  This is potentially misleading; you may wish to consider using the 'by'
    ##   argument.

    ## Conditions used in construction of plot
    ## log_abs_DD: -4.206996
    ## SignDD: NegativeDD
    ## smap_coef_SD: 0.03164694
    ## Function: Prey
    ## Recipient_Donor: Det1_Mp1
    ## Year: 2018

``` r
p2<- plot(visreg(Vulnerability.Ta_JC.TMB, "smap_coef_SD", re.form=NA, plot=FALSE), band=TRUE, partial=TRUE, gg=TRUE, overlay=TRUE, 
          xlab="Interaction temporal variability", ylab="")
```

    ## Warning:   Note that you are attempting to plot a 'main effect' in a model that contains an
    ##   interaction.  This is potentially misleading; you may wish to consider using the 'by'
    ##   argument.

    ## Conditions used in construction of plot
    ## log_abs_DD: -4.206996
    ## SignDD: NegativeDD
    ## smap_coef_mean: 0.05984084
    ## Function: Prey
    ## Recipient_Donor: Det1_Mp1
    ## Year: 2018

``` r
p3<- plot(visreg(Vulnerability.Ta_JC.TMB,  "log_abs_DD", by="SignDD", re.form=NA, plot=FALSE), band=TRUE, partial=TRUE, gg=TRUE, overlay=TRUE, 
          xlab="ln(Interaction density dependence)", ylab="ln(Recipient sensitivity + 1)") + theme(legend.position="none") + ggtitle("Model for I+H vs. C")
p4<- plot(visreg(Vulnerability.Ta_JC.TMB, "Function", re.form=NA, plot=FALSE), band=TRUE, partial=TRUE, gg=TRUE, overlay=TRUE, 
          xlab="Recipient function", ylab="")
```

    ## Warning:   Note that you are attempting to plot a 'main effect' in a model that contains an
    ##   interaction.  This is potentially misleading; you may wish to consider using the 'by'
    ##   argument.

    ## Conditions used in construction of plot
    ## log_abs_DD: -4.206996
    ## SignDD: NegativeDD
    ## smap_coef_SD: 0.03164694
    ## smap_coef_mean: 0.05984084
    ## Recipient_Donor: Det1_Mp1
    ## Year: 2018

``` r
(fig.joint <- p3+p2+p1+p4 + plot_layout(nrow=1, guides="keep"))
```

![](6_Effects_of_interaction_variability_on_population_sensitivity_files/figure-markdown_github/multiple%20regression-3.png)

### Generate Fig. 4

``` r
theme_set(theme_test())
fig5 <- (fig.fipro / fig.pent / fig.joint) + plot_annotation(tag_levels="a") & 
  theme(plot.tag.position=c(0.06, 0.95), plot.tag=element_text(face="bold"), plot.margin=unit(c(8, 1.5, 5.5, 8), "pt"))
#windows(14, 9.5, rescale="fixed")
fig5
```

![](6_Effects_of_interaction_variability_on_population_sensitivity_files/figure-markdown_github/fig%20assemble-1.png)

``` r
ggsave("./figs/fig4.ggplot2.pdf", width=14, height=9, device=cairo_pdf, unit="in")
```

## Correlations between mean interaction strength and interaction temporal variability

### Generate Fig. S3

``` r
p1 <- 
  ggplot(subset(smapcoef_TrEff_SignDD, Treatment_Cont=="Fipro_Cont"), 
         aes(x=abs(smap_coef_mean), y=smap_coef_SD)) +
  geom_point() + geom_smooth(method="lm") + 
  scale_y_continuous(limits=c(0, 0.15)) +
  xlab("Mean interaction strength") + ylab("Interaction temporal variability")
summary(lm(smap_coef_SD ~ abs(smap_coef_mean), subset(smapcoef_TrEff_SignDD, Treatment_Cont=="Fipro_Cont")))
```

    ## 
    ## Call:
    ## lm(formula = smap_coef_SD ~ abs(smap_coef_mean), data = subset(smapcoef_TrEff_SignDD, 
    ##     Treatment_Cont == "Fipro_Cont"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.05353 -0.04029 -0.01915  0.03776  0.08122 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          0.057727   0.007566   7.630 1.12e-10 ***
    ## abs(smap_coef_mean) -0.056861   0.046539  -1.222    0.226    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.04489 on 67 degrees of freedom
    ## Multiple R-squared:  0.02179,    Adjusted R-squared:  0.007194 
    ## F-statistic: 1.493 on 1 and 67 DF,  p-value: 0.2261

``` r
p2 <- 
  ggplot(subset(smapcoef_TrEff_SignDD, Treatment_Cont=="Fipro_Cont"), 
         aes(x=abs(smap_coef_mean), y=smap_coef_SD, color=SignDD)) +
  geom_point() + geom_smooth(method="lm") + 
  scale_y_continuous(limits=c(0, 0.15)) +
  xlab("Mean interaction strength") + ylab("Interaction temporal variability") +
  scale_color_discrete(name="Direction of IDD", labels=c("Negative IDD", "Positive IDD"))
summary(lm(smap_coef_SD ~ abs(smap_coef_mean)*SignDD, subset(smapcoef_TrEff_SignDD, Treatment_Cont=="Fipro_Cont")))
```

    ## 
    ## Call:
    ## lm(formula = smap_coef_SD ~ abs(smap_coef_mean) * SignDD, data = subset(smapcoef_TrEff_SignDD, 
    ##     Treatment_Cont == "Fipro_Cont"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.05458 -0.03335 -0.01236  0.03320  0.08967 
    ## 
    ## Coefficients:
    ##                                       Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                           0.046211   0.009406   4.913 6.37e-06 ***
    ## abs(smap_coef_mean)                  -0.029890   0.049684  -0.602    0.550    
    ## SignDDPositiveDD                      0.015822   0.016458   0.961    0.340    
    ## abs(smap_coef_mean):SignDDPositiveDD  0.238015   0.198985   1.196    0.236    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.04325 on 65 degrees of freedom
    ## Multiple R-squared:  0.1191, Adjusted R-squared:  0.07843 
    ## F-statistic: 2.929 on 3 and 65 DF,  p-value: 0.04012

``` r
p1 + p2 + plot_annotation(tag_levels="a")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](6_Effects_of_interaction_variability_on_population_sensitivity_files/figure-markdown_github/IS%20DD,%20IS%20intrinsic%20variability,%20mean%20IS%202-1.png)

``` r
ggsave("./figs/figS3.ggplot2.pdf", width=9, height=4, device=cairo_pdf, unit="in")
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'
