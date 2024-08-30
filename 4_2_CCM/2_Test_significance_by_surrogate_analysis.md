# Test significance by surrogate analysis

## Set up

### Load required packages

``` r
sessionInfo(); #R v.3.6.3
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
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_3.6.3  magrittr_2.0.3  fastmap_1.1.0   cli_3.4.1      
    ##  [5] tools_3.6.3     htmltools_0.5.2 rstudioapi_0.11 yaml_2.2.1     
    ##  [9] stringi_1.4.6   rmarkdown_2.5   knitr_1.30      stringr_1.4.0  
    ## [13] xfun_0.19       digest_0.6.27   rlang_1.1.0     evaluate_0.14

``` r
library(rEDM); packageVersion("rEDM") #v.0.7.5
```

    ## [1] '0.7.5'

``` r
library(tidyr); packageVersion("tidyr") #v.1.1.2
```

    ## [1] '1.1.2'

``` r
library(dplyr); packageVersion("dplyr") #v.1.0.2
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## [1] '1.0.2'

### Data compilation

``` r
all_Ts <- read.csv("../processed_data/all_Ts240123.csv") 

all_Ts <- 
  all_Ts %>%
  separate(Year_Tank, into=c("Year", "Tank"), sep="_") %>% #year and tank
  mutate(Treatment=substr(Tank, 1, nchar(Tank)-1)) %>% #treatment names
  mutate(Treatment=factor(Treatment, levels=c("Cont", "Fipro", "Pent", "Joint")), #order of treatment names
         Tank=factor(Tank, levels=unique(Tank)), #order of tank names
         Recipient=factor(Recipient, levels=c("Phytopla1", "Roti1", "Zoopla1", "Mp1", 
                                              "Det1", "Herb1", "Pred.P1", "Pred.B1", "Pred.S1", "Moll1"))) 

all_Ts <- as_tibble(all_Ts)
all_Ts$Recipient <- factor(all_Ts$Recipient, levels=unique(all_Ts$Recipient))
Tslist <- split(all_Ts, all_Ts$Recipient)
```

### Params setting

``` r
EPh <- 5; ERo <- 6; EZo <- 6; EMp <- 5; EDet <- 4; EHerb <- 6; EPP <- 5; EPB <- 2; EPS <- 5; EMoll <- 4
E_list <- c(EPh, ERo, EZo, EMp, EDet, EHerb, EPP, EPB, EPS, EMoll)

lib_custom <- matrix(c(c(1,10, 12,21, 23,32, 34,43, 45,54, 56,65, 67,76, 78,87), #define library
                       c(1,10, 12,21, 23,32, 34,43, 45,54, 56,65, 67,76, 78,87)+88, 
                       c(1,10, 12,21, 23,32, 34,43, 45,54, 56,65, 67,76, 78,87)+176), ncol=2, byrow=TRUE)
lib_custom_Mp <- matrix(c(c(1,10, 12,21, 23,32, 34,43, 45,54, 56,65, 67,76, 78,87), #define library of macrophye
                          c(1,10, 12,21, 23,32, 34,43, 45,54, 56,65, 67,76, 78,87)+88, 
                          c(2,10, 13,21, 24,32, 35,43, 46,54, 57,65, 68,76, 79,87)+176), ncol=2, byrow=TRUE)

No.surr <- 1000
No.samp <- 1 #Because only the rho values at the maximum library length are used, one time sampling is sufficient
```

### Define utility functions

``` r
### Make surrogate data

surr_funct <- function(ts, seed=1){
  
  set.seed(seed)
  
  surrogate_season_tss <- make_surrogate_seasonal(ts, num_surr=No.surr, T_period=10)
  surrogate_season_tss <- 
    sapply(seq_len(NCOL(surrogate_season_tss)), 
           function(i) {
             surrogate_ts <- surrogate_season_tss[, i]
             surrogate_ts_Repl <- data.frame(surrogate_ts=surrogate_ts, Repl=gl(24, 10), Week=seq(2, 20, 2))
             surrogate_ts_Repl <- 
               surrogate_ts_Repl %>% 
               spread(key=Repl, value=surrogate_ts) %>%
               rbind(., rep(NaN, ncol(.))) %>% 
               gather(Repl, surrogate_ts, -Week)
             surrogate_ts_Repl$surrogate_ts
           })
  
  return(surrogate_season_tss)
  
}

surrMp_funct <- function(tsMp, seed=1){
  
  make_surrogate_seasonal_custom <- function(ts, num_surr = 100, T_period = 12)
  {
    if (is.data.frame(ts))
    {
      ts <- ts[[1]]
    }
    
    n <- length(ts)
    I_season <- suppressWarnings(matrix(1:T_period, nrow = n, ncol = 1))
    
    lab.na <- which(is.na(ts))
    
    # Calculate seasonal cycle using smooth.spline
    seasonal_F <- smooth.spline(c((I_season - T_period)[-lab.na], I_season[-lab.na], 
                                  (I_season + T_period)[-lab.na]), 
                                c(ts[-lab.na], ts[-lab.na], ts[-lab.na]), )
    seasonal_cyc <- predict(seasonal_F, I_season[-lab.na])$y
    seasonal_resid <- ts[-lab.na] - seasonal_cyc
    
    if (length(lab.na)>0) {
      m <- matrix(unlist(
        lapply(seq(num_surr), function(i) {
          seasonal_cyc + sample(seasonal_resid, n-length(lab.na))
        })
      ), ncol = num_surr)
      M <- m
      for (i in lab.na) {
        M <- apply(M, MARGIN=2, function(x) append(x, NA, after=(i-1)))
      }
      M
    } else {
      matrix(unlist(
        lapply(seq(num_surr), function(i) {
          seasonal_cyc + sample(seasonal_resid, n)
        })
      ), ncol = num_surr)
    }
  }
  
  set.seed(seed)
  
  surrogate_season_Mp <- make_surrogate_seasonal_custom(tsMp, T_period=10, num_surr=No.surr)
  surrogate_season_Mp <- 
    sapply(seq_len(NCOL(surrogate_season_Mp)), 
           function(i) {
             surrogate_ts <- surrogate_season_Mp[, i]
             surrogate_ts_Repl <- data.frame(surrogate_ts=surrogate_ts, Repl=gl(24, 10), Week=seq(2, 20, 2))
             surrogate_ts_Repl <- 
               surrogate_ts_Repl %>% 
               spread(key=Repl, value=surrogate_ts) %>%
               rbind(., rep(NaN, ncol(.))) %>% 
               gather(Repl, surrogate_ts, -Week)
             surrogate_ts_Repl$surrogate_ts
           })
  
  return(surrogate_season_Mp)
}
```

### Define combinations

See the results of testing convergence in CCM

``` r
targetlist <- c("Phytopla1", "Roti1", "Zoopla1", "Mp1", "Det1", "Herb1", "Pred.P1", "Pred.B1", "Pred.S1", "Moll1")

#ph
heterosp_Ph <- data.frame(Roti1=Tslist[["Roti1"]]$scAbundance, 
                          Det1=Tslist[["Det1"]]$scAbundance, 
                          Pred.P1=Tslist[["Pred.P1"]]$scAbundance, 
                          Pred.B1=Tslist[["Pred.B1"]]$scAbundance)

#roti
heterosp_Roti <- data.frame(Phytopla1=Tslist[["Phytopla1"]]$scAbundance, 
                            Zoopla1=Tslist[["Zoopla1"]]$scAbundance, 
                            Herb1=Tslist[["Herb1"]]$scAbundance, 
                            Pred.P1=Tslist[["Pred.P1"]]$scAbundance, 
                            Pred.B1=Tslist[["Pred.B1"]]$scAbundance, 
                            Moll1=Tslist[["Moll1"]]$scAbundance)

#zo
heterosp_Zo <- data.frame(Roti1=Tslist[["Roti1"]]$scAbundance, 
                          Mp1=Tslist[["Mp1"]]$scAbundance, 
                          Herb1=Tslist[["Herb1"]]$scAbundance, 
                          Pred.P1=Tslist[["Pred.P1"]]$scAbundance, 
                          Pred.B1=Tslist[["Pred.B1"]]$scAbundance, 
                          Pred.S1=Tslist[["Pred.S1"]]$scAbundance)

#mp
heterosp_Mp <- data.frame(Phytopla1=Tslist[["Phytopla1"]]$scAbundance, 
                          Roti1=Tslist[["Roti1"]]$scAbundance, 
                          Zoopla1=Tslist[["Zoopla1"]]$scAbundance, 
                          Det1=Tslist[["Det1"]]$scAbundance, 
                          Herb1=Tslist[["Herb1"]]$scAbundance, 
                          Pred.P1=Tslist[["Pred.P1"]]$scAbundance, 
                          Pred.B1=Tslist[["Pred.B1"]]$scAbundance, 
                          Pred.S1=Tslist[["Pred.S1"]]$scAbundance, 
                          Moll1=Tslist[["Moll1"]]$scAbundance)

#det
heterosp_Det <- data.frame(Phytopla1=Tslist[["Phytopla1"]]$scAbundance, 
                           Roti1=Tslist[["Roti1"]]$scAbundance, 
                           Zoopla1=Tslist[["Zoopla1"]]$scAbundance, 
                           Mp1=Tslist[["Mp1"]]$scAbundance, 
                           Herb1=Tslist[["Herb1"]]$scAbundance)

#herb
heterosp_Herb <- data.frame(Phytopla1=Tslist[["Phytopla1"]]$scAbundance, 
                            Mp1=Tslist[["Mp1"]]$scAbundance, 
                            Det1=Tslist[["Det1"]]$scAbundance, 
                            Pred.P1=Tslist[["Pred.P1"]]$scAbundance, 
                            Pred.B1=Tslist[["Pred.B1"]]$scAbundance, 
                            Moll1=Tslist[["Moll1"]]$scAbundance)

#pp
heterosp_PP <- data.frame(Roti1=Tslist[["Roti1"]]$scAbundance, 
                          Zoopla1=Tslist[["Zoopla1"]]$scAbundance)

#pb
heterosp_PB <- data.frame(Herb1=Tslist[["Herb1"]]$scAbundance)

#ps
heterosp_PS <- data.frame(Roti1=Tslist[["Roti1"]]$scAbundance, 
                          Zoopla1=Tslist[["Zoopla1"]]$scAbundance, 
                          Mp1=Tslist[["Mp1"]]$scAbundance, 
                          Pred.B1=Tslist[["Pred.B1"]]$scAbundance)

#moll
heterosp_Moll <- data.frame(Phytopla1=Tslist[["Phytopla1"]]$scAbundance, 
                            Roti1=Tslist[["Roti1"]]$scAbundance, 
                            Zoopla1=Tslist[["Zoopla1"]]$scAbundance, 
                            Mp1=Tslist[["Mp1"]]$scAbundance, 
                            Det1=Tslist[["Det1"]]$scAbundance, 
                            Herb1=Tslist[["Herb1"]]$scAbundance, 
                            Pred.P1=Tslist[["Pred.P1"]]$scAbundance, 
                            Pred.B1=Tslist[["Pred.B1"]]$scAbundance, 
                            Pred.S1=Tslist[["Pred.S1"]]$scAbundance)

heterosp_list <- 
  list("heterosp_Ph"=heterosp_Ph, 
       "heterosp_Roit"=heterosp_Roti, 
       "heterosp_Zo"=heterosp_Zo, 
       "heterosp_Mp"=heterosp_Mp, 
       "heterosp_Det"=heterosp_Det, 
       "heterosp_Herb"=heterosp_Herb, 
       "heterosp_PP"=heterosp_PP, 
       "heterosp_PB"=heterosp_PB, 
       "heterosp_PS"=heterosp_PS, 
       "heterosp_Moll"=heterosp_Moll
  )
```

## Do surrogate analysis in parallel

``` r
surr_list <- lapply(1:10, function(j) {
  if (j==4) {
    surrMp_funct(ts=Tslist[[targetlist[j]]]$scAbundance[-seq(11, 264, 11)], seed=1)
  } else {
    surr_funct(ts=na.omit(Tslist[[targetlist[j]]]$scAbundance), seed=1)
  }
}
)


surr_means_res_list <- lapply(1:10, function(j) {
  if (j==4) {
    lapply(seq_len(NCOL(heterosp_list[[j]])), 
           function(k) {
             set.seed(k)
             surrogate_output <- 
               lapply(seq_len(NCOL(surr_list[[j]])), 
                      function(i) {
                        surrogate_ts <- surr_list[[j]][, i]
                        xmap_output <- ccm(cbind(surrogate_ts, heterosp_list[[j]][k]), E=E_list[[j]], lib=lib_custom_Mp, RNGseed=sample(1:1000000000, 1),
                                           lib_sizes=length(Tslist[[targetlist[j]]]$scAbundance)-(E_list[[j]]-1)*24-24, 
                                           silent=TRUE, replace=FALSE, num_samples=No.samp)
                        xmap_means <- ccm_means(xmap_output)
                      })
             surrogate_means <- do.call(rbind, surrogate_output)
             return(surrogate_means)
           })
  } else {
    lapply(seq_len(NCOL(heterosp_list[[j]])), 
           function(k) {
             set.seed(k)
             surrogate_output <- 
               lapply(seq_len(NCOL(surr_list[[j]])), 
                      function(i) {
                        surrogate_ts <- surr_list[[j]][, i]
                        xmap_output <- ccm(cbind(surrogate_ts, heterosp_list[[j]][k]), E=E_list[[j]], lib=lib_custom, RNGseed=sample(1:1000000000, 1),
                                           lib_sizes=length(Tslist[[targetlist[j]]]$scAbundance)-(E_list[[j]]-1)*24-24, 
                                           silent=TRUE, replace=FALSE, num_samples=No.samp)
                        xmap_means <- ccm_means(xmap_output)
                      })
             surrogate_means <- do.call(rbind, surrogate_output)
             return(surrogate_means)
           })
    
  }
})
```

### Summarize by calculation of P values

``` r
grid.names <- subset(expand.grid(levels(all_Ts$Recipient), levels(all_Ts$Recipient)), Var1!=Var2)
xmap.names <- with(grid.names, paste(Var2, Var1, sep="_xmap_"))

name.list <- lapply(heterosp_list, colnames)
names(surr_means_res_list) <- levels(all_Ts$Recipient)
for (i in 1:10) {
  names(surr_means_res_list[[i]]) <- paste(levels(all_Ts$Recipient)[i], "xmap", name.list[[i]], sep="_")
}

surr_means_res_list2 <- NULL
for (i in 1:10) {
  for (j in 1:length(surr_means_res_list[[i]])){
    surr_means_res_list2 <- c(surr_means_res_list2, list(surr_means_res_list[[i]][[j]]))
    names(surr_means_res_list2)[length(surr_means_res_list2)] <- names(surr_means_res_list[[i]])[j]
  }
}

MCpv <- NULL
for (i in 1:length(xmap.names)){
  ccm_res <- read.csv(paste("../4_2_CCM/ProcessedData/ccm_res_", xmap.names[i], ".csv", sep=""), header=TRUE)
  surrogate_means <- surr_means_res_list2[[xmap.names[i]]]
  if (is.null(surrogate_means)){
    MCpv <- append(MCpv, NA)
    } else {
    xmap <- ccm_means(ccm_res)
    MCpv <- append(MCpv, (sum(xmap[nrow(xmap), "rho"] < subset(surrogate_means, lib_size==max(lib_size))$rho)+1)/
                     (length(subset(surrogate_means, lib_size==max(lib_size))$rho)+1))
    }
  }
MCpv
```

    ##  [1] 0.000999001          NA          NA 0.073926074          NA 0.018981019
    ##  [7] 0.011988012          NA          NA 0.002997003 0.013986014          NA
    ## [13]          NA 0.008991009 0.004995005 0.007992008          NA 0.006993007
    ## [19]          NA 0.001998002 0.170829171          NA 0.001998002 0.014985015
    ## [25] 0.063936064 0.019980020          NA 0.000999001 0.001998002 0.022977023
    ## [31] 0.000999001 0.000999001 0.035964036 0.000999001 0.000999001 0.041958042
    ## [37] 0.006993007 0.094905095 0.183816184 0.003996004 0.061938062          NA
    ## [43]          NA          NA          NA 0.030969031          NA          NA
    ## [49] 0.003996004 0.007992008 0.002997003 0.000999001          NA 0.126873127
    ## [55]          NA 0.010989011 0.184815185          NA          NA          NA
    ## [61]          NA          NA          NA          NA          NA          NA
    ## [67]          NA          NA 0.127872128          NA          NA          NA
    ## [73]          NA 0.102897103 0.304695305 0.000999001          NA          NA
    ## [79]          NA 0.093906094          NA 0.000999001 0.000999001 0.090909091
    ## [85] 0.112887113 0.187812188 0.039960040 0.014985015 0.001998002 0.004995005

``` r
MCpv_mat <- matrix(c(as.vector(rbind(NA, matrix(MCpv, ncol=9))), NA), ncol=10)
rownames(MCpv_mat) <- c("c_Ph", "c_Ro", "c_Zo", "c_Mp", "c_Det", "c_Herb", "c_PP", "c_PB", "c_PS", "c_Moll")
colnames(MCpv_mat) <- c("e_Ph", "e_Ro", "e_Zo", "e_Mp", "e_Det", "e_Herb", "e_PP", "e_PB", "e_PS", "e_Moll")
MCpv_mat <- t(MCpv_mat)
round(MCpv_mat, digits=4)
```

    ##         c_Ph   c_Ro   c_Zo   c_Mp  c_Det c_Herb  c_PP   c_PB  c_PS c_Moll
    ## e_Ph      NA 0.0010     NA     NA 0.0739     NA 0.019 0.0120    NA     NA
    ## e_Ro   0.003     NA 0.0140     NA     NA 0.0090 0.005 0.0080    NA 0.0070
    ## e_Zo      NA 0.0020     NA 0.1708     NA 0.0020 0.015 0.0639 0.020     NA
    ## e_Mp   0.001 0.0020 0.0230     NA 0.0010 0.0010 0.036 0.0010 0.001 0.0420
    ## e_Det  0.007 0.0949 0.1838 0.0040     NA 0.0619    NA     NA    NA     NA
    ## e_Herb 0.031     NA     NA 0.0040 0.0080     NA 0.003 0.0010    NA 0.1269
    ## e_PP      NA 0.0110 0.1848     NA     NA     NA    NA     NA    NA     NA
    ## e_PB      NA     NA     NA     NA     NA 0.1279    NA     NA    NA     NA
    ## e_PS      NA 0.1029 0.3047 0.0010     NA     NA    NA 0.0939    NA     NA
    ## e_Moll 0.001 0.0010 0.0909 0.1129 0.1878 0.0400 0.015 0.0020 0.005     NA

``` r
write.csv(MCpv_mat, "../4_2_CCM/ccm_tables/MCpv_mat.csv")
```
