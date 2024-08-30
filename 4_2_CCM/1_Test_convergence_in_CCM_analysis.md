# Test convergence in CCM analysis (2017-2019 composite)

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

### Data compile

Using unprocessed data (PlaData171819_corrected.csv,
MacrophyteData171819.csv, and
animaldata171819.completegathered.fixed.csv) instead of processed,
exported data (all_Ts240123.csv) to avoid roundoff errors and to
reproduce identical results from the manuscripts. The data compilation
script are copied from “Data_processsing.md”.

``` r
PlanktonData <- read.csv("../raw_data/PlaData171819_corrected.csv", header=TRUE)
MacrophyteData <- read.csv("../raw_data/MacrophyteData171819.csv", header=TRUE)
AnimalData <- read.csv("../raw_data/animaldata171819.completegathered.fixed.csv", header=TRUE)

# Data compilation

data_definefctr <- function(data) {
  data$Insecticide <- as.factor(data$Insecticide)
  data$Herbicide <- as.factor(data$Herbicide)
  data$Treatment <- factor(data$Treatment, levels=c("Cont", "Fipro", "Pent", "Joint")) #Define category order
  data$Tank <- factor(data$Tank, levels=c("Cont1", "Cont2", "Fipro1", "Fipro2", "Pent1", "Pent2", "Joint1", "Joint2"))
  data$Year_Tank <- with(data, paste(Year, Tank, sep="_")) #Define unit of analyses
  data$Year_Tank <- 
    factor(data$Year_Tank, 
           levels=c("2017_Cont1", "2017_Cont2", "2017_Fipro1", "2017_Fipro2", 
                    "2017_Pent1", "2017_Pent2", "2017_Joint1", "2017_Joint2", 
                    "2018_Cont1", "2018_Cont2", "2018_Fipro1", "2018_Fipro2", 
                    "2018_Pent1", "2018_Pent2", "2018_Joint1", "2018_Joint2", 
                    "2019_Cont1", "2019_Cont2", "2019_Fipro1", "2019_Fipro2", 
                    "2019_Pent1", "2019_Pent2", "2019_Joint1", "2019_Joint2")) #Define category order
  return(data)
}

PlanktonData <- data_definefctr(PlanktonData)
MacrophyteData <- data_definefctr(MacrophyteData)
MacrophyteData$SubPlot <- as.factor(MacrophyteData$SubPlot)
AnimalData <- data_definefctr(AnimalData)

# Split by functional groups

Phytopla <- 
  PlanktonData %>%
  filter(Group=="Diatom" | Group=="Flagellate" | Group=="GreenAlgae") %>% # Omit unnecessary taxa
  group_by(Year, Tank, Week, Year_Tank) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup(.) %>% filter(Week>1)
```

    ## `summarise()` regrouping output by 'Year', 'Tank', 'Week' (override with
    ## `.groups` argument)

``` r
Roti <- 
  PlanktonData %>%
  filter((Group=="Rotifer") & Species_each_file!="Nauplius_(E._japonicus)") %>%  #omit a different sampling protocol
  group_by(Year, Tank, Week, Year_Tank) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup(.) %>% filter(Week>1)
```

    ## `summarise()` regrouping output by 'Year', 'Tank', 'Week' (override with
    ## `.groups` argument)

``` r
#Miki Hirao counted all individuals in 1L sample, whereas Ji Cai took 1/10 subsample from 500mL sample and counted individuals in the 1/10 subsample
ZooplaCr <- 
  PlanktonData %>%
  filter((Group=="Cladocera" | Group=="Copepoda" | Group=="Ostracoda") & Species_each_file!="Nauplius_(E._japonicus)") %>%  #omit a different sampling protocol
  group_by(Year, Tank, Week, Year_Tank) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup(.) %>% filter(Week>1)
```

    ## `summarise()` regrouping output by 'Year', 'Tank', 'Week' (override with
    ## `.groups` argument)

``` r
Mp <- 
  MacrophyteData %>%
  group_by(Year, Tank, Week, SubPlot, Year_Tank) %>%
  summarise(Cover=sum(Cover)) %>%
  group_by(Year, Tank, Week, Year_Tank) %>%
  summarise(Abundance=mean(Cover)) %>%
  ungroup(.) %>% filter(Week>1)
```

    ## `summarise()` regrouping output by 'Year', 'Tank', 'Week', 'SubPlot' (override
    ## with `.groups` argument)`summarise()` regrouping output by 'Year', 'Tank',
    ## 'Week' (override with `.groups` argument)

``` r
Det <- 
  AnimalData %>%
  filter(Function=="Detritivore") %>%
  group_by(Year, Tank, Week, Year_Tank) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup(.) %>% filter(Week>1)
```

    ## `summarise()` regrouping output by 'Year', 'Tank', 'Week' (override with
    ## `.groups` argument)

``` r
Herb <- 
  AnimalData %>%
  filter(Function=="Herbivore") %>%
  group_by(Year, Tank, Week, Year_Tank) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup(.) %>% filter(Week>1)
```

    ## `summarise()` regrouping output by 'Year', 'Tank', 'Week' (override with
    ## `.groups` argument)

``` r
Pred <- 
  AnimalData %>%
  filter(Function=="Predator") %>%
  group_by(Year, Tank, Week, Year_Tank) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup(.) %>% filter(Week>1)
```

    ## `summarise()` regrouping output by 'Year', 'Tank', 'Week' (override with
    ## `.groups` argument)

``` r
Pred.P <- 
  AnimalData %>%
  filter(Function=="Predator" & Habitat=="Phytophilous") %>%
  group_by(Year, Tank, Week, Year_Tank) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup(.) %>% filter(Week>1)
```

    ## `summarise()` regrouping output by 'Year', 'Tank', 'Week' (override with
    ## `.groups` argument)

``` r
Pred.B <-   AnimalData %>%
  filter(Function=="Predator" & Habitat=="Benthic") %>%
  group_by(Year, Tank, Week, Year_Tank) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup(.) %>% filter(Week>1)
```

    ## `summarise()` regrouping output by 'Year', 'Tank', 'Week' (override with
    ## `.groups` argument)

``` r
Pred.S <-   AnimalData %>%
  filter(Function=="Predator" & Habitat=="Neustonic") %>%
  group_by(Year, Tank, Week, Year_Tank) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup(.) %>% filter(Week>1)
```

    ## `summarise()` regrouping output by 'Year', 'Tank', 'Week' (override with
    ## `.groups` argument)

``` r
Pred.C <-   AnimalData %>%
  filter(Function=="Predator" & Habitat=="Nektonic") %>%
  group_by(Year, Tank, Week, Year_Tank) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup(.) %>% filter(Week>1)
```

    ## `summarise()` regrouping output by 'Year', 'Tank', 'Week' (override with
    ## `.groups` argument)

``` r
Mollusca <- 
  AnimalData %>%
  filter(Species=="Bivalvia_sp" | Species=="Physa_acuta") %>%
  group_by(Year, Tank, Week, Year_Tank) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup(.) %>% filter(Week>1)
```

    ## `summarise()` regrouping output by 'Year', 'Tank', 'Week' (override with
    ## `.groups` argument)

``` r
Commulist <- list(Phytopla, Roti, ZooplaCr, Mp, Det, Herb, Pred.P, Pred.B, Pred.S, Pred.C, Mollusca)
names(Commulist) <- c("Phytopla", "Roti", "ZooplaCr", "Mp", "Det", "Herb", "Pred.P", "Pred.B", "Pred.S", "Pred.C", "Mollusca")

# Log-transformation and standardization

Tslist <- lapply(1:11, function(i) {
  Commulist[[i]] %>%
    select(Week, Year_Tank, Abundance) %>%
    spread(key=Year_Tank, value=Abundance) %>% # Remove unnecessary columns
    rbind(., rep(NaN, ncol(.))) %>% # Add NaN to separate different tanks, but not necessary if you correctly specify lib parameters
    gather(Year_Tank, Abundance, -Week) %>% # Gather data for EDM
    mutate(Year_Tank=factor(Year_Tank, levels=c("2017_Cont1", "2017_Cont2", "2017_Fipro1", "2017_Fipro2", 
                                                "2017_Pent1", "2017_Pent2", "2017_Joint1", "2017_Joint2", 
                                                "2018_Cont1", "2018_Cont2", "2018_Fipro1", "2018_Fipro2", 
                                                "2018_Pent1", "2018_Pent2", "2018_Joint1", "2018_Joint2", 
                                                "2019_Cont1", "2019_Cont2", "2019_Fipro1", "2019_Fipro2", 
                                                "2019_Pent1", "2019_Pent2", "2019_Joint1", "2019_Joint2")))
} 
)

names(Tslist) <- 
  c("Phytopla1", "Roti1", "Zoopla1", "Mp1", "Det1", "Herb1", "Pred.P1", "Pred.B1", "Pred.S1", "Pred.C1", "Moll1")

lapply(c("Phytopla1", "Roti1", "Det1", "Herb1", "Pred.P1", "Pred.B1", "Pred.S1", "Pred.C1", "Moll1"), 
       function(i) { 
         Tslist[[i]] <<- mutate(Tslist[[i]], scAbundance=as.numeric(scale(log(Abundance+1))))
       }
)
```

    ## [[1]]
    ## # A tibble: 264 x 4
    ##     Week Year_Tank  Abundance scAbundance
    ##    <dbl> <fct>          <dbl>       <dbl>
    ##  1     2 2017_Cont1       799      0.706 
    ##  2     4 2017_Cont1        86     -0.118 
    ##  3     6 2017_Cont1      1873      1.02  
    ##  4     8 2017_Cont1       187      0.168 
    ##  5    10 2017_Cont1       127      0.0258
    ##  6    12 2017_Cont1       305      0.349 
    ##  7    14 2017_Cont1      1001      0.790 
    ##  8    16 2017_Cont1      3170      1.22  
    ##  9    18 2017_Cont1       613      0.608 
    ## 10    20 2017_Cont1      3860      1.29  
    ## # ... with 254 more rows
    ## 
    ## [[2]]
    ## # A tibble: 264 x 4
    ##     Week Year_Tank  Abundance scAbundance
    ##    <dbl> <fct>          <dbl>       <dbl>
    ##  1     2 2017_Cont1         0      -1.19 
    ##  2     4 2017_Cont1        29       0.870
    ##  3     6 2017_Cont1        27       0.828
    ##  4     8 2017_Cont1       124       1.74 
    ##  5    10 2017_Cont1       129       1.76 
    ##  6    12 2017_Cont1        92       1.56 
    ##  7    14 2017_Cont1        75       1.43 
    ##  8    16 2017_Cont1        40       1.06 
    ##  9    18 2017_Cont1        29       0.870
    ## 10    20 2017_Cont1        40       1.06 
    ## # ... with 254 more rows
    ## 
    ## [[3]]
    ## # A tibble: 264 x 4
    ##     Week Year_Tank  Abundance scAbundance
    ##    <dbl> <fct>          <dbl>       <dbl>
    ##  1     2 2017_Cont1        60       1.27 
    ##  2     4 2017_Cont1        17       0.366
    ##  3     6 2017_Cont1        25       0.638
    ##  4     8 2017_Cont1         3      -0.745
    ##  5    10 2017_Cont1         5      -0.446
    ##  6    12 2017_Cont1         8      -0.146
    ##  7    14 2017_Cont1        24       0.609
    ##  8    16 2017_Cont1        36       0.898
    ##  9    18 2017_Cont1        13       0.180
    ## 10    20 2017_Cont1        13       0.180
    ## # ... with 254 more rows
    ## 
    ## [[4]]
    ## # A tibble: 264 x 4
    ##     Week Year_Tank  Abundance scAbundance
    ##    <dbl> <fct>          <dbl>       <dbl>
    ##  1     2 2017_Cont1         2     0.493  
    ##  2     4 2017_Cont1         0    -0.857  
    ##  3     6 2017_Cont1         3     0.847  
    ##  4     8 2017_Cont1         1    -0.00524
    ##  5    10 2017_Cont1         3     0.847  
    ##  6    12 2017_Cont1         1    -0.00524
    ##  7    14 2017_Cont1         2     0.493  
    ##  8    16 2017_Cont1         0    -0.857  
    ##  9    18 2017_Cont1         0    -0.857  
    ## 10    20 2017_Cont1         0    -0.857  
    ## # ... with 254 more rows
    ## 
    ## [[5]]
    ## # A tibble: 264 x 4
    ##     Week Year_Tank  Abundance scAbundance
    ##    <dbl> <fct>          <dbl>       <dbl>
    ##  1     2 2017_Cont1         0      -0.593
    ##  2     4 2017_Cont1         0      -0.593
    ##  3     6 2017_Cont1         5       1.81 
    ##  4     8 2017_Cont1         2       0.882
    ##  5    10 2017_Cont1        13       2.95 
    ##  6    12 2017_Cont1        11       2.74 
    ##  7    14 2017_Cont1         8       2.36 
    ##  8    16 2017_Cont1        14       3.04 
    ##  9    18 2017_Cont1         9       2.50 
    ## 10    20 2017_Cont1        12       2.85 
    ## # ... with 254 more rows
    ## 
    ## [[6]]
    ## # A tibble: 264 x 4
    ##     Week Year_Tank  Abundance scAbundance
    ##    <dbl> <fct>          <dbl>       <dbl>
    ##  1     2 2017_Cont1         2       0.198
    ##  2     4 2017_Cont1         1      -0.170
    ##  3     6 2017_Cont1         8       1.19 
    ##  4     8 2017_Cont1         1      -0.170
    ##  5    10 2017_Cont1         4       0.661
    ##  6    12 2017_Cont1         9       1.29 
    ##  7    14 2017_Cont1         6       0.966
    ##  8    16 2017_Cont1         7       1.09 
    ##  9    18 2017_Cont1         3       0.459
    ## 10    20 2017_Cont1         2       0.198
    ## # ... with 254 more rows
    ## 
    ## [[7]]
    ## # A tibble: 264 x 4
    ##     Week Year_Tank  Abundance scAbundance
    ##    <dbl> <fct>          <dbl>       <dbl>
    ##  1     2 2017_Cont1         0      -1.04 
    ##  2     4 2017_Cont1         3       0.547
    ##  3     6 2017_Cont1         0      -1.04 
    ##  4     8 2017_Cont1         0      -1.04 
    ##  5    10 2017_Cont1         2       0.218
    ##  6    12 2017_Cont1         0      -1.04 
    ##  7    14 2017_Cont1         3       0.547
    ##  8    16 2017_Cont1         4       0.803
    ##  9    18 2017_Cont1         0      -1.04 
    ## 10    20 2017_Cont1         0      -1.04 
    ## # ... with 254 more rows
    ## 
    ## [[8]]
    ## # A tibble: 264 x 4
    ##     Week Year_Tank  Abundance scAbundance
    ##    <dbl> <fct>          <dbl>       <dbl>
    ##  1     2 2017_Cont1         0      -0.263
    ##  2     4 2017_Cont1         0      -0.263
    ##  3     6 2017_Cont1         0      -0.263
    ##  4     8 2017_Cont1         0      -0.263
    ##  5    10 2017_Cont1         0      -0.263
    ##  6    12 2017_Cont1         0      -0.263
    ##  7    14 2017_Cont1         0      -0.263
    ##  8    16 2017_Cont1         0      -0.263
    ##  9    18 2017_Cont1         0      -0.263
    ## 10    20 2017_Cont1         0      -0.263
    ## # ... with 254 more rows
    ## 
    ## [[9]]
    ## # A tibble: 264 x 4
    ##     Week Year_Tank  Abundance scAbundance
    ##    <dbl> <fct>          <dbl>       <dbl>
    ##  1     2 2017_Cont1         0      -1.21 
    ##  2     4 2017_Cont1         0      -1.21 
    ##  3     6 2017_Cont1        33       0.766
    ##  4     8 2017_Cont1         0      -1.21 
    ##  5    10 2017_Cont1        31       0.732
    ##  6    12 2017_Cont1        22       0.546
    ##  7    14 2017_Cont1        15       0.342
    ##  8    16 2017_Cont1         0      -1.21 
    ##  9    18 2017_Cont1         0      -1.21 
    ## 10    20 2017_Cont1         4      -0.311
    ## # ... with 254 more rows

``` r
Tslist[["Mp1"]] <- mutate(Tslist[["Mp1"]], scAbundance=as.numeric(scale(Abundance)))

Tslist[["Zoopla1"]] <- 
  Tslist[["Zoopla1"]] %>%
  separate(Year_Tank, into=c("Year", "Tank")) %>%
  group_by(Year) %>%
  mutate(scAbundance=as.numeric(scale(log(Abundance*ifelse(Year==2017, 1, 20)+1)))) %>%
  unite(Year, Tank, col="Year_Tank")

all_Ts <- 
  lapply(1:11, function(i) {
    mutate(Tslist[[i]], Recipient=names(Tslist)[i])
  }) %>% do.call(rbind, .)

all_Ts <- #define necessary columns
  all_Ts %>%
  separate(Year_Tank, into=c("Year", "Tank"), sep="_") %>% #year and tank
  mutate(Treatment=substr(Tank, 1, nchar(Tank)-1)) %>% #treatment name
  mutate(Treatment=factor(Treatment, levels=c("Cont", "Fipro", "Pent", "Joint")), #order of treatment names
         Tank=factor(Tank, levels=unique(Tank)), #order of tanks
         Recipient=factor(Recipient, levels=c("Phytopla1", "Roti1", "Zoopla1", "Mp1", 
                                              "Det1", "Herb1", "Pred.P1", "Pred.B1", "Pred.S1", "Moll1"))) 

Phytopla1 <- subset(all_Ts, Recipient=="Phytopla1")
Roti1 <- subset(all_Ts, Recipient=="Roti1")
Zoopla1 <- subset(all_Ts, Recipient=="Zoopla1")
Mp1 <- subset(all_Ts, Recipient=="Mp1")
Det1 <- subset(all_Ts, Recipient=="Det1")
Herb1 <- subset(all_Ts, Recipient=="Herb1")
Pred.P1 <- subset(all_Ts, Recipient=="Pred.P1")
Pred.B1 <- subset(all_Ts, Recipient=="Pred.B1")
Pred.S1 <- subset(all_Ts, Recipient=="Pred.S1")
Moll1 <- subset(all_Ts, Recipient=="Moll1")
```

## Save best embedding dimensions and library snd set params

-   See “4_1_Determine_embedding_dimensions”

``` r
EPh <- 5; ERo <- 6; EZo <- 6; EMp <- 5; EDet <- 4; EHerb <- 6; EPP <- 5; EPB <- 2; EPS <- 5; EMoll <- 4
namedlist1 <- function(...){
  setNames(list(...), eval(substitute(alist(...))))
}
E_list <- namedlist1(EPh, ERo, EZo, EMp, EDet, EHerb, EPP, EPB, EPS, EMoll)

lib_custom <- matrix(c(c(1,10, 12,21, 23,32, 34,43, 45,54, 56,65, 67,76, 78,87), #define library
                     c(1,10, 12,21, 23,32, 34,43, 45,54, 56,65, 67,76, 78,87)+88, 
                     c(1,10, 12,21, 23,32, 34,43, 45,54, 56,65, 67,76, 78,87)+176), ncol=2, byrow=TRUE)
lib_custom_Mp <- matrix(c(c(1,10, 12,21, 23,32, 34,43, 45,54, 56,65, 67,76, 78,87), #define library of macrophye
                       c(1,10, 12,21, 23,32, 34,43, 45,54, 56,65, 67,76, 78,87)+88, 
                       c(2,10, 13,21, 24,32, 35,43, 46,54, 57,65, 68,76, 79,87)+176), ncol=2, byrow=TRUE)
No.samp <- 1000 # This takes time
```

## CCM analysis

### Define utility functions

``` r
ccm_conv_test <- function(effect, cause, E) {
  ccm_res <- ccm(cbind(effect, cause), E=E, lib=lib_custom, RNGseed=sample(1:1000000000, 1),
                 lib_sizes=seq(E+1, length(effect)-(E-1)*24-24+2, by=3), silent=TRUE, replace=FALSE, num_samples=No.samp)
  xmap_mean <- ccm_means(na.omit(ccm_res))
  xmap_95 <- ccm_means(na.omit(ccm_res), FUN=quantile, probs=0.95)
  xmap_9 <- ccm_means(na.omit(ccm_res), FUN=quantile, probs=0.9)
  xmap_ucl <- ccm_means(na.omit(ccm_res), FUN=quantile, probs=0.975)
  xmap_lcl <- ccm_means(na.omit(ccm_res), FUN=quantile, probs=0.025)
  xmap_conv <- subset(xmap_mean, lib_size==max(lib_size))$rho - subset(xmap_mean, lib_size==min(lib_size))$rho
  list(ccm_res, xmap_mean, xmap_95, xmap_9, xmap_ucl, xmap_lcl, xmap_conv)
}

ccm_conv_test_Mp <- function(effect, cause, E) {
  ccm_res <- ccm(cbind(effect, cause), E=E, lib=lib_custom_Mp, RNGseed=sample(1:1000000000, 1),
                 lib_sizes=seq(E+1, length(effect)-(E-1)*24-24+2, by=3), silent=TRUE, replace=FALSE, num_samples=No.samp)
  xmap_mean <- ccm_means(na.omit(ccm_res))
  xmap_95 <- ccm_means(na.omit(ccm_res), FUN=quantile, probs=0.95)
  xmap_9 <- ccm_means(na.omit(ccm_res), FUN=quantile, probs=0.9)
  xmap_ucl <- ccm_means(na.omit(ccm_res), FUN=quantile, probs=0.975)
  xmap_lcl <- ccm_means(na.omit(ccm_res), FUN=quantile, probs=0.025)
  xmap_conv <- subset(xmap_mean, lib_size==max(lib_size))$rho - subset(xmap_mean, lib_size==min(lib_size))$rho
  list(ccm_res, xmap_mean, xmap_95, xmap_9, xmap_ucl, xmap_lcl, xmap_conv)
}

ccm_vis <- function(ccm_res_list, main) {
  plot(rho ~ lib_size, ccm_res_list[[2]], type="l", main=main, 
       ylim=c((ccm_res_list[[2]][1, 9]+ccm_res_list[[6]][1, 9]/2), max(ccm_res_list[[5]][, 9])+0.2))
  lines(rho ~ lib_size, ccm_res_list[[5]], type="l", lty=3)
  lines(rho ~ lib_size, ccm_res_list[[6]], type="l", lty=3)
  legend("topright", 
         legend=paste("Diff. rho:", 
                      c(as.character(round(ccm_res_list[[7]], digits=3))), sep=" "), bty="n")
  abline(h=0, lty=2)
}
```

## Do CCM repeatedly

``` r
set.seed(1)
LIST <- NULL
for (i in 1:45) {
  combn.names <- combn(levels(all_Ts$Recipient), m=2)
  combn.E <- combn(names(E_list), m=2)
  dat1 <- get(combn.names[2,i])
  dat2 <- get(combn.names[1,i])
  E1 <- as.numeric(E_list[combn.E[2,i]])
  E2 <- as.numeric(E_list[combn.E[1,i]])
  
  effect <- dat1$scAbundance; cause <- dat2$scAbundance; E <- E1
  res_list1 <- ccm_conv_test(effect, cause, E)
  xmap.name <-paste(combn.names[2, i], combn.names[1, i], sep="_xmap_")
  write.csv(res_list1[[1]], paste("../4_2_CCM/ProcessedData/ccm_res_", xmap.name, ".csv", sep=""), row.names=FALSE)
  LIST <- c(LIST, list(res_list1))
  names(LIST)[length(LIST)] <- xmap.name

  effect <- dat2$scAbundance; cause <- dat1$scAbundance; E <- E2
  res_list2 <- ccm_conv_test(effect, cause, E)
  xmap.name <-paste(combn.names[1, i], combn.names[2, i], sep="_xmap_")
  write.csv(res_list2[[1]], paste("../4_2_CCM/ProcessedData/ccm_res_", xmap.name, ".csv", sep=""), row.names=FALSE)
  LIST <- c(LIST, list(res_list2))
  names(LIST)[length(LIST)] <- xmap.name
}
```

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

    ## Warning in (function (x, y = NULL, use = "everything", method = c("pearson", :
    ## 標準偏差が 0 です

``` r
grid.names <- subset(expand.grid(levels(all_Ts$Recipient), levels(all_Ts$Recipient)), Var1!=Var2)
xmap.names <- with(grid.names, paste(Var2, Var1, sep="_xmap_"))
LIST <- LIST[xmap.names]

# visualize an example
par(mfrow=c(1, 2), las=1, cex=1, family="sans")
ccm_vis(LIST["Mp1_xmap_Det1"][[1]], main="Macrophytes xmap Detritivores")
ccm_vis(LIST["Det1_xmap_Mp1"][[1]], main="Detritivores xmap Macrophytes")
```

![](1_Test_convergence_in_CCM_analysis_files/figure-markdown_github/ccm%20repeat-1.png)

## Summary of CCM results

``` r
z <- NULL
for (i in 1:90){
  z <- append(z, subset(LIST[[i]][[2]], lib_size==max(lib_size))$rho)
}
maxrho_mat <- matrix(c(as.vector(rbind(NA, matrix(z, ncol=9))), NA), ncol=10)
rownames(maxrho_mat) <- c("c_Ph", "c_Ro", "c_Zo", "c_Mp", "c_Det", "c_Herb", "c_PP", "c_PB", "c_PS", "c_Moll")
colnames(maxrho_mat) <- c("e_Ph", "e_Ro", "e_Zo", "e_Mp", "e_Det", "e_Herb", "e_PP", "e_PB", "e_PS", "e_Moll")
maxrho_mat <- t(maxrho_mat)
round(maxrho_mat, digits=4)
```

    ##           c_Ph   c_Ro    c_Zo    c_Mp   c_Det  c_Herb    c_PP   c_PB    c_PS
    ## e_Ph        NA 0.4653  0.0947 -0.1885  0.1559 -0.0662  0.2384 0.2815  0.0144
    ## e_Ro    0.3505     NA  0.3400 -0.2243 -0.0326  0.2721  0.3865 0.3408 -0.0484
    ## e_Zo   -0.1309 0.3481      NA  0.1517  0.0171  0.3089  0.2752 0.1850  0.2804
    ## e_Mp    0.4474 0.3341  0.3586      NA  0.4969  0.4618  0.2334 0.3835  0.4213
    ## e_Det   0.2716 0.1362  0.1572  0.3380      NA  0.1393 -0.1294 0.0754  0.0870
    ## e_Herb  0.2297 0.0418  0.0275  0.3866  0.2955      NA  0.3679 0.5044  0.0034
    ## e_PP    0.2491 0.2320  0.1560 -0.1226 -0.0136 -0.0291      NA 0.3725  0.0022
    ## e_PB    0.0678 0.0211 -0.0116 -0.0087 -0.0467  0.1169  0.2265     NA -0.1347
    ## e_PS    0.0756 0.1243  0.0873  0.4369  0.0339  0.0887  0.0415 0.1442      NA
    ## e_Moll  0.4613 0.3783  0.1276  0.1266  0.0901  0.1825  0.2585 0.3819  0.2742
    ##         c_Moll
    ## e_Ph    0.1952
    ## e_Ro    0.3221
    ## e_Zo   -0.0716
    ## e_Mp    0.2001
    ## e_Det   0.0301
    ## e_Herb  0.1476
    ## e_PP    0.2228
    ## e_PB    0.0916
    ## e_PS   -0.0495
    ## e_Moll      NA

``` r
z <- NULL
for (i in 1:90){
  z <- append(z, LIST[[i]][[7]])
}
conv_mat <- matrix(c(as.vector(rbind(NA, matrix(z, ncol=9))), NA), ncol=10)
rownames(conv_mat) <- c("c_Ph", "c_Ro", "c_Zo", "c_Mp", "c_Det", "c_Herb", "c_PP", "c_PB", "c_PS", "c_Moll")
colnames(conv_mat) <- c("e_Ph", "e_Ro", "e_Zo", "e_Mp", "e_Det", "e_Herb", "e_PP", "e_PB", "e_PS", "e_Moll")
conv_mat <- t(conv_mat)
round(conv_mat, digits=4)
```

    ##           c_Ph   c_Ro   c_Zo    c_Mp   c_Det  c_Herb    c_PP   c_PB    c_PS
    ## e_Ph        NA 0.2767 0.0642 -0.1277  0.1580 -0.0383  0.1295 0.1261  0.0133
    ## e_Ro    0.1867     NA 0.2160 -0.1514  0.0056  0.1658  0.2320 0.2209  0.0021
    ## e_Zo   -0.1080 0.2945     NA  0.1716  0.0755  0.2486  0.2222 0.1936  0.2361
    ## e_Mp    0.3472 0.2825 0.2519      NA  0.2124  0.2655  0.1460 0.3126  0.2114
    ## e_Det   0.2153 0.1020 0.1575  0.1419      NA  0.1280 -0.1186 0.0794  0.0723
    ## e_Herb  0.1406 0.0327 0.0830  0.2041  0.2085      NA  0.1618 0.3325  0.0428
    ## e_PP   -0.0111 0.2079 0.1737 -0.0779  0.0315  0.0150      NA 0.0485  0.0263
    ## e_PB    0.0127 0.0261 0.0127  0.0028 -0.0181  0.1171 -0.0035     NA -0.1252
    ## e_PS    0.0953 0.1241 0.1109  0.1592  0.0307  0.0634  0.0623 0.1389      NA
    ## e_Moll  0.1106 0.1824 0.1129  0.1036  0.1075  0.1321  0.1664 0.1731  0.2664
    ##         c_Moll
    ## e_Ph    0.0659
    ## e_Ro    0.2446
    ## e_Zo   -0.0538
    ## e_Mp    0.1235
    ## e_Det   0.0139
    ## e_Herb  0.1586
    ## e_PP    0.0231
    ## e_PB    0.0485
    ## e_PS   -0.0357
    ## e_Moll      NA

``` r
write.csv(maxrho_mat, "../4_2_CCM/ccm_tables/maxrho_mat.csv")
write.csv(conv_mat, "../4_2_CCM/ccm_tables/conv_mat.csv")
```
