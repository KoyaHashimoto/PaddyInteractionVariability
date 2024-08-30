# Regularized S-map diagnostics

### Load required packages

``` r
sessionInfo() #save session information (R version 3.6.3 (2020))
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
library(ggplot2); packageVersion("ggplot2") #3.3.2
```

    ## [1] '3.3.2'

``` r
library(dplyr); packageVersion("dplyr") #1.0.2
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

``` r
library(tidyr); packageVersion("tidyr") #1.1.2
```

    ## [1] '1.1.2'

``` r
library(patchwork); packageVersion("patchwork") #1.1.1
```

    ## [1] '1.1.1'

### Load data

``` r
all_multsmap_model_outputs <- read.csv("./processed_data/all_regulsmap_model_outputs240123.csv")

head(all_multsmap_model_outputs)
```

    ##   time        obs        pred    target
    ## 1    1 -2.2187036          NA Phytopla1
    ## 2    2  3.0699223 -0.17494158 Phytopla1
    ## 3    3 -2.2993885 -0.69417259 Phytopla1
    ## 4    4 -0.3844117 -0.05691758 Phytopla1
    ## 5    5  0.8715548  0.63859515 Phytopla1
    ## 6    6  1.1861682  0.09632185 Phytopla1

``` r
all_multsmap_model_outputs$target <- 
  factor(all_multsmap_model_outputs$target, 
         levels=unique(all_multsmap_model_outputs$target))
```

## Prediction skills of S-map models

Relationships between predicted vs. observed values. Spearman’s
correlation coefficients rho is shown.

``` r
theme_set(theme_bw())

all_multsmap_model_outputs_list <- all_multsmap_model_outputs %>% split(., f=.$target) 

figs <- lapply(1:9, function(i) {
  ggplot(all_multsmap_model_outputs_list[[i]], aes(x=pred, y=obs)) + 
    geom_point() + geom_abline(slope=1, intercept=0) + 
    labs(title=levels(all_multsmap_model_outputs$target)[i], 
         subtitle=paste("rho=", 
                        round(cor(na.omit(all_multsmap_model_outputs_list[[i]])$obs, 
                                  na.omit(all_multsmap_model_outputs_list[[i]])$pred, method="spearman"), digits=2), sep=""))
})

predobsfigs <- figs[[1]]
for (i in 2:9) {
  predobsfigs <- predobsfigs + figs[[i]]
}
predobsfigs
```

![](4_3_2_Regul_Smap_diagnostics_files/figure-markdown_github/model%20prediction%20skills-1.png)

## Temporal patterns of prediction vs. observation values

Temporal patterns of predicted (red) and observed (black) values. When
prediction is dominated by temporal autocorrelation, shapes of red lines
are sometimes like just a “horizontal shift” from black lines.

``` r
figs <- lapply(1:9, function(i) {
  all_multsmap_model_outputs_list[[i]] %>% 
    mutate(timestep=1:nrow(all_multsmap_model_outputs_list[[i]])) %>%
    ggplot(aes(x=timestep, y=obs)) + geom_line() + 
    geom_line(aes(x=timestep, y=pred, color="red")) + 
    labs(title=levels(all_multsmap_model_outputs$target)[i])
})

for (i in 1:9) {
  print(figs[[i]])
}
```

![](4_3_2_Regul_Smap_diagnostics_files/figure-markdown_github/model%20prediction%20skills%20timeseries-1.png)![](4_3_2_Regul_Smap_diagnostics_files/figure-markdown_github/model%20prediction%20skills%20timeseries-2.png)![](4_3_2_Regul_Smap_diagnostics_files/figure-markdown_github/model%20prediction%20skills%20timeseries-3.png)![](4_3_2_Regul_Smap_diagnostics_files/figure-markdown_github/model%20prediction%20skills%20timeseries-4.png)![](4_3_2_Regul_Smap_diagnostics_files/figure-markdown_github/model%20prediction%20skills%20timeseries-5.png)![](4_3_2_Regul_Smap_diagnostics_files/figure-markdown_github/model%20prediction%20skills%20timeseries-6.png)![](4_3_2_Regul_Smap_diagnostics_files/figure-markdown_github/model%20prediction%20skills%20timeseries-7.png)![](4_3_2_Regul_Smap_diagnostics_files/figure-markdown_github/model%20prediction%20skills%20timeseries-8.png)![](4_3_2_Regul_Smap_diagnostics_files/figure-markdown_github/model%20prediction%20skills%20timeseries-9.png)
