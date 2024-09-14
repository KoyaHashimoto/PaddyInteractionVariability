# Calculation of population sensitiity

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
    ## other attached packages:
    ## [1] rmarkdown_2.5      ggsci_2.9          patchwork_1.1.1    tidyr_1.1.2       
    ## [5] dplyr_1.0.2        ggplot2_3.3.2      RColorBrewer_1.1-2
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] pillar_1.4.6     compiler_3.6.3   tools_3.6.3      digest_0.6.27    evaluate_0.14   
    ##  [6] lifecycle_1.0.3  tibble_3.0.4     gtable_0.3.0     pkgconfig_2.0.3  rlang_1.1.0     
    ## [11] cli_3.4.1        rstudioapi_0.11  yaml_2.2.1       xfun_0.19        fastmap_1.1.0   
    ## [16] stringr_1.4.0    withr_2.5.0      knitr_1.30       generics_0.1.0   vctrs_0.6.1     
    ## [21] grid_3.6.3       tidyselect_1.1.0 glue_1.4.2       R6_2.5.0         fansi_1.0.3     
    ## [26] farver_2.0.3     purrr_0.3.4      magrittr_2.0.3   scales_1.1.1     ellipsis_0.3.2  
    ## [31] htmltools_0.5.2  colorspace_2.0-3 labeling_0.4.2   utf8_1.1.4       stringi_1.4.6   
    ## [36] munsell_0.5.0    crayon_1.3.4

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
library(ggsci); packageVersion("ggsci") #2.9
```

    ## [1] '2.9'

### Data compilation

``` r
# Load standardized time series

all_Ts <- read.csv("./processed_data/all_Ts240123.csv") 

all_Ts <- #add required columns
  all_Ts %>%
  separate(Year_Tank, into=c("Year", "Tank"), sep="_") %>% #year and tank
  mutate(Treatment=substr(Tank, 1, nchar(Tank)-1)) %>% #treatment name
  mutate(Treatment=factor(Treatment, levels=c("Cont", "Fipro", "Pent", "Joint")), #order of treatments
         Tank=factor(Tank, levels=unique(Tank)), #order of tanks
         Recipient=factor(Recipient, levels=c("Phytopla1", "Roti1", "Zoopla1", "Mp1", 
                                              "Det1", "Herb1", "Pred.P1", "Pred.B1", "Pred.S1", "Moll1"))) #order of recipient populations
```

## Calculation of population sensitivity (with Fig. 2b)

``` r
# Averaging density over replicates

all_Ts_smrz <- 
  all_Ts %>%
  group_by(Treatment, Recipient, Year) %>%
  summarise(Abundance_m=mean(Abundance, na.rm=TRUE)) #mean raw density
```

    ## `summarise()` regrouping output by 'Treatment', 'Recipient' (override with `.groups`
    ## argument)

``` r
## Show data sheet

head(all_Ts_smrz)
```

    ## # A tibble: 6 x 4
    ## # Groups:   Treatment, Recipient [2]
    ##   Treatment Recipient Year  Abundance_m
    ##   <fct>     <fct>     <chr>       <dbl>
    ## 1 Cont      Phytopla1 2017        928. 
    ## 2 Cont      Phytopla1 2018       2978. 
    ## 3 Cont      Phytopla1 2019        304. 
    ## 4 Cont      Roti1     2017         40.4
    ## 5 Cont      Roti1     2018         38.2
    ## 6 Cont      Roti1     2019          4.6

``` r
# Calculating log response ratio

Abundance_TrEff_nRR_0.1 <-  # calculate log response ratio following the method from Mratinson and Raupp 2013
  all_Ts_smrz %>%
  ungroup(.) %>%
  spread(value=Abundance_m, key=Treatment) %>%
  mutate(Fipro_Cont=log((Fipro+0.1)/(Cont+0.1)), 
         Pent_Cont=log((Pent+0.1)/(Cont+0.1)),
         Joint_Cont=log((Joint+0.1)/(Cont+0.1))) %>%
  select(-Cont, -Fipro, -Pent, -Joint) %>%
  gather(key=Treatment_Cont, value=Abundance_response_nRR_0.1, -(1:2)) %>%
  mutate(Treatment_Cont=factor(Treatment_Cont, levels=c("Fipro_Cont", "Pent_Cont", "Joint_Cont")))

Abundance_TrEff <- Abundance_TrEff_nRR_0.1
Abundance_TrEff <- filter(Abundance_TrEff, Recipient!="Pred.C1") # omit nektonic predators
Abundance_TrEff$abs_Abundance_response_nRR_0.1 <- abs(Abundance_TrEff$Abundance_response_nRR_0.1) #Absolute values

Abundance_TrEff$Recipient <- factor(Abundance_TrEff$Recipient, 
                                    levels=c("Phytopla1", "Roti1", "Zoopla1", "Mp1", 
                                             "Det1", "Herb1", "Pred.P1", "Pred.B1", "Pred.S1", "Moll1"))

## Show data sheet

head(Abundance_TrEff)
```

    ## # A tibble: 6 x 5
    ##   Recipient Year  Treatment_Cont Abundance_response_nRR_0.1 abs_Abundance_response_nRR_0.1
    ##   <fct>     <chr> <fct>                               <dbl>                          <dbl>
    ## 1 Phytopla1 2017  Fipro_Cont                         0.162                          0.162 
    ## 2 Phytopla1 2018  Fipro_Cont                         0.0631                         0.0631
    ## 3 Phytopla1 2019  Fipro_Cont                         2.72                           2.72  
    ## 4 Roti1     2017  Fipro_Cont                         0.320                          0.320 
    ## 5 Roti1     2018  Fipro_Cont                        -0.547                          0.547 
    ## 6 Roti1     2019  Fipro_Cont                        -0.631                          0.631

``` r
# Visualize sensitivity

smrz_Abundance_TrEff <- 
  Abundance_TrEff %>%
  group_by(Recipient, Treatment_Cont) %>%
  summarise(abs_Abundance_response_nRR_0.1=(exp(mean(log(abs_Abundance_response_nRR_0.1+1)))-1)) %>% #calculating geometric mean
  ungroup(.)
```

    ## `summarise()` regrouping output by 'Recipient' (override with `.groups` argument)

``` r
## Split by treatments

theme_custom <- function() {
  theme_bw() + 
    theme(axis.text=element_text(colour="black"), 
          axis.ticks=element_line(colour="black"),
          panel.border=element_rect(fill=NA, 
                                    colour="black"), 
          panel.spacing=unit(3, "pt"),
          strip.text=element_text(colour="black", lineheight=0.7))
}

palSet1 <- brewer.pal(9, "Set1")
reaction.norm <-
  smrz_Abundance_TrEff %>%
  #  mutate(Treatment_Cont=factor(Treatment_Cont, levels=c("Fipro_Cont", "Joint_Cont", "Pent_Cont"))) %>%
  ggplot(aes(x=Recipient, y=log(abs_Abundance_response_nRR_0.1+1))) +
  geom_point(aes(group=Treatment_Cont, color=Treatment_Cont, shape=Treatment_Cont), size=4) +
  geom_line(aes(group=Treatment_Cont, color=Treatment_Cont)) +
  geom_point(Abundance_TrEff, 
             mapping=aes(x=Recipient, y=log(abs_Abundance_response_nRR_0.1+1), 
                         group=Treatment_Cont, 
                         color=Treatment_Cont, 
                         shape=Treatment_Cont), 
             alpha=0.4, size=2.5) +
  scale_color_manual(values=palSet1[c(1, 2, 3)], labels=c("I vs. C", "H vs. C", "I+H vs. C")) + 
  scale_shape_manual(values=c(16, 17, 15), labels=c("I vs. C", "H vs. C", "I+H vs. C")) + 
  scale_x_discrete(name="Community member", labels=c("Phyt", "Roti", "C.zoop", "Macr", "Detr", "Herb", "P.pred", "B.pred", "N.pred", "Moll")) + 
  guides(color=guide_legend(title=NULL), shape=guide_legend(title=NULL)) +
  ylab("ln(Population sensitivity + 1)") + 
  theme_custom()

reaction.norm #Fig. 2b
```

![](3_Calculation_of_population_sensitivity_files/figure-markdown_github/calc%20sensitiivty-1.png)

``` r
write.csv(Abundance_TrEff, "./processed_data/Abundance_TrEff.csv", row.names=FALSE) # export sensitivity data
```

## Generate Fig. 2

``` r
# Visualize data

labeli <- as_labeller(c(`Phytopla1`="Phytoplankton", #map facet labels
                        `Roti1`="Rotifers",
                        `Zoopla1`="Crustacean\nzooplankton",
                        `Mp1`="Macrophytes",
                        `Det1`="Detritivores",
                        `Herb1`="Herbivores",
                        `Pred.P1`="Phytophilous\npredators",
                        `Pred.B1`="Benthic\npredators",
                        `Pred.S1`="Neustonic\npredators",
                        `Moll1`="Molluscs"))
lab <- data.frame(lab="*", 
                  Recipient=factor(c("Phytopla1", "Mp1", "Mp1", 
                                 "Herb1", "Pred.P1", "Pred.P1", 
                                 "Pred.B1", "Pred.B1", "Moll1"), 
                               levels=c("Phytopla1", "Roti1", "Zoopla1", "Mp1", 
                                        "Det1", "Herb1", "Pred.P1", "Pred.B1", "Pred.S1", "Moll1")))

palSet1 <- brewer.pal(9, "Set1")


set.seed(1)
density.boxplot_gg <- 
  ggplot(na.omit(filter(all_Ts, Recipient!="Pred.C1")), mapping=aes(x=Treatment, y=scAbundance)) +
  geom_boxplot(mapping=aes(fill=Treatment), outlier.shape=NA, size=0.3, colour="black") + 
  scale_fill_manual(values=palSet1[c(9, 1, 2, 3)]) + 
  geom_jitter(width=0.2, size=0.1) + 
  facet_wrap(~Recipient, nrow=2, labeller=labeli) + 
  scale_x_discrete(labels=c("C", "I", "H", "I+H")) +
  theme_custom() + #theme
  theme(strip.background=element_blank(), legend.position="none") +
  geom_text(data=lab, aes(label=lab), x=c(4, 3, 4, 4, 2, 4, 2, 4, 4), y=2.5, size=10) +
  labs(y="Standardized density") #軸ラベル
```

``` r
#windows(8, 7, rescale="fixed")
density.boxplot_gg / reaction.norm + plot_layout(height=c(5, 3)) + plot_annotation(tag_levels="a") & theme(plot.tag=element_text(face="bold"))
```

![](3_Calculation_of_population_sensitivity_files/figure-markdown_github/paper%20fig-1.png)

``` r
ggsave("./figs/fig2.ggplot2.pdf", width=8, height=7, device=cairo_pdf, unit="in")
```
