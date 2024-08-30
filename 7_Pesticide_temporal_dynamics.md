# Pesticide temporal dynamics

## Load data

``` r
PW171819 <- read.csv("./raw_data/PW171819.csv", header=TRUE)
FW171819 <- read.csv("./raw_data/FW171819.csv", header=TRUE)
PS171819 <- read.csv("./raw_data/PS171819.csv", header=TRUE)
FS171819 <- read.csv("./raw_data/FS171819.csv", header=TRUE)
```

## Display pesticide temporal dynamics (Fig. S1)

``` r
PW171819_ <- PW171819
FW171819_ <- FW171819
PS171819_ <- PS171819
FS171819_ <- FS171819
PW171819_$Conc[which(PW171819_$Conc==0)] <- 0.1 #substitute zero by the lower limit finite value
FW171819_$Conc[which(FW171819_$Conc==0)] <- 0.001 #substitute zero by the lower limit finite value
PS171819_$Conc[which(PS171819_$Conc==0)] <- 0.1 #substitute zero by the lower limit finite value
FS171819_$Conc[which(FS171819_$Conc==0)] <- 0.1 #substitute zero by the lower limit finite value

#pdf("./figs/figS1.pdf", width=12, height=7)
par(mfcol=c(2, 2), las=1, family="sans")
plot(log10(Conc) ~ Week2, subset(FW171819_, Tank=="Fipro1"), pch=16, type="o", ylim=c(-3.2, 1.5), 
     xaxt="n", yaxt="n", xlab="", ylab="Water conc. of the insecticide (ug/L)")
lines(log10(Conc) ~ Week2, subset(FW171819_, Tank=="Fipro2"), pch=16, type="o")
lines(log10(Conc) ~ Week2, subset(FW171819_, Tank=="Joint1"), pch=1, type="o", lty=2)
lines(log10(Conc) ~ Week2, subset(FW171819_, Tank=="Joint2"), pch=1, type="o", lty=2)
axis(1, at=seq(0, 70, by=5), 
     labels=c(seq(0, 20, by=5), seq(0, 20, by=5), seq(0, 20, by=5)))
axis(2, lwd=0.5, at=log10(c(0.001, 0.01, 0.1, 1, 10)), lwd.ticks=0.5, 
     labels=c(c(expression(10^{-3}),expression(10^{-2}),expression(10^{-1}),expression(1),expression(10^{1}))))
axis(2, lwd=0.5, lwd.ticks=0.5, at=log10(sapply(-4:2, function(x)(2:9)*10^x)), labels=FALSE, tck=-0.02)
mtext(text=c("Week(2017)", "Week(2018)", "Week(2019)"), 
      side=1, at=c(10, 35, 60), line=2.2)
legend("topright", bty="n", 
       legend=c("Insecticide alone", "Insecticide+Herbicide"), 
       pch=c(16, 1), lty=c(1, 2))
legend("topleft", xpd=NA, legend="a", x.intersp=-1, y.intersp=-3, bty="n", text.font=2, cex=1.2)

plot(log10(Conc) ~ Week2, subset(PW171819_, Tank=="Pent1"), pch=16, type="o", ylim=c(-1.2, 4), 
     xaxt="n", yaxt="n", xlab="", ylab="Water conc. of the herbicide (ug/L)")
lines(log10(Conc) ~ Week2, subset(PW171819_, Tank=="Pent2"), pch=16, type="o")
lines(log10(Conc) ~ Week2, subset(PW171819_, Tank=="Joint1"), pch=1, type="o", lty=2)
lines(log10(Conc) ~ Week2, subset(PW171819_, Tank=="Joint2"), pch=1, type="o", lty=2)
axis(1, at=seq(0, 70, by=5), 
     labels=c(seq(0, 20, by=5), seq(0, 20, by=5), seq(0, 20, by=5)))
axis(2, lwd=0.5, at=log10(c(0.1, 1, 10, 100, 1000, 10000)), lwd.ticks=0.5, 
     labels=c(c(expression(10^{-1}),expression(1),expression(10^{1}),expression(10^{2}),expression(10^{3}),expression(10^{4}))))
axis(2, lwd=0.5, lwd.ticks=0.5, at=log10(sapply(-2:4, function(x)(2:9)*10^x)), labels=FALSE, tck=-0.02)
mtext(text=c("Week(2017)", "Week(2018)", "Week(2019)"), 
      side=1, at=c(10, 35, 60), line=2.2)
legend("topright", bty="n", 
       legend=c("Herbicide alone", "Insecticide+Herbicide"), 
       pch=c(16, 1), lty=c(1, 2))
legend("topleft", xpd=NA, legend="b", x.intersp=-1, y.intersp=-3, bty="n", text.font=2, cex=1.2)

plot(log10(Conc) ~ Week2, subset(FS171819_, Tank=="Fipro1"), pch=16, type="o", ylim=c(-1.2, 3.5), 
     xaxt="n", yaxt="n", xlab="", ylab="Sediment conc. of the insecticide (ug/L)")
lines(log10(Conc) ~ Week2, subset(FS171819_, Tank=="Fipro2"), pch=16, type="o")
lines(log10(Conc) ~ Week2, subset(FS171819_, Tank=="Joint1"), pch=1, type="o", lty=2)
lines(log10(Conc) ~ Week2, subset(FS171819_, Tank=="Joint2"), pch=1, type="o", lty=2)
axis(1, at=seq(0, 70, by=5), 
     labels=c(seq(0, 20, by=5), seq(0, 20, by=5), seq(0, 20, by=5)))
axis(2, lwd=0.5, at=log10(c(0.1, 1, 10, 100, 1000)), lwd.ticks=0.5, 
     labels=c(c(expression(10^{-1}),expression(1),expression(10^{1}),expression(10^{2}),expression(10^{3}))))
axis(2, lwd=0.5, lwd.ticks=0.5, at=log10(sapply(-2:4, function(x)(2:9)*10^x)), labels=FALSE, tck=-0.02)
mtext(text=c("Week(2017)", "Week(2018)", "Week(2019)"), 
      side=1, at=c(10, 35, 60), line=2.2)
legend("topright", bty="n", 
       legend=c("Insecticide alone", "Insecticide+Herbicide"), 
       pch=c(16, 1), lty=c(1, 2))
legend("topleft", xpd=NA, legend="c", x.intersp=-1, y.intersp=-3, bty="n", text.font=2, cex=1.2)

plot(log10(Conc) ~ Week2, subset(PS171819_, Tank=="Pent1"), pch=16, type="o", ylim=c(0, 3.5), 
     xaxt="n", yaxt="n", xlab="", ylab="Sediment conc. of the herbicide (ug/L)")
lines(log10(Conc) ~ Week2, subset(PS171819_, Tank=="Pent2"), pch=16, type="o")
lines(log10(Conc) ~ Week2, subset(PS171819_, Tank=="Joint1"), pch=1, type="o", lty=2)
lines(log10(Conc) ~ Week2, subset(PS171819_, Tank=="Joint2"), pch=1, type="o", lty=2)
axis(1, at=seq(0, 70, by=5), 
     labels=c(seq(0, 20, by=5), seq(0, 20, by=5), seq(0, 20, by=5)))
axis(2, lwd=0.5, at=log10(c(1, 10, 100, 1000)), lwd.ticks=0.5, 
     labels=c(c(expression(1),expression(10^{1}),expression(10^{2}),expression(10^{3}))))
axis(2, lwd=0.5, lwd.ticks=0.5, at=log10(sapply(-1:4, function(x)(2:9)*10^x)), labels=FALSE, tck=-0.02)
mtext(text=c("Week(2017)", "Week(2018)", "Week(2019)"), 
      side=1, at=c(10, 35, 60), line=2.2)
legend("topright", bty="n", 
       legend=c("Herbicide alone", "Insecticide+Herbicide"), 
       pch=c(16, 1), lty=c(1, 2))
legend("topleft", xpd=NA, legend="d", x.intersp=-1, y.intersp=-3, bty="n", text.font=2, cex=1.2)
```

![](7_Pesticide_temporal_dynamics_files/figure-markdown_github/fig-1.png)

``` r
dev.off()
```

    ## null device 
    ##           1