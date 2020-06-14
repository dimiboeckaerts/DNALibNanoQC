########################
# 1. Environment setup #
########################

library("rstudioapi")

setwd(dirname(getActiveDocumentContext()$path))

library("ggplot2")
library("plotly")
library("fitdistrplus")

# Multiple plot function from R cookbook 
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

################################
# 2. Analysis of single blocks #
################################

single.blocks <- read.csv("count_table_single_blocks.csv", header=FALSE, sep="\t", col.names =c("Block", "Position", "Count"))

# subsetting dataset per position in the assembly product
omp.blocks <- single.blocks[single.blocks$Position=="1",]
linker.blocks <- single.blocks[single.blocks$Position=="2",]
cbd.blocks <- single.blocks[single.blocks$Position=="3",]
ead.blocks <- single.blocks[single.blocks$Position=="4",]

# expected average for each position
expected.omp <- sum(omp.blocks$Count)/nrow(omp.blocks)
expected.linker <- sum(linker.blocks$Count)/nrow(linker.blocks)
expected.cbd <- sum(cbd.blocks$Count)/nrow(cbd.blocks)
expected.ead <- sum(ead.blocks$Count)/nrow(ead.blocks)

p1 <- ggplot(data = omp.blocks, aes(x=Position, y=Count, label=Block)) + 
  geom_violin(col="grey", alpha=0.1) +
  geom_jitter(width=0.03, col="darkblue", alpha=0.7) +
  labs(x="OMP blocks", y="Count") +
  geom_boxplot(width=0.2, alpha=0.1) +
  geom_hline(yintercept=expected.omp, linetype="dashed", color="red", alpha=0.9) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")

p2 <- ggplot(data = linker.blocks, aes(x=Position, y=Count, label=Block)) +
 geom_point(width=0.03, col="darkblue", alpha=0.7) +
 labs(x="Linker blocks", y="Count") +
 theme(axis.text.x=element_blank(),
       axis.ticks.x=element_blank(),
       legend.position="none")

p3 <- ggplot(data = cbd.blocks, aes(x=Position, y=Count, label=Block)) +
  geom_violin(col="grey", alpha=0.1) +
  geom_jitter(width=0.03, col="darkblue", alpha=0.7) +
  geom_hline(yintercept=expected.cbd, linetype="dashed", color="red", alpha=0.9) +
  labs(x="CBD blocks", y="Count") +
  geom_boxplot(width=0.2, alpha=0.1) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")

p4 <- ggplot(data = ead.blocks, aes(x=Position, y=Count, label=Block)) +
  geom_violin(col="grey", alpha=0.1) +
  geom_jitter(width=0.03, col="darkblue", alpha=0.7) +
  geom_hline(yintercept=expected.ead, linetype="dashed", color="red", alpha=0.9) +
  labs(x="EAD blocks", y="Count") +
  geom_boxplot(width=0.2, alpha=0.1) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")

# plots side-by-side
multiplot(p1, p2, p3, p4, cols=2) 

# inspecting the graphs interactively with plotly
ggplotly(p1)
ggplotly(p2)
ggplotly(p3)
ggplotly(p4)


########################################################################
# Zooming in on the linker position, which shows the largest imbalance #
########################################################################

# subsetting dataset to linkers only
blocks <- read.csv("reads_pass.csv", header=TRUE, sep=",")
linkers <- subset(blocks, TileName == "Flexiblemedian" | TileName == "Flexibleshort")

ggplot(data=linkers, aes(AlignmentLength)) +
  geom_histogram(aes(y =..density..), fill="blue", col="blue", alpha=.4, binwidth = 2) +
  stat_density(col="darkgreen", geom="line") +
  geom_vline(xintercept = 30, linetype="dashed", color="red")+
  annotate("text", label= "Short linker (30 nt)", x=25, y=0.047) +
  annotate("text", label= "Long linker (54 nt)", x=49, y=0.047) +
  geom_vline(xintercept = 54, linetype="dashed", color= "red")+
  theme_update(plot.title = element_text(hjust = 0.5)) +
  labs(title="Linker size distribution", x="Linker length (nt)", y="Density (count)") 

############################
# 3. Analysis of all tiles #
############################

combinations.blocks <- read.csv("count_table_combinations_blocks.csv", header=FALSE, sep="\t", col.names =c("BlockCombinations", "Count"))

# ploting empirical distribution of the data and exploring skewness-kurtosis plot (Cullen and Frey)
par(mfrow=c(1,1))
plotdist(combinations.blocks$Count, histo=TRUE, discrete=TRUE)
descdist(combinations.blocks$Count, discrete=FALSE, boot=500)

# fitting different distributions with the dataset
fit_p <- fitdist(combinations.blocks$Count, "pois")
fit_u <- fitdist(combinations.blocks$Count, "unif")
fit_e <- fitdist(combinations.blocks$Count, "exp")
fit_ln <- fitdist(combinations.blocks$Count, "lnorm")
fit_w <- fitdist(combinations.blocks$Count, "weibull")
fit_g <- fitdist(combinations.blocks$Count, "gamma")

# plotting side by side all the models
par(mfrow=c(2,2))
plot.legend <- c("poisson", "uniform", "exp", "lognormal", "weibull", "gamma")
denscomp(list(fit_p, fit_u, fit_e, fit_ln, fit_w, fit_g), legendtext = plot.legend)
cdfcomp (list(fit_p, fit_u, fit_e, fit_ln, fit_w, fit_g), legendtext = plot.legend)
qqcomp  (list(fit_p, fit_u, fit_e, fit_ln, fit_w, fit_g), legendtext = plot.legend)
ppcomp  (list(fit_p, fit_u, fit_e, fit_ln, fit_w, fit_g), legendtext = plot.legend)

# parameters of each model
summary(fit_p)
summary(fit_u)
summary(fit_e)
summary(fit_ln)
summary(fit_w)
summary(fit_g)

# goodness of fit
gofstat(list(fit_p, fit_u, fit_e, fit_ln, fit_w, fit_g), fitnames = c("pois", "unif", "exp", "logn", "weibull", "gamma"))

# transforming dataset to have each column being the times combinations occur, ie:
# no combination is missing
# unique combinations appear 766 times (out of 8267) 
# double combinations appear 903 times 
# triple combinations appear 790 times 
# etc.
# This can be obtained with: table(combinations.blocks$Count)
max.combinations <- max(combinations.blocks$Count)

counts.combinations <- data.frame(1:max.combinations)
colnames(counts.combinations)[1] <- "Counts"
tmp.counts <- data.frame(table(combinations.blocks$Count))
colnames(tmp.counts) <- c("Counts", "Freq")
final.counts <- merge(counts.combinations, tmp.counts, by="Counts", all.x=TRUE) # Left outer join
final.counts$Freq[is.na(final.counts$Freq)] <- 0

# generation of datasets using multiple distributions with fitted parameters
simulated.poisson <- dpois(1:max.combinations, lambda=fit_p$estimate["lambda"], log=FALSE)
simulated.uniform <- dunif(1:max.combinations, min=fit_u$estimate["min"], max=fit_u$estimate["max"])
simulated.exp <- dexp(1:max.combinations, rate=fit_e$estimate["rate"])
simulated.ln <- dlnorm(1:max.combinations, meanlog=fit_ln$estimate["meanlog"], sdlog=fit_ln$estimate["sdlog"])
simulated.weibull <- dweibull(1:max.combinations, shape=fit_w$estimate["shape"], scale=fit_w$estimate["scale"])
simulated.gamma <- dgamma(1:max.combinations, shape=fit_g$estimate["shape"], rate=fit_g$estimate["rate"])

# testing empirical distribution vs fitted model using chisquare
chisq.test(x=final.counts$Freq, p=simulated.poisson, simulate.p.value=TRUE, rescale.p = TRUE)
chisq.test(x=final.counts$Freq, p=simulated.uniform, simulate.p.value=TRUE, rescale.p = TRUE)
chisq.test(x=final.counts$Freq, p=simulated.exp, simulate.p.value=TRUE, rescale.p = TRUE)
chisq.test(x=final.counts$Freq, p=simulated.ln, simulate.p.value=TRUE,  rescale.p = TRUE)
chisq.test(x=final.counts$Freq, p=simulated.weibull, simulate.p.value=TRUE, rescale.p = TRUE)
chisq.test(x=final.counts$Freq, p=simulated.gamma, simulate.p.value=TRUE, rescale.p = TRUE)

##################################################
# checking to see if the linker has an influence #
# exploration of linker size distribution        #
##################################################



# considering linkers as a single condition
combinations.blocks.nolinker <- read.csv("count_table_combinations_tiles_no_linkers.csv", header=FALSE, sep="\t", col.names =c("TileCombination", "Count"))

par(mfrow=c(1,1))
plotdist(combinations.blocks.nolinker$Count, histo=TRUE, discrete=TRUE)
descdist(combinations.blocks.nolinker$Count, boot=500, method="sample")

# fitting different distributions with the dataset
fit_ln <- fitdist(combinations.blocks.nolinker$Count, "lnorm")
fit_w <- fitdist(combinations.blocks.nolinker$Count, "weibull")
fit_g <- fitdist(combinations.blocks.nolinker$Count, "gamma")

# parameters of each model
summary(fit_ln)
summary(fit_w)
summary(fit_g)

# plotting side by side all the models
par(mfrow=c(2,2))
plot.legend <- c("lognormal", "weibull", "gamma")
denscomp(list(fit_ln, fit_w, fit_g), legendtext = plot.legend)
cdfcomp (list(fit_ln, fit_w, fit_g), legendtext = plot.legend)
qqcomp  (list(fit_ln, fit_w, fit_g), legendtext = plot.legend)
ppcomp  (list(fit_ln, fit_w, fit_g), legendtext = plot.legend)

# goodness of fit
gofstat(list(fit_ln, fit_w, fit_g), fitnames = c("logn", "weibull", "gamma"))

# transforming dataset to have each column being the times combinations occur, ie:
# no combination is missing
# unique combinations appear 766 times (out of 8267) 
# double combinations appear 903 times 
# triple combinations appear 790 times 
# etc.
# This can be obtained with: table(combinations.blocks$Count)
max.combinations <- max(combinations.blocks.nolinker$Count)

counts.combinations <- data.frame(1:max.combinations)
colnames(counts.combinations)[1] <- "Counts"
tmp.counts <- data.frame(table(combinations.blocks.nolinker$Count))
colnames(tmp.counts) <- c("Counts", "Freq")
final.counts <- merge(counts.combinations, tmp.counts, by="Counts", all.x=TRUE) # Left outer join
final.counts$Freq[is.na(final.counts$Freq)] <- 0

# generation of datasets using lnorm, weibull, and gamma with fitted parameters
simulated.ln <- dlnorm(1:max.combinations, meanlog=fit_ln$estimate["meanlog"], sdlog=fit_ln$estimate["sdlog"])
simulated.weibull <- dweibull(1:max.combinations, shape=fit_w$estimate["shape"], scale=fit_w$estimate["scale"])
simulated.gamma <- dgamma(1:max.combinations, shape=fit_g$estimate["shape"], rate=fit_g$estimate["rate"])

# testing empirical distribution vs fitted model using chisquare
chisq.test(x=final.counts$Freq, p=simulated.ln, simulate.p.value=TRUE,  rescale.p = TRUE)
chisq.test(x=final.counts$Freq, p=simulated.weibull, simulate.p.value=TRUE, rescale.p = TRUE)
chisq.test(x=final.counts$Freq, p=simulated.gamma, simulate.p.value=TRUE, rescale.p = TRUE)

