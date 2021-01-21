###########################################################################
#### Robust compositional analysis of wake-time movement behavior data ####
###########################################################################

library(compositions)
library(ggtern)
library(robCompositions)


# loading data
WMBdata = read.csv("WMBdata.csv", sep=";")
head(WMBdata)
n = nrow(WMBdata)

# WMB composition
WMB = WMBdata[, 1:4] # Choose the appropriate columns
head(WMB)
sum(WMB==0) # Are there any 0? If so, they need to be imputed.
D = ncol(WMB)
cn = colnames(WMB)

zBMI = WMBdata[, 5] # Choose the appropriate column

# underweight (z-score under -2) and normal (z-score between -2 and 1)
# overweight (z-score between 1 and 2) and obese (z-score over 2)
Group_zBMI = ifelse(zBMI>=1, "Overweight/obese", "Underweight/normal")
Group_zBMI = factor(Group_zBMI, levels = c("Underweight/normal", "Overweight/obese"))
summary(Group_zBMI)
indicesNo=which(Group_zBMI=="Underweight/normal")
indicesOv=which(Group_zBMI=="Overweight/obese") 

Age = WMBdata[, 6] # Choose the appropriate column
logAge = log(Age)





###################################################################################################################
##### zBMI~WMB

##############################
### Robust compositional mean

set.seed(1) # for reproducipility

# All
comeanr = mean(acomp(WMB), robust="mcd")
comeanr = 100*as.vector(comeanr)
names(comeanr) = cn
round(comeanr, 2)

# Underweight/normal
comeanrNo = mean(acomp(WMB[indicesNo,]) ,robust="mcd")
comeanrNo = 100*as.vector(comeanrNo)
names(comeanrNo) = cn
round(comeanrNo, 2)

# Overweight/obese
comeanrOv = mean(acomp(WMB[indicesOv,]), robust="mcd")
comeanrOv = 100*as.vector(comeanrOv)
names(comeanrOv) = cn
round(comeanrOv, 2)



#########################################################
### Robust compositional geometric mean barplot by group

set.seed(1)

WMBcr = as.data.frame(acomp(WMB)-acomp(comeanr))
cenr = as.numeric(mean(acomp(WMBcr), robust="mcd"))
cenrNo = as.numeric(mean(acomp(WMBcr[indicesNo, ]), robust="mcd"))
cenrOv = as.numeric(mean(acomp(WMBcr[indicesOv, ]), robust="mcd"))

cenrNovscen = round((cenrNo/cenr-1)*100, 2)
cenrOvvscen = round((cenrOv/cenr-1)*100, 2)

# limits for y-axis
ma = round(max(abs(min(cbind(cenrNovscen, cenrOvvscen))), max(cbind(cenrNovscen, cenrOvvscen))))
mi = -ma

get_legend <- function(myggplot){
  tmp = ggplot_gtable(ggplot_build(myggplot))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)
}

tNo = data.frame(wmb = cn, va = cenrNovscen)
bNo = ggplot(tNo, aes(x = wmb, y = va, fill = wmb, color = wmb)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label = va), vjust = 1.5, color = "black", position = position_dodge(0.9), size = 4.5) +
  scale_x_discrete(limits = cn) + 
  scale_fill_manual(breaks = cn, values = c("firebrick2", "steelblue3", "goldenrod1", "palegreen3")) +   
  scale_color_manual(breaks = cn, values = c("red4", "blue4", "darkorange", "green4")) +  
  ggtitle(levels(Group_zBMI)[1]) +
  ylab("[%]") +
  geom_hline(yintercept = 0, color = "black") +
  scale_y_continuous(limits = c(mi, ma)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey'),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 15, face = "bold"), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 15),
        axis.ticks.x = element_blank(), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.title = element_blank(), 
        legend.position = "top", 
        legend.text = element_text(size = 15),
        legend.spacing.x = unit(0.25,"cm"))
leg = get_legend(bNo)
bNo = bNo + theme(legend.position = "none")

tOv = data.frame(wmb = cn, va = cenrOvvscen)
bOv = ggplot(tOv, aes(x = wmb, y = va, fill = wmb, color = wmb)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label = va), vjust = 1.5, color = "black", position = position_dodge(0.9), size = 4.5) +
  scale_x_discrete(limits = cn) + 
  scale_fill_manual(breaks = cn, values = c("firebrick2", "steelblue3", "goldenrod1", "palegreen3")) +  
  scale_color_manual(breaks = cn, values = c("red4", "blue4", "darkorange", "green4")) +  
  ggtitle(levels(Group_zBMI)[2]) +
  geom_hline(yintercept = 0, color = "black") +
  ylim(mi, ma) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey'),
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.position = "none")

pdf("Barplots.pdf", width = 8)
grid.arrange(bNo, bOv, leg, ncol = 2, nrow=2, layout_matrix = rbind(c(1, 2), c(3, 3)),
             widths = c(2.345, 2), heights = c(2.5, 0.2))
dev.off()



#########################
### Robust MM regression

# with pivot coordinates as explanatory variables
DR = data.frame(zBMI, matrix(rnorm(D*n), n, D))
colnames(DR)[2:(D+1)] = paste("z1_", cn, sep = "")
REG = lmrob(zBMI~., data=DR)
sumREG = summary(REG)

for(i in 1:D) {
  datareg = data.frame(zBMI, pivotCoord(WMB, pivotvar = i))
  reg = lmrob(zBMI~., data=datareg)
  sumreg = summary(reg)
  if (i==1) {
    sumREG$coefficients[1:2,] = sumreg$coefficients[1:2,]
    sumREG$residuals = sumreg$residuals
    sumREG$scale = sumreg$scale
    sumREG$r.squared = sumreg$r.squared
    sumREG$adj.r.squared = sumreg$adj.r.squared
    sumREG$converged = sumreg$converged
    sumREG$iter = sumreg$iter
    sumREG$rweights = sumreg$rweights
    sumREG$fitted.values = sumreg$fitted.values
    sumREG$init.S = sumreg$init.S
    sumREG$init = sumreg$init
    sumREG$rank = sumreg$rank
    sumREG$cov = sumreg$cov
    sumREG$df.residual = sumreg$df.residual
    sumREG$weights = sumreg$weights
    sumREG$control = sumreg$control
    sumREG$call = sumreg$call
    sumREG$control = sumreg$control
  } else {
    sumREG$coefficients[(1+i),] = sumreg$coefficients[2,]
  }
}
sumREG


# with "ordinal" pivot coordinates as explanatory variables
datareg = cbind(zBMI, pivotCoord(WMB[,4:1]))
reg = lmrob(zBMI~., data = datareg)
sumreg = summary(reg)
sumreg

datareg2 = cbind(zBMI, pivotCoord(WMB))
reg2 = lmrob(zBMI~., data = datareg2)
sumreg2 = summary(reg2)
sumreg2


# with "pivoting" balances as explanatory variables
BAL = list()
CODES = list()
VV = list()

codes = matrix(0, D-1, D)
for(i in 1:(D-1)){
  codes[i,] = c(rep(0,i-1),1,rep(-1,D-i))
}   
V = t(gsi.buildilrBase(t(codes)))
cnB = character()
for (i in 1:(D-1)) {
  num = paste(cn[codes[i,]==1], collapse='.')
  denum = paste(cn[codes[i,]==-1], collapse='.')
  cnB[i] = paste(num, denum, sep='_')
}
bal = log(as.matrix(WMB))%*%t(V)
colnames(bal) = cnB
BAL[[1]] = bal
CODES[[1]] = codes
VV[[1]] = V

codes = matrix(c(1,1,-1,-1,
                 1,-1,0,0,
                 0,0,1,-1), D-1, D, byrow=T)
V = t(gsi.buildilrBase(t(codes)))
cnB = character()
for (i in 1:(D-1)) {
  num = paste(cn[codes[i,]==1], collapse='.')
  denum = paste(cn[codes[i,]==-1], collapse='.')
  cnB[i] = paste(num, denum, sep='_')
}
bal = log(as.matrix(WMB))%*%t(V)
colnames(bal) = cnB
BAL[[2]] = bal
CODES[[2]] = codes
VV[[2]] = V

codes = matrix(c(1,1,1,-1,
                 1,-1,-1,0,
                 0,1,-1,0), D-1, D, byrow=T)
V = t(gsi.buildilrBase(t(codes)))
cnB = character()
for (i in 1:(D-1)) {
  num = paste(cn[codes[i,]==1], collapse='.')
  denum = paste(cn[codes[i,]==-1], collapse='.')
  cnB[i] = paste(num, denum, sep='_')
}
bal = log(as.matrix(WMB))%*%t(V)
colnames(bal) = cnB
BAL[[3]] = bal
CODES[[3]] = codes
VV[[3]] = V

lB = length(BAL)
nB = character()
for(i in 1:lB) {
  nB = c(nB, colnames(BAL[[i]])[1])
}
nB

DR = data.frame(zBMI, matrix(rnorm(lB*n), n, lB))
colnames(DR)[2:(lB+1)] = nB
REG = lmrob(zBMI~., data=DR)
sumREG = summary(REG)

for(i in 1:lB) {
  datareg = data.frame(zBMI, BAL[[i]])
  reg = lmrob(zBMI~., data=datareg)
  sumreg = summary(reg)
  if (i==1) {
    sumREG$coefficients[1:2,] = sumreg$coefficients[1:2,]
    sumREG$residuals = sumreg$residuals
    sumREG$scale = sumreg$scale
    sumREG$r.squared = sumreg$r.squared
    sumREG$adj.r.squared = sumreg$adj.r.squared
    sumREG$converged = sumreg$converged
    sumREG$iter = sumreg$iter
    sumREG$rweights = sumreg$rweights
    sumREG$fitted.values = sumreg$fitted.values
    sumREG$init.S = sumreg$init.S
    sumREG$init = sumreg$init
    sumREG$rank = sumreg$rank
    sumREG$cov = sumreg$cov
    sumREG$df.residual = sumreg$df.residual
    sumREG$weights = sumreg$weights
    sumREG$control = sumreg$control
    sumREG$call = sumreg$call
    sumREG$control = sumreg$control
  } else {
    sumREG$coefficients[(1+i),] = sumreg$coefficients[2,]
  }
}
sumREG





###################################################################################################################
##### WMB~Age

#####################
### Ternary diagrams

tab = as.data.frame(cbind(Age, WMB))
head(tab)

t1 = ggtern(data = tab, aes(SB, LPA, MPA)) + 
  geom_point(aes(fill = Age), shape = 21, size = 2) +
  scale_fill_gradient(name = "Age", low = 'black', high = 'gold') +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "grey80"),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.title = element_text(face="bold"),
        tern.axis.title = element_text(size = 15, face = "bold"),
        tern.axis.text = element_text(size = 15, face = "bold")) + 
  theme_legend_position(x = "topright") +
  theme_showarrows()
leg = get_legend(t1)
t1 = t1 + theme(legend.position = "none")
t2 = ggtern(data = tab, aes(SB, LPA, VPA)) + 
  geom_point(aes(fill = Age), shape = 21, size = 2) +
  scale_fill_gradient(name = "Age", low = 'black', high = 'gold') +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "grey80"),
        legend.position = "none",
        axis.title = element_text(face="bold"),
        tern.axis.title = element_text(size = 15, face = "bold"),
        tern.axis.text = element_text(size = 15, face = "bold")) + 
  theme_showarrows()
t3 = ggtern(data = tab, aes(SB, MPA, VPA)) + 
  geom_point(aes(fill = Age), shape = 21, size = 2) +
  scale_fill_gradient(name = "Age", low = 'black', high = 'gold') +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "grey80"),
        legend.position = "none",
        axis.title = element_text(face="bold"),
        tern.axis.title = element_text(size = 15, face = "bold"),
        tern.axis.text = element_text(size = 15, face = "bold")) + 
  theme_showarrows()
t4 = ggtern(data = tab, aes(LPA, MPA, VPA)) + 
  geom_point(aes(fill = Age), shape = 21, size = 2) +
  scale_fill_gradient(name = "Age", low = 'black', high = 'gold') +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "grey80"),
        legend.position = "none",
        axis.title = element_text(face="bold"),
        tern.axis.title = element_text(size = 15, face = "bold"),
        tern.axis.text = element_text(size = 15, face = "bold")) + 
  theme_showarrows()

pdf("Ternary.pdf", height = 10, width = 12)
grid.arrange(t1, t2, t3, t4, leg, ncol = 3, nrow=2, layout_matrix = cbind(c(1, 3), c(2, 4), c(5, NA)),
             widths = c(2.5, 2.5, 0.3), heights = c(2.5, 2.5))
dev.off()


# robust centering
tabR = as.data.frame(cbind(Age, WMBcr))
head(tabR)

t1 = ggtern(data = tabR, aes(SB, LPA, MPA)) + 
  geom_point(aes(fill = Age), shape = 21, size = 2) +
  scale_fill_gradient(name = "Age", low = 'black', high = 'gold') +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        panel.background = element_rect(fill = "white"),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.title = element_text(face="bold"),
        tern.axis.title = element_text(size = 15, face = "bold"),
        tern.axis.text = element_text(size = 15, face = "bold")) + 
  theme_legend_position(x = "topright") +
  theme_showarrows() + 
  theme_nolabels() + 
  theme_nogrid() + 
  theme_noticks()
leg = get_legend(t1)
t1 = t1 + theme(legend.position = "none")
t2 = ggtern(data = tabR, aes(SB, LPA, VPA)) + 
  geom_point(aes(fill = Age), shape = 21, size = 2) +
  scale_fill_gradient(name = "Age", low = 'black', high = 'gold') +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        axis.title = element_text(face="bold"),
        tern.axis.title = element_text(size = 15, face = "bold"),
        tern.axis.text = element_text(size = 15, face = "bold")) + 
  theme_showarrows() + 
  theme_nolabels() + 
  theme_nogrid() + 
  theme_noticks()
t3 = ggtern(data = tabR, aes(SB, MPA, VPA)) + 
  geom_point(aes(fill = Age), shape = 21, size = 2) +
  scale_fill_gradient(name = "Age", low = 'black', high = 'gold') +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        axis.title = element_text(face="bold"),
        tern.axis.title = element_text(size = 15, face = "bold"),
        tern.axis.text = element_text(size = 15, face = "bold")) + 
  theme_showarrows() + 
  theme_nolabels() + 
  theme_nogrid() + 
  theme_noticks()
t4 = ggtern(data = tabR, aes(LPA, MPA, VPA)) + 
  geom_point(aes(fill = Age), shape = 21, size = 2) +
  scale_fill_gradient(name = "Age", low = 'black', high = 'gold') +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        axis.title = element_text(face="bold"),
        tern.axis.title = element_text(size = 15, face = "bold"),
        tern.axis.text = element_text(size = 15, face = "bold")) + 
  theme_showarrows() + 
  theme_nolabels() + 
  theme_nogrid() + 
  theme_noticks()

pdf("TernaryCentered.pdf", height = 10, width = 12)
grid.arrange(t1, t2, t3, t4, leg, ncol = 3, nrow=2, layout_matrix = cbind(c(1, 3), c(2, 4), c(5, NA)),
             widths = c(2.5, 2.5, 0.3), heights = c(2.5, 2.5))
dev.off()



#########################
### Robust MM regression 

# response - selected balance
codes = matrix(c(1,-1,-1,1,
                 1,0,0,-1,
                 0,1,-1,0), D-1, D, byrow=T)
V = t(gsi.buildilrBase(t(codes)))
cnB = character()
for (i in 1:(D-1)) {
  num = paste(cn[codes[i,]==1], collapse='.')
  denum = paste(cn[codes[i,]==-1], collapse='.')
  cnB[i] = paste(num, denum, sep='_')
}
bal = log(as.matrix(WMB))%*%t(V)
colnames(bal) = cnB
head(bal)

reg = lmrob(bal[,1]~logAge)
summary(reg)


# response - each time different first pivot coordinate
sumt = NULL
for (i in 1:D){
  z1 = pivotCoord(WMB, pivotvar = i)[,1]
  reg = lmrob(z1~logAge)
  sumreg = summary(reg)
  sumr = c(sumreg$coefficients[2,], sumreg$r.squared, sumreg$adj.r.squared)
  sumt = rbind (sumt,sumr)
}
rownames(sumt) = paste("z1_", cn, sep = "")
colnames(sumt) = c("Coef. Est.", "Std. Error", "t value", "Pr(>|t|)", "R-sq.", "Adj. R-sq.")
sumt 
