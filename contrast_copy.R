



 

load('../../dds[[sp]][[scheme]][[formula]].RData')



ddsg <- list()
resg <- list()
generationnames <- c('gen1', 'gen2', 'gen3')

for (gi in 1:3) {
  g <- generationnames[gi]
  cat(g) ; cat('\n')
  ddsg[[gi]] <- original[ , colData(original)$Generation == g]
  design(ddsg[[gi]]) <-  ~ Line * Food
  ddsg[[gi]] <- DESeq(ddsg[[gi]], betaPrior = FALSE)
  ## resg[[gi]] <- list()
  ## resultnames <- resultsNames(ddsg[[gi]])
  ## for (ri in 1:length(resultnames)) {
  ##   resultname <- resultnames[ri]
  ##   resg[[gi]][[ri]] <- results(ddsg[[gi]], resultname)
  ## }
}
names(ddsg) <- generationnames
## for (gi in 1:3) {
##   names(resg[[gi]]) <- resultnames
## }

if (!file.exists('D. magna')) { dir.create('D. magna') }


## UP and DOWN
DEgenenames <- list()
for (di in 1:length(ddsg)) { # 1:3, 3 generations
  dirname <- sprintf('D. magna/%s/', names(ddsg)[di])
  cat(dirname) ; cat('\n')
  if (!file.exists(dirname)) { dir.create(dirname) }
  writeDifferntiallyExpressedGenes(ddsg[[di]], dirname)
  DEgenenames[[di]] <-
    read.table(sprintf('%s/%s',  
                       dirname,
                       '04 LineP_line.FoodP gene names all only diff exp by log2FC.txt'),  # I do not understand this step ... I thought we write the results to the directory and not read ...
               as.is = TRUE)[[1]]
}
intersectgenes <-
  intersect(DEgenenames[[1]],
            intersect(DEgenenames[[2]], DEgenenames[[3]])
            )
## uniongenes <-
##   union(DEgenenames[[1]], union(DEgenenames[[2]], DEgenenames[[3]]))
names(DEgenenames) <- names(ddsg)
venn.diagram(DEgenenames, filename = 'D. magna/venndiagram.tiff')

## UP
DEgenenamesU <- list()
for (di in 1:length(ddsg)) { # 1:3, 3 generations
  dirname <- sprintf('D. magna/%s/', names(ddsg)[di])
  DEgenenamesU[[di]] <-
    read.table(sprintf('%s/%s',
                       dirname,
                       '04 LineP_line.FoodP gene names upreg by log2FC.txt'),
               as.is = TRUE)[[1]]
}
intersectgenesU <-
  intersect(DEgenenamesU[[1]],
            intersect(DEgenenamesU[[2]], DEgenenamesU[[3]])
            )
names(DEgenenamesU) <- names(ddsg)
venn.diagram(DEgenenamesU, filename = 'D. magna/venndiagram_up_reg.tiff')

## DOWN
DEgenenamesD <- list()
for (di in 1:length(ddsg)) { # 1:3, 3 generations
  dirname <- sprintf('D. magna/%s/', names(ddsg)[di])
  DEgenenamesD[[di]] <-
    read.table(sprintf('%s/%s',
                       dirname,
                       '04 LineP_line.FoodP gene names downreg by log2FC.txt'),
               as.is = TRUE)[[1]]
}
intersectgenesD <-
  intersect(DEgenenamesD[[1]],
            intersect(DEgenenamesD[[2]], DEgenenamesD[[3]])
            )
names(DEgenenamesD) <- names(ddsg)
venn.diagram(DEgenenamesD, filename = 'D. magna/venndiagram_down_reg.tiff')



# Plots -------------------------------------------------------------------
str(original@assays$data@listData$counts)
counts <- original@assays$data@listData$counts
head(counts)
mu <- original@assays$data@listData$mu # after normalization
head(mu)
round(tail(mu), digits = 7)

png("Contrasts between samples_magna.png", width=480, height=480)
par(mfrow = c(2, 3))
ff <- apply(mu[ , 1:3], 1, mean)
str(ff)
pp <- apply(mu[ , 19:21], 1, mean)
plot(log(pp), log(ff), col = '#11111101', pch = 19, xlim = c(0, 10), ylim = c(0, 10))
abline(b = 1, a = 0)

fff <- apply(mu[ , 4:6], 1, mean)
str(fff)
ppp <- apply(mu[ , 31:33], 1, mean)
plot(log10(ppp), log10(fff), col = '#11111101', pch = 19, xlim = c(0, 10), ylim = c(0, 10))
#model1 <- lm(log10(fff) ~ log10(ppp)-1)
#abline(a = 0, b = coef(model1), col="red") 
abline(b = 1, a = 0)

ffff <- apply(mu[ , 7:9], 1, mean)
pppp <- apply(mu[ , 34:36], 1, mean)
plot(log(pppp), log(ffff), col = '#11111101', pch = 19, xlim = c(0, 10), ylim = c(0, 10))
abline(b = 1, a = 0)


fp <- apply(mu[ , 10:12], 1, mean)
str(fp)
pf <- apply(mu[ , 19:21], 1, mean)
plot(log(pf), log(fp), col = '#11111101', pch = 19, xlim = c(0, 10), ylim = c(0, 10))
abline(b = 1, a = 0)

fpp <- apply(mu[ , 13:15], 1, mean)
str(fp)
pff <- apply(mu[ , 22:24], 1, mean)
plot(log(pff), log(fpp), col = '#11111101', pch = 19, xlim = c(0, 10), ylim = c(0, 10))
abline(b = 1, a = 0)
#log(max(pff))

fppp <- apply(mu[ , 16:18], 1, mean)
str(fp)
pfff <- apply(mu[ , 25:27], 1, mean)
plot(log(pfff), log(fppp), col = '#11111101', pch = 19, xlim = c(0, 10), ylim = c(0, 10))
abline(b = 1, a = 0)
dev.off()



## ## taking contrast [generation2]
## ddsCgen2 <- original[ , colData(original)$Generation == 'gen2']
## design(ddsCgen2) <- ~ Line * Food
## ddsCgen2 <- DESeq(ddsCgen2, betaPrior = FALSE)
## ## > resultsNames(ddsCgen2)
## ## [1] "Intercept"             "Line_P_line_vs_F_line" "Food_P_vs_F"          
## ## [4] "LineP_line.FoodP"

## # resCgen2 <- results(ddsCgen2, contrast = c('Line', 'P_line', 'F_line'))
## resCgen2 <- results(ddsCgen2, 'Line_P_line_vs_F_line')


# Contrast across the generation ------------------------------------------
#17.04.2014_ Maybe Koji could have a look at this part. 
# I think instead of making nested lists I am overwriting one from P_line
# in addion they are not nested 
load('../../dds[[sp]][[scheme]][[formula]].RData')
original <- ddslistlistsbysp[['magna']][['allgen']][[8]]
linenames <- c("P_line", "F_line") # becuase F_line is last F will be produced 
generationnames <- c('gen1', 'gen2', 'gen3')

ddsglist <- list()
ddsllist <- list()
for (li in 1:length(linenames)) {
  l <- linenames[li]
  cat(l); cat('\n')
  
  for (gi in 1:length(generationnames)) {
    g <- generationnames[gi]
    cat(g) ; cat('\n')
    temp <- original[ , colData(original)$Line ==l]
    ddsglist[[gi]] <- temp[ , colData(temp)$Generation == g]  # I needed to split these two lines
    design(ddsglist[[gi]]) <- ~ Food
    ddsglist[[gi]] <- DESeq(ddsglist[[gi]], betaPrior = FALSE)
}
  ddsllist[[li]] <- ddsglist  #  how to nest the generation inside line?
}

### head(ddsglist$gen3@assays$data@listData$mu)
names(ddsllist) <- linenames
names(ddsglist) <- generationnames
str(ddsglist)
str(ddsllist)
str(ddsllist$P_line$gen1)

ddsglist$gen1  # this is only magna  from F line !!! so I am overwriting the P_line somehow
ddsllist$F_line 


# Printing the results to folders and Venn diagrams -------------------------------------

if (!file.exists('D. magna_ffff_fpppp')) { dir.create('D. magna_ffff_fpppp') }

## UP and DOWN
DEgenenames <- list()
for (di in 1:3) { # 1:3, 3 generations #length(ddsglist)
  dirname <- sprintf('D. magna_ffff_fpppp/%s/', names(ddsglist)[di])
  cat(dirname) ; cat('\n')
  if (!file.exists(dirname)) { dir.create(dirname) }
  writeDifferntiallyExpressedGenes(ddsglist[[di]], dirname)
  DEgenenames[[di]] <-
    read.table(sprintf('%s/%s',
                       dirname,
                       '02 Food_P_vs_F gene names all only diff exp by log2FC.txt'),
               as.is = TRUE)[[1]]
}

intersectgenes <-
  intersect(DEgenenames[[1]],
            intersect(DEgenenames[[2]], DEgenenames[[3]])
            )
## uniongenes <-
##   union(DEgenenames[[1]], union(DEgenenames[[2]], DEgenenames[[3]]))
names(DEgenenames) <- names(ddsglist)
venn.diagram(DEgenenames, filename = 'D. magna_ffff_fpppp/venndiagram_ffff_fppp.tiff')


# pppp_pffff for each generation ----------------------------------------------------------
linenames <- c("F_line", "P_line")  # becuase P_line is the end this will be done last
generationnames <- c('gen1', 'gen2', 'gen3')

ddsglist <- list()
ddsllist <- list()
for (li in 1:length(linenames)) {
  l <- linenames[li]
  cat(l); cat('\n')
  
  for (gi in 1:length(generationnames)) {
    g <- generationnames[gi]
    cat(g) ; cat('\n')
    temp <- original[ , colData(original)$Line ==l]
    ddsglist[[gi]] <- temp[ , colData(temp)$Generation == g]  # I needed to split these two lines
    design(ddsglist[[gi]]) <- ~ Food
    ddsglist[[gi]] <- DESeq(ddsglist[[gi]], betaPrior = FALSE)
}
  ddsllist[[li]] <- ddsglist  #  how to nest the generation inside line?
}

### head(ddsglist$gen3@assays$data@listData$mu)
names(ddsllist) <- linenames
names(ddsglist) <- generationnames
str(ddsglist)
str(ddsllist)
str(ddsllist$P_line$gen1)

ddsglist$gen1  # this is only magna  from F line !!! so I am overwriting the P_line somehow
ddsllist$F_line 


if (!file.exists('D. magna_pppp_pffff')) { dir.create('D. magna_pppp_pffff') }

## UP and DOWN
DEgenenames <- list()
for (di in 1:3) { # 1:3, 3 generations #length(ddsglist)
  dirname <- sprintf('D. magna_pppp_pffff/%s/', names(ddsglist)[di])
  cat(dirname) ; cat('\n')
  if (!file.exists(dirname)) { dir.create(dirname) }
  writeDifferntiallyExpressedGenes(ddsglist[[di]], dirname)
  DEgenenames[[di]] <-
    read.table(sprintf('%s/%s',
                       dirname,
                       '02 Food_P_vs_F gene names all only diff exp by log2FC.txt'),
               as.is = TRUE)[[1]]
}

intersectgenes <-
  intersect(DEgenenames[[1]],
            intersect(DEgenenames[[2]], DEgenenames[[3]])
            )
length(intersectgenes)

## uniongenes <-
##   union(DEgenenames[[1]], union(DEgenenames[[2]], DEgenenames[[3]]))
names(DEgenenames) <- names(ddsglist)
venn.diagram(DEgenenames, filename = 'D. magna_pppp_pffff/venndiagram_pppp_pffff.tiff')
head(DEgenenames)


# Now I will look at the intersection accross the generations for F_line and F good -------------
linenames <- c("P_line", "F_line")   # becauase F_line is the end this will be done last 
#generationnames <- c('gen1', 'gen2', 'gen3')
foodnames <- c("P", "F")  # replace later on diet to food!  # F will printed into the list as it is last one foodnames 

ddsllist <- list() # line list
#ddsglist <- list() # generation list
ddsfoodlist <- list() # food list

for (li in 1:length(linenames)) {
  l <- linenames[li]
  cat(l); cat('\n')
  
  for (fi in 1:length(foodnames)) {
    f <- foodnames[fi]
    cat(f) ; cat('\n')
    temp <- original[ , colData(original)$Line ==l]
    ddsfoodlist[[fi]] <- temp[ , colData(temp)$Food == d]  # I needed to split these two lines
    design(ddsfoodlist[[fi]]) <- ~ Generation
    ddsfoodlist[[fi]] <- DESeq(ddsdlist[[fi]], betaPrior = FALSE)
}
  ddsllist[[li]] <- ddsfoodlist  #  how to nest the generation inside line?
}

names(ddsllist) <- linenames
names(ddsglist) <- generationnames
names(ddsfoodlist) <- foodnames

head(ddsdlist$F@assays$data@listData$mu)
head(ddsdlist$P@assays$data@listData$mu)


# Plots for ffvsfffvsffff or pfvspffvspfff-------------------------------------------------
# Plots for ff in Fline-------------------------------------------------------------------

mu <- ddsdlist$F@assays$data@listData$mu
round(tail(mu), digits = 7)

png("Contrasts between ff_fff_ffff_magna.png", width=480, height=480)
par(mfrow = c(1, 3))
title(main="Contrasts between ff_fff_ffff_magna")
ff <- apply(mu[ , 1:3], 1, mean)
fff <- apply(mu[ , 4:6], 1, mean)
plot(log(ff), log(fff), col = '#11111101', pch = 19, xlim = c(0, 10), ylim = c(0, 10))
abline(b = 1, a = 0)

fff <- apply(mu[ , 4:6], 1, mean)
ffff <- apply(mu[ , 7:9], 1, mean)
plot(log(fff), log(ffff), col = '#11111101', pch = 19, xlim = c(0, 10), ylim = c(0, 10))
abline(b = 1, a = 0)

ff <- apply(mu[ , 1:3], 1, mean)
ffff <- apply(mu[ , 7:9], 1, mean)
plot(log(ff), log(ffff), col = '#11111101', pch = 19, xlim = c(0, 10), ylim = c(0, 10))
abline(b = 1, a = 0)
dev.off()

# Plots for pp in Fline-------------------------------------------------------------------

mu <- ddsdlist$P@assays$data@listData$mu
round(tail(mu), digits = 7)

png("Contrasts between fp_fpp_fppp_magna.png", width=480, height=480)
par(mfrow = c(1, 3), title(main="Contrasts between fp_fpp_fppp_magna"))
fp <- apply(mu[ , 1:3], 1, mean)
fpp <- apply(mu[ , 4:6], 1, mean)
plot(log(fp), log(fpp), col = '#11111101', pch = 19, xlim = c(0, 10), ylim = c(0, 10))
abline(b = 1, a = 0, col="salmon")

fpp <- log(apply(mu[ , 4:6], 1, mean))
fppp <- log(apply(mu[ , 7:9], 1, mean))
plot(fpp, fppp, col = '#11111101', pch = 19, xlim = c(0, 12), ylim = c(0, 12))
abline(b = 1, a = 0, col="salmon")

fp <- apply(mu[ , 1:3], 1, mean)
fppp <- apply(mu[ , 7:9], 1, mean)
plot(log(fp), log(fppp), col = '#11111101', pch = 19, xlim = c(0, 10), ylim = c(0, 10))
abline(b = 1, a = 0, col="salmon")
dev.off()











if (!file.exists('D. magna_across_generations_ff_fff_ffff')) { dir.create('D. magna_across_generations_ff_fff_ffff') }

## UP and DOWN
DEgenenames <- list()
for (fi in 1:2) { # 1:2, 2 food quality   #length(ddsfoodlist)
  dirname <- sprintf('D. magna_across_generations_ff_fff_ffff/%s/', names(ddsfoodlist)[fi])
  cat(dirname) ; cat('\n')
  if (!file.exists(dirname)) { dir.create(dirname) }
  writeDifferntiallyExpressedGenes(ddsfoodlist[[fi]], dirname)
  DEgenenames[[fi]] <-
    read.table(sprintf('%s/%s',
                       dirname,
                       '03 Generation_gen2_vs_gen3 gene names all only diff exp by log2FC.txt'),  # this can be also changed to gen1_vs_gen3
               as.is = TRUE)[[1]]
}

tutaj ...
intersectgenes <-
  intersect(DEgenenames[[1]],(DEgenenames[[2]]))
length(intersectgenes) # shared genes 

## uniongenes <-
##   union(DEgenenames[[1]], union(DEgenenames[[2]], DEgenenames[[3]]))
names(DEgenenames) <- names(ddsfoodlist)
venn.diagram(DEgenenames, filename = 'D. magna_across_generations_ff_fff_ffff/venndiagram_gen2_vs_3.tiff', main="gen2_vs_3")
head(DEgenenames)
str(DEgenenames)

# variation for 3 generation at the same time? ----------------------------

if (!file.exists('D. magna_across_generations_ff_fff_ffff_razem')) { dir.create('D. magna_across_generations_ff_fff_ffff_razem') }

## UP and DOWN
DEgenenames1 <- list()
DEgenenames2 <- list()
for (fi in 1:2) { # 1:2, 2 food quality   #length(ddsfoodlist)
  dirname <- sprintf('D. magna_across_generations_ff_fff_ffff_razem/%s/', names(ddsfoodlist)[fi])
  cat(dirname) ; cat('\n')
  if (!file.exists(dirname)) { dir.create(dirname) }
  writeDifferntiallyExpressedGenes(ddsfoodlist[[fi]], dirname)
  DEgenenames1[[fi]] <-
    read.table(sprintf('%s/%s',
                       dirname,
                       '03 Generation_gen2_vs_gen3 gene names all only diff exp by log2FC.txt'),  
               as.is = TRUE)[[1]]
  DEgenenames2[[fi]] <-
  read.table(sprintf('%s/%s',
                       dirname,
                       '02 Generation_gen1_vs_gen3 gene names all only diff exp by log2FC.txt'),
               as.is = TRUE)[[1]]
  
}


tutaj ...
intersectgenes <-
  intersect(DEgenenames[[1]],(DEgenenames[[2]]))
length(intersectgenes) # shared genes 

## uniongenes <-
##   union(DEgenenames[[1]], union(DEgenenames[[2]], DEgenenames[[3]]))
names(DEgenenames) <- names(ddsfoodlist)
venn.diagram(DEgenenames, filename = 'D. magna_across_generations_ff_fff_ffff_razem/venndiagram_razem.tiff', main="gen2_vs_3")
head(DEgenenames)
str(DEgenenames)







