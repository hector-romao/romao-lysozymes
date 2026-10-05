
library(ape)
library(phytools)
library(adephylo)
library(phylobase)
library(picante)
library(geiger)
library(TreeSim)
library(phylosignal)
library(caper)
library(phylolm)
library(here)
library(ggtree)


##Data
traits <- read.table(here("Morans_I","data", "insecta_traits.tsv"), h = T, row.names = 1)
head(traits)
phy <- read.tree(here("Morans_I","data", "insecta_tree.txt"), h = T, row.names = 1)
                     

c_type_count <- as.matrix(traits[, 1])
i_type_count <- as.matrix(traits[, 2])

hist(log(c_type_count) + 1, main = "", xlab = "log(C-type count + 1)")

# Nomeando as linhas com os rótulos da árvore
rownames(c_type_count) <- phy$tip.label
rownames(i_type_count) <- phy$tip.label

# Versões log-transformadas com nomes para uso nos testes
log_c_type_count <- log(c_type_count[, 1] + 1)
names(log_c_type_count) <- rownames(c_type_count)

log_i_type_count <- log(i_type_count[, 1] + 1)
names(log_i_type_count) <- rownames(i_type_count)




###MORAN i###

phy.cor <- vcv(phy, model = "Brownian", cor = T)
diag(phy.cor) <- 0

I_i.type <- Moran.I(log_i_type_count, phy.cor)
I_c.type <- Moran.I(log_c_type_count, phy.cor)




###Correlagram###
dist <- as.matrix(cophenetic(phy))
klim <- c(0, 50, 100, 200, 400, 600, 800, 1030)

Im_c <- numeric()
Ip <- numeric()

for (i in 1:(length(klim) - 1)) {
  W <- ifelse(dist > klim[i] & dist < klim[i + 1], 1, 0)
  diag(W) <- 0
  I_c <- Moran.I(log_c_type_count, W)
  Im_c[i] <- I$observed
  Ip[i] <- I_c$p.value
}

Im_i <- numeric()

for (i in 1:(length(klim) - 1)) {
  W <- ifelse(dist > klim[i] & dist < klim[i + 1], 1, 0)
  diag(W) <- 0
  I_i <- Moran.I(log_i_type_count, W)
  Im_i[i] <- I$observed
  Ip[i] <- I_i$p.value
}

par(mfrow = c(2, 1),
    mar = c(4, 4, 2, 1),
    mgp = c(2.3, 0.8, 0))

plot(klim[2:8], Im_c,
     type = "b",
     pch = 16,
     cex = 1.2,
     lwd = 2,
     xlab = "Cophenetic distance",
     ylab = "Moran's I",
     main = "C-type")


plot(klim[2:8], Im_i,
     type = "b",
     pch = 16,
     cex = 1.2,
     lwd = 2,
     xlab = "Cophenetic distance",
     ylab = "Moran's I",
     main = "I-type")


