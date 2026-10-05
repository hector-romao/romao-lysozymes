library(ape)
library(geiger)
library(here)

## Data ------------------------------------------------------------
traits <- read.table(here("Morans_I", "data", "insecta_traits.tsv"),
                     header = TRUE, row.names = 1)
phy <- read.tree(here("Morans_I", "data", "insecta_tree.txt"))
head(traits)

# Check names and order
name.check(phy, traits)                      # expected: "OK"
traits <- traits[phy$tip.label, ]
identical(rownames(traits), phy$tip.label)   # expected: TRUE
anyNA(traits)                                # expected: FALSE

# Log-transformed counts, named and ordered by tree tip labels
log_c_type_count <- setNames(log(traits[, 1] + 1), phy$tip.label)
log_i_type_count <- setNames(log(traits[, 2] + 1), phy$tip.label)

hist(log_c_type_count, main = "", xlab = "log(C-type count + 1)")
hist(log_i_type_count, main = "", xlab = "log(I-type count + 1)")

## Moran's I (global) ---------------------------------------------
phy.cor <- vcv(phy, model = "Brownian", corr = TRUE)
diag(phy.cor) <- 0

I_c.type <- Moran.I(log_c_type_count, phy.cor)
I_i.type <- Moran.I(log_i_type_count, phy.cor)
I_c.type
I_i.type

## Correlogram -----------------------------------------------------
dist <- cophenetic(phy)
klim <- c(0, 50, 100, 200, 400, 600, 800, Inf)
x_pos <- c(50, 100, 200, 400, 600, 800, 1030)

Im_c <- Ip_c <- Im_i <- Ip_i <- numeric()

for (i in 1:(length(klim) - 1)) {
  W <- ifelse(dist > klim[i] & dist <= klim[i + 1], 1, 0)
  diag(W) <- 0
  Ic <- Moran.I(log_c_type_count, W)
  Ii <- Moran.I(log_i_type_count, W)
  Im_c[i] <- Ic$observed; Ip_c[i] <- Ic$p.value
  Im_i[i] <- Ii$observed; Ip_i[i] <- Ii$p.value
}

padj_c <- p.adjust(Ip_c, "bonferroni")
padj_i <- p.adjust(Ip_i, "bonferroni")

## Table -----------------------------------------------------------
fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))

tab_moran <- data.frame(
  Distance    = x_pos,
  I_Ctype     = round(Im_c, 3),
  p_adj_Ctype = fmt_p(padj_c),
  I_Itype     = round(Im_i, 3),
  p_adj_Itype = fmt_p(padj_i)
)
tab_moran
write.table(tab_moran, here("Morans_I", "moran_correlogram.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## Plot ------------------------------------------------------------
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1), mgp = c(2.3, 0.8, 0))

plot(x_pos, Im_c, type = "b", lwd = 2, cex = 1.2,
     pch = ifelse(padj_c < 0.05, 16, 1),
     xlab = "Cophenetic distance", ylab = "Moran's I", main = "C-type")
abline(h = 0, lty = 2, col = "grey50")

plot(x_pos, Im_i, type = "b", lwd = 2, cex = 1.2,
     pch = ifelse(padj_i < 0.05, 16, 1),
     xlab = "Cophenetic distance", ylab = "Moran's I", main = "I-type")
abline(h = 0, lty = 2, col = "grey50")

