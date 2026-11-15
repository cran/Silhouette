## ----setup, include = FALSE---------------------------------------------------
pkgs <- c("proxy", "cluster", "mclust")
ok <- all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", eval = ok,
  fig.width = 6.5, fig.height = 3.6, dpi = 96, out.width = "100%"
)

## ----not-available, echo = FALSE, results = "asis", eval = !ok----------------
#  cat("This vignette needs the suggested packages",
#      paste0("'", pkgs, "'", collapse = ", "),
#      "to run its code. Install them to see the output.")

## ----library------------------------------------------------------------------
library(Silhouette)
X <- as.matrix(iris[, 1:4])
set.seed(1)
km <- kmeans(X, centers = 3, nstart = 20)
P <- as.matrix(proxy::dist(X, km$centers))   # distances to the 3 centers

## ----classical----------------------------------------------------------------
s_simple <- Silhouette(P)$sil_width
s_classic <- cluster::silhouette(km$cluster, dist(X))[, "sil_width"]
c(correlation = cor(s_simple, s_classic),
  mean_simplified = mean(s_simple),
  mean_simplified_pac = mean(Silhouette(P, method = "pac")$sil_width),
  mean_classical = mean(s_classic))

## ----nonnegative--------------------------------------------------------------
min(s_simple)

## ----squared------------------------------------------------------------------
all.equal(Silhouette(P^2)$sil_width, s_simple * (2 - s_simple))

## ----relabel------------------------------------------------------------------
hc <- cutree(hclust(dist(X), method = "average"), k = 3)
centers <- apply(X, 2, function(v) tapply(v, hc, mean))
P_hc <- proxy::dist(X, centers)
c(reassigned_without_cluster = sum(Silhouette(P_hc)$cluster != hc),
  negative_with_cluster = sum(Silhouette(P_hc, cluster = hc)$sil_width < 0))

## ----negative-----------------------------------------------------------------
suppressPackageStartupMessages(library(mclust))
gm <- Mclust(X, G = 3, verbose = FALSE)
P_gmm <- as.matrix(proxy::dist(X, t(gm$parameters$mean)))
sum(Silhouette(P_gmm, prob_matrix = gm$z)$sil_width < 0)

## ----session, eval = TRUE-----------------------------------------------------
sessionInfo()

