## ----setup, include = FALSE---------------------------------------------------
pkgs <- c("proxy", "ppclust", "mclust")
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
suppressPackageStartupMessages(library(mclust))  # Mclust() needs mclust attached
X <- as.matrix(iris[, 1:4])
set.seed(1)
fm <- ppclust::fcm(X, centers = 3, nstart = 5)   # fuzzy c-means
gm <- Mclust(X, G = 3, verbose = FALSE)          # Gaussian mixture

## ----widths-------------------------------------------------------------------
avg <- function(s) summary(s, print.summary = FALSE)$avg.width
c(PACS  = avg(softSilhouette(gm$z)),
  PPS   = avg(softSilhouette(gm$z, method = "medoid")),
  NLPPS = avg(softSilhouette(gm$z, prob_type = "nlpp", method = "medoid")),
  PDS   = avg(softSilhouette(gm$z, prob_type = "pd", method = "medoid")),
  CeS   = avg(cerSilhouette(gm$z)),
  DBS   = avg(dbSilhouette(gm$z)))

## ----identity-----------------------------------------------------------------
all.equal(softSilhouette(fm$u, method = "medoid")$sil_width,
          Silhouette(fm$d)$sil_width)

## ----exponent-----------------------------------------------------------------
sapply(c(a0 = 0, a0.5 = 0.5, a1 = 1, a2 = 2, a4 = 4), function(a)
  avg(softSilhouette(fm$u, average = "fuzzy", a = a)))

## ----compare------------------------------------------------------------------
P_fcm <- sqrt(fm$d)                                        # Euclidean distances
P_gmm <- as.matrix(proxy::dist(X, t(gm$parameters$mean)))  # to the component means
data.frame(
  algorithm = c("fuzzy c-means", "Gaussian mixture"),
  distance_medoid = c(avg(Silhouette(P_fcm, prob_matrix = fm$u)),
                      avg(Silhouette(P_gmm, prob_matrix = gm$z))),
  membership_pp = c(avg(softSilhouette(fm$u)), avg(softSilhouette(gm$z))),
  ARI_species = c(adjustedRandIndex(max.col(fm$u), iris$Species),
                  adjustedRandIndex(max.col(gm$z), iris$Species)))

## ----choose-k, fig.alt = "Line chart of average width against the number of clusters from 2 to 8 for four simulated groups: the k-means medoid width and the PACS and certainty widths from Gaussian mixtures and from fuzzy c-means."----
set.seed(42)
mu <- rbind(c(0, 0), c(6, 0), c(0, 6), c(6, 6))
sds <- c(0.8, 1.2, 1, 1.5)
id <- rep(1:4, c(120, 80, 100, 60))
Y <- mu[id, ] + matrix(rnorm(2 * length(id)), ncol = 2) * sds[id]

res <- t(sapply(2:8, function(k) {
  g <- Mclust(Y, G = k, modelNames = "VVV", verbose = FALSE)
  set.seed(k)
  km <- kmeans(Y, k, nstart = 20)
  fc <- ppclust::fcm(Y, centers = k)
  c(k = k,
    medoid_kmeans = avg(Silhouette(proxy::dist(Y, km$centers))),
    PACS_mixture = avg(softSilhouette(g$z)),
    CeS_mixture = avg(cerSilhouette(g$z)),
    PACS_fcm = avg(softSilhouette(fc$u)),
    CeS_fcm = avg(cerSilhouette(fc$u)),
    bic = g$bic)
}))
round(res, 3)
cols <- c("#2a78d6", "#eb6834", "#eda100", "#1baf7a", "#8a6fd1")
matplot(res[, "k"], res[, 2:6], type = "b", pch = c(19, 17, 2, 15, 0), lty = 1,
        col = cols, xlab = "Number of clusters", ylab = "Average width")
legend("bottomleft", colnames(res)[2:6], pch = c(19, 17, 2, 15, 0), col = cols,
       bty = "n", cex = 0.8)

## ----choose-k-values, include = FALSE-----------------------------------------
elbow <- function(v) res[which.max(-diff(v)), "k"]
best <- c(medoid = res[which.max(res[, "medoid_kmeans"]), "k"],
          bic = res[which.max(res[, "bic"]), "k"],
          apply(res[, c("PACS_mixture", "CeS_mixture", "PACS_fcm", "CeS_fcm")], 2, elbow))

## ----relabel------------------------------------------------------------------
c(nearest_mean_differs = sum(max.col(gm$z) != max.col(-P_gmm)),
  negative_widths = sum(Silhouette(P_gmm, prob_matrix = gm$z)$sil_width < 0))

## ----pd-----------------------------------------------------------------------
set.seed(4)
z <- c(rnorm(270, 0), rnorm(30, 2.5))
g2 <- Mclust(z, G = 2, modelNames = "E", verbose = FALSE)
pds <- softSilhouette(g2$z, prob_type = "pd", method = "medoid")
pds_map <- softSilhouette(g2$z, prob_type = "pd", method = "medoid",
                          cluster = max.col(g2$z))
c(density_label_differs = sum(pds$cluster != max.col(g2$z)),
  negative_default = sum(pds$sil_width < 0),
  negative_with_map_labels = sum(pds_map$sil_width < 0))

## ----db-zero------------------------------------------------------------------
data(thyroid, package = "mclust")
gt <- Mclust(scale(thyroid[, -1]), G = 3, verbose = FALSE)
db <- dbSilhouette(gt$z, average = "crisp")
avg(db)

## ----session, eval = TRUE-----------------------------------------------------
sessionInfo()

