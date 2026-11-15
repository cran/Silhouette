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

## ----kmeans-------------------------------------------------------------------
X <- as.matrix(iris[, 1:4])
set.seed(1)
km <- kmeans(X, centers = 3, nstart = 20)
P <- proxy::dist(X, km$centers)   # 150 x 3 distances to the centers
sil <- Silhouette(P)
head(sil, 4)

## ----kmeans-summary-----------------------------------------------------------
s <- summary(sil, print.summary = FALSE)
s$avg.width      # average width over all observations
s$sil.sum        # size and average width of each cluster

## ----kmeans-plot, fig.alt = "Silhouette plot of 150 iris flowers in three k-means clusters. One cluster has widths mostly above 0.75; the other two have lower widths, several close to zero, with the average marked by a dashed line near 0.67."----
plot(sil)

## ----pac----------------------------------------------------------------------
sil_pac <- Silhouette(P, method = "pac")
all.equal(sil_pac$sil_width, sil$sil_width / (2 - sil$sil_width))

## ----cluster-arg--------------------------------------------------------------
hc <- cutree(hclust(dist(X), method = "average"), k = 3)
centers <- apply(X, 2, function(v) tapply(v, hc, mean))
sil_hc <- Silhouette(proxy::dist(X, centers), cluster = hc)
summary(sil_hc, print.summary = FALSE)$avg.width
sum(sil_hc$sil_width < 0)   # flowers closer to another group's centroid

## ----fcm-oneline--------------------------------------------------------------
set.seed(1)
sil_fcm <- Silhouette(prox_matrix = "d", prob_matrix = "u", average = "fuzzy",
                      clust_fun = ppclust::fcm, x = X, centers = 3)
summary(sil_fcm, print.summary = FALSE)$avg.width

## ----fcm-sqrt-----------------------------------------------------------------
set.seed(1)
fm <- ppclust::fcm(X, centers = 3)
summary(Silhouette(sqrt(fm$d), prob_matrix = fm$u, average = "fuzzy"),
        print.summary = FALSE)$avg.width

## ----mclust-------------------------------------------------------------------
suppressPackageStartupMessages(library(mclust))  # Mclust() needs mclust attached
gm <- Mclust(X, G = 3, verbose = FALSE)
sil_gm <- softSilhouette(gm$z)
summary(sil_gm, print.summary = FALSE)$avg.width

## ----calsil-------------------------------------------------------------------
calSilhouette(prob_matrix = gm$z)

## ----choose-k, fig.alt = "Line chart of the average simplified silhouette width for k-means on the iris data against the number of clusters from 2 to 6. The average is highest at 2 clusters and lower for every larger number."----
avg_width <- sapply(2:6, function(k) {
  set.seed(k)
  km <- kmeans(X, centers = k, nstart = 20)
  s <- Silhouette(proxy::dist(X, km$centers))
  summary(s, print.summary = FALSE)$avg.width
})
plot(2:6, avg_width, type = "b", pch = 19,
     xlab = "Number of clusters", ylab = "Average silhouette width")

## ----session, eval = TRUE-----------------------------------------------------
sessionInfo()

