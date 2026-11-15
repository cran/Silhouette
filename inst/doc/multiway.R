## ----setup, include = FALSE---------------------------------------------------
pkgs <- c("proxy", "drclust")
ok <- all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", eval = ok,
  fig.width = 7, fig.height = 2.9, dpi = 96, out.width = "100%"
)

## ----not-available, echo = FALSE, results = "asis", eval = !ok----------------
#  cat("This vignette needs the suggested packages",
#      paste0("'", pkgs, "'", collapse = ", "),
#      "to run its code. Install them to see the output.")

## ----library------------------------------------------------------------------
library(Silhouette)
library(drclust)

## ----data---------------------------------------------------------------------
set.seed(2025)
n <- 150; p <- 24
row_group <- sample(1:4, n, replace = TRUE)
col_group <- sample(1:3, p, replace = TRUE)
block_means <- matrix(c( 2, -1,  0,
                        -2,  1,  1,
                         0,  2, -2,
                         1, -2,  2), nrow = 4, byrow = TRUE)
Xb <- scale(1.5 * block_means[row_group, col_group] + matrix(rnorm(n * p), n, p))

## ----block-fun----------------------------------------------------------------
block_sil <- function(Xs, K, Q, seed) {
  set.seed(seed)
  fit <- doublekm(Xs, K = K, Q = Q, Rndstart = 10, prep = 0)
  U <- fit$U
  V <- fit$V
  C <- solve(crossprod(U), t(U) %*% Xs %*% V) %*% solve(crossprod(V))
  list(rows = Silhouette(proxy::dist(Xs, C %*% t(V))),       # n x K distances
       columns = Silhouette(proxy::dist(t(Xs), t(U %*% C))))  # p x Q distances
}
modes <- block_sil(Xb, K = 4, Q = 3, seed = 1)
extSilhouette(modes, dim_names = names(modes))

## ----grid---------------------------------------------------------------------
grid <- expand.grid(K = 2:5, Q = 2:4)
fits <- lapply(seq_len(nrow(grid)),
               function(g) block_sil(Xb, grid$K[g], grid$Q[g], seed = g))
ext <- lapply(fits, extSilhouette)
grid$ExS <- sapply(ext, `[[`, "ext_sil_width")
grid$rows <- sapply(ext, function(e) e$dim_table$avg_sil_width[1])
grid$columns <- sapply(ext, function(e) e$dim_table$avg_sil_width[2])
grid[order(-grid$ExS), ][1:4, ]

## ----grid-plot, fig.alt = "Three heat maps of average width over K from 2 to 5 row clusters and Q from 2 to 4 column clusters, one each for rows, columns and the extended silhouette. The row averages barely depend on Q; the column averages are highest at Q = 3; the extended silhouette follows the row panel closely.", echo = FALSE----
library(ggplot2)
long <- rbind(
  data.frame(grid[, c("K", "Q")], panel = "rows", value = grid$rows),
  data.frame(grid[, c("K", "Q")], panel = "columns", value = grid$columns),
  data.frame(grid[, c("K", "Q")], panel = "extended (ExS)", value = grid$ExS))
long$panel <- factor(long$panel, c("rows", "columns", "extended (ExS)"))
ggplot(long, aes(K, Q, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", value), colour = value > 0.5),
            size = 2.8, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#0b0b0b")) +
  scale_fill_gradient(low = "#cde2fb", high = "#104281", name = "Average width") +
  scale_x_continuous(breaks = 2:5, expand = c(0, 0)) +
  scale_y_continuous(breaks = 2:4, expand = c(0, 0)) +
  facet_wrap(~ panel) +
  labs(x = "Row clusters K", y = "Column clusters Q") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), legend.position = "bottom")

## ----grid-values, include = FALSE---------------------------------------------
top <- grid[which.max(grid$ExS), ]
w_rows <- n / (n + p)
grid$equal <- sapply(fits, function(f) extSilhouette(f, weights = "equal")$ext_sil_width)
top_eq <- grid[which.max(grid$equal), ]

## ----equal-weight-------------------------------------------------------------
grid$equal <- sapply(fits, function(f) extSilhouette(f, weights = "equal")$ext_sil_width)
grid[order(-grid$equal), ][1:3, c("K", "Q", "rows", "columns", "equal")]

## ----session, eval = TRUE-----------------------------------------------------
sessionInfo()

