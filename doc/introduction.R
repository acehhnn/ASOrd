## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ASOrd)

## -----------------------------------------------------------------------------
data("hold_data")
head(hold_data)

## -----------------------------------------------------------------------------
y1 = hold_data[hold_data$group == 1, 1]
y0 = hold_data[hold_data$group == 0, 1]

K = 3

interim_odds(K, k = 1, y1, y0)
interim_rank(K, k = 2, y1, y0, A = 6)

## -----------------------------------------------------------------------------
data("violate_data")
head(violate_data)

## -----------------------------------------------------------------------------
y1 = violate_data[violate_data$group == 1, 1]
y0 = violate_data[violate_data$group == 0, 1]

K = 3

interim_rank(K, k = 2, y1, y0, A = 6)

