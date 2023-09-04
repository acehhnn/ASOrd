# ASOrd-package
An R package for interim analysis with ordinal outcomes. 

- `interim_odds()` returns testing results based on proportional odds model approach
- `interim_rank()` returns testing results based on Wilcoxon rank sum statistics

## Installation
```
if (!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("acehhnn/ASOrd")
```

## Usage
For more datails, please read the vignette (`vignette("introduction", package = "ASOrd")`)
```
# Input data
data('hold_data') # a built-in data set generated from distributions satisfyting the proportional odds assumption
data('violate_data') # a built-in data set generated from distributions violating the proportional odds assumption

y1 = hold_data[hold_data$group == 1, 1]
y0 = hold_data[hold_data$group == 0, 1]

K = 3 # choose the number of interim analysis 

interim_odds(K, k = 1, y1, y0)
interim_rank(K, k = 2, y1, y0, A = 6)
```
