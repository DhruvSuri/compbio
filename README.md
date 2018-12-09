CSE 549.01 Computational Biology - Fall 2018
Project Source Code

-----------------------------------------------

::: The Inferelator :::

Please follow the given steps in RStudio to run the project

1. Load the Data File titled inferelator_data.RData in RStudio

2. install.packages('devtools', dep=T)
3. install.packages('glmnet')
4. install.packages('lars')

5. library(devtools)
6. library(glmnet)
7. library(lars)

8. source('<path/to/project_inferelator.R>')

9. coefficients <- start_inferalator( e, min_choice="min+4se", group_count=999, alpha=0.8, n.boot=1, tau=10,
                r_threshold=Inf, r.filter=Inf, weighted=T, aic.filter=Inf, plot=T )
