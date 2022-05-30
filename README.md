# space-time-associations

## R

The R scripts have been implemented on R version 3.4.4.

### Loading R packages

```{r setup, include=FALSE}

# Load packages

source("R_packages.R")

```

## Data

### Loading bird and geographic data

```{r setup, include=FALSE}

# Load data from the French Breeding Bird Survey (STOC-EPS)

dataprp <- readRDS("raw_data/dataprp.rds")

# Load geographical information on sites

square_centroid <- readRDS("raw_data/square_centroid.rds")

```

## Spatial species associations

Details on the calculation of species associations from spatial data are available in Rigal et al. (2022). The difference here is that spatial associations are not calculated for each year but for the entire period (2001-2017).

### Loading functions

```{r setup, include=FALSE}

# Load functions for spatial associations

source("function_spatial_associations.R")

```


### Calculating species associations

Warning: it might take a while.

```{r}
registerDoParallel(cores = c(detectCores()-1)) # parallelise computation as this might take a while

dataperpoint <- dcast(dataprp,code_point+habit+zonebio~code_sp,fun.aggregate = sum,value.var="abond") # use data from all years

spatial_association <- ddply(dataperpoint, .(habit, zonebio), .fun=g, 
                           N.null=1000, occmin=0, 
                           .parallel = TRUE)

spatial_association$ses_init <- spatial_association$spatial_asso

spatial_association$spatial_asso[spatial_association$pval>0.05] <- NA

```




## References

https://doi.org/10.1111/ecog.01892

Rigal, S., Devictor, V., Gaüzère, P., Kéfi, S., Forsman, J. T., Kajanus, M. H., ... & Dakos, V. (2022). Biotic homogenisation in bird communities leads to large‐scale changes in species associations. Oikos, 2022(3), e08756. https://doi.org/10.1111/oik.08756