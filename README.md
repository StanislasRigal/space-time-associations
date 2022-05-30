# space-time-associations

## R

The R scripts have been implemented on R version 3.4.4.

### Loading R packages

```{r setup, include=FALSE}

# Load packages

source("packages.R")

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

# Parallelise computation as this might take a while

registerDoParallel(cores = c(detectCores()-1))

# Prepare data

dataperpoint <- dcast(dataprp, code_point + habit + zonebio ~ code_sp, fun.aggregate = sum,value.var="abond")

# Run the function

spatial_association <- ddply(dataperpoint, .(habit, zonebio), .fun=main_fun_space, 
                           N.null=1000, occmin=0, 
                           .parallel = TRUE)
                           
# Keep all spatial associations                           

spatial_association$ses_init <- spatial_association$spatial_asso

# Set non-significant spatial associations to 0

spatial_association$spatial_asso[spatial_association$pval>0.05] <- NA

```

## Temporal species associations

### Loading functions

```{r setup, include=FALSE}

# Load functions for spatial associations

source("function_temporal_associations.R")

```

### Transforming data into time-series

```{r}
# Add 0 when site monitored but species not recorded

data_ts_0 <- as.data.frame(dataprp %>% group_by(code_square, zonebio, habit, code_sp, year) %>% summarize(abond=sum(abond), count=n()))

data_ts_0b <- ddply(data_ts_0, .(code_square, zonebio, habit), function(x){
  large_data <- dcast(droplevels(x), code_square + zonebio + habit + code_sp ~ year, fun.aggregate=sum, value.var="abond")
  long_data <- melt(large_data, id.vars=c("code_square","zonebio","habit","code_sp"))
  return(long_data)
})
names(data_ts_0b)[c(5,6)] <- c("year","abond")

data_ts_0b$year <- as.numeric(as.character(data_ts_0b$year))

data_ts_0b <- data_ts_0b[order(data_ts_0b[1],data_ts_0b[2],data_ts_0b[3],data_ts_0b[4],data_ts_0b[5]),]

# Calculate time-series length

data_ts_0b$seq_num <- cumsum(c(1, diff(data_ts_0b$year) != 1))

seq_length <- as.data.frame(data_ts_0b %>% group_by(code_square, zonebio, habit, code_sp, seq_num) %>% summarize(seq=n()))

data_ts_0b <- merge(data_ts_0b, seq_length, by=c("code_square", "zonebio", "habit", "code_sp", "seq_num"), all=T)

# Remove small time-series and species recorded less than 2 times (for using multispatial CCM, Clark et al. 2015)

data_ts_0c <- data_ts_0b[data_ts_0b$seq>5, c("code_square", "zonebio", "habit", "code_sp", "year", "abond","seq")]

to_remove <- data_ts_0c %>% group_by(code_square, zonebio, habit, code_sp) %>% summarize(nb_record = sum(sign(abond)))

data_ts_0c <- merge(data_ts_0c, to_remove, by=c("code_square", "zonebio", "habit", "code_sp"), all.x=T)

data_ts_0c <- data_ts_0c[data_ts_0c$nb_record>1,]

# Add coordinates of sites

data_ts_ready <- merge(data_ts_0c, square_centroid[,c("code_square", "long", "lat", "lon2", "lat2")], by="code_square", all.x=T)

# Detrend time-series if needed

data_ts_ready <- ddply(data_ts_ready, .(code_square, zonebio, habit, code_sp), .fun=detrend_data)
```







## References

Clark, A. T., Ye, H., Isbell, F., Deyle, E. R., Cowles, J., Tilman, G. D., & Sugihara, G. (2015). Spatial convergent cross mapping to detect causal relationships from short time series. Ecology, 96(5), 1174-1181. https://doi.org/10.1890/14-1479.1

Morueta‐Holme, N., Blonder, B., Sandel, B., McGill, B. J., Peet, R. K., Ott, J. E., ... & Svenning, J. C. (2016). A network approach for inferring species associations from co‐occurrence data. Ecography, 39(12), 1139-1150. https://doi.org/10.1111/ecog.01892

Rigal, S., Devictor, V., Gaüzère, P., Kéfi, S., Forsman, J. T., Kajanus, M. H., ... & Dakos, V. (2022). Biotic homogenisation in bird communities leads to large‐scale changes in species associations. Oikos, 2022(3), e08756. https://doi.org/10.1111/oik.08756