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


### Calculating species spatial associations

Warning: it might take a while.

```{r}

# Parallelise computation as this might take a while

registerDoParallel(cores = c(detectCores()-1))

# Prepare data

dataperpoint <- dcast(dataprp, code_point + habit + zonebio ~ code_sp, fun.aggregate = sum,value.var="abond")

# Run the function

spatial_associations <- ddply(dataperpoint, .(habit, zonebio), .fun=main_fun_space, 
                           N.null=1000, occmin=0, 
                           .parallel = TRUE)
                           
# Keep all spatial associations                           

spatial_associations$ses_init <- spatial_associations$spatial_asso

# Set non-significant spatial associations to 0

spatial_associations$spatial_asso[spatial_associations$pval>0.05] <- NA

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

### Assess significance of temporal associations

```{r}

temporal_associations <- ddply(data_ts_ready, .(zonebio, habit), .fun=multisp_CCM, .progress="text")

temporal_associations1 <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="alpine",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")

temporal_associations2a <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="atlantic" & data_ts_ready$habit=="A_1",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")
temporal_associations2b <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="atlantic" & data_ts_ready$habit=="A_3",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")
temporal_associations2c <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="atlantic" & data_ts_ready$habit=="D_1",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")
temporal_associations2d <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="atlantic" & data_ts_ready$habit=="D_2",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")
temporal_associations2e <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="atlantic" & data_ts_ready$habit=="D_3",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")
temporal_associations2f <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="atlantic" & data_ts_ready$habit=="D_4",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")
temporal_associations2g <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="atlantic" & data_ts_ready$habit=="E_1",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")
temporal_associations2h <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="atlantic" & data_ts_ready$habit=="E_2",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")
temporal_associations2i <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="atlantic" & data_ts_ready$habit=="E_3",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")
temporal_associations2j <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="atlantic" & data_ts_ready$habit=="F_1",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")
temporal_associations2k <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="atlantic" & data_ts_ready$habit %in% c("A_2","B_1","B_2","B_3","C_1","C_2","C_4","D_5","G_1"),]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")


temporal_associations3a <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="continental" & data_ts_ready$habit=="A_1",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")
temporal_associations3b <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="continental" & data_ts_ready$habit=="A_3",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")
temporal_associations3c <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="continental" & data_ts_ready$habit=="D_1",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")
temporal_associations3d <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="continental" & data_ts_ready$habit=="D_2",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")
temporal_associations3e <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="continental" & data_ts_ready$habit=="D_3",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")
temporal_associations3f <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="continental" & data_ts_ready$habit=="D_4",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")
temporal_associations3g <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="continental" & data_ts_ready$habit=="E_3",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")
temporal_associations3h <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="continental" & data_ts_ready$habit=="F_1",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")
temporal_associations3i <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="continental" & data_ts_ready$habit %in% c("A_2","B_1","B_2","B_3","C_1","C_2","C_4","D_5","E_1",  "E_2","G_1"),]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")

temporal_associations4 <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="mediterranean",]), .(zonebio, habit), .fun=multisp_CCM, .progress="text")

#

temporal_associations1 <- readRDS("output/temporal_associations1.rds")
temporal_associations2a <- readRDS("output/temporal_associations2a.rds")
temporal_associations2b <- readRDS("output/temporal_associations2b.rds")
temporal_associations2c <- readRDS("output/temporal_associations2c.rds")
temporal_associations2d <- readRDS("output/temporal_associations2d.rds")
temporal_associations2e <- readRDS("output/temporal_associations2e.rds")
temporal_associations2f <- readRDS("output/temporal_associations2f.rds")
temporal_associations2g <- readRDS("output/temporal_associations2g.rds")
temporal_associations2h <- readRDS("output/temporal_associations2h.rds")
temporal_associations2i <- readRDS("output/temporal_associations2i.rds")
temporal_associations2j <- readRDS("output/temporal_associations2j.rds")
temporal_associations2k <- readRDS("output/temporal_associations2k.rds")
temporal_associations3a <- readRDS("output/temporal_associations3a.rds")
temporal_associations3b <- readRDS("output/temporal_associations3b.rds")
temporal_associations3c <- readRDS("output/temporal_associations3c.rds")
temporal_associations3d <- readRDS("output/temporal_associations3d.rds")
temporal_associations3e <- readRDS("output/temporal_associations3e.rds")
temporal_associations3f <- readRDS("output/temporal_associations3f.rds")
temporal_associations3g <- readRDS("output/temporal_associations3g.rds")
temporal_associations3h <- readRDS("output/temporal_associations3h.rds")
temporal_associations3i <- readRDS("output/temporal_associations3i.rds")
temporal_associations4 <- readRDS("output/temporal_associations4.rds")

temporal_associations1$zonebio <- as.character(temporal_associations1$zonebio)
temporal_associations2a$zonebio <- as.character(temporal_associations2a$zonebio)
temporal_associations2b$zonebio <- as.character(temporal_associations2b$zonebio)
temporal_associations2c$zonebio <- as.character(temporal_associations2c$zonebio)
temporal_associations2d$zonebio <- as.character(temporal_associations2d$zonebio)
temporal_associations2e$zonebio <- as.character(temporal_associations2e$zonebio)
temporal_associations2f$zonebio <- as.character(temporal_associations2f$zonebio)
temporal_associations2g$zonebio <- as.character(temporal_associations2g$zonebio)
temporal_associations2h$zonebio <- as.character(temporal_associations2h$zonebio)
temporal_associations2i$zonebio <- as.character(temporal_associations2i$zonebio)
temporal_associations2j$zonebio <- as.character(temporal_associations2j$zonebio)
temporal_associations2k$zonebio <- as.character(temporal_associations2k$zonebio)
temporal_associations3a$zonebio <- as.character(temporal_associations3a$zonebio)
temporal_associations3b$zonebio <- as.character(temporal_associations3b$zonebio)
temporal_associations3c$zonebio <- as.character(temporal_associations3c$zonebio)
temporal_associations3d$zonebio <- as.character(temporal_associations3d$zonebio)
temporal_associations3e$zonebio <- as.character(temporal_associations3e$zonebio)
temporal_associations3f$zonebio <- as.character(temporal_associations3f$zonebio)
temporal_associations3g$zonebio <- as.character(temporal_associations3g$zonebio)
temporal_associations3h$zonebio <- as.character(temporal_associations3h$zonebio)
temporal_associations3i$zonebio <- as.character(temporal_associations3i$zonebio)
temporal_associations4$zonebio <- as.character(temporal_associations4$zonebio)

temporal_associations1$habit <- as.character(temporal_associations1$habit)
temporal_associations2a$habit <- as.character(temporal_associations2a$habit)
temporal_associations2b$habit <- as.character(temporal_associations2b$habit)
temporal_associations2c$habit <- as.character(temporal_associations2c$habit)
temporal_associations2d$habit <- as.character(temporal_associations2d$habit)
temporal_associations2e$habit <- as.character(temporal_associations2e$habit)
temporal_associations2f$habit <- as.character(temporal_associations2f$habit)
temporal_associations2g$habit <- as.character(temporal_associations2g$habit)
temporal_associations2h$habit <- as.character(temporal_associations2h$habit)
temporal_associations2i$habit <- as.character(temporal_associations2i$habit)
temporal_associations2j$habit <- as.character(temporal_associations2j$habit)
temporal_associations2k$habit <- as.character(temporal_associations2k$habit)
temporal_associations3a$habit <- as.character(temporal_associations3a$habit)
temporal_associations3b$habit <- as.character(temporal_associations3b$habit)
temporal_associations3c$habit <- as.character(temporal_associations3c$habit)
temporal_associations3d$habit <- as.character(temporal_associations3d$habit)
temporal_associations3e$habit <- as.character(temporal_associations3e$habit)
temporal_associations3f$habit <- as.character(temporal_associations3f$habit)
temporal_associations3g$habit <- as.character(temporal_associations3g$habit)
temporal_associations3h$habit <- as.character(temporal_associations3h$habit)
temporal_associations3i$habit <- as.character(temporal_associations3i$habit)
temporal_associations4$habit <- as.character(temporal_associations4$habit)

temporal_associations1$species <- as.character(temporal_associations1$species)
temporal_associations2a$species <- as.character(temporal_associations2a$species)
temporal_associations2b$species <- as.character(temporal_associations2b$species)
temporal_associations2c$species <- as.character(temporal_associations2c$species)
temporal_associations2d$species <- as.character(temporal_associations2d$species)
temporal_associations2e$species <- as.character(temporal_associations2e$species)
temporal_associations2f$species <- as.character(temporal_associations2f$species)
temporal_associations2g$species <- as.character(temporal_associations2g$species)
temporal_associations2h$species <- as.character(temporal_associations2h$species)
temporal_associations2i$species <- as.character(temporal_associations2i$species)
temporal_associations2j$species <- as.character(temporal_associations2j$species)
temporal_associations2k$species <- as.character(temporal_associations2k$species)
temporal_associations3a$species <- as.character(temporal_associations3a$species)
temporal_associations3b$species <- as.character(temporal_associations3b$species)
temporal_associations3c$species <- as.character(temporal_associations3c$species)
temporal_associations3d$species <- as.character(temporal_associations3d$species)
temporal_associations3e$species <- as.character(temporal_associations3e$species)
temporal_associations3f$species <- as.character(temporal_associations3f$species)
temporal_associations3g$species <- as.character(temporal_associations3g$species)
temporal_associations3h$species <- as.character(temporal_associations3h$species)
temporal_associations3i$species <- as.character(temporal_associations3i$species)
temporal_associations4$species <- as.character(temporal_associations4$species)

temporal_associations <- rbind(temporal_associations1,temporal_associations2a,temporal_associations2b,
temporal_associations2c,temporal_associations2d,temporal_associations2e,temporal_associations2f,temporal_associations2g,temporal_associations2h,temporal_associations2i,temporal_associations2j,temporal_associations2k,temporal_associations3a,temporal_associations3b,temporal_associations3c,temporal_associations3d,temporal_associations3e,temporal_associations3f,temporal_associations3g,temporal_associations3h,temporal_associations3i,temporal_associations4)

#

data_ts_ready <- readRDS("output/data_ts_ready.rds")
data_ts_ready$code_square <- sprintf("%06d", data_ts_ready$code_square)

spatial_associations <- readRDS("output/spatial_association.rds")

temporal_associations <- readRDS("output/temporal_associations.rds")


# E dim

Edim_temporal_associations1 <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="alpine",]), .(zonebio, habit), .fun=Edim_multisp_CCM, .progress="text")
Edim_temporal_associations2 <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="atlantic",]), .(zonebio, habit), .fun=Edim_multisp_CCM, .progress="text")
Edim_temporal_associations3 <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="continental",]), .(zonebio, habit), .fun=Edim_multisp_CCM, .progress="text")
Edim_temporal_associations4 <- ddply(droplevels(data_ts_ready[data_ts_ready$zonebio=="mediterranean",]), .(zonebio, habit), .fun=Edim_multisp_CCM, .progress="text")

Edim_temporal_associations <- rbind(Edim_temporal_associations1,Edim_temporal_associations2,Edim_temporal_associations3,Edim_temporal_associations4)
names(Edim_temporal_associations)[4] <- "E"
Edim_temporal_associations$species <- substr(Edim_temporal_associations$species,3,21)

temporal_associations_Edim <- merge(temporal_associations, Edim_temporal_associations, by=c("zonebio","habit","species"), all.x=T)

temporal_associations_signif <- droplevels(na.omit(temporal_associations_Edim[temporal_associations_Edim$rho<0.05,]))
temporal_associations_signif$spA <- substr(temporal_associations_signif$species,14,19)
temporal_associations_signif$spB <- substr(temporal_associations_signif$species,1,6)

```

### Quantifying the temporal association between species

```{r}

temporal_associations_quanti <- ddply(temporal_associations_signif, .(zonebio, habit, species), .fun=smap_fun_signif, data_ccm=data_ts_ready, .progress="text")

temporal_associations_Edim <- readRDS("output/temporal_associations_Edim.rds")
temporal_associations_quanti <- readRDS("output/temporal_associations_quanti.rds")

temporal_associations_quanti_sum <- data.frame(temporal_associations_quanti %>% group_by(zonebio,habit,species) %>% summarize(temp_int=mean(spB_smap,na.rm=T), sd_temp_int=sd(spB_smap, na.rm=T), temp_int2=mean(spB_smap2,na.rm=T), sd_temp_int2=sd(spB_smap2, na.rm=T)))

temporal_associations_signif <- merge(temporal_associations_signif,temporal_associations_quanti_sum, by=c("zonebio","habit","species"))

```

### Merge spatial and temporal associations

```{r}

sp_tp_association <- merge(temporal_associations_signif, spatial_associations, by.x=c("habit","zonebio","spA","spB"), by.y=c("habit","zonebio","spi","spj"), all.x=T)

```




# Species distances

### Functional distance

```{r}
# loading the trait matrix from Storchová & Hořák (2018) (https://doi.org/10.1111/geb.12709)

trait <- read.csv("raw_data/life_history_bird2.csv", header = TRUE)

species_name_data <- read.csv("raw_data/species_name_data.csv", header=T)
species_name_data$scientific_name2 <- sub(" ","_",species_name_data$scientific_name)

# removing duplicated traits and traits with multiple NA

table.na <- apply(trait,2,function(x){sum(is.na(x))})

row.names(trait) <- trait$Species

trait <- na.omit(trait[,c("LengthU_MEAN","WingU_MEAN","TailU_MEAN","BillU_MEAN","TarsusU_MEAN","WeightU_MEAN","Sexual.dimorphism","Clutch_MEAN","Broods.per.year","EggL_MEAN","EggW_MEAN","Egg_MASS","Young","Association.during.nesting","Nest.type","Nest.building","Mating.system","Incubation.period","Incubation.sex","Association.outside.the.breeding.season","Territoriality","Sedentary","Facultative.migrant","Short.distance.migrant","Long.distance.migrant","Deciduous.forest","Coniferous.forest","Woodland","Shrub","Savanna","Tundra","Grassland","Mountain.meadows","Reed","Swamps","Desert","Freshwater","Marine","Rocks","Human.settlements","Folivore_Y","Frugivore_Y","Granivore_Y","Arthropods_Y","Other.invertebrates_Y","Fish_Y","Other.vertebrates_Y","Carrion_Y","Omnivore_Y","Folivore_B","Frugivore_B","Granivore_B","Arthropods_B","Other.invertebrates_B","Fish_B","Other.vertebrates_B","Carrion_B","Omnivore_B")])


# transforming into a binomial format the relevant variables

trait[,c("Territoriality","Sedentary","Facultative.migrant","Short.distance.migrant","Long.distance.migrant","Deciduous.forest","Coniferous.forest","Woodland","Shrub","Savanna","Tundra","Grassland","Mountain.meadows","Reed","Swamps","Desert","Freshwater","Marine","Rocks","Human.settlements","Folivore_Y","Frugivore_Y","Granivore_Y","Arthropods_Y","Other.invertebrates_Y","Fish_Y","Other.vertebrates_Y","Carrion_Y","Omnivore_Y","Folivore_B","Frugivore_B","Granivore_B","Arthropods_B","Other.invertebrates_B","Fish_B","Other.vertebrates_B","Carrion_B","Omnivore_B")] <- apply(trait[,c("Territoriality","Sedentary","Facultative.migrant","Short.distance.migrant","Long.distance.migrant","Deciduous.forest","Coniferous.forest","Woodland","Shrub","Savanna","Tundra","Grassland","Mountain.meadows","Reed","Swamps","Desert","Freshwater","Marine","Rocks","Human.settlements","Folivore_Y","Frugivore_Y","Granivore_Y","Arthropods_Y","Other.invertebrates_Y","Fish_Y","Other.vertebrates_Y","Carrion_Y","Omnivore_Y","Folivore_B","Frugivore_B","Granivore_B","Arthropods_B","Other.invertebrates_B","Fish_B","Other.vertebrates_B","Carrion_B","Omnivore_B")], 2, function(x) replace(x,x>0,1))

trait[,c("Sexual.dimorphism","Territoriality","Sedentary","Facultative.migrant","Short.distance.migrant","Long.distance.migrant","Deciduous.forest","Coniferous.forest","Woodland","Shrub","Savanna","Tundra","Grassland","Mountain.meadows","Reed","Swamps","Desert","Freshwater","Marine","Rocks","Human.settlements","Folivore_Y","Frugivore_Y","Granivore_Y","Arthropods_Y","Other.invertebrates_Y","Fish_Y","Other.vertebrates_Y","Carrion_Y","Omnivore_Y","Folivore_B","Frugivore_B","Granivore_B","Arthropods_B","Other.invertebrates_B","Fish_B","Other.vertebrates_B","Carrion_B","Omnivore_B","Young","Association.during.nesting","Nest.type","Nest.building","Mating.system","Incubation.sex","Association.outside.the.breeding.season")] <- lapply(trait[,c("Sexual.dimorphism","Territoriality","Sedentary","Facultative.migrant","Short.distance.migrant","Long.distance.migrant","Deciduous.forest","Coniferous.forest","Woodland","Shrub","Savanna","Tundra","Grassland","Mountain.meadows","Reed","Swamps","Desert","Freshwater","Marine","Rocks","Human.settlements","Folivore_Y","Frugivore_Y","Granivore_Y","Arthropods_Y","Other.invertebrates_Y","Fish_Y","Other.vertebrates_Y","Carrion_Y","Omnivore_Y","Folivore_B","Frugivore_B","Granivore_B","Arthropods_B","Other.invertebrates_B","Fish_B","Other.vertebrates_B","Carrion_B","Omnivore_B","Young","Association.during.nesting","Nest.type","Nest.building","Mating.system","Incubation.sex","Association.outside.the.breeding.season")], factor)

trait[,c("Incubation.period")] <- as.numeric(trait[,c("Incubation.period")])

# Construct a distance matrix using gower metric

mat.trait <- trait
mat.trait$code_sp <- species_name_data$code_sp[match(row.names(trait),species_name_data$scientific_name2)]
mat.trait <- subset(mat.trait, code_sp %in% levels(droplevels(dataprp$code_sp)))
mat.dist <- daisy(mat.trait[,1:58], metric="gower", type=list(asymm=c(22:58)))
attr(mat.dist, "Labels") <- mat.trait$code_sp

# selecting clustering algorithm

library(clue)

clust1 <- hclust(mat.dist, method = "ward.D")
clust2 <- hclust(mat.dist, method = "ward.D2")
clust3 <- hclust(mat.dist, method = "single")
clust4 <- hclust(mat.dist, method = "complete")
clust5 <- hclust(mat.dist, method = "average")
clust6 <- hclust(mat.dist, method = "mcquitty")
#clust7 <- hclust(mat.dist, method = "median") # this two methods are reported to crash
#clust8 <- hclust(mat.dist, method = "centroid")

list.clust <- list(clust1,clust2,clust3,clust4,clust5,clust6)

dissim <- data.frame(ward.D=rep(NA,(2^length(list.clust)-1)),ward.D2=rep(NA,(2^length(list.clust)-1)),single=rep(NA,(2^length(list.clust)-1)),
                   complete=rep(NA,(2^length(list.clust)-1)),average=rep(NA,(2^length(list.clust)-1)),mcquitty=rep(NA,(2^length(list.clust)-1)),
                   median=rep(NA,(2^length(list.clust)-1)),centroid=rep(NA,(2^length(list.clust)-1)),dissimilarity=rep(NA,(2^length(list.clust)-1)))

nb <- 0
for(i in 1:length(list.clust)){
  combinaison <- combn(length(list.clust),i)
  for(j in 1:ncol(combinaison)){
    nb <- nb+1
    print(nb)
    ensemble<-cl_ensemble(list=list.clust[combinaison[,j]])
    if(i==1){consensus <- list.clust[combinaison[,j]]
    }else{consensus <- cl_consensus(ensemble, method = "euclidean")}
    dissim[nb,c(combinaison[,j])] <- combinaison[,j]
    dissim[nb,9] <- cl_dissimilarity(mat.dist, consensus, method = "spectral")
  }
}

dissim[which.min(dissim$dissimilarity),]
dissim[order(dissim$dissimilarity),]
consensus <- cl_consensus(cl_ensemble(list=list.clust[c(5)]))
plot(consensus)

sigma <- var(mat.dist) + var(cophenetic(clust5))
threshold <- 2*sqrt(nrow(as.matrix(mat.dist))*sigma)

# as sigma < threshold we can continue and extract the cophenetic distances

coph.dist <- cophenetic(consensus)
coph.dist2 <- data.frame(t(combn(attr(coph.dist, "Labels"),2)), as.numeric(coph.dist))
names(coph.dist2) <- c("sp1","sp2","distance")

```

### Phylogenetic distance

```{r}

# Data from Ericson (https://birdtree.org/, info on https://doi.org/10.1016/j.cub.2014.03.011) sequenced species, 10 000 trees with 6670OTUS each

tree <- read.nexus("raw_data/output_10000.nex")

tree_cons <- consensus(tree[c(1:2607,2609:8485,8487:10000)], p = 1, check.labels = TRUE)
tree_dist <- cophenetic.phylo(tree_cons)
tree_dist_long <- tabularize(tree_dist, name="distance")
tree_dist_long$distance2 <- tree_dist_long$distance*1e308*1e3 # sometime length issue as some species are very close

tree_dist_long$spi2 <- species_name_data$code_sp[match(tree_dist_long$spi, species_name_data$scientific_name2)]
ch <- levels(tree_dist_long$spi)[which(!(levels(tree_dist_long$spi) %in% levels(as.factor(species_name_data$scientific_name2))))]
for(i in 1:length(ch)){
  tree_dist_long$spi2[which(tree_dist_long$spi %in% ch[i])] <- species_name_data$code_sp[match(str_replace(ch[i], "_", " "), species_name_data$scientific_name)]
}

tree_dist_long$spj2 <- species_name_data$code_sp[match(tree_dist_long$spj, species_name_data$scientific_name2)]
ch <- levels(tree_dist_long$spj)[which(!(levels(tree_dist_long$spj) %in% levels(as.factor(species_name_data$scientific_name2))))]
for(i in 1:length(ch)){
  tree_dist_long$spj2[which(tree_dist_long$spj %in% ch[i])]<-species_name_data$code_sp[match(str_replace(ch[i], "_", " "), species_name_data$scientific_name)]
}

d <- tree_dist_long
d <- data.frame(d, pair=paste(d$spi2,d$spj2,sep="_"))
d <- data.frame(d, pair2=paste(d$spj2,d$spi2,sep="_"))
```

### Niche overlap
```{r}
habitat_sp <- as.data.frame(dataprp %>% group_by(code_sp, habit) %>% summarize(ab=sum(abond)))
habitat_sp <- droplevels(subset(habitat_sp, !(habit=="NA_NA")))
habitat_sp_wide <- dcast(habitat_sp, code_sp~habit)
habitat_sp_wide[is.na(habitat_sp_wide)] <- 0
habitat_sp <- melt(habitat_sp_wide)
names(habitat_sp)[c(2:3)] <- c("habit","ab")
habitat_sp_count <- as.data.frame(habitat_sp %>% group_by(code_sp) %>% summarize(count=sum(ab)))
habitat_sp <- merge(habitat_sp,habitat_sp_count, by="code_sp",all=T)
habitat_sp$freq <- habitat_sp$ab/habitat_sp$count

habitat_sp_wide <- dcast(habitat_sp[,c("code_sp","habit","ab")], code_sp~habit)
row.names(habitat_sp_wide) <- habitat_sp_wide$code_sp 
habitat_sp_wide$code_sp <- NULL

# chi2 distance as we want to use frequency and comparing the distribution profile of each species in the different habitats, independently from the abundance.

habitat_dist_chi <- as.data.frame(habitat_sp_wide) 
habitat_dist_chi <- dudi.coa(habitat_dist_chi,scan=F)
habitat_dist_chi <- dist.dudi(habitat_dist_chi)
habitat_dist_chi <- tabularize(as.matrix(habitat_dist_chi))
habitat_dist_chi$assoc <- 1-habitat_dist_chi$assoc
names(habitat_dist_chi)[3]<-"dist"
```

### Specialisation distance

```{r}
# from Godet et al. (2015) (https://doi.org/10.1111/geb.12266)

sxi <- read.csv("raw_data/SXI_publi.csv")
sxi$code_sp <- as.factor(sxi$code_sp)

specialisation.dist <- function(x){
  sp <- levels(droplevels(x$code_sp))
  d_SGI <- data.frame(spj=sxi$code_sp[!(sxi$code_sp==sp)], dSGI=abs((sxi$SGIo[!(sxi$code_sp==sp)]-x$SGIo)))
  return(d_SGI)
}
specialisation_dist <- ddply(sxi, .(code_sp), .fun=specialisation.dist)

specialisation_dist$pair <- as.factor(paste0(specialisation_dist$code_sp,sep="_",specialisation_dist$spj))
names(specialisation_dist)[1] <- "spi"
```


### Combining species associations with species distances

```{r}

association_distance <- merge(sp_tp_association,coph.dist2[,c("distance","sp1","sp2")], by.x=c("spA","spB"), by.y=c("sp1","sp2"), all.x=T)
association_distance <- merge(association_distance,coph.dist2[,c("distance","sp1","sp2")], by.x=c("spA","spB"), by.y=c("sp2","sp1"), all.x=T)
association_distance$functional_distance <- rep(NA, nrow(association_distance))
association_distance$functional_distance[!is.na(association_distance$distance.x)] <- association_distance$distance.x[!is.na(association_distance$distance.x)]
association_distance$functional_distance[!is.na(association_distance$distance.y)] <- association_distance$distance.y[!is.na(association_distance$distance.y)]
association_distance$distance.x <- association_distance$distance.y <- NULL

association_distance <- merge(association_distance,d[,c("distance","spi2","spj2")], by.x=c("spA","spB"), by.y=c("spi2","spj2"), all.x=T)
names(association_distance)[which(names(association_distance)=="distance")] <- "phylogenetic_distance"
association_distance$phylogenetic_distance[which(is.na(association_distance$phylogenetic_distance))] <- 0

association_distance <- merge(association_distance,specialisation_dist[,c("dSGI","spi","spj")], by.x=c("spA","spB"), by.y=c("spi","spj"), all.x=T)
names(association_distance)[which(names(association_distance)=="dSGI")] <- "specialisation_distance"

association_distance <- merge(association_distance,habitat_dist_chi, by.x=c("spA","spB"), by.y=c("spi","spj"), all.x=T)
association_distance$niche_overlap_distance <- -association_distance$dist
association_distance$dist <- NULL
```

### Combining species association by species pair
```{r}
asso_to_plot <- association_distance
asso_to_plot$spatial_asso_scaled <- scale(asso_to_plot$spatial_asso,center=F)
asso_to_plot$temp_asso_scaled <- scale(asso_to_plot$temp_int2,center=F)

ggplot(droplevels(asso_to_plot), aes(x=spatial_asso_scaled, y=temp_asso_scaled))+
  geom_point(data=asso_to_plot,alpha=0.1,col="black")+
  scale_color_viridis(alpha = 1, direction = 1,discrete = TRUE, option = "D")+
  scale_x_continuous(name="Spatial associations (SES)", limits=c(-5,5))+
  scale_y_continuous(name="Temporal associations", limits=c(-5,5))+
  theme_classic() + theme(text=element_text(size=20),legend.position="none")
```



# Testing the relationship between species associations and niche distances

## For associations significant both in space and time

### Preparing data

```{r}
# remove data from time-series < 500 time steps
data_asso_dist <- droplevels(association_distance)
lev_sp_t <- unique(data_asso_dist$pair[!is.na(data_asso_dist$spatial_asso)])

# save information about combination with association or no association
num_not_zero <- which(!is.na(data_asso_dist$temp_asso2) & data_asso_dist$spatial_asso2!=0)
num_zero <- which(is.na(data_asso_dist$temp_asso2) | data_asso_dist$spatial_asso2==0)
num_not_zero2 <- which(!is.na(data_asso_dist$temp_asso2))
num_zero2 <- which(is.na(data_asso_dist$temp_asso2))
num_not_zero2b <- which(data_asso_dist$spatial_asso2!=0)
num_zero2b <- which(data_asso_dist$spatial_asso2==0)

# add small variation around zero and select only part of the dataset to allow broom::tidy to work
data_asso_dist$temp_asso2_save<-data_asso_dist$temp_asso2
data_asso_dist$spatial_asso2_save<-data_asso_dist$spatial_asso2
data_asso_dist$temp_asso2[num_zero]<-0+runif(length(data_asso_dist$temp_asso2[num_zero]),min=-0.0001,max=0.0001)
data_asso_dist$spatial_asso2[num_zero]<-0+runif(length(data_asso_dist$temp_asso2[num_zero]),min=-0.0001,max=0.0001)
data_asso_dist$temp_asso2_save[num_zero2]<-0+runif(length(data_asso_dist$temp_asso2[num_zero2]),min=-0.0001,max=0.0001)
data_asso_dist$spatial_asso2[num_zero2b]<-0+runif(length(data_asso_dist$temp_asso2[num_zero2b]),min=-0.0001,max=0.0001)
set.seed(1) 
num_zero3<-sample(num_zero,5000)
data_asso_dist2<-data_asso_dist[c(num_not_zero,num_zero3),]
```

### Functional distance (and other niche distances)
```{r}
# quantile regression for spatial_asso2 and functional_distance
Quantr<-rq(data=droplevels(data_asso_dist2),
           tau= seq(0.01,0.99, by=0.01),
           formula = spatial_asso2 ~  functional_distance) # or phylogenetic_distance, specialisation_distance, niche_overlap_distance
sum_rq<-summary(Quantr)
plot(sum_rq, ols=FALSE)
qr_spatial_asso_functional_distance<-broom::tidy(Quantr,se.type = "boot")
qr_spatial_asso_functional_distance<-droplevels(subset(qr_spatial_asso_functional_distance, term=="functional_distance"))

# quantile regression for temp_asso2 and functional_distance
Quantr<-rq(data=droplevels(data_asso_dist2),
           tau= seq(0.01,0.99, by=0.01),
           formula = temp_asso2 ~  functional_distance) # or phylogenetic_distance, specialisation_distance, niche_overlap_distance
sum_rq<-summary(Quantr)
plot(sum_rq, ols=FALSE)
qr_temp_asso_functional_distance<-broom::tidy(Quantr,se.type = "boot")
qr_temp_asso_functional_distance<-droplevels(subset(qr_temp_asso_functional_distance, term=="functional_distance"))

# plot quantile regression

to_plot_functional_distance<-data.frame(functional_distance=data_asso_dist2$functional_distance[], spatial_asso=scale(data_asso_dist2$spatial_asso2, center=F), temp_asso=scale(data_asso_dist2$temp_asso2, center=F))
to_plot_functional_distance<-melt(to_plot_functional_distance, id.vars="functional_distance")
ggplot(droplevels(to_plot_functional_distance), aes(x=functional_distance, y=value, group=variable))+
  geom_point(aes(col=variable),size=1, alpha=0.5) +
  geom_quantile(quantiles=c(0.01,0.99), formula=y ~ x, aes(col=variable)) +
  scale_color_viridis(alpha = 1, direction = 1,discrete = TRUE, option = "D", label=c("Spatial associations","Temporal associations"))+
  labs(y="Associations", x = "Functional distance")+
  geom_smooth(method='lm', formula= y ~ 0, colour="blue", se=FALSE, linetype="11") +
  theme_modern()+theme(legend.position = c(0.8,0.9), legend.title = element_blank())

to_plot_functional_distance2a<-data.frame(tau=as.data.frame(qr_spatial_asso_functional_distance)$tau,
                                  value=scale(as.data.frame(qr_spatial_asso_functional_distance)$estimate, center=F),
                                  sd_value=max(abs(scale(as.data.frame(qr_spatial_asso_functional_distance)$estimate, center=F)))/max(abs(as.data.frame(qr_spatial_asso_functional_distance)$estimate))*as.data.frame(qr_spatial_asso_functional_distance)$std.error,
                                  gr="spatial_asso")
to_plot_functional_distance2b<-data.frame(tau=as.data.frame(qr_temp_asso_functional_distance)$tau,
                                  value=scale(as.data.frame(qr_temp_asso_functional_distance)$estimate, center=F),
                                  sd_value=max(abs(scale(as.data.frame(qr_temp_asso_functional_distance)$estimate, center=F)))/max(abs(as.data.frame(qr_temp_asso_functional_distance)$estimate))*as.data.frame(qr_temp_asso_functional_distance)$std.error,
                                  gr="int")
to_plot_functional_distance2<-rbind(to_plot_functional_distance2a,to_plot_functional_distance2b)
ggplot(droplevels(to_plot_functional_distance2), aes(x=tau, y=value, group=gr))+
  geom_point(aes(col=gr),size=1, alpha=0.5) +
  geom_line(aes(col=gr),size = 1)+ 
  scale_color_viridis(alpha = 1, direction = 1,discrete = TRUE, option = "D", label=c("Spatial associations","Temporal associations"))+
  geom_smooth(method='lm', formula= y ~ 0, colour="blue", se=FALSE, linetype="11")+
  labs(y="Quantile regression")+
  geom_ribbon(aes(ymin=value-sd_value*1.96,ymax=value+sd_value*1.96, fill=gr),alpha=0.25)+
  scale_fill_viridis(alpha = 1, direction = 1,discrete = TRUE, option = "D", label=c("Spatial associations","Temporal associations"))+
  theme_modern()+theme(legend.position = c(0.8,0.9), legend.title = element_blank())
```

## For associations significant either in space or in time

### Functional distance (and other niche distances)

```{r}
# select only part of the dataset to allow broom::tidy to work
set.seed(1)
num_zero4<-sample(num_zero2,5000)
data_asso_functional_distance<-data_asso_dist[c(num_not_zero2,num_zero4),]

set.seed(1)
num_not_zero4b<-sample(num_not_zero2b,5000)
data_asso_dist4<-data_asso_dist[c(num_not_zero4b),]

# quantile regression for spatial_asso2 and functional_distance

Quantr<-rq(data=droplevels(data_asso_dist4),
           tau= seq(0.01,0.99, by=0.01),
           formula = spatial_asso2_save ~  functional_distance) # or phylogenetic_distance, specialisation_distance, niche_overlap_distance
sum_rq<-summary(Quantr)
plot(sum_rq, ols=FALSE)
qr_spatial_asso_functional_distance_tot<-broom::tidy(Quantr,se.type = "boot")
qr_spatial_asso_functional_distance_tot<-droplevels(subset(qr_spatial_asso_functional_distance_tot, term=="functional_distance"))

# quantile regression for temp_asso2 and functional_distance
Quantr<-rq(data=droplevels(data_asso_functional_distance),
           tau= seq(0.01,0.99, by=0.01),
           formula = temp_asso2_save ~  functional_distance) # or phylogenetic_distance, specialisation_distance, niche_overlap_distance
sum_rq<-summary(Quantr)
plot(sum_rq, ols=FALSE)
qr_temp_asso_functional_distance_tot<-broom::tidy(Quantr,se.type = "boot")
qr_temp_asso_functional_distance_tot<-droplevels(subset(qr_temp_asso_functional_distance_tot, term=="functional_distance"))

to_plot_functional_distance_tota<-data.frame(functional_distance=data_asso_dist4$functional_distance[], value=scale(data_asso_dist4$spatial_asso2_save, center=F), variable="ses")
to_plot_functional_distance_totb<-data.frame(functional_distance=data_asso_functional_distance$functional_distance[], value=scale(data_asso_functional_distance$temp_asso2_save, center=F), variable="interaction")
to_plot_functional_distance_tot<-rbind(to_plot_functional_distance_tota,to_plot_functional_distance_totb)

ggplot(droplevels(to_plot_functional_distance_tot), aes(x=functional_distance, y=value, group=variable))+
  geom_point(aes(col=variable),size=1, alpha=0.5) +
  geom_quantile(quantiles=c(0.01,0.99), formula=y ~ x, aes(col=variable)) +
  scale_color_viridis(alpha = 1, direction = 1,discrete = TRUE, option = "D", label=c("Spatial associations","Temporal associations"))+
  labs(y="Associations", x = "Functional distance")+
  geom_smooth(method='lm', formula= y ~ 0, colour="blue", se=FALSE, linetype="11") +
  theme_modern()+theme(legend.position = c(0.8,0.9), legend.title = element_blank())

to_plot_functional_distance_tot2a<-data.frame(tau=as.data.frame(qr_spatial_asso_functional_distance_tot)$tau,
                                      value=scale(as.data.frame(qr_spatial_asso_functional_distance_tot)$estimate, center=F),
                                      sd_value=max(abs(scale(as.data.frame(qr_spatial_asso_functional_distance_tot)$estimate, center=F)))/max(abs(as.data.frame(qr_spatial_asso_functional_distance_tot)$estimate))*as.data.frame(qr_spatial_asso_functional_distance_tot)$std.error,
                                      gr="ses")
to_plot_functional_distance_tot2b<-data.frame(tau=as.data.frame(qr_temp_asso_functional_distance_tot)$tau,
                                      value=scale(as.data.frame(qr_temp_asso_functional_distance_tot)$estimate, center=F),
                                      sd_value=max(abs(scale(as.data.frame(qr_temp_asso_functional_distance_tot)$estimate, center=F)))/max(abs(as.data.frame(qr_temp_asso_functional_distance_tot)$estimate))*as.data.frame(qr_temp_asso_functional_distance_tot)$std.error,
                                      gr="int")
to_plot_functional_distance_tot2<-rbind(to_plot_functional_distance_tot2a,to_plot_functional_distance_tot2b)
ggplot(droplevels(to_plot_functional_distance_tot2), aes(x=tau, y=value, group=gr))+
  geom_point(aes(col=gr),size=1, alpha=0.5) +
  geom_line(aes(col=gr),size = 1)+ 
  scale_color_viridis(alpha = 1, direction = 1,discrete = TRUE, option = "D", label=c("Spatial associations","Temporal associations"))+
  geom_smooth(method='lm', formula= y ~ 0, colour="blue", se=FALSE, linetype="11")+
  labs(y="Quantile regression")+
  geom_ribbon(aes(ymin=value-sd_value*1.96,ymax=value+sd_value*1.96, fill=gr),alpha=0.25)+
  scale_fill_viridis(alpha = 1, direction = 1,discrete = TRUE, option = "D", label=c("Spatial associations","Temporal associations"))+
  theme_modern()+theme(legend.position = c(0.8,0.9), legend.title = element_blank())
```









Simulations from multiccm to build a simulated dataset with make_ccm_data()

```{r}

# data with interaction

n_int  <-  5 # number of species pairs with interaction

data_rand_int_tot  <- data.frame(code_square=NA, year=NA, variable=NA, value=NA)

for(i in 1:n_int){
ccm_data_int  <- make_ccm_data(seednum=i, obs_sd=0.025, Sstr=2)
to_zero <- which(ccm_data_int$Accm<quantile(na.omit(ccm_data_int$Accm),0.3))
ccm_data_int$Accm <- (ccm_data_int$Accm)^10
ccm_data_int$Accm <- round(ccm_data_int$Accm*9/max(na.omit(ccm_data_int$Accm)))+1
ccm_data_int$Accm[to_zero] <- 0
to_zero <- which(ccm_data_int$Bccm<quantile(na.omit(ccm_data_int$Bccm)),0.3)
ccm_data_int$Bccm <- (ccm_data_int$Bccm)^10
ccm_data_int$Bccm <- round(ccm_data_int$Bccm*9/max(na.omit(ccm_data_int$Bccm)))+1
ccm_data_int$Bccm[to_zero] <- 0

data_rand_int <- data.frame(code_square=c(rep("A",10),rep("B",10),rep("C",10),rep("D",10),rep("E",10),
                                  rep("F",10),rep("G",10),rep("H",10),rep("I",10),rep("J",10),
                                  rep("K",10),rep("L",10),rep("M",10),rep("N",10),rep("O",10),
                                  rep("P",10),rep("Q",10),rep("R",10),rep("S",10),rep("T",10)),
                   year=na.omit(ccm_data_int$time_ccm),
                   spa=na.omit(ccm_data_int$Accm), #na.omit(abs(round(10*abs(log(ccm_data_int$Accm)))-2)),
                   spb=na.omit(ccm_data_int$Bccm)) #na.omit(abs(round(10*abs(log(ccm_data_int$Bccm)))-2)))
                   
names(data_rand_int)[3:4] <- c(paste0("sp",i), paste0("sp",(i+100)))

data_rand_int <- melt(data_rand_int, measure.vars= c(paste0("sp",i), paste0("sp",(i+100))))

data_rand_int_tot <- rbind(data_rand_int_tot, data_rand_int)
}

data_rand_int_tot <- data_rand_int_tot[-1,]

# data without interaction

n_no_int  <-  15 # number of species pairs without interaction

data_rand_no_int_tot  <- data.frame(code_square=NA, year=NA, variable=NA, value=NA)

for(j in 1:n_no_int){
ccm_data_no_int  <- make_ccm_data(Sstr = 0, seednum = (1000+j), obs_sd=0.25)

to_zero <- which(ccm_data_no_int$Accm<quantile(na.omit(ccm_data_no_int$Accm),0.3))
ccm_data_no_int$Accm <- (ccm_data_no_int$Accm)^10
ccm_data_no_int$Accm <- round(ccm_data_no_int$Accm*9/max(na.omit(ccm_data_no_int$Accm)))+1
ccm_data_no_int$Accm[to_zero] <- 0
to_zero <- which(ccm_data_no_int$Bccm<quantile(na.omit(ccm_data_no_int$Bccm),0.3))
ccm_data_no_int$Bccm <- (ccm_data_no_int$Bccm)^10
ccm_data_no_int$Bccm <- round(ccm_data_no_int$Bccm*9/max(na.omit(ccm_data_no_int$Bccm)))+1
ccm_data_no_int$Bccm[to_zero] <- 0



data_rand_no_int <- data.frame(code_square=c(rep("A",10),rep("B",10),rep("C",10),rep("D",10),rep("E",10),
                                  rep("F",10),rep("G",10),rep("H",10),rep("I",10),rep("J",10),
                                  rep("K",10),rep("L",10),rep("M",10),rep("N",10),rep("O",10),
                                  rep("P",10),rep("Q",10),rep("R",10),rep("S",10),rep("T",10)),
                   year=na.omit(ccm_data_no_int$time_ccm),
                   spa=na.omit(ccm_data_no_int$Accm), #na.omit(abs(round(10*abs(log(ccm_data_no_int$Accm)))-2)),
                   spb=na.omit(ccm_data_no_int$Bccm)) #na.omit(abs(round(10*abs(log(ccm_data_no_int$Bccm)))-2)))
                   
names(data_rand_no_int)[3:4] <- c(paste0("sp",(j+1000)), paste0("sp",(j+1100)))

data_rand_no_int <- melt(data_rand_no_int, measure.vars= c(paste0("sp",(j+1000)), paste0("sp",(j+1100))))

data_rand_no_int_tot <- rbind(data_rand_no_int_tot, data_rand_no_int)                   
}

data_rand_no_int_tot <- data_rand_no_int_tot[-1,]

# Merge all data with and without interaction

data_rand_tot <- rbind(data_rand_no_int_tot, data_rand_int_tot)

data_rand_tot$sp_sq <- paste0(data_rand_tot$variable,sep="_",data_rand_tot$code_square)

i_ind <- 0
for(i in unique(data_rand_tot$sp_sq)){
i_ind <- i_ind+1
set.seed(i_ind)
piece <- rbinom(1,1,0.5)
if(piece==1){
data_rand_tot$value[which(data_rand_tot$sp_sq==i)] <- 0
}
}

data_rand_tot$code_square <- as.factor(data_rand_tot$code_square)

names(data_rand_tot)[c(3:4)] <- c("code_sp","abond")

data_rand_tot_sp <- data_rand_tot

names(data_rand_tot_sp)[1] <- "code_point"


# Run spatio-temporal analysis

# Prepare data

dataperpoint_rand <- dcast(data_rand_tot_sp, code_point + year ~ code_sp, fun.aggregate = sum,value.var="abond")

# Run the function

spa_tem_associations_rand <- ddply(dataperpoint_rand, .(year), .fun=main_fun_space, 
                           N.null=1000, occmin=1, 
                           .progress = "text")
                           
# Keep all spatial associations                           

spa_tem_associations_rand$ses_init <- spa_tem_associations_rand$spatial_asso

# Set non-significant spatial associations to 0

spa_tem_associations_rand$spatial_asso[spa_tem_associations_rand$pval>0.05] <- NA



# Run spatial analysis

# Prepare data

dataperpoint_rand <- dcast(data_rand_tot_sp, code_point ~ code_sp, fun.aggregate = sum,value.var="abond")

# Run the function

spatial_associations_rand <- main_fun_space(dataperpoint_rand,N.null=1000, occmin=1)
                           
# Keep all spatial associations                           

spatial_associations_rand$ses_init <- spatial_associations_rand$spatial_asso

# Set non-significant spatial associations to 0

spatial_associations_rand$spatial_asso[spatial_associations_rand$pval>0.05] <- NA



# Run temporal analysis

temporal_associations_rand <- multisp_CCM(data_rand_tot)


# E dim

Edim_temporal_associations_rand <- Edim_multisp_CCM(data_rand_tot)

names(Edim_temporal_associations_rand)[2] <- "E"
Edim_temporal_associations_rand$species <- substr(Edim_temporal_associations_rand$species,3,21)

temporal_associations_Edim_rand <- merge(temporal_associations_rand, Edim_temporal_associations_rand, by=c("species"), all.x=T)

temporal_associations_signif_rand <- droplevels(na.omit(temporal_associations_Edim_rand[temporal_associations_Edim_rand$rho<0.05,]))
temporal_associations_signif_rand$spA <- sub(".*_cause_","", temporal_associations_signif_rand$species)
temporal_associations_signif_rand$spB <- sub("_cause_.*","", temporal_associations_signif_rand$species)


# Quantifying the temporal association between species

temporal_associations_quanti_rand <- ddply(temporal_associations_signif_rand, .(species), .fun=smap_fun_signif_rand, data_ccm=data_rand_tot, .progress="text")

#ici

temporal_associations_Edim <- readRDS("output/temporal_associations_Edim.rds")
temporal_associations_quanti <- readRDS("output/temporal_associations_quanti.rds")

temporal_associations_quanti_sum <- data.frame(temporal_associations_quanti %>% group_by(zonebio,habit,species) %>% summarize(temp_int=mean(spB_smap,na.rm=T), sd_temp_int=sd(spB_smap, na.rm=T), temp_int2=mean(spB_smap2,na.rm=T), sd_temp_int2=sd(spB_smap2, na.rm=T)))

temporal_associations_signif <- merge(temporal_associations_signif,temporal_associations_quanti_sum, by=c("zonebio","habit","species"))

```

### Merge spatial and temporal associations obtained for simulation

```{r}

sp_tp_association <- merge(temporal_associations_signif, spatial_associations, by.x=c("habit","zonebio","spA","spB"), by.y=c("habit","zonebio","spi","spj"), all.x=T)


```



## References

Clark, A. T., Ye, H., Isbell, F., Deyle, E. R., Cowles, J., Tilman, G. D., & Sugihara, G. (2015). Spatial convergent cross mapping to detect causal relationships from short time series. Ecology, 96(5), 1174-1181. https://doi.org/10.1890/14-1479.1

Morueta‐Holme, N., Blonder, B., Sandel, B., McGill, B. J., Peet, R. K., Ott, J. E., ... & Svenning, J. C. (2016). A network approach for inferring species associations from co‐occurrence data. Ecography, 39(12), 1139-1150. https://doi.org/10.1111/ecog.01892

Rigal, S., Devictor, V., Gaüzère, P., Kéfi, S., Forsman, J. T., Kajanus, M. H., ... & Dakos, V. (2022). Biotic homogenisation in bird communities leads to large‐scale changes in species associations. Oikos, 2022(3), e08756. https://doi.org/10.1111/oik.08756