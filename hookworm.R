
#' The data is available in: https://cema.shinyapps.io/kenya_ntd/
#' 
# Importing packages ------------------------------------------------------

pacman::p_load(
  dplyr,
  raster,
  terra,
  data.table,
  sf,
  purrr,
  stringi,
  stringr,
  ggplot2,
  caret,
  recipes,
  ggcorrplot,
  patchwork,
  INLA,
  fmesher,
  inlabru,
  ggspatial,
  lme4,
  RColorBrewer,
  pscl
)
source('functions.R')

# Importing data ----------------------------------------------------------

covs <- rast("tif/tif_5KM.tif")
df <- fread('data/hook_schPrevalence.csv')
df2 <- fread("data/schoolPrevalence_coast_west.csv") 
kenya_shapefile <- st_read('shapefiles/gadm41_KEN_0.shp', quiet = T)
county_shapefile1 <- st_read('shapefiles/gadm36_ken_1.shp', quiet = T)
subcounty_shapefile1 <- st_read('shapefiles/gadm36_ken_2.shp', quiet = T)
actualHook <- fread('../../kenya_ntd/data/hookworm_summary.csv')
actualHook_coast_western <- fread('data/subPrevalence_coast_west.csv')
cov <- fread("data/cov.csv")
all_df <- fread("data/all_df.csv")
pred1 <- fread("data/pred.csv")
pp <- fread("data/pp.csv")
covs <- rast("tif/tif_5KM.tif")
dm <- fread("data/dm_cat.csv")
pred <- fread("data/pred_cat.csv")

# Some initial cleaning ----------------------------------------------------

df <- rbind(df, df2)
county_shapefile <- data.table(county_shapefile1) |> 
  _[, .(county = NAME_1, geometry)] |> 
  st_as_sf()
subcounty_shapefile <- data.table(subcounty_shapefile1) |> 
  _[, .(county = NAME_1, subcounty = NAME_2, geometry)] |>
  _[,
    subcounty := case_when(subcounty == "Webute West" ~ "Webuye West",
                           
                           TRUE ~ subcounty)] |>
  st_as_sf()

actualHook_shp <- actualHook |>
  _[ward != "", .(prev = mean(prev)), by = .(county, subcounty)] |>
  rbind(actualHook |>
          _[ward == "", .(county, subcounty, prev)]) |>
  data.table() |>
  rbind(actualHook_coast_western) |> 
  _[, `:=`(
    county = case_when(
      county == 'Murang A' ~ "Murang'a",
      county == 'Nairobi City' ~ 'Nairobi',
      county == 'Elgeyo/Marakwet' ~ 'Elgeyo-Marakwet',
      TRUE ~ county
    ),
    subcounty = case_when(
      subcounty == 'Igambango Ombe' ~ "Chuka/Igambang'Ombe",
      subcounty == 'Mwingi Nort' ~ "Mwingi North",
      subcounty == 'Igembecentral' ~ "Igembe Central",
      subcounty == 'Imenti North' ~ "North Imenti",
      subcounty == 'Lunga Lunga' ~ "Lungalunga",
      subcounty == 'Lunga Lunga' ~ "Lungalunga",
      subcounty == "Rachuonyo South" ~ "Kasipul",
      subcounty == "Cheptais" ~ "Lugari",
      str_detect(subcounty, "^Suba ") ~ "Suba",
      str_detect(subcounty, "^Gem ") ~ "Gem",
      str_detect(subcounty, "^Rachuonyo ") ~ "Kasipul",
      TRUE ~ subcounty
    ),
    prv_ctg =  case_when(
      prev < 2 ~ "<2%",
      prev >= 2 & prev < 10 ~ "2%-<10%",
      prev >= 10 & prev < 20 ~ "10%-<20%",
      prev >= 20 & prev < 50 ~ "20%-<50%",
      prev >= 50 ~ ">=50%",
      TRUE ~ NA_character_
    ) |> factor(levels = c("<2%", "2%-<10%", "10%-<20%", "20%-<50%", ">=50%"))
  )] |>
  merge(data.table(subcounty_shapefile), by = c("county", "subcounty"))

pop_shp <- actualHook_shp[, .(county, subcounty, geometry)]

coords.schs <- all_df[, .(lng, lat)] |>
  st_as_sf(coords = c('lng', 'lat'), crs = crs(kenya_shapefile)) |> 
  st_transform(crs = epsgKM(32737)) |> 
  st_coordinates() |> 
  data.table() |> 
  setNames(c('km_x', 'km_y'))

coords.shp <- all_df[, .(lng, lat)] |>
  st_as_sf(coords = c('lng', 'lat'),
           crs = crs(kenya_shapefile)) |> 
  st_transform(epsgKM(32737))
st_write(coords.shp, "shapefiles/coords.shp.shp", quiet = T)

# School level prevalence categories
kk <- all_df |> 
  mutate(
    prev = prev*1e2,
    prv_ctg =  case_when(
      prev < 2 ~ "<2%",
      prev >= 2 & prev < 10 ~ "2%-<10%",
      prev >= 10 & prev < 20 ~ "10%-<20%",
      prev >= 20 & prev < 50 ~ "20%-<50%",
      prev >= 50 ~ ">=50%",
      TRUE ~ NA_character_
    ) |> factor(levels = c("<2%", "2%-<10%", "10%-<20%", "20%-<50%", ">=50%"))
  )

kk[, .N, by = prv_ctg] |> 
  _[, perc := round(N/ sum(N) *1e2, 2)][]


# Plotting for categories -------------------------------------------------

dm_cat <- dm |> 
  _[, prev := n/N] |> 
  _[, !c("n", "N")] |>
  _[, `:=`(
    soil_temp_cat = factor(dm$soil_temp_cat, 
                           levels = c("[min,19.3]", "(19.3,20.9]", "(20.9,22.6]", "(22.6,max]"),
                           ordered = TRUE)
  )] |> 
  dplyr::select(contains("cat"), prev) |> 
  melt(id.vars = "prev") %>%
  split(.$variable) %>%
  setNames(
    c(
      "Elevation",
      "Improved water source",
      "Day land surface temperature",
      "NDVI",
      "Open defecation",
      'Precipitation',
      'Relative humidity',
      "Soil moisture",
      "Soil pH",
      "Soil Temperature"
    )
  ) %>%
  map2(names(.), \(x, y){
    angle <- ifelse(
      y %in% c(
        "Elevation",
        "Day land surface temperature",
        "Night land surface temperature",
        'Precipitation',
        'Relative humidity',
        "Soil moisture"
      ),
      45,
      0
    )
    hjust <- ifelse(
      y %in% c(
        "Elevation",
        "Day land surface temperature",
        "Night land surface temperature",
        'Precipitation',
        'Relative humidity',
        "Soil moisture"
      ),
      1,
      .5
    )
    x |> 
      _[, !("variable")] |> 
      _[, .(prev = mean(prev)), by = value] |> 
      ggplot(aes(x = value, y = prev)) +
      geom_col(aes(y = prev), fill = "#6baed6", width = .7) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(
          color = 'black',
          angle = angle,
          hjust = hjust
        ),
        axis.text.y = element_text(color = 'black', ),
        axis.ticks.x = element_line(color = "black"),
        axis.title = element_text(color = 'black'),
        plot.title = element_text(color = 'black', hjust = .5)
      ) +
      xlab("") +
      ylab('') +
      ggtitle(y)
  })
  

# Small arrangement of datas ----------------------------------------------

dm2 <- dm |>
  _[, .(
    elevation_cat,
    rel_humidity,
    n,
    N
  )]


# Filtering the shape file and the prediction data  ------------------------

df_c <- df |>
  _[, county := case_when(
    county == "Murang A" ~ "Murang'a",
    county == "Elgeyo/Marakwet" ~ "Elgeyo-Marakwet",
    county == "Nairobi City" ~ "Nairobi",
    county == "Trans  Nzoia" ~ "Trans Nzoia",
    TRUE ~ county
  )]
setdiff(df$county, county_shapefile1$NAME_1)

cc <- county_shapefile1 |>
  dplyr::select(county = NAME_1, geometry) |>
  filter(county %in% c(unique(df_c$county), "Nyeri")) |> 
  st_transform(epsgKM(32737))

cc2 <- county_shapefile1 |>
  dplyr::select(county = NAME_1, geometry) |>
  filter(!county %in% c(unique(df_c$county), "Nyeri")) |> 
  st_transform(epsgKM(32737))

kenya_shapefile2 <- kenya_shapefile |> 
  st_transform(epsgKM(32737))

bb <- cc |>
  st_union() |>   
  st_cast("POLYGON") 

bb2 <- cc2 |>
  st_union() |>   
  st_cast("POLYGON") 
saveRDS(bb, "data/bb.rds")

plot(bb)
points(coords.shp)

covs2 <- brick(covs) %>%
  projectRaster(crs = epsgKM(32737))

grid <- covs2[[1]]
pred.coords <- coordinates(grid)

# Filtering prediction data
pred.dat <- cbind(pred.coords, getValues(covs2))
ind <- apply(pred.dat, 1, \(x) any(is.na(x)))
miss    <- which(ind == TRUE)
nonmiss <- which(ind == FALSE)
coord.p <- pred.dat[nonmiss, 1:2]
pop <- pred.dat[nonmiss,"u15pop"]

# Glm model for the residuals and obtaining the variogram -----------------

dm_logit <- data.table(dm2) |> 
  _[, logit_prev := log((n + .5)/(N - n + .5))]

glm.fit <- glm(
  cbind(n, N - n) ~ .,
  data = dm_logit,
  family = 'binomial'
)

# Obtaining the residuals
resids <- rstandard(glm.fit)
ggvario(
  coords = coords.schs,
  data =  resids,
  bins = 17,
  xlab = "distance (Km)",
  show_nbins = T,
  maxdist = 200
) +
  ggtitle('Empirical variogram for hookworm') +
  theme(
    plot.title = element_text(color = 'black', hjust = .5),
    axis.title = element_text(color = 'black')
  )

ggsave(
  'plots/hookworm_variogram.png',
  bg = NULL,
  dpi = 1e3,
  height = 6,
  width = 10
)

# Creating the mesh and the SPDE model -------------------------------------------------------

mesh <- fm_mesh_2d_inla(
  boundary = kenya_shapefile2,
  loc = coords.schs,
  max.edge = c(20, 40),
  cutoff = 10
)
plot(mesh)

# Matern SPDE model object using inla.pcmatern
r0 <- as.vector(0.05 * (st_bbox(bb)["ymax"] - st_bbox(bb)["ymin"]))
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  alpha = 2,
  prior.range = c(r0, .01),   # P(r0 < .01) = 0.01
  prior.sigma = c(5, .01)   # P(sigma > 5) = 0.01
)

# index set for the SPDE model
indexs <- inla.spde.make.index("s", spde$n.spde)

# Projection matrix
A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(coords.schs))

# Matrix for prediction
Ap <- inla.spde.make.A(mesh = mesh, loc = coord.p)
dim(Ap)

# Creating the stacks --------------------------------------------------------------------

stk.e1 <- inla.stack(
  tag = 'point',
  data = list(y = dm$n, numtrials = dm$N),
  A = list(A, 1, 1),
  effects = list(
    s = indexs,
    rr = 1:length(dm$elevation_cat), 
    data.frame(
      b0 = 1,
      x1 = dm2$elevation_cat,
      x2 = dm2$rel_humidity
    )
  )
)

stk.p1 <- inla.stack(
  tag = "pred",
  data = list(y = NA, numtrials = NA),
  A = list(Ap,1,1),
  effects = list(
    s = indexs,
    rr = (length(pred$elevation_cat) + 1):(length(pred$elevation_cat) + nrow(pred)),
    data.frame(
      b0 = 1,
      x1 = pred$elevation_cat,
      x2 = pred$rel_humidity
    )
  )
)

# Specifying the priors and the formula --------------------------------------------------

hyper.prec = list(theta = list(prior = "pc.prec", param = c(5,0.01)))
formula <- y ~ 0 + b0 + x1 + x2 +
  f(s, model = spde) +
  f(rr, model = 'iid', hyper = hyper.prec)

# The INLA model ----------------------------------------------------------
# stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.e1, stk.p1)

#' The inla() call
mod <- inla(
  formula,
  family = "binomial",
  Ntrials = numtrials,
  data = inla.stack.data(stk.full),
  control.family = control.family(link = 'logit'),
  control.predictor = list(
    compute = TRUE,
    link = 1,
    A = inla.stack.A(stk.full)
  ),
  control.compute = list(config = TRUE, return.marginals.predictor = TRUE, cpo = T),
  control.inla = list(tolerance = 1e-7), 
  num.threads	= 7,
  verbose = F
)
summary(mod)
spde.result <- inla.spde2.result(inla = mod,name = "s",spde = spde)

#Parameters
coeff.reg <- tidy.inla(mod, exp = F) |>
  mutate(
    terms = str_replace_all(terms, "x2", "Open defecation ") |> 
      str_replace_all("x3", "Precipitation ") |> 
      str_replace_all("x1", "Elevation")
  ) |> 
  filter(terms != "b0")
fwrite(coeff.reg, "data/coeff.reg.csv", row.names = F)

# Predictions and rasterizing ------------------------------------------------------------

#' # Mapping the vaccination coverage
index <- inla.stack.index(stack = stk.full, tag = "pred")$data

#> We create vectors with the mean prevalence and lower and upper limits of
#> 95% credible intervals with the values of the columns "mean", "0.025quant"
#> and "0.975quant" and the rows given by index.

prev_mean_init <- mod$summary.fitted.values[index, "mean"]
summary(prev_mean_init)
prev_ll <- mod$summary.fitted.values[index, "0.025quant"]
prev_ul <- mod$summary.fitted.values[index, "0.975quant"]

#' Rasterizing the results
r_prev_mean <- rasterize(x = coord.p,
                         y = grid,
                         field = prev_mean_init)
#' The plot
plot(r_prev_mean)

# Sampling from the posterior distribution ------------------------------------------------

nsamp <- 1000

#Posterior sampling
ps <- inla.posterior.sample(nsamp, mod, num.threads = 7, parallel.configs = T)
contents <- mod$misc$configs$contents

linpred <- inla.posterior.sample.eval("APredictor", ps)[index,]
inv.linpred <- gtools::inv.logit(linpred)
pred.grid <- data.frame(t(apply(inv.linpred, 1, FUN = function(x){ c(mean(x), sd(x), quantile(x, probs = c(0.025,0.5,0.975)))})))
colnames(pred.grid) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
fitted.pred.mean <- pred.grid[["mean"]]
fitted.pred.sd     <- pred.grid[["sd"]]
fitted.pred.median <- pred.grid[["0.5quant"]]
fitted.pred.low    <- pred.grid[["0.025quant"]]
fitted.pred.up     <- pred.grid[["0.975quant"]]

#' Rasterizing the results
hook1.mean <- rasterize(x = coord.p,
                        y = grid,
                        field = fitted.pred.mean)
#crs(hook1.mean) <- epsgKM(32737)
plot(hook1.mean)
#' Rasterizing the low CI
hook1.lower <- rasterize(x = coord.p,
                         y = grid,
                         field = fitted.pred.low)

#' Rasterizing the upper ci
hook1.upper <- rasterize(x = coord.p,
                         y = grid,
                         field = fitted.pred.up)
#crs(hook1.upper) <- epsgKM(32737)

#' The plot
rasters <- list(hook1.lower, hook1.mean, hook1.upper) |>
  setNames(c(
    'Lower CI',
    'Hookworm probability',
    'Upper CI'
  ))
saveRDS(rasters, "../raster results/hookworm_results.RDS")
par(mfrow = c(1,3))
map(rasters, plot)
par(mfrow = c(1,1))

# Obtaining the key statistic's -----------------------------------------------

#Parameters
coeff.reg <- summary(mod)$fixed[,1:5]

#range for spatial RE
range_std <- summary(mod)$hyperpar[1:2,][,1:5]
range_std[2, c(1,3,4,5)] <- (range_std[2, c(1,3,4,5)]) ^ 2

#variance for IID RE
var.ind <- inla.tmarginal(function(x)
  1 / x, mod$marginals.hyperpar[[3]])
var.iid <- inla.zmarginal(var.ind, silent = TRUE)
variance.iid <- c(var.iid$mean,
                  var.iid$sd,
                  var.iid$quant0.025,
                  var.iid$quant0.5,
                  var.iid$quant0.975)
param.all <- rbind(range_std, variance.iid)
param.all <- param.all[, c(1,3,5)]
rownames(param.all) <- c("Spatial range", "Spatial variance", "Nugget")
write.csv(param.all, "data/parameter_output.csv")

# Creating raster considering NA areas ------------------------------------

ll <- 1:length(ind);
ll[nonmiss] <- fitted.pred.mean
ll[miss] <- NA
predicted.df <- cbind(data.table(pred.coords), preds = ll)
hook <- rast(predicted.df)
plot(raster(hook))

par(mfrow = c(1,3))
map(rasters, plot)
par(mfrow = c(1,1))
writeRaster(rasters$`Hookworm probability`, 'rasters/hook_prev.tif', overwrite = T)

# Plotting before categorizing --------------------------------------------

rast_ll <- ggplot() +
  geom_tile(data = as.data.frame(rast(hook1.lower) * 1e2, xy = TRUE, na.rm = TRUE),
            aes(x = x, y = y, fill = layer)) +
  geom_sf(
    data = st_transform(county_shapefile, crs = epsgKM(32737)),
    fill = NA
  ) +
  geom_sf(
    data = st_as_sf(bb2) |>
      st_transform(crs = epsgKM(32737)),
    fill = NA,
    color = "black"
  ) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, 
                       labels = scales::label_number(accuracy = 0.0001)) +
  theme_void() +
  theme(
    plot.title = element_text(colour = "black", face = "bold", hjust = .5, size = 17),
    legend.text = element_text(color = "black", size = 14)
  ) +
  labs(fill = "", title = "Lower CI(%)") 


rast_mean <- ggplot() +
  geom_tile(data = as.data.frame(rast(hook1.mean) * 1e2, xy = TRUE, na.rm = TRUE),
            aes(x = x, y = y, fill = layer)) +
  geom_sf(
    data = st_transform(county_shapefile, crs = epsgKM(32737)),
    fill = NA
  ) +
  geom_sf(
    data = st_as_sf(bb2) |>
      st_transform(crs = epsgKM(32737)),
    fill = NA,
    color = "black"
  ) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  theme_void() +
  theme(
    plot.title = element_text(colour = "black", face = "bold", hjust = .5, size = 17),
    legend.text = element_text(color = "black", size = 14)
  ) +
  labs(fill = "", title = "Mean prevalence(%)")

rast_upper <- ggplot() +
  geom_tile(data = as.data.frame(rast(hook1.upper) * 1e2, xy = TRUE, na.rm = TRUE),
            aes(x = x, y = y, fill = layer)) +
  geom_sf(
    data = st_transform(county_shapefile, crs = epsgKM(32737)),
    fill = NA
  ) +
  geom_sf(
    data = st_as_sf(bb2) |>
      st_transform(crs = epsgKM(32737)),
    fill = NA,
    color = "black"
  ) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  theme_void() +
  theme(
    plot.title = element_text(colour = "black", face = "bold", hjust = .5, size = 17),
    legend.text = element_text(color = "black", size = 14)
  ) +
  labs(fill = "", title = "Upper CI(%)")

rr_all <- rast_ll + rast_mean + rast_upper
ggsave(
  plot = rr_all,
  'plots/rr_all.png',
  dpi = 5e2,
  height = 12,
  width = 20,
  bg = NULL
)

# Categorizing rasters and plotting ---------------------------------------

raster_data <- rast(hook1.mean) * 1e2  
raster_df <- as.data.frame(raster_data, xy = TRUE, na.rm = TRUE)
names(raster_df)[3] <- "value"

breaks <- c(0, 2, 10, 50, Inf)  # Define the breaks
labels <- c("<2%", "2%-<10%", "10%-<50%", ">=50%")  # Labels for the categories

raster_df$category <- cut(raster_df$value, breaks = breaks, labels = labels, include.lowest = TRUE) |> 
  as.factor()

# Step 5: Define colors for the categories
prv_ctg_colors <- c(
  "<2%" = brewer.pal(4, "YlOrRd")[1],
  "2%-<10%" = brewer.pal(4, "YlOrRd")[2],
  "10%-<50%" = brewer.pal(4, "YlOrRd")[3],
  ">=50%" = brewer.pal(4, "YlOrRd")[4]
)

gg_rast <- ggplot() +
  geom_tile(data = raster_df, aes(x = x, y = y, fill = category)) +
  geom_sf(data = st_transform(county_shapefile1, epsgKM(32737)), fill = NA, color = "black", show.legend = T) +
  scale_fill_manual(values = prv_ctg_colors, name = "Prevalence (%)", drop = F) +
  theme_void() +
  annotation_north_arrow(style = north_arrow_fancy_orienteering(), location = "tr") +
  annotation_scale() 

ggsave(
  plot = gg_rast,
  'plots/gg_rast.png',
  dpi = 5e2,
  height = 6,
  width = 10,
  bg = NULL
)

ras_plts <- list(
  hook1.lower,
  hook1.mean,
  hook1.upper
) |> 
  setNames(c(
    'Lower CI',
    'Mean Prevalence',
    'Upper CI'
  )) %>%
  map2(
    names(.), \(x,y) {
      raster_data <- rast(x) * 1e2  
      raster_df <- as.data.frame(raster_data, xy = TRUE, na.rm = TRUE)
      names(raster_df)[3] <- "value"
      
      breaks <- c(0, 2, 10, 50, Inf)  # Define the breaks
      labels <- c("<2%", "2%-<10%", "10%-<50%", ">=50%")  # Labels for the categories
      
      raster_df$category <- cut(raster_df$value, breaks = breaks, labels = labels, include.lowest = TRUE)
      
      # Step 5: Define colors for the categories
      prv_ctg_colors <- c(
        "<2%" = brewer.pal(4, "YlOrRd")[1],
        "2%-<10%" = brewer.pal(4, "YlOrRd")[2],
        "10%-<50%" = brewer.pal(4, "YlOrRd")[3],
        ">=50%" = brewer.pal(4, "YlOrRd")[4]    
      )
      gg_rast <- ggplot() +
        geom_tile(data = raster_df, aes(x = x, y = y, fill = category), show.legend = T) +
        geom_sf(data = st_transform(county_shapefile1, epsgKM(32737)), fill = NA, color = "black", show.legend = T) +
        scale_fill_manual(values = prv_ctg_colors, name = "Prevalence (%)", drop = F) +
        theme_void() +
        theme(
          plot.title = element_text(colour = "black", hjust = .5)
        ) +
        annotation_north_arrow(style = north_arrow_fancy_orienteering(), location = "tr") +
        annotation_scale() +
        ggtitle(y)
    }
  ) |> 
  wrap_plots(guides = "collect")

ggsave(
  plot = ras_plts,
  'plots/ras_plts.png',
  dpi = 5e2,
  height = 12,
  width = 20,
  bg = NULL
)

# Population weighted -----------------------------------------------------

subcounty_shapefile2 <- subcounty_shapefile %>%
  mutate(ID = 1:nrow(.)) |> 
  st_as_sf() |> 
  st_transform(crs = epsgKM(32737))

spol <- as(as(subcounty_shapefile2, 'Spatial'), 'SpatialPolygons')
sp   <- rep(NA, nrow(pred))
crs(spol) <- NA

for (i in 1:length(spol)) {
  sp[as.vector(which(!is.na(over(
    SpatialPoints(coord.p), spol[i]
  ))))] <- i
}

# For each sub-county, it checks which prediction points fall within
# its boundary and assigns the subcounty ID to those points.

# Sub-County estimates
ss <- 1:nrow(subcounty_shapefile2)
ss.un <- unique(sp)
smiss <- which(!ss %in% ss.un) 

# The smiss vector contains the indices of subcounties
# that do not have any prediction points. 

subcounty_shapefile2[subcounty_shapefile2$ID %in% smiss,]

if (length(smiss) > 0)
  ss_num <- ss[-smiss]
if (length(smiss) == 0)
  ss_num <- ss

smat <- matrix(0, length(ss_num), 5)
subthreshNot0.02=subthresh0.02 = subthresh0.1 = 0
for (i in 1:length(ss_num)) {
  if (length(which(sp == ss_num[i])) == 1) {
    temp.pop <- pop[which(sp == ss_num[i])]
    ext <- as.vector(sapply(
      inv.linpred[which(sp == ss_num[i]), ],
      FUN = function(x)
        weighted.mean(x, w = temp.pop, na.rm = TRUE)
    ))
  }
  if (length(which(sp == ss_num[i])) > 1) {
    temp.pop <- pop[which(sp == ss_num[i])]
    ext <- as.vector(apply(
      data.frame(inv.linpred)[which(sp == ss_num[i]),],
      2,
      FUN = function(x)
        weighted.mean(x, w = temp.pop, na.rm = TRUE)
    ))
  }
  
  smat[i, ] <- as.vector(c(mean(ext), sd(ext), quantile(ext, probs = c(0.025, 0.5, 0.975))))
  subthresh0.02[i] <- length(which(ext >= 0.02)) / nsamp
  subthresh0.1[i] <- length(which(ext >= 0.1)) / nsamp
  subthreshNot0.02[i] <- 1 - (length(which(ext >= 0.02)) / nsamp)
}


# The code then calculates population-weighted estimates:
# If there is only one prediction point in the subcounty:
# sapply is used to calculate the weighted mean for each column (each sample) of inv.linpred.
# If there are multiple prediction points in the subcounty:
# apply is used to calculate the weighted mean for each column (each sample) of inv.linpred.

smat <- cbind(ss_num, smat, subthresh0.02, subthresh0.1, subthreshNot0.02)
colnames(smat) <- c("ID", "mean", "sd", "0.025quant", "0.5quant", "0.975quant", "exceed2", "exceed10", "notExceed2")
#fwrite(smat, 'data_insights/subcountyLevel.csv', row.names = F)

# Joining -----------------------------------------------------------------

subcountyLevel <- setDT(merge(subcounty_shapefile2, smat, by = 'ID', all.x = T)) |>
  _[, .(county, subcounty, 
        prev = as.numeric(mean * 1e2),
        lower_ci = as.numeric(`0.025quant` * 1e2),
        upper_ci = as.numeric(`0.975quant` * 1e2),
        exceed2, 
        exceed10,
        notExceed2,
        geometry)] |>
  _[,
    `:=`(
      prv_ctg =  case_when(
        prev < 2 ~ "<2%",
        prev >= 2 & prev < 10 ~ "2%-<10%",
        prev >= 10 & prev < 20 ~ "10%-<20%",
        prev >= 20 & prev < 50 ~ "20%-<50%",
        prev >= 50 ~ ">=50%",
        TRUE ~ NA_character_
      ) |> factor(levels = c("<2%", "2%-<10%", "10%-<20%", "20%-<50%", ">=50%")),
      
      exceed2_cat =  case_when(
        exceed2 < 0.25 ~ "0.0 - <0.25",
        exceed2 >= 0.25 & exceed2 <0.5 ~ "0.25 - <0.5",
        exceed2 >= 0.5 & exceed2 <0.75 ~ "0.5 - <0.75",
        exceed2 >= 0.75 & exceed2 < 0.9 ~ "0.75 - <0.9",
        exceed2 >= 0.9 ~ "0.9 - 1.0"
      ) |> factor(levels = c("0.0 - <0.25", "0.25 - <0.5", "0.5 - <0.75", "0.75 - <0.9", "0.9 - 1.0")),     
      exceed10_cat =  case_when(
        exceed10 < 0.25 ~ "0.0 - <0.25",
        exceed10 >= 0.25 & exceed10 <0.5 ~ "0.25 - <0.5",
        exceed10 >= 0.5 & exceed10 <0.75 ~ "0.5 - <0.75",
        exceed10 >= 0.75 & exceed10 < 0.9 ~ "0.75 - <0.9",
        exceed10 >= 0.9 ~ "0.9 - 1.0"
      ) |> factor(levels = c("0.0 - <0.25", "0.25 - <0.5", "0.5 - <0.75", "0.75 - <0.9", "0.9 - 1.0")),
      notExceed2_cat =  case_when(
        notExceed2 < 0.25 ~ "0.0 - <0.25",
        notExceed2 >= 0.25 & notExceed2 <0.5 ~ "0.25 - <0.5",
        notExceed2 >= 0.5 & notExceed2 <0.75 ~ "0.5 - <0.75",
        notExceed2 >= 0.75 & notExceed2 < 0.9 ~ "0.75 - <0.9",
        notExceed2 >= 0.9 ~ "0.9 - 1.0"
      ) |> factor(levels = c("0.0 - <0.25", "0.25 - <0.5", "0.5 - <0.75", "0.75 - <0.9", "0.9 - 1.0")) 
    )
  ] |>
  st_as_sf()
saveRDS(subcountyLevel, "../raster results/hookworm_subcountyLevel.RDS")

# Plotting after population weighted --------------------------------------

actualHook_shp <- actualHook_shp |> 
  mutate(prv_ctg = factor(prv_ctg, levels = c("<2%", "2%-<10%", "10%-<20%", "20%-<50%", ">=50%"))) |> 
  na.omit()

subcountyLevel <- subcountyLevel |> 
  mutate(prv_ctg = factor(prv_ctg, levels = c("<2%", "2%-<10%", "10%-<20%", "20%-<50%", ">=50%")))

prv_ctg_colors <- c(
  "<2%" = brewer.pal(5, "YlOrRd")[1],
  "2%-<10%" = brewer.pal(5, "YlOrRd")[2],
  "10%-<20%" = brewer.pal(5, "YlOrRd")[3],
  "20%-<50%" = brewer.pal(5, "YlOrRd")[4],
  ">=50%" = brewer.pal(5, "YlOrRd")[5]
)

prv_ctg_colors2 <- c(
  "0.0 - <0.25" = brewer.pal(5, "YlOrRd")[1], 
  "0.25 - <0.5" = brewer.pal(5, "YlOrRd")[2], 
  "0.5 - <0.75" = brewer.pal(5, "YlOrRd")[3], 
  "0.75 - <0.9" = brewer.pal(5, "YlOrRd")[4], 
  "0.9 - 1.0" = brewer.pal(5, "YlOrRd")[5]
)

prv_ctg_colors3 <- c(
  "0.0 - <0.25" = brewer.pal(5, "YlGn")[1], 
  "0.25 - <0.5" = brewer.pal(5, "YlGn")[2], 
  "0.5 - <0.75" = brewer.pal(5, "YlGn")[3], 
  "0.75 - <0.9" = brewer.pal(5, "YlGn")[4], 
  "0.9 - 1.0" = brewer.pal(5, "YlGn")[5]
)

p <- ggplot() +
  geom_sf(data = county_shapefile,
          fill = 'white',
          color = 'black') +
  geom_sf(data = st_as_sf(actualHook_shp),
          aes(
            fill = prv_ctg
          ),
          color = "black") +
  scale_fill_manual(values = prv_ctg_colors, na.value = 'white') +
  labs(fill = "", title = "Actual prevalence") +
  theme_void() +
  theme(legend.position = 'bottom', legend.direction = 'horizontal',
        legend.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(colour = "black", face = "bold", hjust = .5, size = 17)) +
 # annotation_north_arrow(style = north_arrow_fancy_orienteering(), location = "tr") +
  annotation_scale() 

p2 <- ggplot() +
  geom_sf(data = county_shapefile,
          fill = 'white',
          color = 'black') +
  geom_sf(data = st_as_sf(subcountyLevel) |> na.omit(),
          aes(
            fill = prv_ctg
          ),
          color = "black") +
  scale_fill_manual(values = prv_ctg_colors, na.value = 'white') +
  labs(fill = "", title = "Predicted prevalence") +
  theme_void() +
  theme(legend.position = 'bottom', legend.direction = 'horizontal',
        legend.text = element_text(size = 14, colour = "black"),
        plot.title = element_text(colour = "black", face = "bold", hjust = .5, size = 17)) +
# annotation_north_arrow(style = north_arrow_fancy_orienteering(), location = "tr") +
  annotation_scale() 
weighted_plots <- p + p2

ggsave(
  plot = weighted_plots,
  'plots/weighted_plots.png',
  dpi = 5e2,
  height = 12,
  width = 20,
  bg = NULL
)

# Continuous var plotting
pmodel_continous_ll <- ggplot() +
  geom_sf(data = county_shapefile,
          fill = 'white',
          color = 'black') +
  geom_sf(
    data = st_as_sf(subcountyLevel) |> na.omit(),
    aes(fill = lower_ci),
    color = "black"
  ) +
  scale_fill_distiller(palette = 'YlOrRd', direction = 1) +
  labs(fill = "") +
  theme_void() +
  theme(legend.position = 'bottom', legend.direction = 'horizontal') +
  annotation_north_arrow(style = north_arrow_fancy_orienteering(), location = "tr") +
  annotation_scale()

pmodel_continous <- ggplot() +
  geom_sf(data = county_shapefile,
          fill = 'white',
          color = 'black') +
  geom_sf(
    data = st_as_sf(subcountyLevel) |> na.omit(),
    aes(fill = prev),
    color = "black"
  ) +
  scale_fill_distiller(palette = 'YlOrRd', direction = 1) +
  labs(fill = "") +
  theme_void() +
  theme(legend.position = 'bottom', legend.direction = 'horizontal') +
  annotation_north_arrow(style = north_arrow_fancy_orienteering(), location = "tr") +
  annotation_scale()

pmodel_continous_upp <- ggplot() +
  geom_sf(data = county_shapefile,
          fill = 'white',
          color = 'black') +
  geom_sf(
    data = st_as_sf(subcountyLevel) |> na.omit(),
    aes(fill = upper_ci),
    color = "black"
  ) +
  scale_fill_distiller(palette = 'YlOrRd', direction = 1) +
  labs(fill = "") +
  theme_void() +
  theme(legend.position = 'bottom', legend.direction = 'horizontal') +
  annotation_north_arrow(style = north_arrow_fancy_orienteering(), location = "tr") +
  annotation_scale()
model_res <- pmodel_continous_ll +
  pmodel_continous +
  pmodel_continous_upp

ggsave(
  plot = model_res,
  'plots/model_res.png',
  dpi = 5e2,
  height = 12,
  width = 20,
  bg = NULL
)

# Plotting for exceedance 2%
p2 <- ggplot() +
  geom_sf(data = county_shapefile,
          fill = 'white',
          color = 'black') +
  geom_sf(data = st_as_sf(subcountyLevel[!is.na(subcountyLevel$exceed2_cat),]),
          aes(
            fill = exceed2_cat
          ),
          color = "black") +
  scale_fill_manual(values = prv_ctg_colors2, na.value = 'white') +
  labs(fill = "", title = "Probability of exceeding 2%") +
  theme_void() +
  theme(
    legend.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(colour = "black", face = "bold", hjust = .5, size = 17)) +
  #annotation_north_arrow(style = north_arrow_fancy_orienteering(), location = "tr") +
  annotation_scale() 

# Plotting for exceedance 10%
p10 <- ggplot() +
  geom_sf(data = county_shapefile,
          fill = 'white',
          color = 'black') +
  geom_sf(data = st_as_sf(subcountyLevel[!is.na(subcountyLevel$exceed10_cat),]),
          aes(
            fill = exceed10_cat
          ),
          color = "black") +
  scale_fill_manual(values = prv_ctg_colors2, na.value = 'white') +
  labs(fill = "", title = "Probability of exceeding 20%") +
  theme_void() +
  theme(
    legend.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(colour = "black", face = "bold", hjust = .5, size = 17)) +
 # annotation_north_arrow(style = north_arrow_fancy_orienteering(), location = "tr") +
  annotation_scale() 
exceedance <- p2 + p10

p2_not <- ggplot() +
  geom_sf(data = county_shapefile,
          fill = 'white',
          color = 'black') +
  geom_sf(data = st_as_sf(subcountyLevel[!is.na(subcountyLevel$notExceed2_cat),]),
          aes(
            fill = notExceed2_cat
          ),
          color = "black") +
  scale_fill_manual(values = prv_ctg_colors3, na.value = 'white') +
  labs(fill = "", title = "Probability of not exceeding 2%") +
  theme_void() +
  theme(
    legend.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(colour = "black", face = "bold", hjust = .5, size = 17)) +
  # annotation_north_arrow(style = north_arrow_fancy_orienteering(), location = "tr") +
  annotation_scale() 

ggsave(
  plot = exceedance,
  'plots/exceedance.png',
  dpi = 5e2,
  height = 12,
  width = 20,
  bg = NULL
)


ggsave(
  plot = p2,
  'plots/exceedance_two percent.png',
  dpi = 5e2,
  height = 12,
  width = 20,
  bg = NULL
)

ggsave(
  plot = p10,
  'plots/exceedance_twenty percent.png',
  dpi = 5e2,
  height = 12,
  width = 20,
  bg = NULL
)

ggsave(
  plot = p10,
  'plots/exceedance_ten percent.png',
  dpi = 5e2,
  height = 12,
  width = 20,
  bg = NULL
)

ggsave(
  plot = p2_not,
  'plots/exceedance_not two percent.png',
  dpi = 5e2,
  height = 12,
  width = 20,
  bg = NULL
)

# Obtaining some statistics after predicting ------------------------------
# Printing mean prevalence and CI
round(mean(subcountyLevel$prev, na.rm = T), 2)
round(quantile(subcountyLevel$prev, probs = c(0.025, 0.5, 0.975), na.rm = T), 2)
data.table(subcountyLevel) |> 
  na.omit() |> 
  _[, .N, by = prv_ctg]

# Probability of exceedance at the pixel level -----------------------------------------------

#Threshold calculations - 10%, 2%
prob.10 <- apply(inv.linpred, 1, \(.) length(which(. >= 0.1))/nsamp) 
ll=1:length(ind); ll[nonmiss] = prob.10; ll[miss] = NA
prob10 = raster(grid); values(prob10) = ll

prob.2 <- apply(inv.linpred, 1, \(.) length(which(. >= 0.02))/nsamp) 
ll=1:length(ind); ll[nonmiss] = prob.2; ll[miss] = NA
prob2 = raster(grid); values(prob2) = ll

# Forest plot for sub-counties ---------------------------------------------

mm_c <- merge(data.table(subcountyLevel)[, .(county, subcounty, lower_ci, upper_ci,
                                             pred_prev = prev, geometry)], 
              actualHook_shp[, .(county, subcounty, actual_prev = prev)],
              by = c("county", "subcounty")
              ) |> 
  st_as_sf()
nrow(mm_c)

mm_c |>
  ggplot(aes(x = subcounty)) +
  geom_point(aes(y = pred_prev, colour = 'Predicted')) +
  geom_point(aes(y = actual_prev, colour = "Actual")) +
  geom_errorbar(aes(
    ymin = lower_ci,
    ymax = upper_ci,
    y = pred_prev,
    colour = 'Predicted'
  )) +
  coord_flip() +
  scale_colour_manual(values = c("Predicted" = "#1f78b4", "Actual" = "#e31a1c")) +
  theme_classic()