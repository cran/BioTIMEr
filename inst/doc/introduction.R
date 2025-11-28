## -----------------------------------------------------------------------------
#| label: setup
#| include: false
#| cache: false

library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(dggridR)
library(knitr)

set.seed(19)


## ----warning=FALSE, message=FALSE, eval=TRUE----------------------------------
#| label: bt_package_loading
#| cache: false
# Install and load the latest version of BioTIMEr
library(BioTIMEr)


## -----------------------------------------------------------------------------
#| label: data_description
#| cache: false
# you can run the following commands to retrieve more information about the subsets.
?BTsubset_meta
?BTsubset_data


## -----------------------------------------------------------------------------
#| label: "package info"
#| cache: false
# you can also see a full list of BioTIMEr functions and help pages by:
??BioTIMEr


## -----------------------------------------------------------------------------
#| label: tbl-gridding_ex1
grid_samples <- gridding(
  meta = BTsubset_meta,
  btf = BTsubset_data,
  res = 12,
  resByData = FALSE
)

# Get a look at the output
grid_samples |> head() |> kable()


## -----------------------------------------------------------------------------
#| label: tbl-gridding_ex2
check_res_1 <- grid_samples |>
  summarise(
    n_cell = n_distinct(cell),
    n_aID = n_distinct(assemblageID),
    res = "res12",
    .by = c(STUDY_ID, StudyMethod)
  )

check_res_1 |> head(10) |> kable()

## -----------------------------------------------------------------------------
#| label: fig-gridding_ex
# How many samples were there in each study?
ggplot(data = check_res_1) +
  geom_bar(
    mapping = aes(x = as.character(STUDY_ID), y = n_aID, fill = res),
    stat = "identity"
  ) +
  scale_fill_discrete(type = c("#155f49")) +
  xlab("StudyID") +
  ylab("Number of assemblages in a study") +
  theme_bw() +
  theme(legend.position = "none")


## -----------------------------------------------------------------------------
#| label: gridding_ex2
#| results: hold
# define an alternative resolution of 14 (~10.7 km2 cells)
grid_samples_14 <- gridding(
  meta = BTsubset_meta,
  btf = BTsubset_data,
  res = 14,
  resByData = FALSE
)

# allow the spatial extent of the data to define the resolution
grid_samples_auto <- gridding(
  meta = BTsubset_meta,
  btf = BTsubset_data,
  res = 12,
  resByData = TRUE
)
# this option also returns a message with the automatically picked resolution:


## -----------------------------------------------------------------------------
#| label: fig-gridding_ex3
check_res_2 <- grid_samples_14 |>
  summarise(
    n_cell = n_distinct(cell),
    n_aID = n_distinct(assemblageID),
    res = "res14",
    .by = c(StudyMethod, STUDY_ID)
  )

check_res_3 <- grid_samples_auto |>
  summarise(
    n_cell = n_distinct(cell),
    n_aID = n_distinct(assemblageID),
    res = "res15",
    .by = c(StudyMethod, STUDY_ID)
  )

checks <- rbind(check_res_1, check_res_2, check_res_3)

ggplot(data = checks) +
  geom_bar(
    mapping = aes(x = as.character(STUDY_ID), y = n_aID, fill = res),
    stat = "identity",
    position = "dodge"
  ) +
  scale_fill_discrete(type = c("#155f49", "#66c1d1", "#d9d956")) +
  xlab("StudyID") +
  ylab("Number of assemblages in a study") +
  theme_bw()


## -----------------------------------------------------------------------------
#| label: fig-gridding_map
## The following example is built on the demonstrations in
## https://cran.r-project.org/web/packages/dggridR/vignettes/dggridR.html.

# First we build the ~96 km2 global grid
dgg_12 <- dggridR::dgconstruct(res = 12)

# To simplify, we only map the grid cell boundaries for cells which
# have observations.
# NOTE: if you are working with the full BioTIME database, this step may take some time.
grid_12 <- dggridR::dgcellstogrid(dgg_12, as.integer(grid_samples$cell))

# Now let's follow the same steps and build a ~10.7 km2 global grid:
dgg_14 <- dggridR::dgconstruct(res = 14)
grid_14 <- dggridR::dgcellstogrid(dgg_14, as.integer(grid_samples_14$cell))

# And we get some polygons for each country of the world, to create a background:
countries <- ggplot2::map_data("world")

# Now you could map the whole world, but let's just zoom in the UK and have a look at
# STUDY 466 (Marine Fish):
map_uk_locations <- ggplot() +
  geom_polygon(data = countries, aes(x = long, y = lat, group = group), fill = NA, color = "grey") +
  geom_sf(data = grid_12, aes(), color = alpha("blue", 0.4)) +
  coord_sf(xlim = c(-20, 10), ylim = c(50, 60)) +
  geom_rect(aes(xmin = -11, xmax = -0.7, ymin = 57.2, ymax = 59), colour = "red", fill = NA) +
  labs(x = NULL, y = NULL) +
  theme_bw() + theme(text = element_text(size = 8))

zoom_in_map <- ggplot() +
  geom_polygon(data = countries, aes(x = long, y = lat, group = group), fill = NA,
               color = "grey") +
  geom_sf(data = grid_12, aes(), color = alpha("blue", 0.4)) +
  geom_sf(data = grid_14, aes(), color = alpha("red", 0.4)) +
  coord_sf(xlim = c(-11, -0.7), ylim = c(57.2, 59)) +
  theme_bw()

grid::grid.newpage()
main <- grid::viewport(width = 1, height = 1, x = 0.5, y = 0.5)  # the main map
inset <- grid::viewport(width = 0.4, height = 0.4, x = 0.82, y = 0.45)  # the inset in bottom left

# The resulting distribution of different size cells appears as follows:
print(zoom_in_map, vp = main)
print(map_uk_locations, vp = inset)


## -----------------------------------------------------------------------------
#| label: tbl-resampling_ex
# First, if you are not sure you need this step,
# you can always check how many samples there are in every year of the different time series:
check_samples <- grid_samples |>
  summarise(
    n_samples = n_distinct(SAMPLE_DESC),
    n_species = n_distinct(Species),
    .by = c(STUDY_ID, assemblageID, YEAR)
  )

check_samples |> head(10) |> kable()



## -----------------------------------------------------------------------------
#| label: resampling_ex1
# Let's apply resampling() then, using the data frame of the gridded data:
grid_rare <- resampling(
  x = grid_samples,
  measure = "ABUNDANCE",
  resamps = 1,
  conservative = FALSE
)


## -----------------------------------------------------------------------------
#| label: tests1
#| include: false
# Let's apply it then:

#g rid_samples_temp <- subset(grid_samples, !grid_samples$BIOMASS==0) #to be deleted after Faye reviews the data and makes all 0=NAs
# grid_rare <- resampling( x= grid_samples_test, measure ="BIOMASS", resamps = 1, conservative = FALSE)


## -----------------------------------------------------------------------------
#| label: resampling_ex12
# Keep only observations with both abundance and biomass
grid_rare_ab <- resampling(
  x = grid_samples,
  measure = c("ABUNDANCE", "BIOMASS"),
  resamps = 1,
  conservative = FALSE
)

# Keep only sampling events where all observations within had both abundance and biomass to start with
grid_rare_abT <- resampling(
  x = grid_samples,
  measure = c("ABUNDANCE", "BIOMASS"),
  resamps = 1,
  conservative = TRUE)


## -----------------------------------------------------------------------------
#| label: fig-resampling_ex3
# What is the number of samples in the year with the lowest sampling effort?
ggplot(data = check_samples[check_samples$assemblageID == "18_335699",], aes(x = YEAR, y = n_samples)) +
  geom_col(aes(x = YEAR, y = n_samples), fill = "red", alpha = 0.5) +
  geom_segment(aes(x = 1926, y = min(n_samples) + 3,
                   xend = 1927, yend = min(n_samples)),
               arrow = arrow(length = unit(0.2, "cm"))) +
  xlab("Year") +
  ylab("Number of samples") +
  theme_bw()

# In this case,the year 1927 had the lowest sampling effort, with 3 samples (arrow).


## -----------------------------------------------------------------------------
#| label: tbl-resampling_ex4
# Let's implement the sample-based rarefaction by resampling the dataset 10 times.
grid_rare_n10 <- resampling(
  x = grid_samples,
  measure = "ABUNDANCE",
  resamps = 10,
  conservative = FALSE
)

# Note that you may want to resample many more times (e.g. at least 30-100+ times, but up to 199
# if e.g. working with the whole BioTIME data), depending on how many iterations you want a
# subsequent bootstrap analysis to have.
# This may also take some computation time, so if you are working with a big subset of BioTIME
# is advisable to break it down in smaller subsets.

# Each resampling iteration will be identified as resamp = 1:n, in this case 1:10.
# Now we can check if there are differences across the first 3 of these iterations:

check_resamps <- grid_rare_n10[grid_rare_n10$resamp < 4, ] |>
  summarise(
    n_obs = n(),
    n_species = n_distinct(Species),
    n_year = n_distinct(YEAR),
    .by = c(STUDY_ID, assemblageID, resamp)
  )

check_resamps |> head(10) |> kable()


## -----------------------------------------------------------------------------
#| label: metrics
# Get alpha metrics estimates:
alpha_metrics <- getAlphaMetrics(x = grid_rare, measure = "ABUNDANCE")
# see also help("getAlphaMetrics") for more details on the metrics

# Have a quick look at the output
alpha_metrics |> head(6) |> kable()

# Get beta metrics estimates:
beta_metrics <- getBetaMetrics(x = grid_rare, measure = "ABUNDANCE")
#see also help("getBetaMetrics") for more details on the metrics

# Have a quick look at the output
beta_metrics |> head(6) |> kable()
# NOTE the functions used the rarefied data with only one resampling iteration


## -----------------------------------------------------------------------------
#| label: tests2
#| include: false
# # Get alpha metrics estimates:
# alpha_metrics <- getAlphaMetrics(x = grid_rare_ab, measure = "BIOMASS")
# #see also help("getAlphaMetrics") for more details on the metrics
#
# #Have a quick look at the output
# alpha_metrics |> head(6) |> kable()
#
# # Get beta metrics estimates:
# beta_metrics <- getBetaMetrics(x = grid_rare_ab, measure = "BIOMASS")
# #see also help("getBetaMetrics") for more details on the metrics
#
# #Have a quick look at the output
# beta_metrics |> head(6) |> kable()


## -----------------------------------------------------------------------------
#| label: tbl-trends1
# Let's apply it then:
alpha_slopes <- getLinearRegressions(
  x = alpha_metrics,
  pThreshold = 0.05
) #for alpha metrics

# Have a quick look at the output
alpha_slopes |> head(6) |> kable()

## -----------------------------------------------------------------------------
#| label: tbl-trends2
# Let's apply it then:
beta_slopes <- getLinearRegressions(
  x = beta_metrics,
  pThreshold = 0.05
) #for beta metrics

# Have a quick look at the output
beta_slopes |> head(6) |> kable()


## -----------------------------------------------------------------------------
#| label: tbl-trends3
#| fig-width: 10
#| fig-height: 7
# First, how many assemblages in our dataset show a moderate evidence (P < 0.05) of change in alpha diversity?
alpha_slopes |>
  filter(significance == 1) |>   # or use filter(pvalue<0.05)
  summarise(n_sig = n_distinct(assemblageID), mean(slope), .by = metric) |>
  kable()

# We can see that only a few (<40) of the assemblage time series actually show a significant
# trend of change over time, independently of the metric used. This indicates that in most
# time series in the studies we analysed alpha diversity is not really changing through time.


## -----------------------------------------------------------------------------
#| label: fig-trends
#| fig-width: 7.2
#| fig-height: 5
# Get a slope per assemblageID and metric
alpha_slopes_simp <- alpha_slopes |>
  filter(significance == 1) |> #select only the assemblages with significant trends
  summarise(slope = unique(slope), .by = c(assemblageID, metric, pvalue))

# Calculate the mean slope and CIs
stats_alpha <- alpha_slopes_simp |>
  summarise(
    mean_slope = mean(slope), #mean
    ci_slope = qt(0.95, df = length(slope) - 1) *
      (sd(slope, na.rm = TRUE) / sqrt(length(slope))),
    .by = metric
  ) #margin of error

# Let's put it all together
ggplot(data = alpha_slopes_simp) +
  geom_histogram(aes(x = slope, fill = metric), bins = 25) +
  #geom_density(alpha=0.5)+ #in case you want a density plot instead
  geom_vline(
    aes(xintercept = 0),
    linetype = 3,
    colour = "grey",
    linewidth = 0.5
  ) + #add slope=0 line
  geom_vline(
    data = stats_alpha,
    aes(xintercept = mean_slope),
    linetype = 1,
    linewidth = 0.5,
    colour = "black"
  ) + #mean
  geom_vline(
    data = stats_alpha,
    aes(xintercept = mean_slope - ci_slope),
    linetype = 2,
    linewidth = 0.5,
    colour = "black"
  ) + #lower confidence interval
  geom_vline(
    data = stats_alpha,
    aes(xintercept = mean_slope + ci_slope),
    linetype = 2,
    linewidth = 0.5,
    colour = "black"
  ) + #upper confidence interval
  facet_wrap(~metric, scales = "free") +
  scale_fill_biotime() + #using the customize BioTIME colour scale. See help("scale_color_biotime") for more options.
  ggtitle("Alpha diversity change") +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11, hjust = 0.5)
  )



## -----------------------------------------------------------------------------
#| label: fig-trends2
#| fig-width: 7.2
#| fig-height: 1.9

# Get a slope per assemblageID and metric
beta_slopes_simp <- beta_slopes |>
  filter(significance == 1) |> #select only the assemblages with significant trends
  summarise(slope = unique(slope), .by = c(assemblageID, metric, pvalue))

# Calculate the mean slope and CIs
stats_beta <- beta_slopes_simp |>
  summarise(
    mean_slope = mean(slope), #mean
    ci_slope = qt(0.95, df = length(slope) - 1) *
      (sd(slope, na.rm = TRUE) / sqrt(length(slope))),
    .by = metric
  ) #margin of error

# Let's put it all together
ggplot(data = beta_slopes_simp) +
  geom_histogram(aes(x = slope, fill = metric), bins = 25) +
  #geom_density(alpha=0.5)+ #in case you want a density plot instead
  geom_vline(
    aes(xintercept = 0),
    linetype = 3,
    colour = "grey",
    linewidth = 0.5
  ) + #add slope=0 line
  geom_vline(
    data = stats_beta,
    aes(xintercept = mean_slope),
    linetype = 1,
    linewidth = 0.5,
    colour = "black"
  ) + #mean
  geom_vline(
    data = stats_beta,
    aes(xintercept = mean_slope - ci_slope),
    linetype = 2,
    linewidth = 0.5,
    colour = "black"
  ) + #lower confidence interval
  geom_vline(
    data = stats_beta,
    aes(xintercept = mean_slope + ci_slope),
    linetype = 2,
    linewidth = 0.5,
    colour = "black"
  ) + #upper confidence interval
  facet_wrap(~metric, scales = "free") +
  scale_fill_biotime() + #using the customize BioTIME colour scale. See help("scale_color_biotime") for more options.
  ggtitle("Beta diversity change") +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11, hjust = 0.5)
  )



## -----------------------------------------------------------------------------
#| label: trends5
# Hint: If you wish to plot all metrics together, simply merge you alpha_slopes and
# beta_slopes dataframes beforehand.


## -----------------------------------------------------------------------------
#| label: fig-trends6
#| fig-width: 7.2
#| fig-height: 2.3
# First, we need to check the meta file and retrieve the information for the studies of interest
#head(BTsubset_meta)
meta <- select(BTsubset_meta, STUDY_ID, TAXA, REALM, CLIMATE)

# Get a slope per assemblageID and metric
alpha_slopes_simp <- alpha_slopes |>
  summarise(
    slope = unique(slope),
    .by = c(assemblageID, metric, pvalue, significance)
  )

# Get back the Study ID by separating our assemblageID column into multiple other columns.
alpha_slopes_simp <- alpha_slopes_simp |>
  separate_wider_delim(
    cols = assemblageID,
    names = c("STUDY_ID", "cell"),
    delim = "_",
    cols_remove = FALSE
  ) |>
  as.data.frame()

# Merge it all
alpha_slopes_meta <- merge(alpha_slopes_simp, meta, by = "STUDY_ID")

# Select relevant data for plotting
alpha_slopes_set1 <- subset(
  alpha_slopes_meta,
  alpha_slopes_meta$metric == "Shannon" &
    (alpha_slopes_meta$TAXA == "Fish" | alpha_slopes_meta$TAXA == "Birds")
)

# Let's put it all together
ggplot(data = alpha_slopes_set1, aes(x = slope, fill = TAXA)) +
  geom_histogram(fill = "grey", bins = 25) + #all assemblages
  geom_histogram(
    data = alpha_slopes_set1[alpha_slopes_set1$significance == 1, ],
    bins = 25
  ) + #only significant change assemblages
  geom_vline(
    aes(xintercept = 0),
    linetype = 3,
    colour = "black",
    linewidth = 0.5
  ) + #add slope=0 line
  facet_wrap(~TAXA, scales = "free") +
  scale_fill_biotime() +
  ggtitle("Shannon-Weiner Species Diversity Index") +
  theme_bw() +
  theme(plot.title = element_text(size = 11, hjust = 0.5))



## -----------------------------------------------------------------------------
#| label: trends7
# Now we see the full distribution of slopes of change for the Shannon index for the birds
# and fish time series in our subset: all slopes are shown in grey, while in color we show
# the subset of assemblages for which evidence (P < 0.05) of change was detected.


## -----------------------------------------------------------------------------
#| label: fig-trends8
#| fig-width: 7.2
#| fig-height: 5

# Let's go back to the data frame with the yearly values and select our relevant data for plotting
alpha_metrics_set <- subset(alpha_metrics, assemblageID == "18_335699")

# Transform data
alpha_metrics_set_long <- pivot_longer(alpha_metrics_set,
                                   cols = c("S", "N", "maxN", "Shannon", "expShannon",
                                   "Simpson", "invSimpson", "PIE", "DomMc"),
                                   names_to = "metric",
                                   names_transform = as.factor)

# Get assemblage ID slope of change
alpha_slopes_set2 <- subset(alpha_slopes, assemblageID == "18_335699")

# Get assemblage ID
name_assemblage <- unique(alpha_slopes_set2$assemblageID)

# Merge info
alpha_metr_slop<- left_join(
  x = alpha_slopes_set2,
  y = alpha_metrics_set_long,
  by = join_by(assemblageID, metric))

# Let's put it all together
ggplot(data = alpha_metr_slop, aes(x = YEAR, y = value)) +
    geom_point(colour = "#155f49", size = 1) +   #plot the yearly estimates
    stat_smooth(method = "lm", se = FALSE, linetype = 2, colour = "grey") + #draw all regression line
    stat_smooth(data = subset(alpha_metr_slop, alpha_metr_slop$significance == 1),
                aes(x = YEAR, y = value), method = "lm", se = FALSE,
                linetype = 2, colour = "#155f49") + #color only trends that are significant
    facet_wrap(~metric, scales = "free") +
    scale_fill_biotime() +
    ggtitle(paste("Assemblage", name_assemblage)) + ylab("Diversity") +
    theme_bw() +
    theme(plot.title = element_text(size = 11, hjust = 0.5))



## -----------------------------------------------------------------------------
#| label: trends9
#| include: false

# Hint: If you want to draw the p-value on the plot, you can try using the #ggpmisc::stat_fit_glance which takes anything passed through lm() in R and
#allows it to processed and printed using ggplot2.


## -----------------------------------------------------------------------------
#| label: citation
#| cache: false

citation("BioTIMEr")

