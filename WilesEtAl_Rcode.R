devtools::install_github('NeotomaDB/neotoma2', force = TRUE)
install.packages("tidyverse")
install.packages("readxl")

library(neotoma2)
library(rioja)
library(vegan)
library(tidyverse)
library(dplyr)
library(RRatepol)
library(viridis)
library(readxl)
library(gridExtra)
library(cowplot)
library(pacman)

# Set bounding coordinates to capture sites only within Lower MI
lower.mich <- '{"type": "Polygon",
              "coordinates": [[
            [
              -84.72185848189304,
              45.759888069635394
            ],
            [
              -85.09491447007726,
              45.58610573744363
            ],
            [
              -85.21926646613841,
              45.30984667381506
            ],
            [
              -85.34361846219996,
              45.0395574329832
            ],
            [
              -85.46797045826149,
              44.76798529984478
            ],
            [
              -85.64413578601486,
              45.13467010950731
            ],
            [
              -86.01719177419906,
              44.87823673234129
            ],
            [
              -86.3487970970296,
              44.21356571056529
            ],
            [
              -86.4938744257682,
              44.05738387790626
            ],
            [
              -86.43169842773744,
              43.84849894028574
            ],
            [
              -86.5353250911219,
              43.63887990545999
            ],
            [
              -86.15190643659916,
              42.762643642841425
            ],
            [
              -86.34888439575836,
              42.206070858452364
            ],
            [
              -86.49686136408835,
              42.10246409149872
            ],
            [
              -86.80925718611714,
              41.77223572184076
            ],
            [
              -84.80334717098384,
              41.76610425881796
            ],
            [
              -84.77868434292873,
              41.71089472211128
            ],
            [
              -83.47155445601769,
              41.72316760970082
            ],
            [
              -83.17560051935897,
              42.03533405287635
            ],
            [
              -83.1344958059338,
              42.2608528272018
            ],
            [
              -82.88786752538465,
              42.39456409234691
            ],
            [
              -82.88786752538465,
              42.49163010361863
            ],
            [
              -82.83032092659003,
              42.57643948561662
            ],
            [
              -82.76455338511039,
              42.66717886139358
            ],
            [
              -82.65768113020556,
              42.68531086633996
            ],
            [
              -82.52614604724563,
              42.61880098458735
            ],
            [
              -82.50148321919117,
              42.800024024572934
            ],
            [
              -82.42749473502647,
              42.98071782805525
            ],
            [
              -82.4850413338211,
              43.08887997324277
            ],
            [
              -82.54258793261572,
              43.44207450002557
            ],
            [
              -82.67412301557563,
              43.917707802152194
            ],
            [
              -82.89608846806969,
              44.05375439623339
            ],
            [
              -83.23314711815357,
              43.976896806800255
            ],
            [
              -83.49621728407278,
              43.69225212969249
            ],
            [
              -83.67707802314229,
              43.5970688889478
            ],
            [
              -83.91548536100639,
              43.68036247245078
            ],
            [
              -83.94836913174653,
              43.816951238717564
            ],
            [
              -83.86615970489683,
              44.00646920672585
            ],
            [
              -83.70174085119739,
              44.01238191844956
            ],
            [
              -83.54554294018237,
              44.106905140719164
            ],
            [
              -83.53732199749795,
              44.266070869531006
            ],
            [
              -83.3153565450033,
              44.35431055872073
            ],
            [
              -83.29891465963321,
              44.857657119628186
            ],
            [
              -83.43044975241611,
              44.92754475132142
            ],
            [
              -83.43044975241611,
              45.084481407949056
            ],
            [
              -83.32357749751128,
              45.142496986458895
            ],
            [
              -83.46333352315561,
              45.33930880775543
            ],
            [
              -84.1374508233234,
              45.59874679503446
            ],
            [
              -84.72185848189304,
              45.759888069635394
            ]
          ]
        ]
      }'

# download datasets within the defined geography, narrowing to pollen records
lowermich_sites <- neotoma2::get_sites(loc = lower.mich, all_data = TRUE)
lowermich_datasets <- neotoma2::get_datasets(lowermich_sites, all_data = TRUE)

lowermich_pollen <- lowermich_datasets %>%
  neotoma2::filter(datasettype == "pollen" & !is.na(age_range_young))

lowermich_dl <- lowermich_pollen %>% get_downloads(all_data = TRUE) # download sample data
# based on features we've filtered


allSamp <- samples(lowermich_dl) # extract all samples we've collected; note that
# this will not yet include Sunrise Lake record so those data will be added manually

# constrain taxa by overall abundance throughout record
plottingTaxa <- taxa(lowermich_dl) %>%
  filter(ecologicalgroup %in% c("UPHE", "TRSH")) %>%
  arrange(desc(samples)) %>% 
  head(n = 10)

# edit/harmonize and filter variables of extracted taxa data
allSamp <- allSamp %>%
   mutate(variablename = replace(variablename,
                                 stringr::str_detect(variablename, "Pinus.*"),
                                 "Pinus"),
          variablename = replace(variablename,
                                 stringr::str_detect(variablename, "Acer.*"),
                                 "Acer"),
          variablename = replace(variablename,
                                 stringr::str_detect(variablename, "Alnus.*"),
                                 "Alnus"),
          variablename = replace(variablename,
                                 stringr::str_detect(variablename, "Fraxinus*"),
                                 "Fraxinus"),
          variablename = replace(variablename,
                                 stringr::str_detect(variablename, "Picea.*"),
                                 "Picea")) %>%
  filter(variablename %in% c(plottingTaxa$variablename, "Acer", "Alnus", "Ambrosia", "Artemisia", "Fraxinus", "Pinus")) %>%
  group_by(siteid, sitename,
           sampleid, variablename, units, age,
           agetype, depth, datasetid,
           long, lat) %>%
  summarise(value = sum(value), .groups='keep')

# sort entries by sample age, taxa, and site ID respectively
allSites <- allSamp %>%
  group_by(age) %>%
  mutate(pollencount = sum(value, na.rm = TRUE)) %>%
  group_by(variablename) %>%
  mutate(prop = value / pollencount) %>%
  group_by(siteid) %>%
  arrange(desc(age), .by_group = TRUE)

# Spread the data to a "wide" table, with taxa as column headings.
widetable <- allSites %>%
  dplyr::select(siteid, sitename, sampleid, long, lat, age, variablename, prop) %>%
  mutate(prop = as.numeric(prop))
  #filter(variablename %in% c("Pinus", "Betula", "Quercus",
   #                          "Tsuga", "Ulmus", "Picea"))

write.csv(x = widetable, file = "LowerMichPollen.csv")

# Count the abundance (as a proportion) of each taxon at a given depth
counts <- tidyr::pivot_wider(widetable,
                             id_cols = c(sampleid, siteid, sitename, long, lat, age),
                             names_from = variablename,
                             values_from = prop,
                             values_fill = 0)

# Convert species counts to percent abundances
counts_spp <- counts[,7:21]
counts_spp_prop <- 100*(counts_spp / rowSums(counts_spp))

# add in Sunrise pollen data (full pollen timeseries code below)
rise.poll <- readxl::read_excel("RISE_Poll_2.xlsx")
rise.sum <- colSums(rise.poll[,3:126]); rise.sum
rise.avg <- mean(col(rise.poll))
rise.poll2 <- rise.poll[, which(rise.sum > 25) + 2]
rise.poll2 <- cbind(rise.poll[,1], rise.poll2[,-c(1:4, 9, 10, 12, 19, 20)])
colnames(rise.poll2)[c(5, 11:15)] <- c("Fagus", "Ambrosia", "Artemisia", "Pinus", "Acer", "Fraxinus")
rise.poll2 <- cbind("long" = -85.329383, "lat" = 44.032833, "sitename" = "Sunrise Lake", rise.poll2, rise.poll$`Picea undiff`)
rise_spp <- rise.poll2[,5:19]
rise_spp <- rise_spp %>%
  select("Acer", "Alnus", "Ambrosia", "Artemisia", "Betula", "Carya", "Fagus",
         "Fraxinus", "Picea", "Pinus", "Quercus", "Tilia", "Tsuga", "Ulmus", "Ostrya/Carpinus")
rise_spp_prop <- rise_spp / 100
rise_spp_prop <- 100*(rise_spp_prop / rowSums(rise_spp_prop))
counts_spp_prop <- rbind(counts_spp_prop, rise_spp_prop)
counts_spp_prop_mat <- as.matrix(counts_spp_prop)

# Run NMDS analysis; for first biplot, use default choices=c(1,2), and for second change to changes=c(3,4)
pollen.nmds <- vegan::metaMDS(counts_spp_prop_mat, k=4, distance = "chisq", autotransform = FALSE)
plot(pollen.nmds, choices=c(1,2), type = 'n', cex.axis = 1.5, cex.lab = 1.5)
stems <- colSums(counts_spp_prop)
cl <- hclust(vegdist(counts_spp_prop, method='chisq'))
vegan::ordihull(pollen.nmds, cutree(cl, 3), choices=c(1,2), col = cutree(cl, 3))
vegan::ordicluster(pollen.nmds, cl, choices=c(1,2), prune=3, col = cutree(cl, 3), lty = 3, lwd = .7)
vegan::orditorp(pollen.nmds,display="species", priority = stems, pch = '+', choices=c(1,2), font = 2, cex = 1.31, col="red", air = 0.53)

# Test stress of NMDS, i.e. goodness of fit
stressplot(pollen.nmds)

rise.meta <- rise.poll2[,1:5]
rise.meta <- cbind(rise.meta, "siteid" = 1)
rise.meta[,1] <- as.character(rise.meta[,1])
rise.meta[,6] <- as.character(rise.meta[,6])
rise.meta <- cbind(rise.meta, "sampleid" = 1:63)
rise.meta <- rise.meta %>%
  select('sampleid', 'siteid', 'sitename', 'long', 'lat', 'age')

# create and format dataframe for appending site data to NMDS scores
lower.mich.nmds <- as.data.frame(pollen.nmds$points)
lower.mich.nmds <- cbind(lower.mich.nmds, rbind(counts[,1:6], rise.meta))
data.scores <- as.data.frame(lower.mich.nmds)
data.scores$MDS1 <- as.numeric(data.scores$MDS1)
data.scores$MDS2 <- as.numeric(data.scores$MDS2)
data.scores$MDS3 <- as.numeric(data.scores$MDS3)
data.scores$MDS4 <- as.numeric(data.scores$MDS4)
data.scores$age <- as.numeric(data.scores$age)

# subset each site values from the NMDS results to create individual NMDS timeseries
chip.data <- data.scores[data.scores$siteid == 360, ]
chip.clust <- rioja::chclust(dist(sqrt(chip.data[1:5])), method = 'coniss')

demont.data <- data.scores[data.scores$siteid == 661, ]
demont.clust <- rioja::chclust(dist(sqrt(demont.data[1:5])), method = 'coniss')

minden.data <- data.scores[data.scores$siteid == 13582, ]
minden.clust <- rioja::chclust(dist(sqrt(minden.data[1:5])), method = 'coniss')

irwin.data <- data.scores[data.scores$siteid == 8559, ]
irwin.clust <- rioja::chclust(dist(sqrt(irwin.data[1:5])), method = 'coniss')

wintergreen.data <- data.scores[data.scores$siteid == 2864, ]
wintergreen.clust <- rioja::chclust(dist(sqrt(wintergreen.data[1:5])), method = 'coniss')

cub.data <- data.scores[data.scores$siteid == 518, ]
cub.clust <- rioja::chclust(dist(sqrt(cub.data[1:5])), method = 'coniss')

rise.data <- data.scores[data.scores$sitename == "Sunrise Lake", ]
rise.clust <- rioja::chclust(dist(sqrt(rise.data[1:5])), method = 'coniss')

hicks.data <- data.scores[data.scores$sitename == "Hicks Lake", ]
hicks.clust <- rioja::chclust(dist(sqrt(hicks.data[1:5])), method = 'coniss')

vest.data <- data.scores[data.scores$sitename == "Vestaburg Bog", ]
vest.clust <- rioja::chclust(dist(sqrt(vest.data[1:5])), method = 'coniss')

otter.data <- data.scores[data.scores$sitename == "Otter Lake", ]
otter.clust <- rioja::chclust(dist(sqrt(otter.data[1:5])), method = 'coniss')

cowden.data <- data.scores[data.scores$sitename == "Cowden Lake", ]
cowden.clust <- rioja::chclust(dist(sqrt(cowden.data[1:5])), method = 'coniss')

frains.data <- data.scores[data.scores$sitename == "Frains Lake", ]
frains.clust <- rioja::chclust(dist(sqrt(frains.data[1:5])), method = 'coniss')

green.data <- data.scores[data.scores$sitename == "Green Lake", ]
green.clust <- rioja::chclust(dist(sqrt(green.data[1:5])), method = 'coniss')

horse.data <- data.scores[data.scores$sitename == "Horseshoe Lake", ]
horse.clust <- rioja::chclust(dist(sqrt(horse.data[1:5])), method = 'coniss')

obrien.data <- data.scores[data.scores$sitename == "O'Brien Lake", ]
obrien.clust <- rioja::chclust(dist(sqrt(obrien.data[1:5])), method = 'coniss')

cliff.data <- data.scores[data.scores$sitename == "Clifford Lake", ]
cliff.clust <- rioja::chclust(dist(sqrt(cliff.data[1:5])), method = 'coniss')

green.data <- data.scores[data.scores$sitename == "Green Lake", ]
green.clust <- rioja::chclust(dist(sqrt(green.data[1:5])), method = 'coniss')

morrison.data <- data.scores[data.scores$sitename == "Morrison Lake", ]
morrison.clust <- rioja::chclust(dist(sqrt(morrison.data[1:5])), method = 'coniss')

huff.data <- data.scores[data.scores$sitename == "Huffman Lake", ]
huff.clust <- rioja::chclust(dist(sqrt(huff.data[1:5])), method = 'coniss')

lake16.data <- data.scores[data.scores$sitename == "Lake Sixteen", ]
lake16.clust <- rioja::chclust(dist(sqrt(lake16.data[1:5])), method = 'coniss')

lake27.data <- data.scores[data.scores$sitename == "Lake 27", ]
lake27.clust <- rioja::chclust(dist(sqrt(lake27.data[1:5])), method = 'coniss')

# plot the NMDS timeseries for each site
chip.plot <- rioja::strat.plot(chip.data[1:4], yvar = chip.data$age, y.rev = TRUE, title = chip.data$sitename[1], ylabel = "Calibrated Years BP", ylabPos = 2.7, xlabel = "NMDS Scores", clust = chip.clust, plot.line = TRUE, col.line = "blue", col.bar = "blue")
demont.plot <- rioja::strat.plot(demont.data[1:4], yvar = demont.data$age, y.rev = TRUE, title = demont.data$sitename[1], ylabel = "Calibrated Years BP", ylim = c(min(demont.data$age), max(demont.data$age)), ylabPos = 2.7, xlabel = "NMDS Scores", clust = demont.clust, plot.line = TRUE, col.line = "darkgreen", col.bar="darkgreen")
minden.plot <- rioja::strat.plot(minden.data[1:4], yvar = minden.data$age, y.rev = TRUE, title = minden.data$sitename[1], ylabel = "Calibrated Years BP", ylabPos = 2.7, xlabel = "NMDS Scores", clust = minden.clust, plot.line = TRUE, col.line = "darkred", col.bar = "darkred")
irwin.plot <- rioja::strat.plot(irwin.data[1:4], yvar = irwin.data$age, y.rev = TRUE, title = irwin.data$sitename[1], ylabel = "Calibrated Years BP", ylabPos = 2.7, xlabel = "NMDS Scores", clust = irwin.clust, plot.line = TRUE, col.line = "gold", col.bar = "gold")
wintergreen.plot <- rioja::strat.plot(wintergreen.data[1:4], yvar = wintergreen.data$age, y.rev = TRUE, title = wintergreen.data$sitename[1], ylabel = "Calibrated Years BP", ylabPos = 2.7, xlabel = "NMDS Scores", clust = wintergreen.clust, plot.line = TRUE, col.line = "pink", col.bar = "pink")
cub.plot <- rioja::strat.plot(cub.data[1:4], yvar = cub.data$age, y.rev = TRUE, title = cub.data$sitename[1], ylabel = "Calibrated Years BP", ylabPos = 2.7, xlabel = "NMDS Scores", clust = cub.clust, plot.line = TRUE, col.line = "lightgreen", col.bar = "lightgreen")
rise.plot <- rioja::strat.plot(rise.data[,c(2,4)], clust.width = 0, yvar = rise.data$age, y.rev = TRUE, title = rise.data$sitename[1], ylabel = "Calibrated Years BP", ylabPos = 2.7, xlabel = "NMDS Scores", clust = rise.clust, plot.line = TRUE, col.line = "red", col.bar = "red")
hicks.plot <- rioja::strat.plot(hicks.data[1:4], yvar = hicks.data$age, y.rev = TRUE, title = hicks.data$sitename[1], ylabel = "Calibrated Years BP", ylabPos = 2.7, xlabel = "NMDS Scores", clust = hicks.clust, plot.line = TRUE, col.line = "lightblue", col.bar = "lightblue")
vest.plot <- rioja::strat.plot(vest.data[1:4], yvar = vest.data$age, y.rev = TRUE, title = vest.data$sitename[1], ylabel = "Calibrated Years BP", ylabPos = 2.7, xlabel = "NMDS Scores", clust = vest.clust, plot.line = TRUE, col.line = "orange", col.bar = "orange")
otter.plot <- rioja::strat.plot(otter.data[1:4], yvar = otter.data$age, y.rev = TRUE, title = otter.data$sitename[1], ylabel = "Calibrated Years BP", ylabPos = 2.7, xlabel = "NMDS Scores", clust = otter.clust, plot.line = TRUE, col.line = "purple", col.bar = "purple")
cowden.plot <- rioja::strat.plot(cowden.data[1:4], yvar = cowden.data$age, y.rev = TRUE, title = cowden.data$sitename[1], ylabel = "Calibrated Years BP", ylabPos = 2.7, xlabel = "NMDS Scores", clust = cowden.clust, plot.line = TRUE, col.line = "cyan", col.bar = "cyan")
frains.plot <- rioja::strat.plot(frains.data[1:4], yvar = frains.data$age, y.rev = TRUE, title = frains.data$sitename[1], ylabel = "Calibrated Years BP", ylabPos = 2.7, xlabel = "NMDS Scores", clust = frains.clust, plot.line = TRUE, col.line = "black", col.bar = "black")
green.plot <- rioja::strat.plot(green.data[1:4], yvar = green.data$age, y.rev = TRUE, title = green.data$sitename[1], ylabel = "Calibrated Years BP", ylabPos = 2.7, xlabel = "NMDS Scores", clust = green.clust, plot.line = TRUE, col.line = "gray", col.bar = "gray")
horse.plot <- rioja::strat.plot(horse.data[1:4], yvar = horse.data$age, y.rev = TRUE, title = horse.data$sitename[1], ylabel = "Calibrated Years BP", ylabPos = 2.7, xlabel = "NMDS Scores", clust = horse.clust, plot.line = TRUE, col.line = "maroon", col.bar = "maroon")
obrien.plot <- rioja::strat.plot(obrien.data[1:4], yvar = obrien.data$age, y.rev = TRUE, title = obrien.data$sitename[1], ylabel = "Calibrated Years BP", ylabPos = 2.7, xlabel = "NMDS Scores", clust = obrien.clust, plot.line = TRUE, col.line = "darkblue", col.bar = "darkblue")
cliff.plot <- rioja::strat.plot(cliff.data[1:4], yvar = cliff.data$age, y.rev = TRUE, title = cliff.data$sitename[1], ylabel = "Calibrated Years BP", ylabPos = 2.7, xlabel = "NMDS Scores", clust = cliff.clust, plot.line = TRUE, col.line = "hotpink", col.bar = "hotpink")
morrison.plot <- rioja::strat.plot(morrison.data[1:4], yvar = morrison.data$age, y.rev = TRUE, title = morrison.data$sitename[1], ylabel = "Calibrated Years BP", ylabPos = 2.7, xlabel = "NMDS Scores", clust = morrison.clust, plot.line = TRUE, col.line = "plum", col.bar = "plum")
huff.plot <- rioja::strat.plot(huff.data[1:4], yvar = huff.data$age, y.rev = TRUE, title = huff.data$sitename[1], ylabel = "Calibrated Years BP", ylabPos = 2.7, xlabel = "NMDS Scores", clust = huff.clust, plot.line = TRUE, col.line = "yellow", col.bar = "yellow")
lake16.plot <- rioja::strat.plot(lake16.data[1:4], yvar = lake16.data$age, y.rev = TRUE, title = lake16.data$sitename[1], ylabel = "Calibrated Years BP", ylabPos = 2.7, xlabel = "NMDS Scores", clust = lake16.clust, plot.line = TRUE, col.line = "green", col.bar = "green")
lake27.plot <- rioja::strat.plot(lake27.data[1:4], yvar = lake27.data$age, y.rev = TRUE, title = lake27.data$sitename[1], ylabel = "Calibrated Years BP", ylabPos = 2.7, xlabel = "NMDS Scores", clust = lake27.clust, plot.line = TRUE, col.line = "magenta", col.bar = "magenta")


rioja::addClustZone(minden.plot, minden.clust, 4, col = "orange")

subset.nmds <- as.data.frame(chip.data)
subset.nmds <- rbind(chip.data,demont.data, minden.data, irwin.data, wintergreen.data, cub.data)

library(ggplot2)
ggplot(subset.nmds, aes(x = age, y = MDS2, fill = siteid)) +
  geom_area() +
  #geom_area(aes(y = MDS1 * 2), alpha = 0.3, show.legend = FALSE) +
  scale_x_reverse(expand = c(0.02, 0)) +
  facet_grid(vars(sitename)) +
  labs(y = "NMDS Scores", x = "Calibrated Years BP")


# Pull ecotone data from ArcGIS Pro via Excel sheet to extract position and distance values
avg_nmds_positions <- read.csv(file="Ecotone_position_full.csv")
colnames(avg_nmds_positions)[4:14] <- c("0-1ka", "1-2ka", "2-3ka", "3-4ka", "4-5ka", "5-6ka", "6-7ka","7-8ka", "8-9ka", "9-10ka", "10-11ka")
nmds_positions.long <- avg_nmds_positions %>%
  dplyr::select("Latitude", "0-1ka", "1-2ka", "2-3ka", "3-4ka", "4-5ka", "5-6ka", "6-7ka","7-8ka", "8-9ka", "9-10ka", "10-11ka") %>%
  tidyr::pivot_longer(-Latitude, names_to = "Timebins", values_to = "values")

site_eco_dists <- read_excel("Ecotone_site_distances.xlsx")
site_eco_dists <- cbind(site_eco_dists, "Timebins" = c("0-1ka", "1-2ka", "2-3ka", "3-4ka", "4-5ka", "5-6ka", "6-7ka","7-8ka", "8-9ka", "9-10ka", "10-11ka"))
init <- 1; site_eco_dists <- cbind(site_eco_dists, "Average Distance")
dist.list <- c()
for (x in site_eco_dists$Sitename) {
  avg.dist <- mean(site_eco_dists$AbsDistkm[site_eco_dists$Sitename == x])
  site_eco_dists$`"Average Distance"`[site_eco_dists$Sitename == x] <- avg.dist
  dist.list[init] <- avg.dist
  init <- init + 1
}
dist.list <- unique(dist.list)
dist.list <- cbind(dist.list, "Sitename" = unique(site_eco_dists$Sitename))
dist.ecotone.mtrx <- matrix(data=dist.list, nrow = 20, ncol = 2)
dist.ecotone.mtrx <- as.data.frame(dist.ecotone.mtrx)
dist.ecotone.mtrx$V2[17] <- "Lake 27"; dist.ecotone.mtrx$V2[19] <- "Lake Sixteen"; dist.ecotone.mtrx$V2[15] <- "Irwin Smith"
dist.ecotone.mtrx <- dist.ecotone.mtrx %>%
  arrange(V2)
lower.mich.roc <- lower.mich.roc %>%
  arrange(Sitename)
lower.mich.roc <- cbind(lower.mich.roc, "Distance" = dist.ecotone.mtrx$V1)

# Plot NMDS scores, mean RoC, and max RoC as functions of ecotonal proximity
mdstimeplot <- ggplot(nmds_positions.long, aes(Latitude, values, colour = Timebins)) + 
  geom_point() + geom_smooth(method = lm, se = FALSE, linetype = 2, linewidth = 0.4) +
  labs(y="MDS2 Averages", fill = "Timebins") + scale_color_discrete(breaks=c("0-1ka", "1-2ka", "2-3ka", "3-4ka", "4-5ka", "5-6ka", "6-7ka","7-8ka", "8-9ka", "9-10ka", "10-11ka"))

meanrocdistplot <- ggplot(lower.mich.roc, aes(as.numeric(Distance), roc.site.scores, colour = Latitude)) + 
  geom_point() + geom_smooth(method = "lm", se = TRUE, linetype = 2, linewidth = 0.4) + stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"))) + scale_color_viridis(discrete = FALSE) +
  labs(y="Mean Rate of Change", x = "Average Absolute Distance to Ecotone (km)", fill = 'Latitude') + 
  scale_x_continuous(limits=c(0, 300), breaks=seq(0, 300, 50))

maxrocdistplot <- ggplot(lower.mich.roc, aes(as.numeric(Distance), roc.max, colour = Latitude)) + 
  geom_point() + geom_smooth(method = "lm", se = TRUE, linetype = 2, linewidth = 0.4) + stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"))) + scale_color_viridis(discrete = FALSE) +
  labs(y="Maximum Rate of Change", x = "Average Absolute Distance to Ecotone (km)", fill = 'Latitude') + 
  scale_x_continuous(limits=c(0, 300), breaks=seq(0, 300, 50))

legend <- get_legend(meanrocdistplot)
maxrocdistplot <- maxrocdistplot + theme_classic() + theme(axis.title = element_text(size = 16), axis.text = element_text(size = 12)) + theme(legend.position = "none")
meanrocdistplot <- meanrocdistplot + theme_classic() +  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 12)) + theme(legend.position = "none")
#blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
 # cowplot::theme_nothing()
grid.arrange(meanrocdistplot, legend, maxrocdistplot, legend, nrow=2, ncol = 2, widths=c(5, 1), heights=c(3.5, 3.5))

plot_grid(mdstimeplot, meanrocdistplot, maxrocdistplot, legend, nrow=3, ncol=1)


# Run Rate-of-Change analyses on each lake
rise.poll2 <- cbind("sampleid" = 1:63, rise.poll2)
rise.community <- as.data.frame(cbind("sample_id" = as.character(rise.poll2$sampleid), rise_spp_prop))
rise.age <- as.data.frame(cbind("sample_id" = as.character(rise.poll2$sampleid), "age" = as.numeric(rise.poll2$age)))
rise.age[,2] <- as.numeric(rise.age[,2])

rise.roc.sequence <-
  estimate_roc(
    data_source_community = rise.community,
    data_source_age = rise.age,
    age_uncertainty = NULL,
    smooth_method = "shep",
    working_units = "MW",
    rand = 1e3,
    dissimilarity_coefficient = "chisq"
  )

rise.seq.peak <-
  detect_peak_points(
    rise.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  rise.seq.peak,
  age_threshold = 11.8e3,
  roc_threshold = 2,
  peaks = TRUE,
  trend = "trend_non_linear"
) 

mean(rise.seq.peak$ROC)

cub.community <- as.data.frame(cbind("sample_id" = 
                                       as.character(counts$sampleid[counts$sitename 
                                       == "Cub Lake"]), counts[counts$sitename 
                                       == "Cub Lake", 7:23]))
cub.age <- as.data.frame(cbind("sample_id" = 
                               as.character(counts$sampleid[counts$sitename 
                               == "Cub Lake"]), age = counts$age[counts$sitename
                               == "Cub Lake"]))
cub.age[,2] <- as.numeric(cub.age[,2])

cub.roc.sequence <- estimate_roc(
  data_source_community = cub.community, 
  data_source_age = cub.age,
  age_uncertainty = NULL,
  smooth_method = "shep",
  working_units = "MW",
  rand = 1e3,
  dissimilarity_coefficient = "chisq")

cub.seq.peak <-
  detect_peak_points(
    cub.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  cub.seq.peak,
  age_threshold = 11.8e3,
  roc_threshold = 2,
  peaks = TRUE,
  trend = "trend_non_linear"
) 

cowden.community <- as.data.frame(cbind("sample_id" = 
                                       as.character(counts$sampleid[counts$sitename 
                                                                    == "Cowden Lake"]), counts[counts$sitename 
                                                                                            == "Cowden Lake", 7:23]))
cowden.age <- as.data.frame(cbind("sample_id" = 
                                 as.character(counts$sampleid[counts$sitename 
                                                              == "Cowden Lake"]), age = counts$age[counts$sitename
                                                                                                == "Cowden Lake"]))
cowden.age[,2] <- as.numeric(cowden.age[,2])

cowden.roc.sequence <- estimate_roc(
  data_source_community = cowden.community, 
  data_source_age = cowden.age,
  age_uncertainty = NULL,
  smooth_method = "shep",
  working_units = "MW",
  rand = 1e3,
  dissimilarity_coefficient = "chisq")

cowden.seq.peak <-
  detect_peak_points(
    cowden.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  cowden.seq.peak,
  age_threshold = 3e3,
  roc_threshold = 1,
  peaks = TRUE,
  trend = "trend_non_linear"
) 
mean(cowden.roc.sequence2$ROC)

chip.community <- as.data.frame(cbind("sample_id" = 
                                       as.character(counts$sampleid[counts$sitename 
                                                                    == "Chippewa Bog"]), counts[counts$sitename 
                                                                                            == "Chippewa Bog", 7:23]))
chip.age <- as.data.frame(cbind("sample_id" = 
                                 as.character(counts$sampleid[counts$sitename 
                                                              == "Chippewa Bog"]), age = counts$age[counts$sitename
                                                                                                == "Chippewa Bog"]))
chip.age[,2] <- as.numeric(chip.age[,2])

chip.roc.sequence <- estimate_roc(
  data_source_community = chip.community, 
  data_source_age = chip.age,
  age_uncertainty = NULL,
  smooth_method = "shep",
  working_units = "MW",
  rand = 1e3,
  dissimilarity_coefficient = "chisq")

chip.seq.peak <-
  detect_peak_points(
    chip.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  chip.seq.peak,
  age_threshold = 11.8e3,
  roc_threshold = 2,
  peaks = TRUE,
  trend = "trend_non_linear"
) 

cliff.community <- as.data.frame(cbind("sample_id" = 
                                       as.character(counts$sampleid[counts$sitename 
                                                                    == "Clifford Lake"]), counts[counts$sitename 
                                                                                            == "Clifford Lake", 7:23]))
cliff.age <- as.data.frame(cbind("sample_id" = 
                                 as.character(counts$sampleid[counts$sitename 
                                                              == "Clifford Lake"]), age = counts$age[counts$sitename
                                                                                                == "Clifford Lake"]))
cliff.age[,2] <- as.numeric(cliff.age[,2])

cliff.roc.sequence <- estimate_roc(
  data_source_community = cliff.community, 
  data_source_age = cliff.age,
  age_uncertainty = NULL,
  smooth_method = "shep",
  working_units = "MW",
  rand = 1e3,
  dissimilarity_coefficient = "chisq")

cliff.seq.peak <-
  detect_peak_points(
    cliff.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  cliff.seq.peak,
  age_threshold = 11.8e3,
  roc_threshold = 2,
  peaks = TRUE,
  trend = "trend_non_linear"
) 

demont.community <- as.data.frame(cbind("sample_id" = 
                                       as.character(counts$sampleid[counts$sitename 
                                                                    == "Demont Lake"]), counts[counts$sitename 
                                                                                            == "Demont Lake", 7:23]))
demont.age <- as.data.frame(cbind("sample_id" = 
                                 as.character(counts$sampleid[counts$sitename 
                                                              == "Demont Lake"]), age = counts$age[counts$sitename
                                                                                                == "Demont Lake"]))
demont.age[,2] <- as.numeric(demont.age[,2])

demont.roc.sequence <- estimate_roc(
  data_source_community = demont.community, 
  data_source_age = demont.age,
  age_uncertainty = NULL,
  smooth_method = "shep",
  working_units = "MW",
  rand = 1e3,
  dissimilarity_coefficient = "chisq")

demont.seq.peak <-
  detect_peak_points(
    demont.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  demont.seq.peak,
  age_threshold = 11.7e3,
  roc_threshold = 2,
  peaks = TRUE,
  trend = "trend_non_linear"
) 

mean(demont.seq.peak$ROC[demont.seq.peak$Age > 300 & demont.seq.peak$Age < 11700])
frains.community <- as.data.frame(cbind("sample_id" = 
                                       as.character(counts$sampleid[counts$sitename 
                                                                    == "Frains Lake"]), counts[counts$sitename 
                                                                                            == "Frains Lake", 7:23]))
frains.age <- as.data.frame(cbind("sample_id" = 
                                 as.character(counts$sampleid[counts$sitename 
                                                              == "Frains Lake"]), age = counts$age[counts$sitename
                                                                                                == "Frains Lake"]))
frains.age[,2] <- as.numeric(frains.age[,2])

frains.roc.sequence <- estimate_roc(
  data_source_community = frains.community, 
  data_source_age = frains.age,
  age_uncertainty = NULL,
  smooth_method = "shep",
  working_units = "MW",
  rand = 1e3,
  dissimilarity_coefficient = "chisq")

frains.seq.peak <-
  detect_peak_points(
    frains.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  frains.seq.peak,
  age_threshold = 18.9e3,
  roc_threshold = 2,
  peaks = TRUE,
  trend = "trend_non_linear"
) 

green.community <- as.data.frame(cbind("sample_id" = 
                                       as.character(counts$sampleid[counts$sitename 
                                                                    == "Green Lake"]), counts[counts$sitename 
                                                                                            == "Green Lake", 7:23]))
green.age <- as.data.frame(cbind("sample_id" = 
                                 as.character(counts$sampleid[counts$sitename 
                                                              == "Green Lake"]), age = counts$age[counts$sitename
                                                                                                == "Green Lake"]))
green.age[,2] <- as.numeric(green.age[,2])

green.roc.sequence <- estimate_roc(
  data_source_community = green.community, 
  data_source_age = green.age,
  age_uncertainty = NULL,
  smooth_method = "shep",
  working_units = "MW",
  rand = 1e3,
  dissimilarity_coefficient = "chisq")

green.seq.peak <-
  detect_peak_points(
    green.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  green.seq.peak,
  age_threshold = 14.2e3,
  roc_threshold = 2,
  peaks = TRUE,
  trend = "trend_non_linear"
) 

hicks.community <- as.data.frame(cbind("sample_id" = 
                                       as.character(counts$sampleid[counts$sitename 
                                                                    == "Hicks Lake"]), counts[counts$sitename 
                                                                                            == "Hicks Lake", 7:23]))
hicks.age <- as.data.frame(cbind("sample_id" = 
                                 as.character(counts$sampleid[counts$sitename 
                                                              == "Hicks Lake"]), age = counts$age[counts$sitename
                                                                                                == "Hicks Lake"]))
hicks.age[,2] <- as.numeric(hicks.age[,2])

hicks.roc.sequence <- estimate_roc(
  data_source_community = hicks.community, 
  data_source_age = hicks.age,
  age_uncertainty = NULL,
  smooth_method = "shep",
  working_units = "MW",
  rand = 1e3,
  dissimilarity_coefficient = "chisq")

hicks.seq.peak <-
  detect_peak_points(
    hicks.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  hicks.seq.peak,
  age_threshold = 3.2e3,
  roc_threshold = 1,
  peaks = TRUE,
  trend = "trend_non_linear"
) 

horse.community <- as.data.frame(cbind("sample_id" = 
                                       as.character(counts$sampleid[counts$sitename 
                                                                    == "Horseshoe Lake"]), counts[counts$sitename 
                                                                                            == "Horseshoe Lake", 7:23]))
horse.age <- as.data.frame(cbind("sample_id" = 
                                 as.character(counts$sampleid[counts$sitename 
                                                              == "Horseshoe Lake"]), age = counts$age[counts$sitename
                                                                                                == "Horseshoe Lake"]))
horse.age[,2] <- as.numeric(horse.age[,2])

horse.roc.sequence <- estimate_roc(
  data_source_community = horse.community, 
  data_source_age = horse.age,
  age_uncertainty = NULL,
  smooth_method = "shep",
  working_units = "MW",
  rand = 1e3,
  dissimilarity_coefficient = "chisq")

horse.seq.peak <-
  detect_peak_points(
    horse.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  horse.seq.peak,
  age_threshold = 11.8e3,
  roc_threshold = 2,
  peaks = TRUE,
  trend = "trend_non_linear"
) 

huff.community <- as.data.frame(cbind("sample_id" = 
                                       as.character(counts$sampleid[counts$sitename 
                                                                    == "Huffman Lake"]), counts[counts$sitename 
                                                                                            == "Huffman Lake", 7:23]))
huff.age <- as.data.frame(cbind("sample_id" = 
                                 as.character(counts$sampleid[counts$sitename 
                                                              == "Huffman Lake"]), age = counts$age[counts$sitename
                                                                                                == "Huffman Lake"]))
huff.age[,2] <- as.numeric(huff.age[,2])

huff.roc.sequence <- estimate_roc(
  data_source_community = huff.community, 
  data_source_age = huff.age,
  age_uncertainty = NULL,
  smooth_method = "shep",
  working_units = "MW",
  rand = 1e3,
  dissimilarity_coefficient = "chisq")

huff.seq.peak <-
  detect_peak_points(
    huff.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  huff.seq.peak,
  age_threshold = 3e3,
  roc_threshold = 2,
  peaks = TRUE,
  trend = "trend_non_linear"
) 
mean(huff.seq.peak$ROC[huff.seq.peak$Age < 3100])

irwin.community <- as.data.frame(cbind("sample_id" = 
                                       as.character(counts$sampleid[counts$sitename 
                                                                    == "Irwin Smith Bog"]), counts[counts$sitename 
                                                                                            == "Irwin Smith Bog", 7:23]))
irwin.age <- as.data.frame(cbind("sample_id" = 
                                 as.character(counts$sampleid[counts$sitename 
                                                              == "Irwin Smith Bog"]), age = counts$age[counts$sitename
                                                                                                == "Irwin Smith Bog"]))
irwin.age[,2] <- as.numeric(irwin.age[,2])

irwin.roc.sequence <- estimate_roc(
  data_source_community = irwin.community, 
  data_source_age = irwin.age,
  age_uncertainty = NULL,
  smooth_method = "shep",
  working_units = "MW",
  rand = 1e3,
  dissimilarity_coefficient = "chisq")

irwin.seq.peak <-
  detect_peak_points(
    irwin.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  irwin.seq.peak,
  age_threshold = 7.4e3,
  roc_threshold = 2,
  peaks = TRUE,
  trend = "trend_non_linear"
) 

lake27.community <- as.data.frame(cbind("sample_id" = 
                                       as.character(counts$sampleid[counts$sitename 
                                                                    == "Lake 27"]), counts[counts$sitename 
                                                                                            == "Lake 27", 7:23]))
lake27.age <- as.data.frame(cbind("sample_id" = 
                                 as.character(counts$sampleid[counts$sitename 
                                                              == "Lake 27"]), age = counts$age[counts$sitename
                                                                                                == "Lake 27"]))
lake27.age[,2] <- as.numeric(lake27.age[,2])

lake27.roc.sequence <- estimate_roc(
  data_source_community = lake27.community, 
  data_source_age = lake27.age,
  age_uncertainty = NULL,
  smooth_method = "shep",
  working_units = "MW",
  rand = 1e3,
  dissimilarity_coefficient = "chisq")

lake27.seq.peak <-
  detect_peak_points(
    lake27.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  lake27.seq.peak,
  age_threshold = 3.3e3,
  roc_threshold = 2,
  peaks = TRUE,
  trend = "trend_non_linear"
) 

lake16.community <- as.data.frame(cbind("sample_id" = 
                                       as.character(counts$sampleid[counts$sitename 
                                                                    == "Lake Sixteen"]), counts[counts$sitename 
                                                                                            == "Lake Sixteen", 7:23]))
lake16.age <- as.data.frame(cbind("sample_id" = 
                                 as.character(counts$sampleid[counts$sitename 
                                                              == "Lake Sixteen"]), age = counts$age[counts$sitename
                                                                                                == "Lake Sixteen"]))
lake16.age[,2] <- as.numeric(lake16.age[,2])

lake16.roc.sequence <- estimate_roc(
  data_source_community = lake16.community, 
  data_source_age = lake16.age,
  age_uncertainty = NULL,
  smooth_method = "shep",
  working_units = "MW",
  rand = 1e3,
  dissimilarity_coefficient = "chisq")

lake16.seq.peak <-
  detect_peak_points(
    lake16.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  lake16.seq.peak,
  age_threshold = 12.9e3,
  roc_threshold = 2,
  peaks = TRUE,
  trend = "trend_non_linear"
) 

minden.community <- as.data.frame(cbind("sample_id" = 
                                       as.character(counts$sampleid[counts$sitename 
                                                                    == "Minden Bog"]), counts[counts$sitename 
                                                                                            == "Minden Bog", 7:23]))
minden.age <- as.data.frame(cbind("sample_id" = 
                                 as.character(counts$sampleid[counts$sitename 
                                                              == "Minden Bog"]), age = counts$age[counts$sitename
                                                                                                == "Minden Bog"]))
minden.age[,2] <- as.numeric(minden.age[,2])

minden.roc.sequence <- estimate_roc(
  data_source_community = minden.community, 
  data_source_age = minden.age,
  age_uncertainty = NULL,
  smooth_method = "shep",
  working_units = "MW",
  rand = 1e3,
  dissimilarity_coefficient = "chisq")

minden.seq.peak <-
  detect_peak_points(
    minden.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  minden.seq.peak,
  age_threshold = 3.3e3,
  roc_threshold = 1,
  peaks = TRUE,
  trend = "trend_non_linear"
) 

morrison.community <- as.data.frame(cbind("sample_id" = 
                                       as.character(counts$sampleid[counts$sitename 
                                                                    == "Morrison Lake"]), counts[counts$sitename 
                                                                                            == "Morrison Lake", 7:23]))
morrison.age <- as.data.frame(cbind("sample_id" = 
                                 as.character(counts$sampleid[counts$sitename 
                                                              == "Morrison Lake"]), age = counts$age[counts$sitename
                                                                                                == "Morrison Lake"]))
morrison.age[,2] <- as.numeric(morrison.age[,2])

morrison.roc.sequence <- estimate_roc(
  data_source_community = morrison.community, 
  data_source_age = morrison.age,
  age_uncertainty = NULL,
  smooth_method = "shep",
  working_units = "MW",
  rand = 1e3,
  dissimilarity_coefficient = "chisq")

morrison.seq.peak <-
  detect_peak_points(
    morrison.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  morrison.seq.peak,
  age_threshold = 2.8e3,
  roc_threshold = 2,
  peaks = TRUE,
  trend = "trend_non_linear"
) 

obrien.community <- as.data.frame(cbind("sample_id" = 
                                       as.character(counts$sampleid[counts$sitename 
                                                                    == "O'Brien Lake"]), counts[counts$sitename 
                                                                                            == "O'Brien Lake", 7:23]))
obrien.age <- as.data.frame(cbind("sample_id" = 
                                 as.character(counts$sampleid[counts$sitename 
                                                              == "O'Brien Lake"]), age = counts$age[counts$sitename
                                                                                                == "O'Brien Lake"]))
obrien.age[,2] <- as.numeric(obrien.age[,2])

obrien.roc.sequence <- estimate_roc(
  data_source_community = obrien.community, 
  data_source_age = obrien.age,
  age_uncertainty = NULL,
  smooth_method = "shep",
  working_units = "MW",
  rand = 1e3,
  dissimilarity_coefficient = "chisq")

obrien.seq.peak <-
  detect_peak_points(
    obrien.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  obrien.seq.peak,
  age_threshold = 14.0e3,
  roc_threshold = 2,
  peaks = TRUE,
  trend = "trend_non_linear"
) 

otter.community <- as.data.frame(cbind("sample_id" = 
                                       as.character(counts$sampleid[counts$sitename 
                                                                    == "Otter Lake"]), counts[counts$sitename 
                                                                                            == "Otter Lake", 7:23]))
otter.age <- as.data.frame(cbind("sample_id" = 
                                 as.character(counts$sampleid[counts$sitename 
                                                              == "Otter Lake"]), age = counts$age[counts$sitename
                                                                                                == "Otter Lake"]))
otter.age[,2] <- as.numeric(otter.age[,2])

otter.roc.sequence <- estimate_roc(
  data_source_community = otter.community, 
  data_source_age = otter.age,
  age_uncertainty = NULL,
  smooth_method = "shep",
  working_units = "MW",
  rand = 1e3,
  dissimilarity_coefficient = "chisq")

otter.seq.peak <-
  detect_peak_points(
    otter.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  otter.seq.peak,
  age_threshold = 3.8e3,
  roc_threshold = 2,
  peaks = TRUE,
  trend = "trend_non_linear"
) 

vest.community <- as.data.frame(cbind("sample_id" = 
                                       as.character(counts$sampleid[counts$sitename 
                                                                    == "Vestaburg Bog"]), counts[counts$sitename 
                                                                                            == "Vestaburg Bog", 7:23]))
vest.age <- as.data.frame(cbind("sample_id" = 
                                 as.character(counts$sampleid[counts$sitename 
                                                              == "Vestaburg Bog"]), age = counts$age[counts$sitename
                                                                                                == "Vestaburg Bog"]))
vest.age[,2] <- as.numeric(vest.age[,2])

vest.roc.sequence <- estimate_roc(
  data_source_community = vest.community, 
  data_source_age = vest.age,
  age_uncertainty = NULL,
  smooth_method = "shep",
  working_units = "MW",
  rand = 1e3,
  dissimilarity_coefficient = "chisq")

vest.seq.peak <-
  detect_peak_points(
    vest.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  vest.seq.peak,
  age_threshold = 20.2e3,
  roc_threshold = 2,
  peaks = TRUE,
  trend = "trend_non_linear"
) 

winter.community <- as.data.frame(cbind("sample_id" = 
                                       as.character(counts$sampleid[counts$sitename 
                                                                    == "Wintergreen Lake"]), counts[counts$sitename 
                                                                                            == "Wintergreen Lake", 7:23]))
winter.age <- as.data.frame(cbind("sample_id" = 
                                 as.character(counts$sampleid[counts$sitename 
                                                              == "Wintergreen Lake"]), age = counts$age[counts$sitename
                                                                                                == "Wintergreen Lake"]))
winter.age[,2] <- as.numeric(winter.age[,2])

winter.roc.sequence <- estimate_roc(
  data_source_community = winter.community, 
  data_source_age = winter.age,
  age_uncertainty = NULL,
  smooth_method = "shep",
  working_units = "MW",
  rand = 1e3,
  dissimilarity_coefficient = "chisq")

winter.seq.peak <-
  detect_peak_points(
    winter.roc.sequence,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  winter.seq.peak,
  age_threshold = 15.2e3,
  roc_threshold = 2,
  peaks = TRUE,
  trend = "trend_non_linear"
) 

# Create lists/dataframes for plotting RoC data, including means and maxima
roc.site.scores <- c(mean(chip.seq.peak$ROC[chip.seq.peak$Age > 300 & chip.seq.peak$Age < 11700]), mean(cliff.seq.peak$ROC[cliff.seq.peak$Age > 300 & cliff.seq.peak$Age < 11700]), mean(cowden.seq.peak$ROC[cowden.seq.peak$Age > 300 & cowden.seq.peak$Age < 11700]),
                     mean(cub.seq.peak$ROC[cub.seq.peak$Age > 300 & cub.seq.peak$Age < 11700]), mean(demont.seq.peak$ROC[demont.seq.peak$Age > 300 & demont.seq.peak$Age < 11700]), mean(frains.seq.peak$ROC[frains.seq.peak$Age > 300 & frains.seq.peak$Age < 11700]),
                     mean(green.seq.peak$ROC[green.seq.peak$Age > 300]), mean(hicks.seq.peak$ROC[hicks.seq.peak$Age > 300]), mean(horse.seq.peak$ROC[horse.seq.peak$Age > 300 & horse.seq.peak$Age < 11700]),
                     mean(huff.seq.peak$ROC[huff.seq.peak$Age > 300 & huff.seq.peak$Age < 11700]), mean(irwin.seq.peak$ROC[irwin.seq.peak$Age > 300 & irwin.seq.peak$Age < 11700]), mean(lake16.seq.peak$ROC[lake16.seq.peak$Age > 300 & lake16.seq.peak$Age < 11700]),
                     mean(lake27.seq.peak$ROC[lake27.seq.peak$Age > 300 & lake27.seq.peak$Age < 11700]), mean(minden.seq.peak$ROC[minden.seq.peak$Age > 300 & minden.seq.peak$Age < 11700]), mean(morrison.seq.peak$ROC[morrison.seq.peak$Age > 300 & morrison.seq.peak$Age < 11700]),
                     mean(obrien.seq.peak$ROC[obrien.seq.peak$Age > 300 & obrien.seq.peak$Age < 11700]), mean(otter.seq.peak$ROC[otter.seq.peak$Age > 300 & otter.seq.peak$Age < 11700]), mean(vest.seq.peak$ROC[vest.seq.peak$Age > 300 & vest.seq.peak$Age < 11700]),
                     mean(winter.seq.peak$ROC[winter.seq.peak$Age > 300 & winter.seq.peak$Age < 11700]), mean(rise.seq.peak$ROC[rise.seq.peak$Age > 300 & rise.seq.peak$Age < 11700]))

roc.sites <- data.frame("Sitename" = c("Chippewa Bog", "Clifford Lake", "Cowden Lake",
                                       "Cub Lake", "Demont Lake", "Frains Lake",
                                       "Green Lake", "Hicks Lake", "Horseshoe Lake",
                                       "Huffman Lake", "Irwin Smith Bog", "Lake Sixteen",
                                       "Lake 27", "Minden Bog", "Morrison Lake", "O'Brien Lake",
                                       "Otter Lake", "Vestaburg Bog", "Wintergreen Lake", "Sunrise Lake"), 
                        "Latitude" = c(counts$lat[counts$siteid == 360][1], counts$lat[counts$siteid == 28100][1],
                                       counts$lat[counts$siteid == 10239][1], counts$lat[counts$siteid == 518][1],
                                       counts$lat[counts$siteid == 661][1], counts$lat[counts$siteid == 821][1],
                                       counts$lat[counts$siteid == 952][1], counts$lat[counts$siteid == 10242][1],
                                       counts$lat[counts$siteid == 28113][1], counts$lat[counts$siteid == 28059][1],
                                       counts$lat[counts$siteid == 8559][1], counts$lat[counts$siteid == 1541][1],
                                       counts$lat[counts$siteid == 1540][1], counts$lat[counts$siteid == 13582][1],
                                       counts$lat[counts$siteid == 10251][1], counts$lat[counts$siteid == 28101][1],
                                       counts$lat[counts$siteid == 16313][1], counts$lat[counts$siteid == 2814][1],
                                       counts$lat[counts$siteid == 2864][1], 44.032833), "Longitude" = c(counts$long[counts$siteid == 360][1], counts$long[counts$siteid == 28100][1],
                                                                                             counts$long[counts$siteid == 10239][1], counts$long[counts$siteid == 518][1],
                                                                                             counts$long[counts$siteid == 661][1], counts$long[counts$siteid == 821][1],
                                                                                             counts$long[counts$siteid == 952][1], counts$long[counts$siteid == 10242][1],
                                                                                             counts$long[counts$siteid == 28113][1], counts$long[counts$siteid == 28059][1],
                                                                                             counts$long[counts$siteid == 8559][1], counts$long[counts$siteid == 1541][1],
                                                                                             counts$long[counts$siteid == 1540][1], counts$long[counts$siteid == 13582][1],
                                                                                             counts$long[counts$siteid == 10251][1], counts$long[counts$siteid == 28101][1],
                                                                                             counts$long[counts$siteid == 16313][1], counts$long[counts$siteid == 2814][1],
                                                                                             counts$long[counts$siteid == 2864][1], -85.329383))

roc.medians <- c(median(chip.seq.peak$ROC[chip.seq.peak$Age > 300]), median(cliff.seq.peak$ROC[cliff.seq.peak$Age > 300]), median(cowden.seq.peak$ROC[cowden.seq.peak$Age > 300]),
                 median(cub.seq.peak$ROC[cub.seq.peak$Age > 300]), median(demont.seq.peak$ROC[demont.seq.peak$Age > 300]), median(frains.seq.peak$ROC[frains.seq.peak$Age > 300]),
                 median(green.seq.peak$ROC[green.seq.peak$Age > 300]), median(hicks.seq.peak$ROC[hicks.seq.peak$Age > 300]), median(horse.seq.peak$ROC[horse.seq.peak$Age > 300]),
                 median(huff.seq.peak$ROC[huff.seq.peak$Age > 300]), median(irwin.seq.peak$ROC[irwin.seq.peak$Age > 300]), median(lake16.seq.peak$ROC[lake16.seq.peak$Age > 300]),
                 median(lake27.seq.peak$ROC[lake27.seq.peak$Age > 300]), median(minden.seq.peak$ROC[minden.seq.peak$Age > 300]), median(morrison.seq.peak$ROC[morrison.seq.peak$Age > 300]),
                 median(obrien.seq.peak$ROC[obrien.seq.peak$Age > 300]), median(otter.seq.peak$ROC[otter.seq.peak$Age > 300]), median(vest.seq.peak$ROC[vest.seq.peak$Age > 300]),
                 median(winter.seq.peak$ROC[winter.seq.peak$Age > 300]), median(rise.seq.peak$ROC[rise.seq.peak$Age > 300]))

roc.max <- c(max(chip.seq.peak$ROC[chip.seq.peak$Age > 300 & chip.seq.peak$Age < 11700]), max(cliff.seq.peak$ROC[cliff.seq.peak$Age > 300 & cliff.seq.peak$Age < 11700]), max(cowden.seq.peak$ROC[cowden.seq.peak$Age > 300 & cowden.seq.peak$Age < 11700]),
             max(cub.seq.peak$ROC[cub.seq.peak$Age > 300 & cub.seq.peak$Age < 11700]), max(demont.seq.peak$ROC[demont.seq.peak$Age > 300 & demont.seq.peak$Age < 11700]), max(frains.seq.peak$ROC[frains.seq.peak$Age > 300 & frains.seq.peak$Age < 11700]),
             max(green.seq.peak$ROC[green.seq.peak$Age > 300]), max(hicks.seq.peak$ROC[hicks.seq.peak$Age > 300]), max(horse.seq.peak$ROC[horse.seq.peak$Age > 300 & horse.seq.peak$Age < 11700]),
             max(huff.seq.peak$ROC[huff.seq.peak$Age > 300 & huff.seq.peak$Age < 11700]), max(irwin.seq.peak$ROC[irwin.seq.peak$Age > 300 & irwin.seq.peak$Age < 11700]), max(lake16.seq.peak$ROC[lake16.seq.peak$Age > 300 & lake16.seq.peak$Age < 11700]),
             max(lake27.seq.peak$ROC[lake27.seq.peak$Age > 300 & lake27.seq.peak$Age < 11700]), max(minden.seq.peak$ROC[minden.seq.peak$Age > 300 & minden.seq.peak$Age < 11700]), max(morrison.seq.peak$ROC[morrison.seq.peak$Age > 300 & morrison.seq.peak$Age < 11700]),
             max(obrien.seq.peak$ROC[obrien.seq.peak$Age > 300 & obrien.seq.peak$Age < 11700]), max(otter.seq.peak$ROC[otter.seq.peak$Age > 300 & otter.seq.peak$Age < 11700]), max(vest.seq.peak$ROC[vest.seq.peak$Age > 300 & vest.seq.peak$Age < 11700]),
             max(winter.seq.peak$ROC[winter.seq.peak$Age > 300 & winter.seq.peak$Age < 11700]), max(rise.seq.peak$ROC[rise.seq.peak$Age > 300 & rise.seq.peak$Age < 11700]))

lower.mich.ages <- c(diff(range(chip.age$age)), diff(range(cliff.age$age)), diff(range(cowden.age$age)), diff(range(cub.age$age)), diff(range(demont.age$age)), diff(range(frains.age$age)), diff(range(green.age$age)), diff(range(hicks.age$age)), diff(range(horse.age$age)), diff(range(huff.age$age)), diff(range(irwin.age$age)), diff(range(lake16.age$age)), diff(range(lake27.age$age)), diff(range(minden.age$age)), diff(range(morrison.age$age)), diff(range(obrien.age$age)), diff(range(otter.age$age)), diff(range(rise.age$age)), diff(range(vest.age$age)), diff(range(winter.age$age)))
lower.mich.coms <- c(length(rownames(chip.community)), length(rownames(cliff.community)), length(rownames(cowden.community)), length(rownames(cub.community)), length(rownames(demont.community)), length(rownames(frains.community)), length(rownames(green.community)), length(rownames(hicks.community)), length(rownames(horse.community)), length(rownames(huff.community)), length(rownames(irwin.community)), length(rownames(lake16.community)), length(rownames(lake27.community)), length(rownames(minden.community)), length(rownames(morrison.community)), length(rownames(obrien.community)), length(rownames(otter.community)), length(rownames(rise.community)), length(rownames(vest.community)), length(rownames(winter.community)))

roc.sites <- roc.sites %>%
  arrange(Sitename)
colnames(roc.sites)[1] <- "Sitename"
lower.mich.meta <- as.data.frame(cbind(roc.sites, "Age Range" = lower.mich.ages, "RecordLength" = lower.mich.coms, "AverageRoC" = lower.mich.roc$roc.site.scores, "MaxRoC" = lower.mich.roc$roc.max))
write.csv(lower.mich.meta, "Lower Michigan Site Table.csv")

plot(lower.mich.meta$`Record Length`, lower.mich.roc$roc.site.scores)
plot(lower.mich.meta$`Record Length`, lower.mich.roc$roc.max)

ggplot(lower.mich.meta, aes(RecordLength, MaxRoC, colour = Sitename)) + 
  geom_point() + geom_smooth(method = 'gam', se = FALSE, linetype = 2, linewidth = 0.4) +
  labs(y="Max RoC Scores", x="Record Length", fill = "Site Name")# + scale_color_discrete(breaks=c("0-1ka", "1-2ka", "2-3ka", "3-4ka", "4-5ka", "5-6ka", "6-7ka","7-8ka", "8-9ka", "9-10ka", "10-11ka"))

# Add RoC scores by site to metadata
lower.mich.roc <- as.data.frame(cbind(roc.sites, roc.site.scores, roc.max))
write.csv(lower.mich.roc, file='RoC_by_Site5.csv')
mean(lower.mich.roc$roc.site.scores)

# Compile full site metadata into an Excel sheet
site.length = c()
oldest <- c(); youngest <- c(); lat <- c(); long <- c()
siteid.list <- c("360", "518", "661", "821", "952", "1540", "1541", "2864", "2814", "8559", "10239", 
                 "10242", "10251", "13582", "16313", "28059", "28100", "28101", "28113")
counts$age[length(counts$sampleid[counts$siteid == "360"])]
counts$age[counts$siteid == "360"][1]
site.inventory <- data.frame(siteid.list, lat, long, site.length, youngest, oldest)
write.csv(site.inventory, file = "siteinventory.csv")

# --- Below: attached code from separate file generating the Sunrise pollen diagram ---

rise.poll <- readxl::read_excel("RISE_Poll_2.xlsx")
test.poll <- rise.poll[, -c(10:15, 25:26, 124)]
test.sum <- colSums(test.poll[,3:117])
rise.poll23 <- test.poll[, which(test.sum > 20) + 2]
rise.poll23 <- cbind(rise.poll[,1], rise.poll23[,-c(16, 20)])
colnames(rise.poll23)[c(4,8,17,18,20,21)] <- c("Pinus undiff.", "Cupressaceae", "Ambrosia", "Artemisia", "Acer", "Fraxinus")

rise.poll23 <- rise.poll23 %>% 
  dplyr::select(age,'Pinus banksiana/resinosa', 'Pinus strobus', 'Pinus undiff.', 'Tsuga', 'Acer', 'Alnus',
                'Betula', 'Carya', 'Cupressaceae', 'Fagus grandifolia', 'Fraxinus',
                'Ostrya/Carpinus', 'Platanus', 'Quercus', 'Salix', 'Tilia',
                'Ulmus', 'Ambrosia', 'Artemisia', 'Poaceae')

clust2 <- rioja::chclust(dist(sqrt(rise.poll23)),
                         method = "coniss")

plot <- rioja::strat.plot(rise.poll23[,-1], yvar = rise.poll23$age,
                          #                    title = "Sunrise Lake",
                          ylabel = "Calibrated Years BP",
                          y.tks = seq(0, 10000, 2000),
                          ylim = c(0, 11000),
                          xlabel = "Pollen (%)",
                          srt.xlabel = 57,
                          y.rev = TRUE,
                          #clust = clust,
                          plot.poly = TRUE, col.poly = c(rep("darkgreen", 17), rep("gold", 3)),
                          col.bar = c(rep("darkgreen", 17), rep('gold', 3)), ylabPos = 2.7,
                          scale.percent = TRUE, exag = c(rep(FALSE, 5), TRUE, FALSE, TRUE, TRUE, FALSE, rep(TRUE, 3), FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
                          col.exag = "auto", exag.alpha = 0.5, exag.mult = 5)
rioja::addClustZone(plot, clust2, 5, col = "red")

