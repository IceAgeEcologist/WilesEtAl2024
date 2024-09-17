# README: Analytical Workflow Documentation
This repository holds the R code and related input datasets for analyses and figures for a manuscript submitted to Global Ecology and Biogeography under double-blind review.  Analyses for this manuscript used a combination of R and ESRI ArcGIS software.  The R scripts are publicly archived here; the ESRI ArcGIS software employed a graphical user interface and so no code workflow exists.  ESRI analytical steps and functions are described in the manuscript.

The R scripts included here achieve the following functions:
1. Find and download fossil pollen data for Lower Michigan from the Neotoma Paleoecology Database (Neotoma) and add in the locally generated dataset for Sunrise Lake
2. Taxonomically harmonize datasets and prepare standard percentages.
3. Run non-metric multidimensional scaling (NMDS) on the Lower Michigan pollen sites to ordinate the data.
4. Calculate rates of change (RoC) for the pollen data, as an indicator of climate-driven community variability
5. Run GAM and LOESS models for NMDS and distance-to-site data.
4. Produce figures that are not map-based 
    1. Fig. 2:  NMDS ordination biplots
    2. Fig. 3:  Pollen diagram for Sunrise Lake 
    3. Fig. 5B,C:  Rate of Change vs. distance to ecotone
    4. Fig. S3: NMDS stress plot
    5. Fig. S5: NMDS ordinations for individual sites

ArcGIS was used to
1. Spatially interpolate the NMDS scores using empirical Bayesian kriging (EBK)
2. Calculate the spatial distance between sites and the inferred position of the Michigan Tension Zone (MTZ) ecotone
3. Produce map-based figures in the ms. (Figs. 1, 4, 5A, 6, S4)


Note that, at the time of code development, the Sunrise pollen dataset and its age-depth model had not yet been uploaded to Neotoma and so this code uses a local version.  Since code development, the Sunrise pollen data has been added to Neotoma. 
