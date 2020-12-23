This is the data and code in R associated with the manuscript: Medeiros, L. P., Song, C., and Saavedra, S. (2021). Merging dynamical and structural indicators to measure resilience in multispecies systems. Journal of Animal Ecology.

This repository contains 4 folders: code, data, figs, and results.

Inside code we provide the scripts to reproduce Figures 1-5 (and SI Figures). To generate the theoretical systems one has to run generate_theoretical_systems.R. The theoretical systems are stored in data/random_matrices. Then, one has to run theoretical_systems_eigenvalues_distances.R to compute full and partial recovery and resistance for all theoretical systems. The results from these analyzes are stored in results/random_matrices. Then, one can make Figures 1, 2, and 3 using the codes in fig1.R and fig2_and_3.R.

The experimental data are stored in data/friedman2017_data. To compute full and partial recovery and resistance for the experimental systems, one has to run experimental_systems_eigenvalues_distances.R. The results from these analyzes are stored in results/friedman2017_data. Then, one can make Figures 4 and 5 using the codes in fig4.R and fig5.R. Additional code files can be used to generate the SI Figures.

All codes are commented and use different R functions located in the folder code/functions. The folder figs contains all the figures generated.
