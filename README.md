# Coral_Light_modelcode
Model code and data files used to perform simulations for study on a corals dependence on light and temperature, part of the study of Kruijt, A. L., et al. 2024 (under review).

Paper abstract:
The latitudinal range of shallow-water tropical corals is controlled by temperature, and presently limited to waters warmer than 16-18 °C yearround. However, even during Cenozoic climates with such temperatures in polar regions, coral reefs are not found beyond >50° latitude. Here, we test the hypothesis that daily available solar radiation limited poleward expansion of coral reefs during warm climates, using a new box model of shallow marine coral calcification. Our results show that calcification rates start to decline beyond 40° and more quickly beyond 50°, suggesting that winter light intensity and day length prohibits further poleward expansion. This implies that fossil coral reef distribution is not a robust proxy for water temperatures and that poleward expansion of reefs is not an expected carbon cycle feedback of climate warming.

Code info:


R code works with:
R version 4.2.3 (2023-03-15 ucrt) -- "Shortstop Beagle"
Copyright (C) 2023 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64 (64-bit)

Files are in r-markdown and r format.

Installation steps:
--> r version: https://cran.r-project.org/bin/windows/base/
--> r markdown: https://rmarkdown.rstudio.com/authoring_quick_tour.html
--> r studio: https://rstudio-education.github.io/hopr/starting.html

Running the code:
To run the code on your own device, dowload the entire repository. To perform the simulations as they were done for the paper, follow the steps in the MASTER_script in the model_code folder. The script reads-in the required temperature data from the data_files folder. Model output is stored in the model_output folder.

(This is not necessary for the simulations but) if you want to download the modern ocean temperature data from NOAA, follow the steps in the ocean_atlas_modernTemp_extraction.R file.

Making the figures:
To re-do the figures as presented in the paper, the plotting scripts in the model_code folder can be run. These read-in the original model output from the simulations presented in the paper, which are stored in the output_for_plotting folder. 
