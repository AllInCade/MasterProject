# MastersThesis_Rcode
Repository for R scripts and models for the Master's Thesis

---

# README: Exploring Spline Regression Models in glmmTMB

## Project Overview

This repository contains the R scripts and related materials for me and my colleagues Master's thesis titled " Exploring Spline Based Models in glmmTMB." The main goal of this thesis is to implement and test smooth term functionality in the `glmmTMB` package using the existing machinery from `mgcv`. This work includes assessing the compatibility of `glmmTMB` with `mgcv` model objects, identifying challenges in integration and presentation, and evaluating the effectiveness of spline based model in the 'glmmTMB' framework. The thesis also explores the ease of use of these models and tests them on real data sets, including financial time series data, insurance data, and weather data.

## Repository Structure

1. Project Document: Contains the files from the Overleaf project.
2. Spline models: Contains all R scripts for models developed for this project. Each script is named according to its purpose and data set.
3. Ridge models: Contains the R scripts for the models and methods we implement manually with a ridge penalty.
4. Misc R code: Other relevant scripts of R code for the project.
5. Visuals: Plots, images etc used in the document. May upload more visuals that are not present in the delievered document.

## Getting Started

To use the contents of this repository:
1. Clone or download the repository to your local machine.
2. Install the required R packages listed in `requirements.R`.
3. Run the R scripts in the `/Spline models` and `/Ridge models` directories to replicate the analysis.
4. You may need to install glmmTMB by `remotes::install_github("https://github.com/glmmTMB/glmmTMB/tree/smooth_fixes")`, but
   `ìnstall.packages("glmmTMB", type = "source")` should work and have the latest updates.

## Contributions and Contact

While this repository primarily serves as a part of our Master's thesis at the University of Bergen, contributions, suggestions, or questions are welcome. Please feel free to submit an issue or pull request.

## Acknowledgements

Special thanks to our thesis advisor Hans J. Skaug, colleagues at the University of Bergen, and the open-source community for their invaluable support and guidance.

## License

MIT
