# MHWdetection

This GitHub repository contains all of the code used in the _Detecting marine heatwaves with sub-optimal data_ paper.   
It is set up using the package `workflowr` to ensure openness and reproducibility standards.

Please follow this link to the `workflowr` documentation website: https://robwschlegel.github.io/MHWdetection/. This site provides a more interactive and annotated walkthrough of the data analysis pipeline for the paper.

Should one want to view the code in it's raw state the following scripts will be of highest interest:
  
* `code/workflow.R`: Contains the step-by-step process to perform the full analysis  
* `code/figures.R`: Contains the code used to create the figures in the manuscript and supplementary material
* `code/functions.R`: Contains the functions written to perform the analyses in `code/workflow.R` and create most of the figures in `code/figures.R`
