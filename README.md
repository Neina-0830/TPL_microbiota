# TPL_microbiota
Code for "Experimental Warming Amplifies Temporal Variance Scaling Across Diverse Soil Microbiota"

## System Requirements

R (code tested on v 4.3.2)

packages: vegan (v 2.6.4), psych (v 2.4.6.26), lme4 (v 3.1.3), MuMIn (v 1.48.4), ropls (v 1.34.0), iCAMP (v 1.8.1), ieggr (4.17), castor (v 1.8.4), ape (v 5.8.1), picante (v 1.8.2), ggplot2 (v 4.0.0), dplyr (v 1.1.4), reshape2 (v 1.4.4), scatterplot3d (v 0.3.44), hrbrthemes (v 0.8.7)


## Installation

R https://www.r-project.org/

rstudio https://rstudio.com/

(time required to install software and packages <2hr)

## Instruction to Use

data/
##this folder contains the data (in .csv format) used for the analysis.


code/
##this folder contains code to replicate results. 

Taylor_law_power_tTTPL_test.R ##Script for testing Type IV Taylor's power law of different microbiota groups.

Taylor_law_power_tTTPL_test_phylum.R  ##Script for testing Type IV Taylor's power law of different microbiota groups at Phylum level.

Taylor_law_power_pTTPL_test.R  ##Script for testing phylogenetic-based Taylor's power law of different microbiota groups.

Taylor_law_power_pTTPL_test_phylum.R  ##Script for testing phylogenetic-based Taylor's power law of different microbiota groups at Phylum level.
 
Taylor_law_power_pTTPL_test_null_model.R  ##Null model analysis using randomized phylogenetic trees with fixed species abundances.

Calculate_MPDik.R  ##Script for calculating MPDik of different microbiota groups.

Calculate_null_MPDik.R  ##Script for calculating MPDik in null model of different microbiota groups.

Fig. 2&S5.R ##Script for ploting figures.
