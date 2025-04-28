# R For Movement Ecologists: Species Distribution Modelling using Acoustic Data
**2025 Ocean Tracking Network Early Career Researcher Workshop, April 2025**

*By Reid Steele and Esteban Salazar*

![alt text](https://github.com/Future-of-Marine-Ecosystems/fome_otn_sdm/otn-dal-logo-whitev2.png "test")
![alt text](https://github.com/Future-of-Marine-Ecosystems/fome_otn_sdm/FOME_Primary_WhiteGradient.png "test")

## Welcome!
Hello, and thanks for your interest in our workshop! This repository contains all the code necessary to build species distribution models (SDMs) using acoustic telememtry, as well as example data from the Nova Scotia Blue Shark tagging program.
It also contains the slide deck we used during the workshop.

## Setup
Before we start, please go to https://ocean-tracking-network.github.io/2025-mrf-ecr-workshop/setup.html and copy and run the SDM workshop requirements code block. This will download all of the packages you'll need to follow along during the workshop. It's also always a good idea to update R and RStudio if you have not done so recently.

## Workshop Scripts
The workflow contained in this repository that we'll be using for the workshop is as follows:

### *1. preprocess.R*
This script loads in the example acoustic data, performs some basic data explorations, and cleans and reformats the data for later use in building SDMs.

### *2. download_rasters.R*
This script downloads environmental data layers from Bio-ORACLE (https://www.bio-oracle.org/index.php), plots them, and reformats them for use in building SDMs.

### *3. sdm_blue.R*
This script runs and evaluates SDMs using three statistical methods (GLM, GAM, and Random Forest), as well as an ensemble model between the three. It also plots projections of all four SDMs.

## Bonus Scripts!
We've also included two bonus scripts for you to try out on your own time in case you're intested:

### *sdm_blue_regional.R*
This script downloads blue shark occurrence data from the Ocean Biodiversity Information System (OBIS, https://obis.org/) and runs an SDM within the study area of the acoustic data we used for the workshop. 
This provides a good example of how to download occurrence data from OBIS, and how it can potentially be used to generate regional SDMs.

### *sdm_blue_global.R*
This script downloads all global records of blue shark occurrence data from OBIS and generates a global SDM. Note that this script takes a **long** time to run (think hours, not minutes). 
This is an example of how to run global SDMs such as those found at AquaMaps (https://www.aquamaps.org/search.php)

## Acknowledgements
Thanks to everyone who made this workshop possible! 
That includes our supervisor Dr. Derek Tittensor, who made most of the slides we used in the workshop, as well as the other members of the Future of Marine Ecosystems Lab (https://www.fomelab.org/) who helped us practice and refine this workshop.
A special thanks as well to the Ocean Tracking Network for giving us the opportunity to put on this workshop, with a special thanks to the OTN data team members who helped us refine the workshop content, identify technical issues, created the setup script, and helped all of you with and coding issues you may have had throughout the workshop.
Finally, thanks to all of you for your interest and attention! This is our first time running a workshop, and we're humbled that you chose to take an hour and half of your day to listen to us blather about SDMs. We hope you got something worthwhile out of it!

