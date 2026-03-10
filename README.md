# healthdis

This repository holds all the scripts for the manuscript: *Social connection can mitigate depression risk*.

March, 2026

**Katherine N Thompson**, Shawn Bauldry, Elisabeth S Noland, Evelina T Akimova, Erin C Dunn, Yeongmi Jeong, Geyu Zhou, Moritz Herle, and Robbee Wedow. 

All code created by [Katherine N Thompson](https://scholar.google.co.uk/citations?user=xD4dn1IAAAAJ&hl=en). 

Pre-registration is available on the OSF: [https://osf.io/akzhd/](https://osf.io/zhsxr/)

***

Analyses for this project were conducted in R (Version 4.0.3) on the Rossmann cluster hosted by Purdue Univeristy. 

Below provides the file structure with a brief description of each script. 

***

**Data Prep**

1. `/prep/dataprep.Rmd` All data cleaning and variable creation
2. `/prep/efa.R` Script to conduct exploratory factor analysis for social connection across W2 and W4 datasets 

**Main analyses: Health Disparity Framework** 

1. `pgs.R` Function to create mitigated and total effects for PGS as the exposure: Monte Carlo simulation on a 1000-fold expanded dataset and with 1000 bootstrap samples
2. `sex.R` Function to create mitigated and total effects for Sex as the exposure: Monte Carlo simulation on a 1000-fold expanded dataset and with 1000 bootstrap samples
3. `race.R` Function to create mitigated and total effects for Race as the exposure: Monte Carlo simulation on a 1000-fold expanded dataset and with 1000 bootstrap samples
4. `ses.R` Function to create mitigated and total effects for SES as the exposure: Monte Carlo simulation on a 1000-fold expanded dataset and with 1000 bootstrap samples
                
***

**Robustness checks:** 

*Item-level*

1. `/robust/pgs_item.R` PGS as the exposure and *each item (31 items)* as individual mediators
2. `/robust/sex_item.R` Sex as the exposure and *each item (31 items)* as individual mediators
3. `/robust/race_item.R` Race as the exposure and *each item (31 items)* as individual mediators
4. `/robust/ses_item.R` SES as the exposure and *each item (31 items)* as individual mediators

*By-sex*

1. `/robust/pgs_bysex.R` PGS as the exposure and stratified by sex
2. `/robust/sex_bysex.R` Sex as the exposure and stratified by sex
3. `/robust/race_bysex.R` Race as the exposure and stratified by sex
4. `/robust/ses_bysex.R` SES as the exposure and stratified by sex

*By-race*

1. `/robust/pgs_byrace.R` PGS as the exposure and stratified by race
2. `/robust/sex_byrace.R` Sex as the exposure and stratified by race
3. `/robust/race_byrace.R` Race as the exposure and stratified by race
4. `/robust/ses_byrace.R` SES as the exposure and stratified by race
                
***
