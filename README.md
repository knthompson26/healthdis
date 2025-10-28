# healthdis

This repository holds all the scripts for the manuscript: *Are social support interventions an avenue to reduce disparities in depression symptoms? Applying the health disparity framework in a representative US longitudinal cohort*.

September, 2025

**Katherine N Thompson**, Evelina Akimova, Geyu Zhou, Yeongmi Jeong, Erin C. Dunn, Elisabeth Noland, Moritz Herle, and Robbee Wedow. 

All code created by [Katherine N Thompson](https://scholar.google.co.uk/citations?user=xD4dn1IAAAAJ&hl=en). 

Pre-registration is available on the OSF: [https://osf.io/akzhd/](https://osf.io/zhsxr/)

***

Analyses for this project were conducted in R (Version 4.0.3) on the Rossmann cluster hosted by Purdue Univeristy. Below lists all scripts with a brief explanation. 

***

**Main analyses: Health Disparity Framework** 

1. bootstrap_pgs_ALL.R Function to create Interventional Disparity Effect and Total Effect caluclations for PGS as the exposure: Monte Carlo simulation on a 1000-fold expanded dataset and with 1000 bootstrap samples.
2. bootstrap_race.R. Function to create Interventional Disparity Effect and Total Effect caluclations for race as the exposure: Monte Carlo simulation on a 1000-fold expanded dataset and with 1000 bootstrap samples.
3. bootstrap_ses.R. Function to create Interventional Disparity Effect and Total Effect caluclations for SES as the exposure: Monte Carlo simulation on a 1000-fold expanded dataset and with 1000 bootstrap samples.
4. bootstrap_sex.R. Function to create Interventional Disparity Effect and Total Effect caluclations for sex as the exposure: Monte Carlo simulation on a 1000-fold expanded dataset and with 1000 bootstrap samples.
                
***

**Sensitivity analysis: Item-level** 

1. bootstrap_pgs_ALL_item.R. Function to create Interventional Disparity Effect and Total Effect caluclations for PGS as the exposure and *each item (31 items)* as the mediator: Monte Carlo simulation on a 1000-fold expanded dataset and with 1000 bootstrap samples.
2. bootstrap_race_item.R. Function to create Interventional Disparity Effect and Total Effect caluclations for race as the exposure and *each item (31 items)* as the mediator: Monte Carlo simulation on a 1000-fold expanded dataset and with 1000 bootstrap samples.
3. bootstrap_ses_item.R. Function to create Interventional Disparity Effect and Total Effect caluclations for SES as the exposure and *each item (31 items)* as the mediator: Monte Carlo simulation on a 1000-fold expanded dataset and with 1000 bootstrap samples.
4. bootstrap_sex_item.R. Function to create Interventional Disparity Effect and Total Effect caluclations for sex as the exposure and *each item (31 items)* as the mediator: Monte Carlo simulation on a 1000-fold expanded dataset and with 1000 bootstrap samples.
                
***

**Sensitivity analysis: Within-race** 

1. bootstrap_pgs_ALL_byrace.R. Function to create Interventional Disparity Effect and Total Effect caluclations for PGS as the exposure seperately *within each racial group (White, Black, Asian, Native American)*: Monte Carlo simulation on a 1000-fold expanded dataset and with 1000 bootstrap samples.
2. bootstrap_ses_byrace.R. Function to create Interventional Disparity Effect and Total Effect caluclations for SES as the exposure seperately *within each racial group (White, Black, Asian, Native American)*: Monte Carlo simulation on a 1000-fold expanded dataset and with 1000 bootstrap samples.
3. bootstrap_sex_byrace.R. Function to create Interventional Disparity Effect and Total Effect caluclations for sex as the exposure seperately *within each racial group (White, Black, Asian, Native American)*: Monte Carlo simulation on a 1000-fold expanded dataset and with 1000 bootstrap samples.

***
