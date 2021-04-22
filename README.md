# *Spatial Mapping Reveals Adipocyte Heterogeneity as a Determinant of Insulin Response (preliminary title)*   

**Jesper Bäckdahl§, Lovisa Franzén§, Lucas Massier, Qian Li, Jutta Jalkanen, Hui Gao, Alma Andersson, Nayanika Bhalla, Anders Thorell, Mikael Rydén§, Patrik L. Ståhl§, Niklas Mejhert§**

§ These authors contributed equally to the work.


## Description  

This repo contains all R code and smaller tables related to the publication xxx (TBD) by Bäckdahl J & Franzén L, et al. 

In this project, we performed spatial transcriptomics, using the 10x Genomic's Visium platform, on abdominal subcutaneous white adipose tissue (scWAT) to generate data from a total of ten human subjects. The subjects are of a range of different ages and with different body mass index (BMI). Four of the subjects were subjected to euglycemic hyperinsulinemic clamp with samples collected before and after – the data from these samples were used in the "insulin" analyses.  


Publication: *TBD*


![Manuscript figure](/doc/manus_fig_overview.png)


## Content  

* `data`
  * `visium`  
    * Download data from Mendeley Data (DOI: 10.17632/3bs5f8mvbs.1) and place in this folder  
  * `gene_annotation`: two tables with gene annotation information used for gene filtering and ID conversions  
  * `CAGE`: Data containing FANTOM5 CAGE bulk data used for comparison (Rydén et al., 2016, DOI: 10.1016/j.celrep.2016.07.070)  
* `scripts`: All main R scripts used for processing and analysing the Visium data  
* `image_analysis_lm`: R  code, ImageJ macro, and final data output related to the adipocyte size determination using image analysis
* `doc`: Contains pdf with data analysis overview flowchart  


## Contact  

For questions related to the Visium data and related code, please contact Lovisa Franzén (lovisa.franzen@scilifelab.se)  

For questions related to the image analysis and adipocyte size determination, please contact Lucas Massier (lucas.massier@ki.se)  


