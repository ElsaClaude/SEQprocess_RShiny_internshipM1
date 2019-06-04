# How to use SEQprocess package from R and use it to create an user interface with RShiny package to analyse cancer data.
-----------------

## SEQprocess
### How to install

Github : https://github.com/omicsCore/SEQprocess  
Supported installation : Linux  
**/!\ need at least the 3.5 version of R /!\**  
Open R and create a script with the following commands :
  ```
    source("https://bioconductor.org/biocLite.R")
    biocLite("GenomeInfoDbData")
    biocLite("DelayedArray")

    install.packages("devtools")
    library(devtools)

    install_github("omicsCore/SEQprocess")
  ```
Then install extern tools that are necessary for the pipelines (see list on the github page).  
**/!\ MuTect2 seems to be already install when GATK is downloaded /!\**

04/06 11:58 : reste ANNOVAR , HT-Seq et Cufflinks à installer
