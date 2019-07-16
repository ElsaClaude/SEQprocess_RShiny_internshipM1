# How to install and use SEQprocess package from R.
-----------------

## SEQprocess
### How to install SEQprocess and its tools

This manual complete the one provided by the research team which developped SEQprocess in order to give access to the pipelines to beginners in computer science.

SEQprocess Github : https://github.com/omicsCore/SEQprocess  
Supported installation : Linux  
Supported files format : paired-end fastq  

**/!\ need the 3.5 version of R /!\**  
R can be downloaded with a shell command in your terminal or directly on the website of the cran-project. It is higly recommanded to download a 3.5 version of it.
Some packages can be heavy. Set R folder in a place where you have enough memory space.  
Then, open Rstudio and create a script with the following commands :
  ```
  #Change lib path if necessary with : .libPaths()
  #example :
  .libPaths("/net/travail/elclaude/R/x86_64-pc-linux-gnu-library/3.5")

  #source codes
  source("https://bioconductor.org/biocLite.R")

  #GenomeInfoDbData package
  biocLite("GenomeInfoDbData")

  #DelayedArray package
  biocLite("DelayedArray")

  #devtools packages
  install.packages("devtools")

  #call devtools library
  library(devtools)

  #install packages through github
  install_github("omicsCore/SEQprocess")

  #load packages and documentation
  library(SEQprocess) # Loads the package
  library(help="SEQprocess") # Lists package info

  #install libraries
  #some packages sometimes need options to be installed such as : dependencies=TRUE or INSTALL_opts = c('--no-lock')
  install.packages("sequenza")
  install.packages("data.table")
  install.packages("fastqcr")
  install.packages("pander", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
  install.packages("knitr")
  install.packages("gridExtra")
  install.packages("ggplot2")
  install.packages("reshape2")
  ```
Then you have to install external tools that are necessary for the pipelines. For that there are multiple ways. Most of the libraries are available with an Conda installation, others have to be installed manually in the system.  
If you do not know how to install packages with Conda, be sure to read the next part.

*The version mentioned is the minimal one*

Package | Installation
:-: | -:
 ** FASTQC v0.11.5 ** | Conda
 ** Trim Galore v0.4.2 ** | Conda
 ** Cutadapt v1.11** | Manually
 ** BWA v0.7.15 ** | Manually
 ** SAMtools v0.1.18** | Conda
 ** STAR v2.5.2b** | Conda
 ** Tophat2 v2.1.1** | Conda
 ** Bowtie2 v2.2.9** | Conda
 ** Picard v2.17.4** | Conda
 ** GATK v3.7** | Manually
 ** SomaticSniper v1.0.5.0** | Conda
 ** Varscan2 v2.4.3** | Conda
 ** MuSE v1 Orc** | Conda
 ** Variant Effect Predictor v91** | Conda
 ** ANNOVAR ** | Manually
 ** HT-Seq v0.6.1** | Conda
 ** Cufflinks v2.2.1** | Conda

**/!\ MuTect2 seems to be already installed when GATK is downloaded /!\**

#### How to install packages with Conda

You need to have a version of conda. The fatest way is to install Miniconda. For that, follow the guide provided by Conda : https://docs.conda.io/en/latest/miniconda.html
Now you can create a virtual environment dedicated to SEQprocess usage if you want. Here is an example of command you could use :

  ```
  #Creation of a new virtual environment
  conda create --name MY_ENVIRONMENT_NAME

  #activation of this environment
  #It has to be used in the terminal from which you will launch Rstudio and work on the pipelines
  conda activate MY_ENVIRONMENT_NAME
  ```
You can find more details on the conda documentation site : https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

Once you have installed Conda and activate your environment you can now install packages on it with the specific command. Here is an example :

  ```
  conda install fastqc samtools
  ```
You can enter many packages name in this command. To check what is the right name of a package (or command if you have an error after enter this one) you can check on the Anaconda cloud by entering the name of the library you need : https://anaconda.org/  

### How to configure SEQprocess pipelines

When you have downloaded SEQprocess through Github you can access the folder by searching in the R libraries folder. There is a config.R file which contains all details about the 6 pipelines of SEQprocess.

Example of path to find config.R file :
  ```
  /net/travail/elclaude/R/x86_64-pc-linux-gnu-library/3.5/SEQprocess/data/config.R
  ```

This script need to be modified in ordered to fit your parameters and your tools and references files paths.
You should pay attention to some important modifications :
* Your fastqc files need to have specific names. It must fit the regular expressions determined by the configuration script to correctly find your data or you can modify it so your files can be found. This is the line you need to consider :
```
fq1.idx=".1.fastq$|_R1.fastq$|.1_val_1.fq$|_1.fq$|.1_val_1.fq$|.R1_val_1.fq$"
fq2.idx=".2.fastq$|_R2.fastq$|.2_val_2.fq$|_2.fq$|.2_val_2.fq$|.R2_val_2.fq$"
```
For example if one of your files is named "R1_paired.fq" it will not be recognized by the script unless you add : **|R1_paired.fq$** to the first line.
*
