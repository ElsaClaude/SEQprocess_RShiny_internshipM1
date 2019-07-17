# How to install and use SEQprocess package from R.
-----------------

## SEQprocess
### 1. How to install SEQprocess and its tools
#### 1.1 General information

This manual completes the one provided by the research team which developped SEQprocess in order to give access to the pipelines to beginners in computer science.

SEQprocess Github : https://github.com/omicsCore/SEQprocess  
Supported installation : Linux  
Supported files format : paired-end fastq  

**/!\ need the 3.5 version of R /!\**  
R can be downloaded with a shell command in your terminal or directly on the website of the cran-project. It is higly recommanded to download a 3.5 version of it.
Some packages can be heavy. Set R folder in a place where you have enough memory space.  

In the next parts you will have to install packages on your system. You can do it with shell commands or with Conda installation. If you do not know how to install packages with Conda, be sure to read the "How to install packages with Conda".

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
The devtools package maybe need some extra package to be installed correctly. You may have to add the **openssl-dev**, **libcurl** and the **libxml2-dev** libraries to your system. If you work on a virtual environment you should use Conda to do it. Or if you are on the base environment (the one per default if you do not activate another) you can also use a shell command to have it directly on your system :
  ```
  #Try to search for the last version available
  sudo apt-get install libssl-dev
  sudo apt-get install libcurl4-gnuls-dev
  sudo apt get install libxml2-dev
  ```

Be sure to read every line of the Rstudio console to check the installation of every packages.

Then you have to install external tools that are necessary for the pipelines. For that there are multiple ways. Most of the libraries are available with a Conda installation, others have to be installed manually in the system. If you install it manually you should create a folder dedicated to this.

*The version mentioned is the minimal one.*

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

#### 1.2 Small guide to install packages with Conda

You need to have a version of conda. The fatest way is to install Miniconda. For that, follow the guide provided by Conda : https://docs.conda.io/en/latest/miniconda.html.  
Now you can create a virtual environment dedicated to SEQprocess usage if you want. Here is an example of command you could use :

  ```
  #Creation of a new virtual environment
  conda create --name MY_ENVIRONMENT_NAME

  #Activation of this environment
  #It has to be used in the terminal from which you will launch Rstudio and work on the pipelines
  conda activate MY_ENVIRONMENT_NAME
  ```
You can find more details on the conda documentation site : https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

Once you have installed Conda and activate your environment you can now install packages on it with the specific command. Here is an example :

  ```
  conda install fastqc samtools
  ```
You can enter many packages name in this command. To check what is the right name of a package (or command if you have an error after enter this one) you can check on the Anaconda cloud by entering the name of the library you need : https://anaconda.org/.  

### 2. How to configure SEQprocess pipelines

When you have downloaded SEQprocess through Github you can access the folder by searching in the R libraries folder. There is a config.R file which contains all details about the 6 pipelines of SEQprocess.

Example of path to find config.R file :
  ```
  /net/travail/elclaude/R/x86_64-pc-linux-gnu-library/3.5/SEQprocess/data/config.R
  ```

This script need to be modified in ordered to fit your parameters and your tools and references files paths.
You should pay attention to some important modifications.

In the current Github repository you can find a config.R file which is an example of customization of the files paths etc.

#### 2.1 Fastqc names

Your fastqc files need to have specific names. It must fit the regular expressions determined by the configuration script to correctly find your data or you can modify it so your files can be found. This is the line you need to consider :
  ```
  fq1.idx=".1.fastq$|_R1.fastq$|.1_val_1.fq$|_1.fq$|.1_val_1.fq$|.R1_val_1.fq$"
  fq2.idx=".2.fastq$|_R2.fastq$|.2_val_2.fq$|_2.fq$|.2_val_2.fq$|.R2_val_2.fq$"
  ```
For example if one of your files is named "R1_paired.fq" it will not be recognized by the script unless you add : **|R1_paired.fq$** to the first line.

#### 2.2 Tools paths

Now that you have the external tools installed in your system or on your Conda virtual environment you have to set their path in the configuration file in order to let the pipeline find the various programs.
For example, tools installed with Conda can be found in a path such as :
  ```
  "/home/elsa/soft/miniconda3/envs/SEQprocess_internship/bin"
  ```
In the previous example, **SEQprocess_internship** is the name of the Conda virtual environment.


If this is in your system, it could be in the following path :

  ```
  ""/usr/bin"
  #example
  ""/usr/bin/cutadapt"
  ```

#### 2.2 Reference files

As you can see in the pre-formatted config.R script, SEQprocess needs various sort of references files.
Some can be found on the internet.

### 3. How to use SEQprocess pipelines

donner droit accès fichiers chmod
