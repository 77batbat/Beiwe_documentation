# Running the Beiwe-Analysis code locally

This documentation will demonstrate a guideline of how to run the Beiwe-Analysis code locally and will also provide introduction for each function.

## Table of Contents
- Language and Environment Setup
- Introduction to the functions involved in Beiwe-Analysis
    - Preprocessing
    - Processing
    - Outputs
    - Utility 

## Language, Packages, and Environment Setup

Some of the Beiwe-Analysis code references **`Python`** or **`C`**-related code for speed, but most of this code is in **`R`**, so we suggest running the code in the **`Rstudio`**. If you already have **`Rstudio`** installed, we suggest updating it to the latest version.

The following steps instruct how to install necessary packages including **Rcpp** that facilitate interfacing **`C`** code in R Packages. 

- **Rcpp requires gcc and gfortran installed**
    - **For macOS users**:
        - Open the terminal console, and check if you already installed Xcode Command Line Tools by typing:
            ```
            $ xcode-select -p
            ```
         
         If you see the following returns then the full Xcode package is already installed
            
                ```
                /Applications/Xcode.app/Contents/Developer
                ```
         
         If not, enter the following to install
            
            ```
            $xcode-select --install
            ```
         
         The following can be used to verify **`gcc`** is installed:
            
            ```
            $gcc --version
            ```
        - The next step is to install **`gfortran`**. Follow the steps in [How to install gfortran on Mac OS X](http://skipperkongen.dk/2012/04/27/how-to-install-gfortran-on-mac-os-x/).
    
    - **For windows users**, ...

- **Install/Load Required R Packages and Set Directories in Rstudio**


    - Open **Rstudio** and run the following code to install/load all the R packages needed for the Beiwe-Analysis:
        ```
        list.of.packages = c(
          "Rcpp",
          "RcppArmadillo",
          "mvtnorm",
          "Matrix",
          "MCMCpack",
          "VGAM",
          "stringr",
          "plyr",
          "pryr",
          "dplyr",
          "tidyr",
          "purrr",
          "tibble",
          "lme4",
          'lmerTest',
          'glmmLasso'
        )

        new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
        if(length(new.packages)) install.packages(new.packages)
        lapply(list.of.packages, require, character.only = TRUE)


        detach("package:plyr", unload=TRUE)

        ```

