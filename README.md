# Running the Beiwe-Analysis code locally

This documentation will demonstrate a guideline of how to run the Beiwe-Analysis code locally and will also provide introduction for each function.

## Table of Contents
- [Environment Setup](#environment-setup) 
- [Introduction to the functions involved in Beiwe-Analysis](#introduction-to-the-functions-involved-in-beiwe-analysis)
    - [Preprocessing](#preprocessing)
    - [Processing](#processing)
    - [Outputs](#outputs)
    - [Utility](#utility) 

## Environment Setup

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
    - Set filepaths for code source, patient data and analysis output in Rstudio. Below is an example to set these filepaths:
    	
	   ```
	   source_filepath    = "/Users/OnnelaLab/Beiwe-Analysis"
	   data_filepath      = "/Users/OnnelaLab/Sample-Data"
	   output_filepath    = "/Users/OnnelaLab/output"
	   root_filepath      = "/Users/OnnelaLab/"
	   ```

## Introduction to the functions involved in Beiwe-Analysis
### Preprocessing

#### surveys_preprocessing
```
function(patient_name, ...){
  ...
          specific_survey_data[,"survey_id"]      = survey_name
          specific_survey_data["date"]            = date
          specific_survey_data[,"patient_name"]   = patient_name
          specific_survey_data[,"timestamp"]      = timestamp
          specific_survey_data[,"date"]           = as.Date(date)
          survey_data[[specific_survey_filepath]] = specific_survey_data
        }
  ...
```

survey preprocessing:
in this preprocessing, we call each patient name by loop and read survey data  from “specific survey filepath”. If there exists information in the survey file, we read the information by such order that: survey name, date, patient name, timestamp, etc. We store this information into several lists and that save then as RDS file.


#### text_preprocessing
```
function(patient_name, ...){
  ...
		textmat[,"timestamp"] = textmat[,"timestamp"] / 1000
		textmat = textmat[,-2]
		textmat[,c("hours","days")] = hours(textmat[,"timestamp"])
    }
  }
  saveRDS(textmat, paste(output_filepath, "/Preprocessed_Data/Individual/",patient_name, "/text_data.rds",sep=""))
}
```

calling the text preprocessing function, we create textmat by adding information from texts file. To be specific, we first convert textmat by dividing it by 1000. Then we created the variable the variables "hours","days" by splitting the variable  "timestamp”. The reason to do this is since timestamp combine the information of hours and days and we want to extract them from this variable. 

#### calls_preprocessing
```
function(patient_name, ...){
 ...
			callmat = rbind(callmat, data = read.csv(paste(calls_filename,"/",call_file,sep=""),header=T))
		callmat[,"timestamp"] = callmat[,"timestamp"] / 1000
		callmat = callmat[,-2]
		callmat[,c("hours","days")] = hours(callmat[,"timestamp"])
    }
  }
  saveRDS(callmat, paste(output_filepath, "/Preprocessed_Data/Individual/",patient_name, "/call_data.rds",sep=""))
}
```
to call the calls preprocessing, we read files from our working directory and basically did the same thing as the text preprocessing. We also want to combine information from each file and split the variable "timestamp” into “hours” and “days” and convert them information and store them in the format of call.


#### powerstate_preprocessing
```
function(patient_name, ...){
	...
		  statemat = do.call(rbind, statemat)
		  statemat[,1] = statemat[,1] / 1000
		  statemat = statemat[,-2]
		  statemat[,c("hours","days")] = hours(statemat[,"timestamp"])
		  saveRDS(statemat, power_state_filename)
		}
	}
}
```

to call the powerstate preprocessing, we read files from our working directory and basically did the same thing as the text preprocessing. We also want to combine information from each file and split the variable "timestamp” into “hours” and “days” and convert them information and store them in the format of statemat and store them in RDS file.


#### accelerometer_preprocessing
```
function(patient_name, minutes, verbose = TRUE, ...){
...
		accmat[,1] = accmat[,1] / 1000
		accmat = accmat[,-2]
		accmat[,c("hours","days")] = hours(accmat[,"timestamp"])
		saveRDS(accmat, patient_data_filename_RDS)
		file.remove(patient_data_filename_TXT)
	  }
  }
}
```

to call the accelerometer_preprocessing, we read files from our working directory and basically did the same thing as the text preprocessing. We also want to combine information from each file and split the variable "timestamp” into “hours” and “days” and convert them information and store them in the format of accmat and store them in RDS file.

}}


### Processing

### Outputs

#### Individual Results
The functions below generate plots for individual patient analysis. All plots of analyses for an individual patient are saved in a folder named as the patient’s id in `output/Results/Individual/`.

-`ContinuousDataCollectionTracks()` 

Plots daily data collection status for **Surveys**, **Screen on/off**, **GPS**, and **Accelerometer**. 

<p align="center"> 
<img src="./example_plots/DataCollectionTracks.png" width="500">
</p>

#### Group Results

The functions below generate plots for group analysis results and save pdf files in `output/Results/Group`.

-`daily_adherence_grid()` 

Plots grid graphs demonstrating the daily adherence status of all patients. The columns of each grid graph are days, and the rows are survey responses, # of missed calls, call duration, total length of texts received, # of texts received, total length of texts sent, circadian routine, # of significant locations visited, max distance from home, distance travelling, time at home, and GPS amount recorded. A blank cell will be displayed on a given day if the data is 0 or was not collected for that category on that day. If the data is available, the greater the value is, the deeper color the cell will show. An example plot for one subject is shown below. 

<p align="center"> 
<img src="./example_plots/adherence.png">
</p>

-` plot_data_quality()` 

Generate 9 individual scatter plots to show daily quality for accelerometer or GPS data for all subjects in the study. When plotting accelerometer data quality, run ` plot_data_quality(stream = "accelerometer", acc_frequency, acc_burst_duration, acc_break_duration,legend=FALSE)`. When plotting GPS data quality, run ` plot_data_quality(stream = "gps",  gps_frequency, gps_burst_duration, gps_break_duration,legend=FALSE)`. 


Among the 9 plots, 4 are plotting **Number of Busts Per Day**, **Average Frequency Per Burst**, **Average Duration per Burst Over Time**, and **Average Duration Between Bursts Over Time** over **Unique Daily Measurements**. Points of each subject’s records are shown in a unique color. One example of Number of Busts Per Day over Unique Daily Measurements with 3 patients’ GPS data is shown below: 

<p align="center"> 
<img src="./example_plots/quality_unique_day_measurement.png" width="500">
</p>



The x-axes of the other 5 plots are **Day**, and the y-axes are **Number of Bursts**, **Average Frequency Per Burst**, **Average Duration per Burst**, **Average Duration per Burst**, and **Overall Coverage**. One example of Number of Number of Busts Per Day over Day with 2 patients’ GPS data is shown below:

<p align="center"> 
<img src="./example_plots/quality_day.png" width="500">
</p>

-`plot_survey_responsiveness()` 

Generates 4 individual scatter plots to show daily survey responsiveness (measured by **Time to First Response** and **Time to Complete After First Response**) for subjects in the study. **Time to First Response** and **Time to Complete After First Response** are plotted over **Unique Daily Measurements** and over **Day**.


-`plot_survey_completion()` 

Generates 2 scatter plots to visualize weekly survey completion. 

-`plot_accelerometer()` 

Generates individual plot of **Daily Accelerometer Data** for each patient. The x-axis for each plot is **Time of day** (in 24 hours format), and the y-axis is **Day** over the entire study period.  COLOR ??????????????

<p align="center"> 
<img src="./example_plots/slat_plot.png" width="500">
</p>

### Utility
