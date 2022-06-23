# Spectral Super Resolution Analyzer

Data analysis software that takes localisations (usually from Super Resolution Microscopy) and creates pairs from them.

It does this based on the assumption that a grating was fixed in front of the camera sensor and thus each feature produced multiple localisations, with a distance based on the wavelength of the original signal.

An ImageJ (JAVA-based) plugin is available, which is recommended for general use. MATLAB scripts as provided as well, with mainly identical functionality.

Please cite our work when using this software: Martens et al., "Enabling spectrally resolved single-molecule localization microscopy at high emitter densities" (2022) (DOI available soon)

# Installation
One can get the Plugin by using the download link below, or by subscribing to the Hohlbein lab update site. You can subscribe by going to Help > Update... > Manage Update Sites and checking the Hohlbein lab site. It will automatically be downloaded and updated.

This plugin requires two files: sSMLMA.jar and jblas.jar:  
The latest version of sSMLA.jar (latest = 0.11.1) can be found [Here](https://github.com/HohlbeinLab/sSMLMAnalyzer/target)  
The latest version of jblas.jar can be found [Here](https://repo1.maven.org/maven2/org/jblas/jblas/1.2.5/jblas-1.2.5.jar)


Once downloaded it can be installed by launching ImageJ > Plugins > Install... Selecting the downloaded jar and restarting ImageJ (for each jar file).  
The plugin should then show up in the Plugins menu under "Spectral Analyzer".

# Usage

This plugin takes as an input one csv file. The ones generated by ThunderSTORM will work without any changes.  
This csv file should have a header and can have either a tab, semicolon or comma as seperator.  
The required columns are: frame, x, y, intensity. 

The ImageJ plugin features a module that tries to determine the angle of the grating, as well as the wavelength range used to detect pairs. These values can also be determined manually by following the instructions in AngleDistance.md.

The settings and their effects are as follows:  
The start and end settings will calculate any value that is unset(=0).
* Angles start and end - The range between which the angle must be (rad)
* Distance start and end - The range between the distance between features must me
* Number of Orders - The maximum number of orders to search for
  
* Restrict Delta Z - A pair is not allowed to have more than delta z difference  
* max delta z [nm] - The maximum delta z.  
* Intensity Order Required - Require that each next order has less intensity
* Ratio between orders - The amount by which the next order has less intensity


* Flip angle/Mirror Angle - The angle could be calculated incorrectly due to quirks with FFT. These settings allow you to manipulate which angle is found
* (Deep) Search for Angles - Search for the first (or best with the most pairs) option that gives good results. Select these options if you are not getting good results.


* Remove lone points - Select this options if some sparse noise appears. It will filter out pairs that do not fulfil the next settings
* Required Neighbours - This will require points to have at least N neighbours within the Required Distance
* Required Distance - This sets the distance within which the N neighbours must be found
  

* Visualise Results - Set this to show a variety of graphs and histograms displaying the results
* Visualise with ZOLA - Set this to include visualisation with ZOLA. The values are prefilled but you need to click OK.
* Histogram binwidth - Sets the width of the distance histograms bins. This value is calibrated for a sample in nm.
* LUT - This shows all LUT's available for your ImageJ and allows you to select one
* LUT start and end - Allows you to display a custom range (same unit as distance) for the LUT visualization



# Running from a Macro

This plugin can also be run from a macro.
An example: run("Analyze Pairs", "csv_in='F:\ThesisData\output\combined_drift.csv' angle_start=-0.094 angle_end=0.22 distance_start=1500 distance_end=2200 visualisation=true")
All keywords must be provided in the format: keyword=value, seperated by spaces.
The keywords available are:

* csv_in - Filepath to the csv input file. Ensure there are no ' in the path
* csv_out - Directory Path to the folder where the output is written to
* angle_start - lower angle boundary (rad)
* angle_end - upper angle boundary (rad)
* distance_start - lower distance boundary
* distance_end - upper distance boundary
* order_number - The maximum number of orders to search for  
* check_z - Require that a pair has a maximum delta z
* check_z_margin - The maximum delta z a pair can have (nm)
* check_order_intensity - Require that each next order has less intensity
* check_order_ratio - The amount by which the next order has less intensity
* angle_flip - Flip the Angle
* angle_mirror - Mirror the Angle
* angle_search - Search for the first best result
* angle_deep_search - Search all possibilities for the best results
* lone_pair_remove - Select this options if some sparse noise appears. It will filter out pairs that do not fulfil the next settings
* lone_pair_neighbours - This will require points to have at least N neighbours within the Required Distance
* lone_pair_distance -  This sets the distance within which the N neighbours must be found
* visualisation - Set this to show a variety of graphs and histograms displaying the results
* visualisationZOLA - Set this to visualise using ZOLA. This does not properly render the 3D version due to issues with ZOLA
* hist_binwidth - Sets the width of the distance histograms bins. This value is calibrated for a sample in nm.
* LUT - The LUT to select
* LUT_start - LUT start (same unit as distance)
* LUT_end - LUT end (same unit as distance)

# JBLAS

The Linear Algebra Plugin is large because it has different libraries depending on the platform.  
I have not tested this fully, but in theory ImageJ on linux should be the fastest version.  
If JBLAS ever gets updated with improved functionality I will try and integrate this into this plugin.

Developement can be tracked [here](https://github.com/jblas-project/jblas).
