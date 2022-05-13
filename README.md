# Anisotropy
Overview of code for analyzing global and local spatiotemporal anisotropy
--------------------------
This repository contains Matlab code for processing SPT trajectory data and for analyzing global and local anisotropy at multiple spatiotemporal scales. 

## Global anisotropy
The global anisotropy is calculated from the angles of all free segments (i.e., free segments from both free trajectories and free-bound mixed trajectories). The code we used to calculate global anisotropy is deposited under the folder "Merge2FrameRate", and they are completely adapted from Anders Sejr Hansen's paper "Guided nuclear exploration increases CTCF target search efficiency" published in Nature Chemical Biology, 2020. For detailed explanations on the calculation steps about global anisotropy, we encourage readers to move forward to "https://gitlab.com/anders.sejr.hansen/anisotropy".

## Local anisotropy
Similar to the definition of free and bound segments in the calculation of local anisotropy, each segment from SPT trajectories is first classified as either bound or free state using 2-state Hidden Markov Model (HMM). The local anisotropy is then calculated from the free segments with defined lengths occured right before or after a bound segments. It should be noted that we merged the quality control step and the HMM classification step from the calculation of global anisotropy to get the state assigned trajectory data before the local anisotropy calculation. Below, we will describe each step we used to calculate local anisotropy in more details.

### Detailed description of each step when calculating local anisotropy

#### Step 1 - Obtain HMM-classified trajectory data
Here, we will use the provided example SPT data of U2OS XLONE-FOXA2-Halo in directory ./ExampleData/MTT_trajectory. The mat files are the output of tracked SPT data using the MTT tracking scripts （"https://gitlab.com/tjian-darzacq-lab/SPT_LocAndTrack"). The tracking data are captured at 2 frame rates: 223Hz and 133 Hz. We follow the global anisotropy code step 1 to 3 to obtain HMM-classified trajectories (see README.md under the folder ./GlobalAnisotropy/). The results of HMM-classified trajectories will be saved in the folder ./ExampleData/HMM_first_QC_data.  
#### Step 2 - Calculate local anisotropy
Run the script 'script_batch_f180_FreeAfBfBoundSegment.m' under the folder ./LocalAnisotropy. The script will call the main function 'FreeSegmentAsymCalculation.m' to first find the maximum allowed free segments before and after one bound segment and then calculate angles of these free segments and their corresponding 'f(180/0)'. For example, if you max_FreSeg_size to be 3, then only the first 3 free segments before or after one bound segments will be considered. The variable 'MaxJump' is a maximum jump length threshold to avoid unreasonable long free segments during the angle calculations, and the number is chosen the same as the Anders's global anisotropy. The calculated results will be saved in the variable 'FinalResults' in a mat file under the folder ./ExampleData/LocalAnisotropy.
#### Step3 - Plot the local anisotropy at different free segment lengths
Run the script 'script_draw_GlobalandFreSegAnisotropy.m' under the folder ./LocalAnisotropy to plot the 'f(180/0)' under different maximum allowed free segments defined in Step 2. Also, the polar histogram of the angle distribution will be plotted. The figures will be saved under the folder ./ExampleData/LocalAnisotropy/Figures.
#### Step4 - Visualize the confinement of the random selected bound-free mixed trajectories
Run the script 'script_RandomDrawTrajectories.m' under the folder ./LocalAnisotropy to randomly select trajectory segments in the following way. The script will follow the same threshold of the segment length that used in Step2 when calculating the angles. Then, either one bound segment followed by the defined number of free segments or the defined number of free segments followed by one bound segment will be randomly selected and visualized. When rendering bound to free segments, blue line indicates the first bound segment and black line indicates the following free segments. When rendering free to bound segments, black line indicates the  free segments, and red line indicates the last bound segment. A red circle with the defined radius centered around the end localization of the bound segments will be overlaid to emphasize the spatial scale of the confinment. The figures will be saved under the folder ./ExampleData/LocalAnisotropy/Figures.

## Acknowledgement
This project makes heavy use of below two code sources:

1) Anisotropy (https://gitlab.com/anders.sejr.hansen/anisotropy)
Hansen, Anders S., Assaf Amitai, Claudia Cattoglio, Robert Tjian, and Xavier Darzacq. "Guided Nuclear Exploration Increases CTCF Target Search Efficiency." Nature Chemical Biology 16, no. 3 (March 2020): 25766. https://doi.org/10.1038/s41589-019-0422-3.

2) vbSPT for HMM-classification of trajectories into *BOUND* and *FREE* segments. 
Please see the full vbSPT paper below for details (see also SourceForge for the latest version https://sourceforge.net/projects/vbspt/). Persson, Fredrik, Martin Lindén, Cecilia Unoson, and Johan Elf. "Extracting Intracellular Diffusive States and Transition Rates from Single-Molecule Tracking Data." Nature Methods 10, no. 3 (March 2013): 265–69. https://doi.org/10.1038/nmeth.2367.
