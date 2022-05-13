Global Anisotropy
--------------------------

# Overview of code for data processing and spatiotemporal anisotropy calculations
This folder contains Matlab code for processing SPT trajectory data
and for calculating global anisotropy at multiple spatiotemporal
scales.

## Quick tutorial: going through each step

1. **Step 1**: Obtain SPT data at multiple temporal scales
   1. Here we will use the provided example data for U2OS XLONE-FOXA2-Halo, which is provided in the directory
      ./ExampleData/MTT_trajectory and contains data at 2 frame rates:
      ~133 Hz, ~74 Hz. 
2. **Step 2**: Merge and QC the SPT data from many different single
cells
	1. Open script `Step1_MergeQC_SPT_data_merge2time.m` and click run (~1-2 min).
	2. Use the script to merge data from multiple cells. Dependent
       function: `RemoveAmbigiousTracks.m`
	3. Adjust `ClosestDist` to set the threshold in micrometers for
       when particles are too close and trajectories should be
       aborted.
    4. At the end of this step, the directory `QC_data` should contain
       a single MAT file for each frame rate,
       e.g. `U2OS_XLONE-FOXA2-Halo_74Hz_pooled_QC_CD2.mat`
    5. Step 3-5 will be automatically run in this script.
3. **Step 3** Classify all trajectories using a Hidden-Markov Model (HMM)
   1. This step is using the script 'Step2_Batch_vbSPT_classify_merge2time.m'.
   2. This script uses the HMM vbSPT to classify trajectories into
      *BOUND* and *FREE* segments using a 2-state model. Please see
	  acknowlegdements for a full citation of Persson et al. for vbSPT.
   3. Since vbSPT cannot handle gaps, `Step2_Batch_vbSPT_classify_merge2time.m` goes
         through trajectories with gaps splits them into
         subtrajectories. vbSPT runs in parallel by default. 
   4. At the end, HMM-classified trajectories are saved to
     `HMM_first_QC_data` containing two variables: `CellTracks` is a
     cell array with the xy data. `CellTrackViterbiClass` is the
     HMM-classification, with `1`=*BOUND* and `2`=*FREE*. 
4. **Step 4** Temporally subsample the HMM-classified SPT data
	1. This step is using the script 'Step3_CompileTemporalSubSamplesOfHMM_merge2time.m'.
	2. This script subsamples the data to generate trajectories at
    longer lag times (e.g. 100 Hz --> 50 Hz).
	3. It also carries over the HMM-classification from the faster
       frame rates. Dependent function: `TemporallyReSampleCellTracks.m`
5. **Step 5** Perform full SpatioTemporal analysis of anisotropy
   1. This step is using the script 'Step4_Process_SpatioTemporal_AngleAnalysis_v2_merge2time.m'.
   2. This will do a full analysis of anistropy at multiple
      spatiotemporal scales. The relevant parameters for the analysis
      should be specificied in
      `Process_SpatioTemporal_AngleAnalysis_v2.m`.
   3. The bulk of the analysis is performed in
      `angleFWHM_Amp_HMM_analyzer_v5.m`. Additional dependent functions
      `AngleMatrix_analyzer.m` and `ComputeAmpFWHM.m`.
   4. At the end of the run, the analysis results are saved as MAT
      file: `PA646_SpatioTemporalAnalysis_2022-05-13_17_00_26.mat` 
6. **Step 6** Plot the results
    Depending on how many samples you run, if you run the script on a single sample, e.g., U2OS XLONE-FOXA2-Halo, you can run the script 'Step5_Worked_PLOT_SpatioTemporalAnalysisResults_v5.m' to visualize the results as the figure saved as ".ExampleData\AngularAnalysis\QC_Plots\U2OS XLONE FOXA2-Halo_HMMfirst_Plot1.pdf".
    If you run multiple samples, you can then run the script "Step5_Merge_FinalResults.m" to first merge the results from multiple samples, and then run the script "AfterMerge_Worked_PLOT_SpatioTemporalAnalysisResults_v5.m" to visualize the anisotropy from each sample and their overlays.

## Detailed description of each step
We strongly encourage readers to move forward to "https://gitlab.com/anders.sejr.hansen/anisotropy" for detailed explanation of each step.

#### Compatibility
This code was tested with Matlab 2021b on a Windows 10 and Matlab 2021a on a Ubuntu 20.04 LTS.

## License
These programs are released under the GNU General Public License version 3 or upper (GPLv3+).

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


## Acknowledgements

This project is forked from "https://gitlab.com/anders.sejr.hansen/anisotropy", please cite their work if you use them.
Hansen, Anders S., Assaf Amitai, Claudia Cattoglio, Robert Tjian, and Xavier Darzacq. "Guided Nuclear Exploration Increases CTCF Target Search Efficiency." Nature Chemical Biology 16, no. 3 (March 2020): 25766. https://doi.org/10.1038/s41589-019-0422-3.

This project makes heavy use of vbSPT for HMM-classification of
trajectories into *BOUND* and *FREE* segments. Please see the full
vbSPT paper below for details (see also SourceForge for the latest
version https://sourceforge.net/projects/vbspt/ ):

    Extracting intracellular reaction rates from single molecule tracking data
    Person F,Lindén M, Unoson C, Elf J
    Nature Methods 10, 265–269 (2013). doi:10.1038/nmeth.2367







