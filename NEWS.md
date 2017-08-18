## DLMtool 
The current version of the DLMtool package is available for download from [CRAN](https://CRAN.R-project.org/package=DLMtool).


### DLMtool v4.4.1

#### Minor Changes 

- fixed a typo in V4.4 that affected effort controls
- removed redundant code to speed up sampling of fleet parameters


## Previous Versions

### DLMtool v4.4   

#### Changes 

- can now specify different retention and vulnerability at size curves to evaluate impact of discarding for sub-legal fish
- can now include general discarding across all size/age classes
- The structure of Input control MPs has been simplified. See the userguide for more details.
- The userguide can now be directly accessed from the console with: userguide()
- Additional operating models and other objects can be loaded from the DLMdata package, use: GetMoreData()
- Age or size-specific M has been updated (was in previous version) with plotting functions to visualize M-at-age
- All MPs have been thoroughly tested and modified to improve robustness and avoid crashing the model
- plotting functions have been improved 
- SRA has been converted to RCpp for significant speed gain 

### DLMtool v4.1

#### Major Changes

- Objects names have been shortened and modified for consistency: 
  - Observation has been renamed to Obs 
  - DLM_fease has been renamed to Fease 
  - DLM_output has been renamed Output
  - DLM_input has been renamed Input 
  - DLM_data has been renamed Data 

- Implementation error and a dedicated implementation error object called 'Imp', (like Stock, Fleet and Obs)
   that can account for variability and overages/underages in both effort, catch and size limit advice (discarding
   rate and post release mortality rate are coming soon)

- At least 10 new example operating models from real DLMtool applications in the US and Canada

- Stochastic SRA operating model functions (i.e. build a full operating model accounting for correlations in
   parameters using historical catch trends and some composition data)

- Full plots of operating model conditions

- Stock-synthesis assessment to DLM function for specifying operating models (currently MLE only but adding
   MCMC and hessian options)

- iSCAM assessment to DLM support for specifying operating models (MPD only but adding MCMC and hessian options soon)

- Tracking of historical versus current simulated age composition in catches and population over simulations

- New function 'Replace' for copying parts of an operating model to another operating model (e.g. Robin Hood
  approach,where you may want to borrow say, the fleet characteristics from another operating model)

- Estimation of a new reference point 'Blow' for calculation of conservation-related performance metrics
   (Blow is the spawning biomass where it takes HZN number of mean generation times to reach Bfrac 
   fraction of SSBMSY given zero catches, where the user can specify HZN and Bfrac). 

- Canadian DFO performance plots (DFO_hist, DFO_proj, DFO_plot2)

- Biomass recovery relative to Blow plots (COSEWIC_plot)

- A dedicated custom parameters slot in the operating model object. This is a list where the user can place custom 
   parameters samples (from any distribution / correlation structure you wish) that are named as they appear in 
   the operating model. 

#### Minor Changes 

- Fsd slot in Fleet object has been re-named Esd 

### DLMtool V3.2.3

#### Major Changes
- all objects previously in `DLMdat` have now been added as separate data objects. This means that it is no longer neccessary to unpack data objects at the beginning of an R session.
- parallel processing can be initialized with a new function `setup()`
- `runMSE` and `runMSErobust` functions have a new logical argument `Hist`. When `Hist=TRUE` the model returns the historical simulations only. 
- `CheckConverg` has been deprecated - use `Converge` now to check MSE for convergence
- several optimization functions in the historical simulations have been re-coded in C++ using the Rcpp package (about 10x faster)

#### Minor Changes
- package now uses Roxygen2 for all documentation. 
- exported several previously internal functions that make it easier for users to examine code
- removed LBSPR functions and replaced with dependency on LBSPR package 
- moved `custompars` section in `runMSE.r` and renamed parameters to match OM object
- modified MSE code so that `R0` can be a vector `nsim` long
- `custompars` has been moved in the MSE code so that more parameters can be passed in to the MSE. A message alerts users which valid parameters have been found in `custompars` and any invalid parameters are displayed in a warning message 
- changed variable name `BMSY_B0` to `SSBMSY_SSB0` to avoid confusion
- added projected spawning stock biomass `SSB` and vulnerable biomass `VB` to MSE object
- added `SpAbun` slot to data object for abundance of spawning stock. The `Abun` slot relates to vulnerable biomass
- added `ntrials` and `fracD` control arguments to `runMSE`. Allows user to modify the maximum number of re-samples in the optimization to achieve the sampled depletion. The model will stop if a proprtion greater than `fracD` of the simulated depletion values are not the same as the sampled values after more than `ntrials` re-samples


#### Bug Fixes 
- fixed typo in `getr` function  
- fixed indexing error in Ricker SRR fixed
- fixed bug in catch-at-length and catch-at-age generation code where years were being indexed incorrectly
- fixed bug in `ChooseEffort` function where nyears was being set to 0 in some cases 
- fixed bug in calculation of Lbar
- fixed colors in `Pplot2` for F/FMSY
- fixed calculation of bias in steepness when `hcv` is 0 in the runMSE code
- LFC and LFS in future projections weren't being updated correctly
- fixed F calculations at end of MSE so that `FMSYref` results in exactly F/FSMY = 1 
- `maxF` now also applies to the catch that is taken from the population


### New Additions to Version 3.2.2
A number of additional plotting functions, and a few new MPs have been added in this version.  Also a few minor changes to improve performance and reliability of the model.

For improved stability, especially with large files, the `runMSErobust` function has been changed so that it now uses the `saveRDS` function to write the MSE objects to disk. MSE objects saved with this version of the function need to be loaded with `readRDS`.  

The `plotFun` function can be used to print out all available plotting functions for objects of class `MSE` or `DLM_data`

**MSE Object**

Pplot and Kplot functions have been modified for extra control of various features of the plots.

Two functions, `barplot.MSE` and `boxplot.MSE`, have been added to plot the MSE object. Call them using the generic `barplot` or `boxplot` functions, and see `?barplot.MSE` and `?boxplot.MSE` for information on the arguments for the function.

New functions: `Jplot`, `Splot`, `Cplot`, and `VOIplot` have been added.

**DLM_data Object**

`boxplot.DLM_data` has been added and can be used to plot boxplots of the TACs recommendations from different methods. Call with `boxplot(DLM_data)`

**New MPs**

Three new input control methods, developed by Helena Geremont, have been added to the package: `EtargetLopt`, `L95target`, and `minlenLopt1`.

### New Additions to Version 3.2.1
A new slot (Effort) has been added to the MSE object. This stores the fishing effort for each year, simulation, and MP in the projection years. The addition of the new slot may cause a warning message to be thrown up if an MSE object from a previous version of the DLMtool is loaded.

You can update the old MSE object by adding an empty `Effort` slot:
`MSEobj <- updateMSE(MSEobj)`.

There were some issues with a couple of the input control MPs, which have now been addressed (thanks, Helena, for identifying these). There was also a problem with effort and selectivity being calculated for the input controls, which has also been fixed.

### New Additions to Version 3.2
A number of small but important bugs have been fixed, with special thanks to Liz Brooks, Helena Geromont, and Bill Harford for alerting us to some of these issues.

Quang Huynh has recoded the mean length methods in C++, and they now run much faster, so should pass the time-limit constraint.

A new function (`runMSErobust`) has been added which is a wrapper for the `runMSE` function.  In time, this may replace `runMSE` as the primary function to use when running a MSE.  `runMSErobust` splits large simulations into a series of smaller packets and stitches them together to return an MSE object. This has the benefit of increasing speed and efficiency, particularly for runs with large number of simulations.  The function also checks for errors and re-starts the MSE if the model crashes.  

A set of functions `OM_xl` and `Fease_xl` have been added. These are used to read in operating model and feasibility parameters from an Excel spreadsheet rather than a CSV file.  These are essentially wrappers for the `new` function, but allow you to store all operating model tables in a single spreadsheet rather than a number of CSV files.  This is mainly useful if you are working on multiple species/stocks.

The size limit feature has been updated to include an upper slot limit. See `slotlim` for an example MP.  The slot limit is specified as the last element in the input control vector.  Similar to the lower size limit, all individuals above the slot limit experience no fishing mortality. 

A number of new MPs have been added. There are now 63 output and 22 input control MPs in the DLMtool.

A new function `makePerf` has been added. This function takes an OM object, and returns the same OM object with no process or observation error. This is useful for testing the performance of methods under perfect conditions, to see if they work as expected. And for debugging!

Two new plotting functions have been added: `wormplot` which creates worm plots of the likelihood of meeting biomass targets in future years, and `VOIplot` which is another value of information plot, similar to the `VOI` function, and shows how observation and operating model parameter values affect trends in long-term yield and biomass.

Coming soon: bag limit MPs for recreational fisheries

### Notes from Version 3.1
In order to simulate fisheries that have experienced important shifts in historical length selectivity, this can now be user specified using a graphical user interface (the 'ChooseSelect' function) or by manually editing a series of new slots in the Fleet object (`SelYears`, `AbsSelYears`, `L5Upper`, `L5Lower`, `LFSUpper`, `LFSLower`, `VmaxUpper`, `VmaxLower`).  

Persistent shifts in stock productivity are a particular concern for fishery management. These can now be generated in the toolkit using a new function `SetRecruitCycle` that generates cyclical pattern in recruitment strength. 

Length-based spawning potential ratio (SPR) MPs have been added.  Currently these methods are slow and often don't pass the time constraint. 

Two features have been added to allow MPs to return additional information for future reference. (1) The DLM_data object that MPs operate on now has a miscellaneous slot `Misc`. (2) MPs can now return a list. The first position is the management recommendation (e.g. TAC) the second is information that is stored in the `Misc` slot that can be used by the MP in the next iteration. This can be useful for storing information that is time-consuming to calculate each year. 

A new generic trade-off performance plot `TradePlot` has also been added. 

### A Note on Version 2.11 
Operating model effort is now simulated by a time-series of year vertices and relative magnitude of effort at each vertex. It follows that the slot `Fleet@Fgrad`, and has been replaced by three slots with vectors of equal length: `Fleet@EffYear`, `Fleet@EffUpper` and `Fleet@EffLower`. These effort trajectories can now be specified by a new graphical interface (function `ChooseEffort()`) which uses points to determine the three slots described above.

Operating model fleet selectivity has been robustified to prevent users from specifying length at first capture (`Fleet@L5`) and length at full selection (`Fleet@LFS`) that are unrealistically high. According to our view of reality these now have upper limits of L50 and maximum length, respectively. 

A function `DOM()` has been added that evaluates how often one MP outperforms another across simulations. It is possible that an MP could have higher average performance but perform worse on higher fraction of simulations. The `DOM()` function provides a diagnostic to analyze this.

An additional function `Sub()` has been added which allows users to subset an MSE object according to either (or both) a vector of MPs and simulations. This means you no longer have to rerun everything to provide results for a smaller number of MPs or particular simulations. 

### A Note on Bug Fixes in 2.1.1 and 2.1.2
A bug was found in which length at first capture was being sampled from a uniform distribution `U(LB,UB*2)` rather than `U(LB,UB)`. When depletion could not be simulated by even very high fishery catchabilities an error could occur after more than 10 attempts to find a suitable value of depletion. Length composition simulation in 2.1.1 was not correctly implemented leading to minor biases. 

### A Note on Version 2.1
In response to popular demand, simulation and data are entirely length-based now. It follows that many objects that worked with 2.0 will no longer be compatible. In most cases it is very quick to make files/objects compatible with version 2.1, but nonetheless we apologize if this is frustrating!

The DLMtool package is stochastic, so if you run into problems with the code, please report them (along with a random seed). In the meantime, simply try running it again; the problem may be attributable to a rare combination of sampled parameters.

Be warned that if you abort a parallel process (e.g. `runMSE()`) half-way through, you are in the lap of the gods! It will often be necessary to restart the cluster `sfInit()` or even restart R.

The package is not designed for very short lived stocks (that live for less than 5 years) due to the problems with approximating fine-scale temporal dynamics with an annual model. Technically, you could just divide all your parameters by a sub-year resolution, but the TAC would be set by sub-year and the data would also be available at this fine-scale, which is highly unlikely in a data-limited setting.

### New to Version 2.1
(1) The DLMtool has moved to a length-based simulator (maturity, fisheries selectivity by length)

(2)	The spatial targeting function has been removed for the moment as its implementation was flawed , so could not distribute fishing correctly with respect to both density and the amount of the resource among the two areas.

(3) `Tplot2` adds a different set of trade-offs including long-term and short-term probability of achieving 50% of FMSY yield and average annual variability in yields.
 
(4) Version 2.0 did not include observation error in estimates of current stock abundance and depletion (only biases were simulated). Many thanks to Helena Geromont for spotting this. This has now been corrected.

(5) `DLM_data` objects now have a slot `LHYear` which is a numeric value corresponding with the last historical year. This is needed for some MPs that want to run off only the past data rather than the updated (projected, closed-loop simulation) data. 

(6) Post-MSE, you can now run a Convergence function `CheckConverg()` to see if performance metrics are stable.

(7) The package now contains `CSRA`, a tool for calculating very rough estimates of current depletion and fishing mortality rate from mean catch data.

(8) `getAFC` now can be used for converting length estimates to age estimates through a stochastic growth model.

(9)	The value of information function (`VOI`) contained bugs in version 2.0. This now has been fixed.

(10) Users can now send their own parameter values to the `runMSE` function allowing outputs from stock assessments or correlated parameters (e.g. *K* and age at maturity) values. 

(11) After deliberation, Pope's approximation has been used to account for intra-year mortality (i.e., TACs are taken from biomass at the start of the year subject to half of natural mortality rate). This is probably a reasonable approximation in a data-limited setting: alternative structural assumptions for *M* are eclipsed by uncertainty in *M* itself and other operating model parameters such as selectivity and bias in observation of data such as annual catches. 

(12) The simulation of length composition data was bugged in version 2.0. The variability in length at age was taken from the observation model. Using the perfect information observation model therefore led to no variability in length at age and hence very odd length composition data. This has been solved; now a fixed 10% CV in length-at-age is assumed (normally distributed).

(13) A bug with Delay-Difference MPs has been fixed (`DD` and `DD4010`) in which stochastic TACs were sampled when reps =1. This should just be the mean estimate. The result is that DD is much less variable between years but comes with less contrast in the data. In addition to the much less variable catch recommendations, long-term mean performance of the MP is reduced while medium-term performance has been improved. 

(14) In the move to length-based inputs it is possible to prescribe wild biases for maximum length and length at maturity. In this version these sampled biases are not correlated so it is possible to create simulated data sets where maximum length is lower than length-at-95% maturity and length-at-50% maturity. We put a hard ceiling on this such that length at 95 percent maturity must be below 90 percent of maximum length and length-at-50% maturity must be below 90% of length-at-95% maturity. This isn't great and this will be improved for v2.11. 

(15) The package now works without initiating a cluster `sfInit()`. 

(16) A simple modification to `DCAC` has been added `EDCAC` (Harford and Carruthers, 2015) that better accounts for absolute stock depletion.

(17) Three new slots are available to run MPs on that related to mean length of catches (`ML`), modal length of captures (`Lc`), and the mean length of catches of fish over `Lc` (`Lbar`)

### New to Version 2.0 
(1) Much has changed in the DLMtool terminology to make it more generally applicable. For example, OFL (overfishing limits, FMSY x current biomass), now belongs to a larger class of TACs (Total Allowable Catches). 

(2) There are now just two classes of DLM MPs, DLM\_output (MPs linked to output controls e.g. TACs) and DLM\_input (MPs linked to input controls such as time-area closures, age selectivity and effort). The new DLM\_input function classes have four components, fractional reallocation of spatial effort, fraction of effort in final historical year prescribed in the current year, spatial limits on fishing mortality and a user-defined age-selectivity curve. For example, given an hypothetical stock with 8 age classes a DLM\_input method might return a vector c(0.5,  0.8,  0,1,  0,0,0,0,1,1,1,1). This is interpreted as a 50% reallocation (Allocation = 0.5) of spatial effort, with a total effort that is 80% of historical levels (Effort = 0.8) with a closure in area 1 and full fishing in area 2 (Spatial = c(0,1)) and knife-edge selectivity at age class 5 (Selectivity = c(0,0,0,0,1,1,1,1)) [note that Selectivity has changed in newer versions of the package]. To demonstrate this new feature there are four new input controls, current effort (curE), 75% of current effort (curE75), age selectivity that matches the maturity ogive (matagelim) and a marine reserve in area 1 (area1MR) [note that matagelim has changed to matlenlim in recent versions].

(3) A 'dumb' MP has been added: Mean Catch Depletion (MCD) that simply calculates a TAC based on mean catches and depletion. This is to demonstrate the (theoretically) very high information content of a reliable estimate of current stock depletion. 

(4) A better length composition simulator has been added. Note that this still renews the normal length structure between ages and does not properly simulate the higher mortality rate of larger, faster growing fish (a growth type group simulator is on its way). 

(5) Help documentation has been much improved including complete guides for `Fleet`, `Stock`, `Observation` and `MSE` objects. Eg `class?MSE`.

(6) Minor bugs have been found with the help of Helena Geromont including a problem with update intervals of 1 and low simulated steepness values. 

(7) Reliability is much improved following a full combinatorial test of all Fleet, Stock, Observation objects against all MPs. 

(8) A dedicated Value of information function is now available for MSE objects: `VOI(MSEobject)` which is smarter than the former version which was included in plot(MSE object class). 

(9) Plotting functions have been improved, particularly `Tplot`, `Kplot`, `Pplot` and `plot(DLM_data object class)`.

(10) `SPmod` has been robustified to stop strongly negative surplus production estimates from leading to erratic behavior.

(11) The butterfish stock type now has less variable recruitment and slightly lower natural mortality rate as previous values were rather extreme and led to data generation errors (with natural mortality rate as high as 0.9, butterfish is right at the limit of what can be simulated reasonably with an annual age-structured operating model). 


