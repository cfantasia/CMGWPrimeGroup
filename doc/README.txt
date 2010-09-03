=========================================================================
(A) Step 1: Producing ROOT-tuples
    (Code runs on RECO and produces custom ROOT-tuples with MC, RECO and
    HLT info) 
=========================================================================

(1) Copy and run the installation script: 
UserCode/CMGWPrimeGroup/macros/setup_Wprime.sh

It will:

o Setup a new working area (Default name: V140, release:
CMSSW_3_6_3)

o Check out the UserCode/CMGWPrimeGroup package (default version:
V00-01-40)

o Compile the code

o Copy example config file to the top directory.


(2) You can now modify the parameters in configuration file, change the
location of the input file and run.



NB The following ustom ROOT-tuples produced with 356 (Spring10 STARTUP) 

- filtered (muon pt>100 GeV) MC bacgkround 
- Wprime->munu MC signal for masses: 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 
  1.5 and 2.0 TeV 
- collision data  for a muon skim with pt > 10 GeV


can be found in

(a) /castor/cern.ch/user/c/cleonido/wprime/v121/ (made with V00-01-21, and
1.32 pb^-1 of real data; use
UserCode/CMGWprimeGroup/config/samples_cross_sections_v4.txt) 

(b) /castor/cern.ch/user/c/cleonido/wprime/v105/ (made with V00-01-05, and
255 nb^-1 of real data; use
UserCode/CMGWprimeGroup/config/samples_cross_sections_v2.txt) 



=========================================================================
(B) Step 2: Analyzing custom ROOT-tuples
    (Code runs on output of previous step; applies analysis cuts and
produces histograms (muon-pt, charge asymmetries, etc)
=========================================================================

(1) Example #1 (very simple)

ln -s UserCode/CMGWPrimeGroup/root_macros/Make_example.C
root -b Make_example.C

It will run read_wprime.C. The code demonstrates how to access muon
candidates in the Event and prints out some kinematic info.

-----------------------------------------------------------------


(2) Example #2 (more involved)

ln -s UserCode/CMGWPrimeGroup/root_macros/Make.C
root -b Make.C

The macro will compile and load the following files 
(directory: UserCode/CMGWPrimeGroup/root_macros):

(a) run.C

The driver that runs the rest of the functions.

"option":
Set to 1 (default) for muon-pt distribution
Set to 2 for charge-asymmetry distribution


(b) loadInputFiles.C 

Comment: loads input files (W' and W, QCD, Z/DY, ttbar, data). Location of
top-level directory with input files to be specified in text file 
top_directory.txt (located in UserCode/CMGWPrimeGroup/config/). The
parameters for all input MC files (description, cross-section, # of events
produced etc), can be found in text file samples_cross_sections.txt (in 
the same directory). 

(c) loadCutsAndThresholds.C (.h)

Comment: this is the place to define the cuts and the actual threshold
values. 

(d) GetDistributionGeneric.C

Comment: Wrapper that calls 
o GetMuonPtDistribution (option = 1), or
o GetChargePtDistribution (option = 2)

Details:

(d1) GetMuonPtDistribution

Comment: Makes distributions for muon-pt for various cuts (default
values). To see plots with nice colors, legends, etc. run afterwards
  o plotMuPt.C to see signal + individual bgd contributions
  o plotMuPt2.C to see signal +  total bgd contributions

[NB: To fit the muon-pt spectrum see separate instructions at
UserCode/CMGWPrimeGroup/fitting/README_mupt_fits.txt]


(d2) GetChargePtDistribution

Comment: Makes distributions for charge asymmetry for various cuts 
(default values). To see nice plots, run macro plotChargeAsym.C
 
