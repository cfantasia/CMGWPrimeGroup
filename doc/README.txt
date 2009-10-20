========================================
(A) For producing ROOT files ("ntuples"):
========================================

(1) Copy and run the installation script: 
UserCode/CMGWPrimeGroup/macros/setup_Wprime.sh

It will:

o Setup a new working area (Default name: WPrime_work, release:
CMSSW_2_1_19)

o Check out the UserCode/CMGWPrimeGroup package (default version:
V00-00-60)

o Compile the code

o Copy example config file to the top directory.


(2) You can now modify the parameters in configuration file, change the
location of the input file and run.





============================================
(B) For reading back ROOT files ("ntuples"):
============================================

(1) Example #1 (very simple)

ln -s UserCode/CMGWPrimeGroup/root_macros/Make_example.C
root -b Make_example.C

It will run read_wprime.C. The code demonstrates how to access muon
candidates in the Event and prints out some kinematic info.

==================================================================


(2) Example #2 (more involved)

ln -s UserCode/CMGWPrimeGroup/root_macros/Make.C
root -b Make.C

The macro will compile and load the following files:
(a) loadInputFiles.C 

Comment: loads input files. Location of top-level directory with input
files to be specified with string top_level_dir at top of file.

(b) loadCrossSections.C

Comment: must update # of produced events if using a different MC
production. This is the place to modify the cross-section if need to study
systematic effects.

(c) loadCuts.C

Comment: this is the place to define the cuts (NB: but not the actual
values!). The "cuts" should be called with the actual thresholds. 

(d) util.C

Comment: this is the place where the actual values for the cuts are
defined. 

(e) GetDistributionGeneric.C

Comment: Wrapper that calls 
o GetMuonPtDistribution_JetIso.C (option = 1), or
o GetMuonPtDistribution.C (option = 2), or
o GetChargePtDistribution.C. (option = 3)
Make sure to preserve the number and the types of arguments when
modifying/adding new functions! 

(f) GetMuonPtDistribution_JetIso.C

Comment: Makes distributions for muon-pt and # of muons in event with and
without cuts on jet-activity and isolation. 

(g) GetMuonPtDistribution.C

Comment: Makes distributions for muon-pt for various cuts (default
values). To see plots with nice colors, legends, etc. run afterwards
  o plotMuPt.C to see signal + individual bgd contributions
  o plotMuPt2.C to see signal +  total bgd contributions


(h) GetChargePtDistribution.C

Comment: Makes distributions for charge asymmetry for various cuts 
(default values). To see nice plots, run macro plotChargeAsym.C
 
