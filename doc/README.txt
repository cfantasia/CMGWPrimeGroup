=========================================================================
(A) Step 1: Producing ROOT-tuples
    (Code runs on RECO and produces custom ROOT-tuples with MC, RECO and
    HLT info) 
=========================================================================

(1) Copy and run the installation script: 
UserCode/CMGWPrimeGroup/macros/setup_Wprime.sh

It will:

o Setup a new working area (Default name: V102, release:
CMSSW_3_5_8)

o Check out the UserCode/CMGWPrimeGroup package (default version:
V00-01-03)

o Compile the code

o Copy example config file to the top directory.


(2) You can now modify the parameters in configuration file, change the
location of the input file and run.



NB: Custom ROOT-tuples produced with 

(a) 21x IDEAL RECO for W' and bacgkround samples with V00-00-50 can be
found at /castor/cern.ch/user/c/cleonido/v55/

(b) 22x 50pb^-1 RECO for W' and W with V00-00-67 can be found at
 /castor/cern.ch/user/c/cleonido/v66_50pb-1_hack/

The "hack" refers to the reversion of this format change:
http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/CMGWPrimeGroup/interface/wprimeEvent.h?r1=1.2&r2=1.3
so that we can combine the mis-aligned root-tuples with the ideal
non-W background in (a) above




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

(pre-a) run.C

The driver that runs the rest of the functions.

Flag #1: "detector_conditions"; 
Set to 1 (IDEAL) or 2 (50 pb^-1 for W and W',IDEAL for top, QCD and Z/DY) 

Flag #2: "option":
Set to 1 (default) for muon-pt distribution
Set to 2 for charge-asymmetry distribution


(a) loadInputFiles.C 

Comment: loads input files (W' and W, QCD, Z/DY, ttbar). Location of
top-level directory with input files to be specified with string
top_level_dir_flag1 (_flag2) at top of file. 

(b) loadCrossSections.C

Comment: must update # of produced events if using a different MC
production. This is the place to modify the cross-section if need to study
systematic effects.

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
 
