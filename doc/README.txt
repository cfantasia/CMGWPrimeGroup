=========================================================================
(A) Step 1: Producing PAT-tuples
    (Code runs on GEN-SIM-RECO/AOD and produces custom PAT-tuples 
    with MC, RECO and HLT info) 
=========================================================================

(1) Copy and run the installation script: 
UserCode/CMGWPrimeGroup/macros/setup_Wprime.sh

It will:

o Setup a new working area (Default name: V390, release:
CMSSW_4_2_8_patch7)

o Check out the UserCode/CMGWPrimeGroup package (default version:
V00-03-90)

o Check out the recommended versions of DataFormats/PatCandidates,
PhysicsTools/PatAlgos, PhysicsTools/Utilities, RecoLuminosity/LumiDB and
UserCode/SHarper/HEEPAnalyzer

o Hack the content of PAT Muon collections in order to include the 
dedicated high-pt reconstructors (DYT, cocktail, etc)

o Compile the code

o Copy example config file for producing PAT-tuple to the top directory.



(2) You can now modify the parameters in configuration file, change the
location of the input file and run.

*********************************************************************
Notes: documentation on all available PAT-tuples can be found in
doc/data_and_MC_samples.txt 
*********************************************************************


(3) Structure of python files for PAT-tuple making:

(a) In directory python (low-level files): 
 o patTuple_common_cfg.py contains common definitions (e.g. addition of
pfMET, event-content for trigger, pileup info,  etc); it is not called
directly from user/high-level cfg file 
 
 o patTuple_mumet_cfg.py and patTuple_elmet_cfg.py: files containing info
specific to mu+MET and el+MET analyses (e.g. additional event-content to
keep; for the muon channel, the subset of pfParticles that needs to be
stored for the correction of pfMET when using TeV muon reconstructors);
this is where the call to patTuple_common is made. 
 
(b) In directory test (high-level files):
 o patTuple_[X]_[Y].py: high-level files for running muon+MET/electron+MET
analyses on MC/data. Depending on flags, MC info is in(ex)cluded, different
filters are applied on lepton-pt, etc.  
 [X] = elmet or mumet
 [Y] = data or MC_cfg_lowPtSkim or MC_cfg_highPtSkim



=========================================================================
(B) Step 2: Analyzing custom PAT-tuples
    (Code runs on output of previous step; applies analysis cuts and
produces histograms (Mt distributions, etc)
=========================================================================

To run:
 WPrimeAnalyzer UserCode/CMGWPrimeGroup/bin/analyzeWprimeMuNu_cfg.py
(mu+MET)

 WPrimeAnalyzer UserCode/CMGWPrimeGroup/bin/analyzeWprimeElNu_cfg.py
(ele+MET)

 WPrimeAnalyzer UserCode/CMGWPrimeGroup/bin/analyzeWprimeWZ_cfg.py
(WZ)

 WPrimeAnalyzer UserCode/CMGWPrimeGroup/bin/analyzeWprimeWgamma_cfg.py
(W[mu]+gamma)

 WPrimeAnalyzer UserCode/CMGWPrimeGroup/bin/analyzeWprimeHadVZ_cfg.py
(HadVZ)

 WPrimeAnalyzer UserCode/CMGWPrimeGroup/bin/analyzeWprimeHadWZ_cfg.py
(HadWZ)

 WPrimeAnalyzer UserCode/CMGWPrimeGroup/bin/analyzeWprimeTB_cfg.py
(TB)


Options:
- enable/disable analyzers for Mu+MET (available), W(mu)+gamma (available),
El+MET (available), WZ (available), hadV+Z (available), hadV + W
(available), tb (available)
- specify input collections 
- specify parameters for specific analysis



Code structure:

(a) In directory UserCode/CMGWPrimeGroup/bin/
- WPrimeAnalyzer.cc (dummy wrapper)

(b) In directory UserCode/CMGWPrimeGroup/interface/ (and src/)

- class AnalyzerBase: Base class that individual channels inherit from.  
Expected analyzer methods that get called by this
class: 
  o defineHistos
  o fillHistos
  o eventLoop (run analysis, ie. selection cuts and fill in histograms)

- class WPrimeUtil: helper class to carry out work that is common across
channels. Examples: 
  o parse description of input files, extract integrated luminosity and
weight 
  o calculate transverse mass for lepton and MET
  o calculate hadronic MET correction from Z->mumu events
for MC W->mu samples


- class MuMETAnalyzer: Mu+MET analysis (applies selection cuts and creates
kinematic distributions)

- class EleMETAnalyzer: Ele+MET analysis (applies selection cuts and creates
kinematic distributions)

- class WZAnalyzer.cc: W(mu)+Z analysis (applies selection cuts and
creates kinematic distributions)

- class WgammaAnalyzer.cc: W(mu)+gamma analysis (applies selection cuts and
creates kinematic distributions)

- etc

Location of top-level directory with input files to be specified in text
file samples_cross_sections_*txt (located in
UserCode/CMGWPrimeGroup/config/). The parameters for all input MC files
(description, cross-section, # of events produced etc), can be found in the
same file. 

If the filename given for a sample ends in .txt, it will be read in as a list of root files to be analyzed.

A subdirectory can be specified for all root files for a given sample.
