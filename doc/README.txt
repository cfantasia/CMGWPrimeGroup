=========================================================================
(A) Step 1: Producing PAT-tuples
    (Code runs on GEN-SIM-RECO/AOD and produces custom PAT-tuples 
    with MC, RECO and HLT info) 
=========================================================================

(1) Copy and run the installation script: 
UserCode/CMGWPrimeGroup/macros/setup_Wprime.sh

It will:

o Setup a new working area (Default name: V239, release:
CMSSW_4_1_4)

o Check out the UserCode/CMGWPrimeGroup package (default version:
V00-02-39)

o Check out the default versions of DataFormats/PatCandidates and
PhysicsTools/PatAlgos and V00-05-00 UserCode/SHarper/HEEPAnalyzer 

o Hack the content of PAT Muon collections in order to include the 
dedicated high-pt reconstructors (picky, cocktail, etc)

o Hack the BuildFiles of the HEEP code to make it 41x-compatible

o Compile the code

o Copy example config file for producing PAT-tuple to the top directory.



(2) You can now modify the parameters in configuration file, change the
location of the input file and run.

NB The following PAT-tuples can be found in 
/castor/cern.ch/user/c/cleonido/wprime/V235/Data (a), and
/castor/cern.ch/user/c/cleonido/wprime/V235/MC_MuMET (b)

(a) 
o SingleMu skims for datasets:
/SingleMu/Run2011A-PromptReco-v1/AOD (8.3M events)
/SingleMu/Run2011A-PromptReco-v2/AOD (11.2 M events)
JSON: Cert_160404-162917_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt
67.7 ipb in total

o SingleElectron skims for datasets
/SingleElectron/Run2011A-PromptReco-v1/AOD (1.3M events)
/SingleElectron/Run2011A-PromptReco-v2/AOD (4.2M events)
/SingleElectron/Run2011A-PromptReco-v3/AOD (1.6M events)
JSON: Cert_160404-163369_7TeV_PromptReco_Collisions11_JSON.txt
43.4 ipb in total

NB: 
o v1, v2 can be combined, but suffer from bug in ECAL (see 
https://hypernews.cern.ch/HyperNews/CMS/get/datasets/920.html)
o v3 contains a bugfix, but it is a subset of v1+v2

(b) Spring 2011 MC samples for Mu+MET analysis. 
NB: still missing: QCD-muEnriched-high-pt and W- -> tau samples



(3) NB: Structure of python files for PAT-tuple making:

(a) In directory python (low-level files): 
 o patTuple_common_cfg.py contains common defitions (e.g. addition of
pfMET, event-content for trigger, pileup info,  etc); it is not called
directly from user/high-level cfg file 
 
 o patTuple_mumet_cfg.py and patTuple_elmet_cfg.py: files containing info
specific to mu+MET and el+MET analyses (e.g. aditional event-content to
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

 WPrimeAnalyzer UserCode/CMGWPrimeGroup/bin/analyzeWprimeWgamma_cfg.py
(W[mu]+gamma)



Options:
- enable/disable analyzers for Mu+MET (available), W(mu)+gamma (available),
El+MET (available), WZ (missing), tb (missing)
- specify input collections 
- specify parameters for Mu+MET/El+MET/W+gamma analysis



Code structure:

(a) In directory UserCode/CMGWPrimeGroup/bin/
- WPrimeAnalyzer.cc (dummy wrapper)

(b) In directory UserCode/CMGWPrimeGroup/interface/ (and src/)

- class WPrimeFinder: high-level class/wrapper invoking analyzers for
different channels. Expected analyzer methods that get called by this
class: 
  o beginFile (operations to be done when changing input file, e.g create
new histograms) 
  o endFile (operations to be done when closing input file, e.g. save
histograms, print summary)
  o endAnalysis (e.g. print summary of expected events for all samples)
  o eventLoop (run analysis, ie. selection cuts and fill in histograms)

- class WPrimeUtil: helper class to carry out work that is common across
channels. Examples: 
  o parse description of input files, extract integrated luminosity and
weight 
  o calculate tranverse mass for lepton and MET
  o calculate hadronic MET correction from Z->mumu events
for MC W->mu samples


- class MuMETAnalyzer: Mu+MET analysis (applies selection cuts and creates
kinematic distributions)

- class EleMETAnalyzer: Ele+MET analysis (applies selection cuts and creates
kinematic distributions)

- class WgammaAnalyzer.cc: W(mu)+gamma analysis (applies selection cuts and
creates kinematic distributions)


Location of top-level directory with input files to be specified in text
file  top_directory.txt (located in UserCode/CMGWPrimeGroup/config/). The
parameters for all input MC files (description, cross-section, # of events
produced etc), can be found in text file samples_cross_sections_*.txt (in 
the same directory). 

