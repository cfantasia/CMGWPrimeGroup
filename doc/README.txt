=========================================================================
(A) Step 1: Producing PAT-tuples
    (Code runs on GEN-SIM-RECO/AOD and produces custom PAT-tuples 
    with MC, RECO and HLT info) 
=========================================================================

(1) Copy and run the installation script: 
UserCode/CMGWPrimeGroup/macros/setup_Wprime.sh

It will:

o Setup a new working area (Default name: V280, release:
CMSSW_4_2_3)

o Check out the UserCode/CMGWPrimeGroup package (default version:
V00-02-80)

o Check out the default versions of DataFormats/PatCandidates and
PhysicsTools/PatAlgos and V00-05-00 UserCode/SHarper/HEEPAnalyzer 

o Hack the content of PAT Muon collections in order to include the 
dedicated high-pt reconstructors (DYT, cocktail, etc)

o Hack the BuildFiles of the HEEP code to make it 42x-compatible

o Compile the code

o Copy example config file for producing PAT-tuple to the top directory.



(2) You can now modify the parameters in configuration file, change the
location of the input file and run.

Notes:
The following 42x-compatible PAT-tuples are available in:
/castor/cern.ch/user/c/cleonido/wprime/V270/Data

The directory contains skims (lepton pt > 25 GeV) for SingleMu and
SingleElectron datasets. 

(a) May10ReRecoV1: pat-tuples created (with
JSON/json_160404-165620_DCSonly.txt)   from
o /SingleMu/Run2011A-May10ReReco-v1/AOD (21.6 M events, or 221.43 ipb with
MuonPhys JSON)
o /SingleElectron/Run2011A-May10ReReco-v1/AOD (9.1 M events, or 204.16 ipb
with golden JSON)


(b) PromptRecoV4_SingleMuon and PromptRecoV4_SingleElectron: pat-tuples
created (with JSON/json_160404-167151_DCSonly.txt) from 
/SingleMu/Run2011A-PromptReco-v4/AOD (28.0M events)
/SingleElectron/Run2011A-PromptReco-v4/AOD (16.1M events)

which are futher grouped into these directories by run range:
(b1) PromptRecoV4_SingleMuon/165071_166763
(b2) PromptRecoV4_SingleMuon/166764_167043
(c1) PromptRecoV4_SingleElectron/165071_166763
(c2) PromptRecoV4_SingleElectron/166764_167043

o For Mu+MET analysis:
(a)+(b1)+(b2) correspond to 825.73 ipb (if used with
JSON/json_160404-167151_DCSonly.txt), or 756 ipb (if used with
JSON/Cert_160404-166861_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt) 

For El+MET analysis:
(a)+(c1)+(c2) correspond to 817.92 ipb (if used with
JSON/json_160404-167151_DCSonly.txt), or 715 ipb (if used with 
JSON/Cert_160404-166861_7TeV_PromptReco_Collisions11_JSON.txt) 

NB: all these files must be transferred into directory specified in
config/top_directory.txt (which should be edited and reflect the actual
location!)



(3) Structure of python files for PAT-tuple making:

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

 WPrimeAnalyzer UserCode/CMGWPrimeGroup/bin/analyzeWprimeWZ_cfg.py
(W[mu]Z)

 WPrimeAnalyzer UserCode/CMGWPrimeGroup/bin/analyzeWprimeWgamma_cfg.py
(W[mu]+gamma)



Options:
- enable/disable analyzers for Mu+MET (available), W(mu)+gamma (available),
El+MET (available), WZ (available), tb (missing)
- specify input collections 
- specify parameters for Mu+MET/El+MET/WZ/W+gamma analysis



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

- class WZAnalyzer.cc: W(mu)+Z analysis (applies selection cuts and
creates kinematic distributions)

- class WgammaAnalyzer.cc: W(mu)+gamma analysis (applies selection cuts and
creates kinematic distributions)


Location of top-level directory with input files to be specified in text
file  top_directory.txt (located in UserCode/CMGWPrimeGroup/config/). The
parameters for all input MC files (description, cross-section, # of events
produced etc), can be found in text file samples_cross_sections_*.txt (in 
the same directory). 

