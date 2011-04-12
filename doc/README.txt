=========================================================================
(A) Step 1: Producing PAT-tuples
    (Code runs on GEN-SIM-RECO/AOD and produces custom PAT-tuples 
    with MC, RECO and HLT info) 
=========================================================================

(1) Copy and run the installation script: 
UserCode/CMGWPrimeGroup/macros/setup_Wprime.sh

It will:

o Setup a new working area (Default name: V220, release:
CMSSW_4_1_4)

o Check out the UserCode/CMGWPrimeGroup package (default version:
V00-02-20)

o Check out the default versions of DataFormats/PatCandidates and
PhysicsTools/PatAlgos 

o Hack the content of PAT Muon collections in order to include the 
dedicated high-pt reconstructors (picky, cocktail, etc)

o Compile the code

o Copy example config file for producing PAT-tuple to the top directory.



(2) You can now modify the parameters in configuration file, change the
location of the input file and run.

NB The following Mu+MET PAT-tuples produced with V21x can be found in 
/castor/cern.ch/user/c/cleonido/wprime/V210

(a) Data_Run2011A_160404_161365_19.30ipb.root (/SingleMu/Run2011A-PromptReco-v1/AOD),
processed ~8.3M events with 4_1_4
JSON: Cert_160404-161216_7TeV_PromptReco_Collisions11_JSON_MuonPhys.txt and
Cert_161079-161352_7TeV_PromptReco_Collisions11_JSON_noESpbl_v2.txt

(b) (MC)  ttbar (600 events) produced with 4_1_3 (START311_V2_v1) 
Input file used: /castor/cern.ch/cms/store/relval/CMSSW_4_1_3/RelValTTbar/GEN-SIM-RECO/START311_V2-v1/0038/12763BEE-5A52-E011-8988-003048679048.root 


=========================================================================
(B) Step 2: Analyzing custom PAT-tuples
    (Code runs on output of previous step; applies analysis cuts and
produces histograms (Mt distributions, etc)
=========================================================================

To run:
 WPrimeAnalyzer UserCode/CMGWPrimeGroup/bin/analyzeWprimeMuNu_cfg.py
(mu+MET)

 WPrimeAnalyzer UserCode/CMGWPrimeGroup/bin/analyzeWprimeWgamma_cfg.py
(W[mu]+gamma)



Options:
- enable/disable analyzers for Mu+MET (available), W(mu)+gamma (available),
El+MET (missing), WZ (missing), tb (missing)
- specify input collections 
- specify parameters for Mu+MET/W+gamma analysis (may want to move to
separate cfg file) 



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

- class WgammaAnalyzer.cc: W(mu)+gamma analysis (applies selection cuts and
creates kinematic distributions)


Location of top-level directory with input files to be specified in text
file  top_directory.txt (located in UserCode/CMGWPrimeGroup/config/). The
parameters for all input MC files (description, cross-section, # of events
produced etc), can be found in text file samples_cross_sections.txt (in 
the same directory). 

