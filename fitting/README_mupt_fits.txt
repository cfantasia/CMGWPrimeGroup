Prerequisites:
--------------
For producing the input file to fits, see 
UserCode/CMGWPrimeGroup/doc/README.txt

Must run Step-2, Example #2 with option = 1 set in file
GetDistributionGeneric.C 


=======================================================
Before you start, you may consider making these changes:

- cp gsmear_IDEAL_22x.root gsmear.root for running on IDEAL 22x wprime samples
- cp gsmear_50PB-1_22x.root gsmear.root for running on 50 pb^-1 22x wprime samples
- cp gsmear_START3X_356.root gsmear.root for running on START3x 356 wprime
samples (default)

=======================================================





Description of files used for fits of muon-pt spectrum:
=======================================================

(1) fit_wprime.h, cpp  
Fitting functions for bacgkround & signal muon-pt distributions. The
functions are normalized, and one of their parameters returns the # of
events in the fitting region (typically: p[0]). 

For background: using landau, exponential or relativistic Breit-Wigner
(myRBW). myRBW (ie. Relativistic Breit-Wigner) works better for IDEAL 22x
W' samples, whereas myLandau (i.e. Landau) works better for 50 pb^-1 22x 
W'samples. 

Default background function: myLandau (set in myBgd)

For signal: using smeared_sig, performing a convolution of mySig with
muon-pt resolution. mySig is performing a convolution over the RBW space
(see CMS AN-2009/157 for details).

For the combined (Sig+Bgd) spectrum: using mySigBgd = smeared_sig + myBgd



(2) fitSigBgd.C

High-level macro that 
o reads reference muon-pt histograms for signal and background
sources from input file (default: Wprime_analysis.root) 
o reads muon-pt resolution histograms (default: from gsmear.root)
--> make sure this is a copy of either gsmear_IDEAL_22x.root,
gsmear_50PB-1_22x.root or gsmear_START3x_356.root) 
o sets up output ROOT file storing fit results
o runs fitSigBgd_eventLoop macro that does the real fitting work




(3) fitSigBgd_eventLoop.C

This is where the fit is engineered. Code gets called with
reference histograms (W', W, QCD, etc) from which "running" histograms
for given pseudo-experiment are to be derived. Event-loop
- creates histograms for signal, background, total distributions
- creates function to be used by fit
- calls fit
- extracts fit results and stores in output root file



(4) Results.h, C

Class for storing fit results. Useful when running fits over large number
of pseudo-experiments. 



(5) common_fit.h: contains fit settings for (experiences) users
Fitting range (fXMIN, fXMAX) can be set in file fit_wprime.h



(6) runFit.C

Driver for running the program (i.e. "root runFit.C")
Options for user: # of pseudo-experiments, mass value, whether to fit with
H0 or H1 hypothesis


(7) plot_sigbgd.C

Accesses ROOT file created by fitSigBgd_eventLoop.C, and plots data
distribution, along with fitting function and fit results

