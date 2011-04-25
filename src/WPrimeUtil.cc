#include "UserCode/CMGWPrimeGroup/interface/WPrimeUtil.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

using std::cout; using std::cerr; using std::endl;
using std::string;

WPrimeUtil::WPrimeUtil(const char * out_filename, edm::InputTag genParticles, string cross_sections)
{
  fs = new fwlite::TFileService(out_filename);

  hRecoilPerp = 0;
  hRecoilParalvsVBPt = 0;
  histRecoilParal = NULL;
  lumi_ipb = -1;
  genParticles_ = genParticles;
  sample_cross_sections = cross_sections;
  setupZMETcorrection();

}
WPrimeUtil::~WPrimeUtil()
{
  delete fs; delete [] histRecoilParal;
}

void WPrimeUtil::setupZMETcorrection()
{
  string filename = "ZMET_data.root";
  // open the Z data file with info about recoil
  TFile* File = new TFile(filename.c_str(), "READONLY");
  if(!File || File->IsZombie())
    {
      cerr << " **** Oops! Missing file " << filename << endl;
      cerr << " Exiting... " << endl;
      abort();
    }

  TDirectoryFile * demo=(TDirectoryFile*)File->Get("pflow");
  hRecoilParalvsVBPt = (TH2D*)demo->Get("hMETParalvsVBPt");
  hRecoilParalvsVBPt->SetName("hRecoilParalvsVBPt");
  hRecoilPerp = (TH1D*)demo->Get("hMETPerp");
  hRecoilPerp->SetName("hRecoilPerp");
  
  
  setRecoilProjections();

}

void WPrimeUtil::setRecoilProjections()
{
  assert(hRecoilParalvsVBPt != 0);
  assert(histRecoilParal == 0);
  int N = hRecoilParalvsVBPt->GetXaxis()->GetNbins();
  histRecoilParal = new TH1D *[N];

  for(int bin_no = 1; bin_no <= N; ++bin_no)
    {
      // Get projection in the W pt bin
      histRecoilParal[bin_no-1] = new TH1D(*(hRecoilParalvsVBPt->ProjectionY("_pbinWpt_paral", bin_no, bin_no)));

      if(histRecoilParal[bin_no-1]->Integral()==0) {
	cout << " *** Couldn't correct for recoil for bin_no = " << bin_no
	     << endl;

      }
    }

}

// get input files (to be retrieved from samples_cross_sections.txt)
void WPrimeUtil::getInputFiles(std::vector<wprime::InputFile> & inputFiles)
{
  ifstream topdir_file("UserCode/CMGWPrimeGroup/config/top_directory.txt");
  getline(topdir_file, top_level_dir);
  if(top_level_dir.empty())
    {
      cerr << " *** Failed to load top level directory! " << endl;
      return;
    }

  cout << "\n Reading PAT-tuples from directory " << top_level_dir << endl;

  string txt_file = "UserCode/CMGWPrimeGroup/config/" + sample_cross_sections;
  ifstream in(txt_file.c_str());
  string new_line; wprime::InputFile * new_file = 0;
  while(getline(in, new_line))
    {
      // if DONE, we are done!
      if(new_line == "DONE")break;

      if(new_line.find("samplename = ") != string::npos)
	// new file found! create structure to put in info
	new_file = new wprime::InputFile();
      
      if(!new_line.empty())
	parseLine(new_line, new_file);
      else
	{
	  if(new_file->samplename != "data")
	    new_file->weight = lumi_ipb*(new_file->x_sect)
	      /(new_file->Nprod_evt);
	  else
	    new_file->weight = 1;
	  cout << " Weight to be applied on " << new_file->samplename
	       << " sample = " << new_file->weight << endl << endl;
	  // all info should now be logged in; check!
	  new_file->checkFile();
	  // if we made it here, everything looks good: 
	  // add to vector of input files
	  inputFiles.push_back(*new_file);
	  // release memory
	  delete new_file;
	}
    }
}


// Calculate efficiencies
//------------------------------------------------------------------------
void WPrimeUtil::getEff(float & eff, float & deff,float Num,float Denom)
{
  //------------------------------------------------------------------------
  eff = Num/Denom;
  deff = TMath::Sqrt(eff * (1-eff)/Denom);
}//---------------getEff


// used for the parsing of samples_cross_sections.txt
void WPrimeUtil::parseLine(const string & new_line, wprime::InputFile * in_file)
{
  size_t i = 0;
  i = new_line.find("samplename = ");
  if(i != string::npos)
    {
      in_file->samplename = new_line.substr(13, new_line.length() - 13);
      return;
    }

  i = new_line.find("description = ");
  if(i != string::npos)
    {
      in_file->description = new_line.substr(14, new_line.length() - 14);

      // this is where we extract the integrated luminosity
      // from the description of the data sample; expect something like
      // "blah-blah-blah random description, XYZ pb^-1" ==> XYZ is what we want
      // NB: the data sample must appear first in samples_cross_sections.txt !
      if(lumi_ipb < 0)
	{
	  unsigned int begin = in_file->description.find(",") + 2;
	  unsigned int end = in_file->description.find("pb^-1");
	  string integ_lumi = in_file->description.substr(begin, end-begin);
	  lumi_ipb = atof(integ_lumi.c_str());

	  cout << "\n Will apply weights to MC samples to get distributions for "
	       << lumi_ipb << " pb^-1" << endl << endl;
	  if(lumi_ipb <= 0)
	    {
	      cout << " *** Oops, read in \"" << integ_lumi 
		   << "\" which I failed to parse! Please seek help... " 
		   << endl;
	      abort();
	    }
	}
      return;
    }

  i = new_line.find("pathname = ");
  if(i != string::npos)
    {
      string pathname = top_level_dir + new_line.substr(11, new_line.length() - 11);
      in_file->pathname = pathname;
      
      cout << " Input file: " << in_file->samplename;
      cout << " (" << in_file->description << ") " << endl;

      return;
    }

  i = new_line.find("x-section = ");
  if(i != string::npos)
    {
      in_file->x_sect = atof(new_line.substr(12, new_line.length() - 12).c_str());
      return;
    }

  i = new_line.find("Nprod_evt = ");
  if(i != string::npos)
    {
      in_file->Nprod_evt = atoi(new_line.substr(12, new_line.length() - 12).c_str());
      return;
    }
  
}

// get hadronic MET component (that needs to be corrected 
// if applyMETCorrection=true)from Z data; this will be done according to hadronic 
// activity from Z->mumu reconstructed events
TVector2 WPrimeUtil::getHadronicMET(edm::EventBase const & event)
{
  if(hadronicMETcalculated_)
    return hadronicMETcached;

  assert(event.isRealData() == false);

  int W_index = -1;
  event.getByLabel(genParticles_, genParticles);
  for(size_t i = 0; i != genParticles->size(); ++i) {
    const reco::GenParticle & W_p = (*genParticles)[i];
    if (W_p.pdgId() == 24 && W_p.status() == 3)
      {
	W_index = i;
	break;
      }
  } // loop over genParticles

  assert(W_index >= 0);
  const reco::GenParticle & W_p = (*genParticles)[W_index];

  // correct hadronic MET by taking into account recoil for given W pt
  
  // Find bin corresponding to W pt
  int binWpt = (hRecoilParalvsVBPt->GetXaxis())->FindBin(W_p.pt());

  // protect against outliers: EITHER use the last bin OR return MET w/o correction
  if(binWpt > hRecoilParalvsVBPt->GetXaxis()->GetNbins())
    binWpt = hRecoilParalvsVBPt->GetXaxis()->GetNbins();

  // Shoot random MET (parallel) from the Z histograms
  double dataSampledMETParal= histRecoilParal[binWpt-1]->GetRandom();
  // Shoot perpendicular component
  double dataSampledMETPerp = hRecoilPerp->GetRandom();
  
  // Rotate back from system in which the boson is in the x axis
  double cosW = W_p.px()/W_p.pt();
  double sinW = W_p.py()/W_p.pt();
  
  double dataSampledMEx = 
    cosW*dataSampledMETParal - sinW*dataSampledMETPerp;
  double dataSampledMEy = 
    sinW*dataSampledMETParal + cosW*dataSampledMETPerp;
  
  hadronicMETcached.Set(dataSampledMEx, dataSampledMEy);
  setHadronicMETCalculated(true);
  return hadronicMETcached;    

}
