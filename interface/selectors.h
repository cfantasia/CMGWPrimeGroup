#ifndef _selectors_h
#define _selectors_h

#include "UserCode/CMGWPrimeGroup/interface/util.h"
#include "PhysicsTools/SelectorUtils/interface/Selector.h"

////Selectors//////////

/// Class derived from Selector that implements some convenience functions
template<class T>
class MinMaxSelector : public Selector<T> {
 public:
  bool cutOnMin(std::string param) {
    bool useMin = param.find("min") != std::string::npos;
    bool useMax = param.find("max") != std::string::npos;
    if (!useMin && !useMax) {
      std::cout << param << " has neither 'min' nor 'max' in its name!" << std::endl;
      exit(1);
    }
    return useMin;
  }
  template<class C>
    void loadFromPset(Pset params, std::string param, bool shouldset = true) {
    C defaultValue = 0;
    if (!this->cutOnMin(param))
      defaultValue = std::numeric_limits<C>::max();
    C val = params.getUntrackedParameter<C>(param, defaultValue);
    if(!params.exists(param)){
      shouldset = false;//If you don't list it, I'll ignore it
      std::cout<<" You didn't specify param "<<param<<", so it will be ignored\n";
    }else{
      std::cout<<" Param: "<<param;
      if(this->cutOnMin(param)) std::cout<<" >= ";
      else                      std::cout<<" <= ";
      std::cout<<val<<std::endl;
    }
    this->push_back(param, val);
    this->set(param, shouldset);
  }
  template<class C>
    void setpassCut(std::string param, C value, pat::strbitset & ret) {
    bool useMin = this->cutOnMin(param);
    if (this->ignoreCut(param) || 
        ( useMin && value >= this->cut(param, C())) ||
        (!useMin && value <= this->cut(param, C())) )
      this->passCut(ret, param);
  }

  pat::strbitset bitmask;
};

/// Selector for electrons (either barrel or endcap) based on params
class ElectronSelectorBase : public MinMaxSelector<pat::Electron> {
public:
  ElectronSelectorBase();
  ElectronSelectorBase(Pset const params);

  virtual bool operator()(const pat::Electron & p, pat::strbitset & ret);
  virtual bool operator()(const pat::Electron & p, const float pu, const reco::Vertex & vtx=reco::Vertex());
  virtual bool operator()(const heep::Ele & p, pat::strbitset & ret);
  virtual bool operator()(const heep::Ele & p, const float pu);
  
  static float pfIso(const pat::Electron & p, const float & pu);
};

/// A wrapper to handle the barrel/endcap split for electrons
class ElectronSelector {
 public:
  ElectronSelector();
  ElectronSelector(Pset pset, std::string selectorName);
  bool operator()(const pat::Electron & p, const float pu=0., const reco::Vertex & vtx=reco::Vertex());
  bool operator()(const heep::Ele & p, const float pu=0.);
  pat::strbitset getBitTemplate();
 private:
  ElectronSelectorBase barrelSelector_;
  ElectronSelectorBase endcapSelector_;
};


/// Selector for muons based on params
class MuonSelector : public MinMaxSelector<TeVMuon> {
public:
  MuonSelector();
  MuonSelector(Pset pset, std::string selectorName);
  virtual bool operator()(const TeVMuon & p, pat::strbitset & ret);
  bool operator()(const TeVMuon & p, const float pu=0, const reco::Vertex & vtx=reco::Vertex());
};


/// Selector for jets based on params
class JetSelector : public MinMaxSelector<pat::Jet> {
public:
  JetSelector();
  JetSelector(Pset pset, std::string selectorName);
  virtual bool operator()(const pat::Jet & p, pat::strbitset & ret);
  virtual bool operator()(const pat::Jet & p);
  
};

class PhotonSelector : public MinMaxSelector<pat::Photon> {
public:
  PhotonSelector();
  PhotonSelector(Pset pset, std::string selectorName);

  virtual bool operator()(const pat::Photon & p, pat::strbitset & ret);
  virtual bool operator()(const pat::Photon & p);
 
};

#endif // _selectors_h
