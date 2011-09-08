from UserCode.CMGWPrimeGroup.patTuple_common_cfg import *
from UserCode.CMGWPrimeGroup.patTuple_mu_cfg import *
from UserCode.CMGWPrimeGroup.patTuple_met_cfg import *

def mumet_config(process, reportEveryNum=100, maxEvents=-1) :
    common_config(process, reportEveryNum, maxEvents)
    mu_config(process)
    met_config(process)



