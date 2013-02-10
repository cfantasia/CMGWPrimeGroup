#! /usr/bin/env python
import sys,os

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-M", "--masses", dest="Masses", default=-1,
                  help="What masses to run on", metavar="MASSES")
parser.add_option("-o", "--observed", action='store_true',
                  dest="observed", default=False,
                  help="Do Observed Limit")
parser.add_option("-p", "--plot", action='store_true',
                  dest="plot", default=False,
                  help="Print Plots")
parser.add_option("-c", "--condor", action='store_true',
                  dest="condor", default=False,
                  help="Submit Via Condor")
parser.add_option("-t", "--ntoys",
                  dest="ntoys", default="10",
                  help="Number of Toys")
parser.add_option("-l", "--LtShift",
                  dest="LtShift", default="0",
                  help="Lt Shift")
parser.add_option("-w", "--WindShift",
                  dest="WindShift", default="0",
                  help="Wind Shift")
parser.add_option("-s", "--setup", metavar="DIR",
                  dest="setup", default="",
                  help="Directory to run in")

(options, args) = parser.parse_args()
print options

def run(Mass):
    NToys=options.ntoys
    Seed="-1"

    cardFile="card_"+Mode+Suffix+"_M"+Mass+".txt"
    cardCoreFile="card_Core_"+Mode+Suffix+"_M"+Mass+".txt"
    
    if not os.path.exists(cardFile):
        print cardFile+" does not exist, creating."
        os.system("root -b -l -q 'makeLimitCard.C+(\"../../../WprimeWZ.root\", \""+cardCoreFile+"\", "+Mass+", "+options.LtShift+", "+options.WindShift+")'")
        os.system("cat card_Head.txt "+cardCoreFile+" card_Tail.txt > "+cardFile)
        os.system("rm "+cardCoreFile)
        
    CombineCMD="combine -H "+HintAlgo+" -M "+Algo+" -s "+Seed+" -n "+Mode+Suffix+" -m "+Mass+" --rMax 0.1 "+cardFile  #MarkovChainMC Line
    #CombineCMD="combine -H "+HintAlgo+" -M "+Algo+" -s "+Seed+" -n "+Mode+Suffix+" -m "+Mass+" "+cardFile+" --tries 200 --iteration 30000 --burnInSteps 100 --rMax 0.1"  #MarkovChainMC Line

    #CombineCMD="combine -H "+HintAlgo+" -M "+Algo+" -s "+Seed+" -n "+Mode+" -m "+Mass+" "+cardFile+""
    #CombineCMD="combine -H "+HintAlgo+" --rMax 1 -M "+Algo+" -s "+Seed+" -n "+Mode+" -m "+Mass+" "+cardFile+""
    #CombineCMD="combine -H "+HintAlgo+" -V -v 2 -M "+Algo+" -s "+Seed+" -n "+Mode+" -m "+Mass+" "+cardFile+""
    print CombineCMD
    if options.observed:
        os.system(CombineCMD + "              ") # >& /dev/null   #Observed Limit
    else:
        os.system(CombineCMD + " --toys "+NToys) # >& /dev/null   #Expected Limit

    name="higgsCombine"+Mode+Suffix+"."+Algo+".mH"+Mass+".root"
    os.system("hadd -f "+name+" higgsCombine"+Mode+Suffix+"."+Algo+".mH"+Mass+".*.root >& /dev/null")


def makePlots():
    #Now run code to extract limts and make plots
    extractCMD="root -b -l -q 'extractLimits.C+(\"../../../WprimeWZ.root\",\""+Mode+Suffix+"\",\""+Algo+"\")'"
    print extractCMD
    os.system(extractCMD)
    os.system("cat nLimit_"+Mode+Suffix+"_"+Algo+".txt")

    os.chdir("../Limits")
    os.system("root -b -l -q '../Limits/PlotLimit.C+(\"WprimeWZ\", \"../combined_limits/nLimit_"+Mode+Suffix+"_"+Algo+".txt\")'")
    os.system("scp  ../Limits/limitVsMass_WZ.pdf buphy.bu.edu:~/public_html")

def submitCondor(mass):
    condorFile="submit_M"+mass+"_LtShift"+options.LtShift+"_WindShift"+options.WindShift+".condor"
    f = open(condorFile, 'w')
    f.write("universe = vanilla\n")
    f.write("Executable = condorLimits.sh\n")
    #f.write("Executable = doCombinedLimits.py\n")
    f.write('Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000\n')
    f.write("Should_Transfer_Files = YES\n")
    f.write("WhenToTransferOutput = ON_EXIT\n")
    label="M"+mass+"LtShift"+options.LtShift+"_WindShift"+options.WindShift
    f.write("Output = out_"+label+".$(Process)_$(cluster)\n")
    f.write("Error = err_"+label+".$(Process)_$(cluster)\n")
    f.write("Log = log_"+label+".$(Process)_$(cluster)\n")
    #f.write("notify_user = fantasia@bu.edu\n")
    f.write("Arguments = "+mass+" "+options.LtShift+" "+options.WindShift+" "+options.ntoys+" "+os.getcwd()+"\n")
    #f.write("Arguments = -M "+mass+" --LtShift="+options.LtShift+" --WindShift="+options.WindShift+" --setup="+os.getcwd()+"\n")
    f.write("Queue 1\n")
    f.close()
    #print condorFile
    os.system("/opt/condor/bin/condor_submit "+ condorFile)


    
#########################
#if len(sys.argv) != 1:
#    print "usage ./runCoresCron.py"
#    sys.exit(1)
##########################
if __name__ == '__main__':
    Mode="WprimeWZ"
    Suffix=""
    if options.LtShift!= "0" or options.WindShift!= "0":
        Suffix="_LtShift"+options.LtShift+"_WindShift"+options.WindShift
    print "Suffix is "+Suffix

    #HintAlgo=ProfileLikelihood
    HintAlgo="Asymptotic"
    #Algo=BayesianToyMC
    Algo="MarkovChainMC"      #Bayesian
    #Algo=ProfileLikelihood  #?
    #Algo=Asymptotic          #Asymptotic CLs

    if options.setup != "":
        os.system("source /uscmst1/prod/sw/cms/setup/shrc prod")
        os.chdir(options.setup)
        os.system("eval `scramv1 runtime -sh`")
        print "Changed to "+os.getcwd()

    if options.Masses != -1:
        Masses=options.Masses.split(",")
        for M in Masses:
            print "Running mass "+M
            if not options.condor:
                run(M)
            else:
                submitCondor(M)


    if options.plot:
        makePlots()
        sys.exit(0)



##########################
#for line in txtfile:
#    addressCore = line.split(":")[0].rstrip()
#    currCore = line.split(":")[1].rstrip()
#    currCorePID = line.split("tmp/")[1].rstrip()
#############

 