#! /usr/bin/env python
import sys,os

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-M", "--masses", dest="Masses", default="",
                  help="What masses to run on", metavar="MASSES")
parser.add_option("-o", "--observed", action='store_true',
                  dest="observed", default=False,
                  help="Do Observed Limit")
parser.add_option("-p", "--plot", action='store_true',
                  dest="plot", default=False,
                  help="Print Plots")
parser.add_option("-d", "--doDuplicates", action='store_true',
                  dest="doDuplicates", default=False,
                  help="Skip Limit if Card file exists")
parser.add_option("--combine", action='store_true',
                  dest="combine", default=True,
                  help="Combine multiple root file")
parser.add_option("-c", "--condor", action='store_true',
                  dest="condor", default=False,
                  help="Submit Via Condor")
parser.add_option("--doSepCh", action='store_true',
                  dest="doSepCh", default=False,
                  help="Calculate limit for separately")
parser.add_option("-t", "--ntoys",
                  dest="ntoys", default="10",
                  help="Number of Toys")
parser.add_option("-l", "--LtShift",
                  dest="LtShifts", default="0",
                  help="Lt Shift")
parser.add_option("-w", "--WindShift",
                  dest="WindShifts", default="0",
                  help="Wind Shift")
parser.add_option("-i", "--input", metavar="FILE",
                  dest="input", default="../../../WprimeWZ.root",
                  help="Input Root File to Use")

(options, args) = parser.parse_args()
print options

def runLimit(Mode, Mass, LtShift, WindShift):
    if not os.path.exists(cardFile):
        print cardFile+" does not exist, creating."
        os.system("root -b -l -q 'makeLimitCard.C+(\""+options.input+"\", \""+cardFile+"\", "+Mass+", "+LtShift+", "+WindShift+")'")
    else:
        if not options.doDuplicates: #Card file exists
            return

    if not os.path.exists(cardFile):
        print "Card file "+cardFile+" STILL does not exist, check for errors!"
        sys.exit(1)


    if options.condor:
        return

    NToys=options.ntoys
    Seed="-1"
    
    CombineCMD="combine -H "+HintAlgo+" -M "+Algo+" -s "+Seed+" -n "+Mode+" -m "+Mass+" --rMax 0.1 --tries 100 "+cardFile  #MarkovChainMC Line
    print CombineCMD
    
    if options.observed:
        os.system(CombineCMD                   ) # >& /dev/null   #Observed Limit
    else:
        os.system(CombineCMD + " --toys "+NToys) # >& /dev/null   #Expected Limit

    #Combine root files
    if options.combine:
        name="higgsCombine"+Mode+"."+Algo+".mH"+Mass+".root"
        os.system("hadd -f "+name+" higgsCombine"+Mode+"."+Algo+".mH"+Mass+".*.root >& /dev/null")

            
def makePlots(Mode, Factor):
    if options.combine:
        print "Merging root files now ..."
        import glob
        for Mass in range(200, 2001): #Combine root files
            name="higgsCombine"+Mode+"."+Algo+".mH"+str(Mass)+".root"
            if len(glob.glob("higgsCombine"+Mode+"."+Algo+".mH"+str(Mass)+".*.root")):
                os.system("hadd -f "+name+" higgsCombine"+Mode+"."+Algo+".mH"+str(Mass)+".*.root >& /dev/null")

    #Now run code to extract limts and make plots
    extractCMD="root -b -l -q 'extractLimits.C+(\""+options.input+"\",\""+Mode+"\",\""+Algo+"\", "+Factor+")'"
    print extractCMD
    os.system(extractCMD)
    os.system("cat nLimit_"+Mode+"_"+Algo+".txt")

    os.chdir("../Limits")
    os.system("root -b -l -q '../Limits/PlotLimit.C+(\"WprimeWZ\", \"../combined_limits/nLimit_"+Mode+"_"+Algo+".txt\", \"limitVsMass_"+Mode+".pdf\", "+Factor+")'")
    os.system("scp  ../Limits/limitVsMass_"+Mode+".pdf buphy.bu.edu:~/public_html")
    os.chdir("../combined_limits")

def submitCondor(Mass):
    condorFile="submit_"+cardFile
    inputFiles = "doCombinedLimits.py"
    for X in Channels:
        inputFiles = inputFiles+","+cardFile.replace(Search,Search+X)

    f = open(condorFile, 'w')
    f.write("universe = vanilla\n")
    f.write("Executable = condorLimits.sh\n")
    f.write('Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000\n')
    f.write("Transfer_Input_Files = "+inputFiles+"\n")
    f.write("Should_Transfer_Files = YES\n")
    f.write("WhenToTransferOutput = ON_EXIT\n")
    f.write("Output = out_"+cardFile+".$(Process)_$(cluster)\n")
    f.write("Error = err_"+cardFile+".$(Process)_$(cluster)\n")
    f.write("Log = log_"+cardFile+".$(Process)_$(cluster)\n")
    f.write("notify_user = @bu.edu\n")
    f.write("notification = Error\n")
    f.write("Arguments = "+Mass+" "+options.LtShifts+" "+options.WindShifts+" "+options.ntoys+" "+os.getcwd()+"\n")
    f.write("Queue 1\n")
    f.close()
    #print condorFile
    os.system("/opt/condor/bin/condor_submit "+ condorFile)
    

# Main Function
if __name__ == '__main__':
    Search="WprimeWZ"

    #HintAlgo=ProfileLikelihood
    HintAlgo="Asymptotic"
    #Algo=BayesianToyMC
    Algo="MarkovChainMC"      #Bayesian
    #Algo=ProfileLikelihood  #?
    #Algo=Asymptotic          #Asymptotic CLs

    
    Channels = [""]
    if options.doSepCh:
        Channels += ["Ch0", "Ch1", "Ch2", "Ch3"]
    print "Channels are "
    print Channels

    Masses=options.Masses.split(",")
    Masses=filter(None, Masses)
    print "Masses are "
    print Masses

    LtShifts=options.LtShifts.split(",")
    print "LtShifts are "
    print LtShifts

    WindShifts=options.WindShifts.split(",")
    print "WindShifts are "
    print WindShifts


    #Loop over  Lt,WindShift, Mass, ch 
    for Ch in Channels:
        for LtShift in LtShifts:
            for WindShift in WindShifts:
                Suffix=""
                if options.LtShifts != "0" or options.WindShifts != "0":
                    Suffix="_LtShift"+LtShift+"_WindShift"+WindShift
                    print "Suffix is "+Suffix
                Mode = Search + Ch + Suffix
                for M in Masses:

                    print "Running mass "+M+" Ch "+Ch+" LtShift "+LtShift+" WindShift "+WindShift
                    
                    cardFile="card_"+Mode+"_M"+M+".txt"
                    runLimit(Mode,M,LtShift,WindShift)

                    if options.condor and Ch == "":
                        submitCondor(M)


                if options.plot:
                    if Ch == "":
                        Scale="1."
                    else:
                        Scale="0.25"
                    makePlots(Mode,Scale)




    sys.exit(0)

