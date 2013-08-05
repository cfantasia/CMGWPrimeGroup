import re, os, sys

# Function that returns the value of coupling corresponding to limit
def findIntersection(sigmas, limit):
    intersection = 0
    y1 = 0 # coupling value for sigma < limit
    y2 = 0 # coupling value for sigma > limit
    x1 = 0
    x2 = 0
    smaller = 1
    listOfCouplings = list(sigmas.keys())
    listOfCouplings.sort()
    for coupling in listOfCouplings:
        if sigmas[coupling] < limit:
            y1 = sigmas[coupling]
            x1 = coupling
            #print "Smaller!", y1, "compared with", limit
        if sigmas[coupling] > limit:
            y2 = sigmas[coupling]
            x2 = coupling
            #print "Larger!", y2, "compared with", limit            
            break
    # Now a simple interpolation
    slope = (y1 - y2)/(x1 - x2)
    offset = y1 - slope*x1
    result = (limit - offset)/slope
    return result
            
leftValues = [0, 5, 10]
massValues = []
couplingValues = []

# Read in Cory's numbers
file = open("../../combined_limits/nLimit_WprimeWZ_MarkovChainMC.txt", "r")
#file = open("input2.txt", "r")
lines = file.readlines()
file.close()

# Dictionaries (to ease lookup values by mass)
obsLimit = {}
expLimit = {}
expP1Limit = {}
expM1Limit = {}
expP2Limit = {}
expM2Limit = {}

for line in lines:
    if ":" in line or "#" in line: #Ignore that first line
        continue
    notZero = []
    for i in re.split(" |\t", line):
        if i <> '' and i <> '\t':
            notZero += [i,]
    currentMass = float(notZero[0])
    massValues += [notZero[0],]
    obsLimit[currentMass]   = notZero[4]
    expLimit[currentMass]   = notZero[5]
    expP1Limit[currentMass] = notZero[6]
    expM1Limit[currentMass] = notZero[7]
    expP2Limit[currentMass] = notZero[8]
    expM2Limit[currentMass] = notZero[9][:-1]

for i in range(2, 102, 2):
    couplingValues += [i,]


#kfactor = {200: 1.357, 250: 1.357, 300: 1.357, 400: 1.357, 500: 1.357, 600: 1.351, 700: 1.352, 800: 1.347, 900: 1.341, 1000:1.341, 1100: 1.341, 1200:1.341, 1300:1.341, 1400:1.341, 1500:1.341, 1600:1., 1700:1., 1800:1., 1900:1., 2000:1.}#Fix k-facotrs
kfactor = {                           200: 1.347,  300: 1.347,  400: 1.355,
            500: 1.363,  600: 1.357,  700: 1.351,  800: 1.349,  900: 1.347,
           1000: 1.339, 1100: 1.331, 1200: 1.324, 1300: 1.317, 1400: 1.305,
           1500: 1.293, 1600: 1.275, 1700: 1.257, 1800: 1.244, 1900: 1.230,
           2000: 1.214}#Fix k-facotrs

expC = {}
obsC = {}
expP1C = {}
expM1C = {}
expP2C = {}
expM2C = {}

# consider 1000 GeV point for now
for l in leftValues:
    if l <> 0:
        continue
    for m in massValues:
        m = int(m)
        if m%100 is not 0:
            continue
        sigmas = {}
        
        for c in couplingValues:
            left = str(l)
            mass = str(m)
            coupling = str(c)

            name = "LOGS/8TeV/wprime_m_" + mass + "_l_" + left + "_c_" + coupling + "_pythia6_cff_py_GEN.py.log"
            if not os.path.isfile(name):
                continue
            sigma = []
            file = open(name, "r")
            lines = file.readlines()
            file.close()
            for line in lines:
                elements = re.split("I 142", line)
                if len(elements) > 1:
                    sigma += [line[:-1],]
            line = sigma[0]
            elements = re.split("I", line)
            xs = elements[3]
            xs1 = re.split("D", xs)[0]
            xs2 = re.split("D", xs)[1]
            xs = xs1 + "E" + xs2
            #print str(c)+" "+str(m)+" "+str(xs)
            #print len(sigmas)
            sigmas[float(c)/10] = float(xs)*1E9*kfactor[m]
            #print float(c)/100, float(xs)*1E9, expLimit[m]
            # Find the intersection
            
        expC[m]   = findIntersection(sigmas, float(expLimit[m]))
        obsC[m]   = findIntersection(sigmas, float(obsLimit[m]))
        expP1C[m] = findIntersection(sigmas, float(expP1Limit[m]))
        expM1C[m] = findIntersection(sigmas, float(expM1Limit[m]))
        expP2C[m] = findIntersection(sigmas, float(expP2Limit[m]))
        expM2C[m] = findIntersection(sigmas, float(expM2Limit[m]))

masses = expC.keys()
masses.sort()

# Now producing the macro
header = """void plotExclusion() {
   gStyle->SetOptStat(0);
   setTDRStyle();
   gROOT->ForceStyle();
   TH2F* frame = new TH2F("frame", "", 100, 150, 2050, 100, 0.05, 10.0);
   frame->GetYaxis()->SetTitle("W'WZ coupling");
   frame->GetXaxis()->SetTitle("M(W') (GeV)");
   TCanvas* c1 = new TCanvas("c1");
   c1->SetLogy(1);

"""
file = open("plotExclusion.C", "w")
file.write("#include \"../../root_macros/CMSStyle.C\"\n");
file.write("#include \"../../root_macros/setTDR_modified.C\"\n");
file.write(header);
# Mass ===================================================
file.write("   float masses[" + str(len(masses)) + "];\n")
i = 0
for mass in masses:
    file.write("   masses[" + str(i) + "] = " + str(mass) + ";\n")
    i += 1
file.write("\n")
# Mass for 2D contours
file.write("   float masses2D[" + str(len(masses)*2) + "];\n")
i = 0
for mass in masses:
    file.write("   masses2D[" + str(i) + "] = " + str(mass) + ";\n")
    file.write("   masses2D[" + str(2*len(masses) - i - 1) + "] = " + str(mass) + ";\n")
    i += 1

# Exp/Obs ===================================================
file.write("   float exp[" + str(len(masses)) + "];\n")
file.write("   float obs[" + str(len(masses)) + "];\n")
file.write("   float exp1[" + str(2*len(masses)) + "];\n")
file.write("   float exp2[" + str(2*len(masses)) + "];\n")
i = 0
for mass in masses:
    file.write("   exp[" + str(i) + "] = " + str(expC[mass]) + ";\n")
    file.write("   obs[" + str(i) + "] = " + str(obsC[mass]) + ";\n")
    # now +/- 1 sigma and 2 sigma contours
    file.write("   exp1[" + str(i) + "] = " + str(expP1C[mass]) + ";\n")
    file.write("   exp1[" + str(2*len(masses) - i - 1) + "] = " + str(expM1C[mass]) + ";\n")
    file.write("   exp2[" + str(i) + "] = " + str(expP2C[mass]) + ";\n")
    file.write("   exp2[" + str(2*len(masses) - i - 1) + "] = " + str(expM2C[mass]) + ";\n")
    i += 1

file.write("   TGraph* grExp = new TGraph(" + str(len(masses)) + ", masses, exp);\n")
file.write("   TGraph* grObs = new TGraph(" + str(len(masses)) + ", masses, obs);\n")
file.write("   TGraph* grExp1 = new TGraph(" + str(2*len(masses)) + ", masses2D, exp1);\n")
file.write("   TGraph* grExp2 = new TGraph(" + str(2*len(masses)) + ", masses2D, exp2);\n")

tail = """
   grExp->SetLineColor(kBlack);
   grExp->SetMarkerColor(kBlack);
   grExp->SetLineStyle(kDashed);
   grObs->SetLineColor(kBlack);
   grObs->SetMarkerColor(kBlack);
   grExp1->SetFillStyle(1001);
   grExp1->SetFillColor(kGreen);
   //grExp1->SetLineColor(kGreen);
   grExp2->SetFillStyle(1001);
   grExp2->SetFillColor(kYellow);
   //grExp2->SetLineColor(kYellow);
   grExp->SetLineWidth(3);
   grObs->SetLineWidth(3);
"""
file.write(tail)

file.write("   frame->Draw();\n")
file.write("   frame->GetXaxis()->SetNdivisions(505);\n")
file.write("   grExp2->Draw(\"F\");\n")
file.write("   grExp1->Draw(\"F\");\n")
file.write("   grExp->Draw(\"C\");\n")
file.write("   grObs->Draw(\"C\");\n")
#file.write("   c1->RedrawAxis();\n")
file.write("   TLine* line = new TLine(150, 1, 2050, 1);\n")
file.write("   line->SetLineWidth(3);\n")
file.write("   line->SetLineStyle(1);\n")
file.write("   line->SetLineColor(kBlue);\n")
file.write("   line->Draw();\n")

file.write("  TLatex latexLabel;\n")
file.write("  latexLabel.SetNDC();\n")
file.write("  latexLabel.SetTextSize(0.05);\n")
file.write("  latexLabel.SetTextFont(42);\n")
#file.write("  latexLabel.DrawLatex(0.42, 0.96, \"CMS 2012\");\n")
file.write("  latexLabel.DrawLatex(0.33, 0.96, \"CMS Preliminary 2012\");\n")
file.write("  latexLabel.DrawLatex(0.59, 0.30, \"#sqrt{s} = 8 TeV\");\n")
file.write("  latexLabel.DrawLatex(0.57, 0.20, Form(\"#intL dt = %.1f fb^{-1}\",19.6));\n")

tail = """
   TLegend* leg = new TLegend(0.17, 0.62, 0.58, 0.92);
   leg->AddEntry(grObs, "Obs. 95% C.L.", "l");
   leg->AddEntry(grExp, "Exp. 95% C.L.", "l");
   leg->AddEntry(grExp1, "Exp. #pm 1#sigma", "f");
   leg->AddEntry(grExp2, "Exp. #pm 2#sigma", "f");
   leg->AddEntry(line, "SSM W'WZ coupling", "l");
   leg->SetTextSize(0.05);
   leg->SetTextFont(42);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->Draw();



   c1->RedrawAxis();
   c1->Print("coupling.pdf");
"""
file.write(tail)
file.write("}\n")
file.close()


for mass in masses:
    print mass, obsC[mass], expC[mass], expP2C[mass], expP1C[mass], expM1C[mass], expM2C[mass]

#Run the program!
os.system("root -b -l -q plotExclusion.C")
