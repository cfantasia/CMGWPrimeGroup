#!/bin/bash

#######  Eff Rates   ####################

rm -f EffRates.dat

date >> EffRates.dat
echo "\n#######Now on Data T #######"  >> EffRates.dat

#Elec Barrel Data
echo "\n#######Now on Elec Barrel #######" >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 1, 0,  0, 20)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 1, 0, 20, 25)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 1, 0, 25, 30)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 1, 0, 30, 50)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 1, 0, 50, 9e9)' >> EffRates.dat

#Elec EC Data
echo "\n#######Now on Elec EndCap #######" >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 1, 1,  0, 20)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 1, 1, 20, 25)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 1, 1, 25, 30)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 1, 1, 30, 50)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 1, 1, 50, 9e9)' >> EffRates.dat

#Muons Data
echo "\n#######Now on Muon #######" >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 0, 1, 0,  0, 20)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 0, 1, 0, 20, 25)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 0, 1, 0, 25, 30)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 0, 1, 0, 30, 50)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 0, 1, 0, 50, 9e9)' >> EffRates.dat

echo "\n#######Now on Data TT #######"  >> EffRates.dat

#Elec Barrel Data
echo "\n#######Now on Elec Barrel #######" >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 1, 0,  0, 20)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 1, 0, 20, 25)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 1, 0, 25, 30)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 1, 0, 30, 50)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 1, 0, 50, 9e9)' >> EffRates.dat

#Elec EC Data
echo "\n#######Now on Elec EndCap #######" >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 1, 1,  0, 20)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 1, 1, 20, 25)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 1, 1, 25, 30)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 1, 1, 30, 50)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 1, 1, 50, 9e9)' >> EffRates.dat

#Muons Data
echo "\n#######Now on Muon #######" >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 0, 1, 0,  0, 20)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 0, 1, 0, 20, 25)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 0, 1, 0, 25, 30)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 0, 1, 0, 30, 50)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 0, 1, 0, 50, 9e9)' >> EffRates.dat

#############MC
echo "\n#######Now on MC T #######" >> EffRates.dat

#Elec Barrel MC
echo "\n#######Now on Elec Barrel #######" >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 0, 0,  0, 20)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 0, 0, 20, 25)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 0, 0, 25, 30)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 0, 0, 30, 50)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 0, 0, 50, 9e9)' >> EffRates.dat

#Elec EC MC
echo "\n#######Now on Elec EndCap #######" >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 0, 1,  0, 20)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 0, 1, 20, 25)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 0, 1, 25, 30)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 0, 1, 30, 50)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 1, 0, 1, 50, 9e9)' >> EffRates.dat

#Muons MC
echo "\n#######Now on Muon #######" >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 0, 0, 0,  0, 20)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 0, 0, 0, 20, 25)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 0, 0, 0, 25, 30)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 0, 0, 0, 30, 50)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate.root", 0, 0, 0, 50, 9e9)' >> EffRates.dat

echo "\n#######Now on MC TT #######" >> EffRates.dat

#Elec Barrel MC
echo "\n#######Now on Elec Barrel #######" >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 0, 0,  0, 20)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 0, 0, 20, 25)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 0, 0, 25, 30)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 0, 0, 30, 50)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 0, 0, 50, 9e9)' >> EffRates.dat

#Elec EC MC
echo "\n#######Now on Elec EndCap #######" >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 0, 1,  0, 20)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 0, 1, 20, 25)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 0, 1, 25, 30)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 0, 1, 30, 50)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 1, 0, 1, 50, 9e9)' >> EffRates.dat

#Muons MC
echo "\n#######Now on Muon #######" >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 0, 0, 0,  0, 20)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 0, 0, 0, 20, 25)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 0, 0, 0, 25, 30)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 0, 0, 0, 30, 50)' >> EffRates.dat
root -b -l -q 'Est_Zjets.C+("../../../WZEffRate-TT.root", 0, 0, 0, 50, 9e9)' >> EffRates.dat


echo "\n ####  Summary ####\n" >> EffRates.dat
grep "Now\|e_tight" EffRates.dat >> EffRates.dat

#######  Fake Rates   ####################
rm -f FakeRates.dat

date >> FakeRates.dat
echo "\n#######Now on Data #######" >> FakeRates.dat
root -b -l -q 'CalcFakeRate.C+("../../../WZFakeRateElec.root", 1)' >> FakeRates.dat
root -b -l -q 'CalcFakeRate.C+("../../../WZFakeRateMuon.root", 1)' >> FakeRates.dat
root -b -l -q 'CalcFakeRate.C+("../../../WZFakeRateElec-TT.root", 1)' >> FakeRates.dat
root -b -l -q 'CalcFakeRate.C+("../../../WZFakeRateMuon-TT.root", 1)' >> FakeRates.dat

#############MC
echo "\n#######Now on MC #######" >> FakeRates.dat

root -b -l -q 'CalcFakeRate.C+("../../../WZFakeRateElec.root", 0)' >> FakeRates.dat
root -b -l -q 'CalcFakeRate.C+("../../../WZFakeRateMuon.root", 0)' >> FakeRates.dat
root -b -l -q 'CalcFakeRate.C+("../../../WZFakeRateElec-TT.root", 0)' >> FakeRates.dat
root -b -l -q 'CalcFakeRate.C+("../../../WZFakeRateMuon-TT.root", 0)' >> FakeRates.dat
