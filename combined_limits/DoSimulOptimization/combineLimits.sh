#!/bin/bash

for Code in 200 250 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500
  do

  Frac=$(($Code % 100))
  Frac=$(($Frac / 10))

  Base=$(($Code / 100))

    grep "SignalCode"  nLimits_${Code}_part0.txt > nLimits_${Code}.txt
    cat nLimits_${Code}_part?.txt | grep "^${Base}.${Frac}" >> nLimits_${Code}.txt

done

