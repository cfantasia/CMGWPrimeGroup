#!/bin/bash

for FILE in `ls *.txt`
  do
  FILEPREFIX=${FILE%%.txt}_part
  
  NLINES_TOT=$(cat $FILE | wc -l)
  NFILES=$((($NLINES_TOT / 300) + 1))
  NLINES=$((($NLINES_TOT / $NFILES) + 1))
  
  if [ $NFILES -gt 1 ] 
      then
      echo Splitting $FILE into $NFILES files
  
      #echo "Count is greater than 1"
      #echo $FILEPREFIX
      #echo NLines_TOT $NLINES_TOT
      #echo NFILES $NFILES
      #echo NLines $NLINES
      
      split --suffix-length=1 --numeric-suffixes -l $NLINES $FILE $FILEPREFIX
      
      for PART_I in `ls ${FILEPREFIX}*` 
        do
        mv $PART_I{,.txt}  
      done
      
      rm $FILE
  fi
done