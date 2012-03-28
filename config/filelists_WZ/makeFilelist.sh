#!/bin/bash
Ver="W06-06-05A"
#Base=/hdfs/store/user/jklukas
Base=/pnfs/cms/WAX/11/store/user/fantasia/42X
#Base=~/nobackup/42X/filelists

for Directory in `ls ${Base}/ | grep ${Ver} | grep ${Ver}`
  do
  #echo Directory is $Directory
  for Sample in `ls ${Base}/${Directory}/`
    do
    #echo Sample is $Sample
    for List in `ls ${Base}/${Directory}/${Sample}/`
      do
      #echo List is $List
      if [ -d $Base/$Directory/$Sample/$List ]; then
          find ${Base}/${Directory}/$Sample/$List | grep .root >& ${Ver}-${List}.txt
      fi
    done
  done
done

#Strip off prefix
/usr/bin/perl -p -i -e "s/\/pnfs\/cms\/WAX\/11//g" *.txt
#/usr/bin/perl -p -i -e "s/\/hdfs//g" *.txt
#./mergeFilelists.sh
