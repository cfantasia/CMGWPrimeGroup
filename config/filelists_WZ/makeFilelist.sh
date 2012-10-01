#!/bin/bash
Ver="W06-09-08"
#Base=/hdfs/store/user/jklukas
Base=/pnfs/cms/WAX/11/store/user/fantasia/42X
#Base=~/nobackup/42X/filelists

for Directory in `ls ${Base}/ | grep ${Ver} | grep ${Ver}`
  do
  #echo Directory is $Directory
  for Sample in `ls ${Base}/${Directory}/`
    do
    echo Sample is $Sample
      if [ -d $Base/$Directory/$Sample ]; then
          find ${Base}/${Directory}/$Sample | grep .root >& ${Sample}.txt
      fi
      #for List in `ls ${Base}/${Directory}/${Sample}/`
      #do
      #echo List is $List
     #done
  done
done

#Strip off prefix
/usr/bin/perl -p -i -e "s/\/pnfs\/cms\/WAX\/11//g" *.txt
#/usr/bin/perl -p -i -e "s/\/hdfs//g" *.txt
#./mergeFilelists.sh
