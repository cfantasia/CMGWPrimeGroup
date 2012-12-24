#!/bin/bash
Ver="07-02-01"
#Base=/hdfs/store/user/jklukas
Base=/pnfs/cms/WAX/11/store/user/fantasia/53X
#Base=~/nobackup/42X/filelists

for Directory in `ls ${Base}/ | grep ${Ver}`
  do
  echo Directory is $Directory
  find ${Base}/${Directory} | grep .root >& ${Directory}.txt
done

#Strip off prefix
/usr/bin/perl -p -i -e "s/\/pnfs\/cms\/WAX\/11//g" *.txt
#/usr/bin/perl -p -i -e "s/\/hdfs//g" *.txt
