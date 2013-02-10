#!/bin/bash

if [ "$#" -ge 1 ]; then
    Ver=$1
else
    Ver="07-03-0"
fi

echo Checking Directory $Ver
#Base=/hdfs/store/user/jklukas
Base=/pnfs/cms/WAX/11/store/user/clint/53X
#Base=/pnfs/cms/WAX/11/store/user/fantasia/53X
#Base=~/nobackup/42X/filelists

for Directory in `ls ${Base}/ | grep ${Ver}`
  do
  echo Directory is $Directory
  find ${Base}/${Directory} | grep .root >& ${Directory}.txt
done

#Strip off prefix
/usr/bin/perl -p -i -e "s/\/pnfs\/cms\/WAX\/11//g" *.txt
#/usr/bin/perl -p -i -e "s/\/hdfs//g" *.txt
