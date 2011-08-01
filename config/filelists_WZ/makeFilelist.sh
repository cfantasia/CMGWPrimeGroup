#!/bin/bash
Ver="W05-08-02"
Dir=/hdfs/store/user/jklukas
#Dir=/pnfs/cms/WAX/11/store/user/fantasia/42X
#Dir=~/nobackup/42X/filelists

for Directory in `ls ${Dir} | grep ${Ver} | grep PatTuple`
  do
  if [ -d $Dir/$Directory ]; then
      find ${Dir}/${Directory}/ | grep .root >& ${Directory}.txt
  fi
done

./mergeFilelists.sh