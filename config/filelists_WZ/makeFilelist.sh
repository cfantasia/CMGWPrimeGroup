#!/bin/bash
Ver="05-07-01"
#Dir=/pnfs/cms/WAX/11/store/user/fantasia/42X
Dir=~/nobackup/42X/filelists

for Directory in `ls ${Dir} | grep ${Ver}`
  do
  if [ -d $Dir/$Directory ]; then
      find ${Dir}/${Directory}/ | grep .root >& ${Directory}.txt
  fi
done

./mergeFilelists.sh