#!/bin/bash
Ver="B05-06-07"
Dir=/pnfs/cms/WAX/11/store/user/fantasia/42X

for Directory in `ls ${Dir} | grep ${Ver}`
  do
  find ${Dir}/${Directory}/ | grep .root >& ${Directory}.txt

done

./mergeFilelists.sh