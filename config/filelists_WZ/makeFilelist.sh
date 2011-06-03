#!/bin/bash
Ver="B05-06-03"
Dir=/pnfs/cms/WAX/11/store/user/fantasia/42X

for Directory in `ls ${Dir} | grep ${Ver}`
  do
  find ${Dir}/${Directory}/ | grep .root >& ${Directory}.txt

done
