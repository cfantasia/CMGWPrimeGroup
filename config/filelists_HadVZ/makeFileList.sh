#!/bin/bash


for i in `cat ListOfData.txt`;do
    srmls srm://osg-se.sprace.org.br:8443/pnfs/sprace.org.br/data/cms/store/user/fladias/${i} | grep root | awk '{print $2}' > ${i}.txt
done
