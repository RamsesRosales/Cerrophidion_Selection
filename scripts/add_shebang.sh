#!/bin/bash

###small script to add the correct shebang to my previous scripts that dont have it, this is because I did not do that for my rscripts
##first argument is the name of the script
##second argument is the shebang


echo adding $2 to $1

echo $2 > tmp
cat $1>>tmp
cat tmp>$1
rm tmp
perl -p -i -e 's/\#\\\!/\#\!/g' $1



