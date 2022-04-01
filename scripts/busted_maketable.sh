#!/bin/bash

echo -e gen,LRT,p-value,CW1,CW2,CW3,UCW1,UCW2,UCW3 > Busted_Results
for i in *.BUSTED.json
do echo $i
grep -A29 -e 'Constrained model' ${i} > tmp.file
grep -A29 -e 'Unconstrained model' ${i} > tmp1.file
grep -n -e ".*" tmp.file > tmp1.1.file
grep -n -e ".*" tmp1.file > tmp1.2.file
echo 0`grep -e '^21' tmp1.1.file`,0`grep -e '^25' tmp1.1.file`,0`grep -e '^29' tmp1.1.file`,0`grep -e '^21' tmp1.2.file`,0`grep -e '^25' tmp1.2.file`,0`grep -e '^29' tmp1.2.file` > tmp2.file
perl -p -i -e "s/21\: //g" tmp2.file
perl -p -i -e "s/25\: //g" tmp2.file
perl -p -i -e "s/29\: //g" tmp2.file
perl -p -i -e "s/\"omega\"\://g" tmp2.file
perl -p -i -e "s/,\n/\n/g" tmp2.file
perl -p -i -e "s/,00\./,0\./g" tmp2.file
perl -p -i -e "s/^00\./0\./g" tmp2.file
perl -p -i -e 's/,0([0-9]+)\./,$1\./g' tmp2.file
perl -p -i -e 's/0([0-9]+),,/$1,,/g' tmp2.file
echo $i,`grep -e '\"LRT\"\:' $i`,`grep -e '\"p-value\"\:' $i`, `cat tmp2.file` >> Busted_Results
done
perl -p -i -e "s/^Cgodm_//g" Busted_Results
perl -p -i -e "s/_noStop\.fasta\.BUSTED\.json//g" Busted_Results
perl -p -i -e "s/\"LRT\"\://g" Busted_Results
perl -p -i -e "s/\"p-value\"\://g" Busted_Results
perl -p -i -e "s/,,/,/g" Busted_Results
