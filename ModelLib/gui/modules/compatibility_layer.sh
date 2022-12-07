#usr/bin/bash
tmp="temptemptemp.temp"
cat $1 > $tmp

sed -i "s/INDEX =/NUMBER =/g" $tmp

if [ $(diff $1 $tmp |wc -l) -gt 0 ]
then
  echo InitFile was updated to current version. Changes are:
  diff $1 $tmp
  mv $tmp $1
else
  rm $tmp
fi
