#!/bin/bash
SRC=src
GRAPH=graphs
MOD=module
MOM=$(pwd)

mkdir -p $GRAPH
mkdir -p $MOD

#Generate the  .html images by stripping the svg part from the src files
cd $SRC
pwd 
for file in *;do
    echo $file
    sed -n "/<svg id/,/<\/svg>/p" $file > $MOM/$GRAPH/$file
done
cd $MOM


cd $GRAPH

#satisfy my ocd
for i in *.html; do
  sed -e "s/white/\#fcfcfc/g" $i > $i.new
  mv $i.new $i
done


#Generete a list of the actual .html files used in the src images, sorted and uniq-ed
grep  -ri "../module/" * | awk -F/ '/"/{ print $3 }' | awk '{print $1}' | sed -e "s/\"//g" |sort -u > $MOM/list2
cd $MOM
sed -e "s/html/rst/g" list2 > list
rm list2


#generate .rst
for file in $(cat list); do
  name=$(echo $file | sed -e "s/\.rst//g")
  name_upper=$(echo $name | tr '[:lower:]' '[:upper:]' | sed -e "s/HXV/HxV/g")
  relativepath=$(find ../../ -type f -name "*${name_upper}.*90" | awk -Fsrc '{print $2}')
  githubpath="https://github.com/aamaricci/EDIpack2.0/tree/master/src${relativepath}"
  echo $name_upper > module/$file
  echo "=====================================" >> module/$file
  echo " " >> module/$file
  if [ -f "graphs/${name}.html" ]; then
    echo "Found graph for $name"
    echo ".. raw:: html" >> module/$file
    echo "   :file:  ../graphs/${name}.html" >> module/$file
  echo " " >> module/$file
    echo "|" >> module/$file
  else
    echo "Not found graph for $name"
  fi
  echo " " >> module/$file
  echo "\`Open source file <${githubpath}>\`_ on GitHub" >> module/$file
  echo " " >> module/$file
done

# rm list
