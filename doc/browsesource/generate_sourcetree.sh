#!/bin/bash

cd graphs
grep  -ri "../module" * | awk -F/ '/"/{ print $3 }' | awk '{print $1}' | sed -e "s/\"//g" > ../list2
cd ..
sort -u list > list2
sed -e "s/html/rst/g" list2 > list
rm list2


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
