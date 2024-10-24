#!/bin/bash

cd graphs
grep  -ri "../module" * | awk -F/ '/"/{ print $3 }' | awk '{print $1}' | sed -e "s/\"//g" > ../list
cd ..
sed -e "s/html/rst/g" list2 > list
rm list2


for file in $(cat list); do
  name=$(echo $file | sed -e "s/\.rst//g")
  echo $name > module/$file
  echo "=====================================" >> module/$file
  echo " " >> module/$file
  if [ -f "graphs/${name}.html" ]; then
    echo $name
    echo ".. raw:: html" >> module/$file
    echo "   :file:  ../graphs/${name}.html" >> module/$file
  fi
done
