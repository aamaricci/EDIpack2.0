#!/bin/bash

for file in $(cat list2); do
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
