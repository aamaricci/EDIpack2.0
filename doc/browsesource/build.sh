#!/bin/bash


# #-----------------------------
# #UNCOMMENT THIS PART IF YOU
# #NEED TO SETUP AND RUN FORD
# #-----------------------------
# cat <<EOF > ford_doc.md
# ---
# project: EDIpack2.0
# preprocess:false
# display:none
# hide_undoc:true
# print_creation_date:true
# src_dir:../../src
# output_dir:./ford_doc
# extensions: f90
# quiet:false
# parallel:0
# graph:true
# graph_maxdepth:5
# graph_maxnodes:20
# ---
# This is my Fortran project!

# EOF

# #Run FORD (not checking actual presence)
# ford ford_doc.md
# #Sync module/*.html files in local src
# rsync -avPhHO --del ford_doc/module/*.html src/
#-----------------------------
#-----------------------------




SRC=src
GRAPH=graphs
MOD=module
MOM=$(pwd)

mkdir -p $GRAPH
mkdir -p $MOD
rm $GRAPH/*
rm $MOD/*


#Generate the  .html images by stripping the svg part from the src files
cd $SRC
pwd 
for file in *.html;do
    echo $file
    sed -n "/<svg id/,/<\/svg>/p" $file > $MOM/$GRAPH/$file
done
cd $MOM


cd $GRAPH

#color tweaks:

#satisfy my ocd
for i in *.html; do
  sed -e "s/white/\#fcfcfc/g" $i > $i.new
  mv $i.new $i
done

#generic blocks in purple
for ifile in *.html; do
  #if there's link
  awk '/<title>/ {print; next_line=1; next} next_line {gsub(/#337ab7/, "#c061cb"); next_line=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
  #if there's no link
  awk '/<title>/ {flag=2; print; next} flag == 2 {flag--; print; next} flag == 1 {gsub(/#337ab7/, "#c061cb"); flag=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
done

#edipack2.0 in consistent blue. Use both regex
for ifile in *.html; do
  #if there's link
  awk '/title="ED/ {print; next_line=1; next} next_line {gsub(/#c061cb/, "#2980b9"); next_line=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
  #if there's no link
  awk '/title="ED/ {flag=2; print; next} flag == 2 {flag--; print; next} flag == 1 {gsub(/#c061cb/, "#2980b9"); flag=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
  #if there's link
  awk '/<title>ED/ {print; next_line=1; next} next_line {gsub(/#c061cb/, "#2980b9"); next_line=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
  #if there's no link
  awk '/<title>ED/ {flag=2; print; next} flag == 2 {flag--; print; next} flag == 1 {gsub(/#c061cb/, "#2980b9"); flag=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
done


#scifortran blocks in green 
for ifile in *.html; do
  #if there's link
  awk '/<title>SF_/ {print; next_line=1; next} next_line {gsub(/#c061cb/, "#26a269"); next_line=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
  #if there's no link
  awk '/<title>SF_/ {flag=2; print; next} flag == 2 {flag--; print; next} flag == 1 {gsub(/#c061cb/, "#26a269"); flag=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
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
  echo "\`Open source file <${githubpath}>\`_ for :f:mod:\`${name}\` on GitHub" >> module/$file
  echo " " >> module/$file
done

# rm list
