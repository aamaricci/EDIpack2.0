#Building EDIpack2
#Errors
set -e

cd edipack2
mkdir build
cd build

echo "cmake .."
cmake ..

echo "make"
make

echo "make install"
make install

echo "source ~/opt/edipack2/gnu/*/bin/edipack2_config_user.sh" >> ~/.edipack2_config_user
echo -e "\e[32m EDIpack2 installed and sourced \e[0m"

