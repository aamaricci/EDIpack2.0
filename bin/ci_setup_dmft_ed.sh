#Building LIB_DMFT_ED
#Errors
set -e

cd EDIpack2.0
mkdir build
cd build

echo "cmake .."
cmake ..

echo "make"
make

echo "make install"
make install

echo "source ~/opt/dmft_ed/gnu/*/bin/dmft_ed_config_user.sh" >> ~/.dmft_ed_config_user
echo -e "\e[32m EDIpack2.0 installed and sourced \e[0m"

