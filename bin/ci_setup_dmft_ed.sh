#Building LIB_DMFT_ED
#Errors
set -e

cd LIB_DMFT_ED
mkdir build
cd build

echo "cmake .."
cmake ..

echo "make"
make

echo "make install"
make install

echo "source ~/opt/dmft_ed/gnu/*/bin/dmft_ed_config_user.sh" >> ~/.dmft_ed_config_user
echo -e "\e[32m LIB_DMFT_ED installed and sourced \e[0m"

