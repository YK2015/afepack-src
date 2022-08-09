#!/bin/bash

echo ==========================
echo = 请不要用 sudo 运行本程序 =
echo = 请将 dell.II 和 AFEPack 压缩包放在 $HOME 目录下 =
echo = 我们将 deall.II 和 AFEPack 安装在 $HOME/local 这个目录中 =
echo =================
sleep 2

echo ================
echo = 编译easymesh =
echo ================
cd $HOME
gcc -o easymesh -O3 easymesh.c -lm
sudo mv easymesh /usr/local/bin
sleep 2

echo =================
echo = 先开始更新机器 =
echo =================
sudo apt update
sudo apt upgrade -y
sleep 2

echo ===================
echo = 安装必须的支持包 =
echo ===================
sudo apt install -y cmake g++ libboost-all-dev libtbb-dev automake liblapack-dev libmumps-dev libptscotch-dev libsuitesparse-dev libarpack2-dev trilinos-all-dev
sleep 2

echo ====================
echo = 编译并安装dell.II =
echo ====================
cd $HOME
tar zxvf deal.II.tar.gz
cd deal.II
mkdir -p build
cd build
#cmake -DCMAKE_INSTALL_PREFIX=$HOME/local -DDEAL_II_WITH_MPI=ON -DDEAL_II_WITH_UMFPACK=OFF -DDEAL_II_COMPONENT_MESH_CONVERTER=OFF -DDEAL_II_WITH_PETSC=OFF -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_FC_COMPILER="mpif90" ..
cmake -DCMAKE_INSTALL_PREFIX=$HOME/local  -DDEAL_II_COMPONENT_MESH_CONVERTER=OFF -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_FC_COMPILER="mpif90" ..
make
#make test
make install
sleep 2

echo ====================
echo = 编译并安装AFEPack =
echo ====================
cd $HOME
tar zxvf AFEPack.tar.gz
cd AFEPack
CC="mpicc" CXX="mpicxx" FC="mpifort" CXXFLAGS="-O2 -fno-delete-null-pointer-checks -I/usr/include/trilinos" LDFLAGS="-L$HOME/local/lib" CPPFLAGS="-I$HOME/local/include -I$HOME/local/include/deal.II -I/usr/include/trilinos" CFLAGS="-I/usr/include/trilinos" ./configure --prefix="$HOME/local"
make
make install

#gcc -o showmesh -O3 showmesh.c -lX11

echo ====================
echo = 成功完成 AFEPack 安装！ =
echo = 你可以使用 AFEPack 样例中的 poisson equation 来测试 =
echo ====================

