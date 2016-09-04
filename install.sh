#!/bin/bash

#LAMPVER="3.0.0"
#LAMPURL="https://github.com/a-terada/lamp/archive/lamp-${LAMPVER}.tar.gz"
#LAMPDIR="lamp"
#LAMPFILE="lamp-${LAMPVER}.tar.gz"

PLINKURL="http://pngu.mgh.harvard.edu/~purcell/plink/dist/plink-1.07-src.zip"
PLINKDIR="src"
PLINKFILE="plink-1.07-src.zip"

## Install LAMP
#echo "Downloading LAMP library..."
#echo
#pushd ./lib
#if [ ! -e ${LAMPFILE} ];then
#    curl -O ${LAMPURL} -L
#else
#    echo "Skip downloading LAMP"
#fi
#if [ ! -d ${LAMPDIR} ];then
#    tar zxf ${LAMPFILE}
#    mv lamp-${LAMPVER} ${LAMPDIR}
#fi
#pushd ${LAMPDIR}
#echo "Installing LAMP library..."
#echo
#make
#popd
## Add test for LAMP
#if [ ! -e ${LAMPDIR}/lib-lamp.py ];then
#    mv lib-lamp.py ${LAMPDIR}
#fi
#popd

# Download PLINK source code
echo "Downloading PLINK..."
echo
if [ ! -e ${PLINKFILE} ];then
    curl -O ${PLINKURL} -L
else
    echo "Skip downloading PLINK"
fi
if [ ! -d ${PLINKDIR} ];then
    unzip ${PLINKFILE}
    mv plink-1.07-src ${PLINKDIR}
fi

# Install LAMPLINK
echo "Create symbolic link of LAMPLINK files in PLINK directory"
echo
if [ -e lamp.cpp ];then
    ln -sf "$PWD/"*.cpp ${PLINKDIR}/
    ln -sf "$PWD/"*.h ${PLINKDIR}/
    ln -sf "$PWD/"Makefile ${PLINKDIR}/
    ln -sf "$PWD/"lcm53 ${PLINKDIR}/
    ln -sf "$PWD/"functions ${PLINKDIR}/
    ln -sf "$PWD/"frepattern ${PLINKDIR}/
fi
pushd ${PLINKDIR}
echo "Installing LAMPLINK..."
make
popd
