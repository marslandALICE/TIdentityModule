#!/bin/bash
if [ ! $ROOTSYS ]
  then
    echo '************'
    echo '$ROOTSYS is not found, set your ROOT enviroment'
    echo '************'
exit 1
fi
if [ ! $TIdentityDIR ]
  then
    echo '************'
    echo 'export TIdentityDIR=PATHtoTIdentity'
    echo '************'
else
echo $TIdentityDIR
cd $TIdentityDIR
rm -rf build
rm -rf lib
rm -rf bin
mkdir build
cd build
cmake ..
#make
cmake --build .
cd $TIdentityDIR/test
fi
