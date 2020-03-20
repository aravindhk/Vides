#! /bin/tcsh -f
echo "Compiling NANOTCAD ViDES module for python"
cd ./src
python ./configure.py
make
make install
cd ..
