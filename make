#!/bin/bash

currentDir=$(pwd)
cd $(dirname $0)
pathToSuite=$(pwd)

# writing the path to the root of the repo to be read by Magma
echo "pathToSuite:=\"${pathToSuite}/\";" > $pathToSuite/src/pathToSuite.mag
echo "pathToSuite=\"${pathToSuite}/\";" > $pathToSuite/src/sage_integrator/pathToSuite.py

# prepares folders for storage and for temporary files
mkdir -p ${pathToSuite}"/""fermat_data";
mkdir -p ${pathToSuite}"/""ode_storage/incinerator";
touch ${pathToSuite}"/""ode_storage/incinerator/.PERIODSUITE-this-directory-is-safe-to-rm-fr";

cd $currentDir
