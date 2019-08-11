#!/bin/bash

pathToSuite=$(pwd)
suite="suite.mag"
integrator="integrator.sage"

mkdir -p ${pathToSuite}"/""ode_storage/incinerator";
mkdir -p ${pathToSuite}"/""fermat_data";

sed -i.backup "1d" $suite
sed -i.backup "1i\\
pathToSuite:=\"${pathToSuite}"/"\";
" $suite

sed -i.backup "1d" $integrator
sed -i.backup "1i\\
pathToSuite=\"${pathToSuite}"/"\";
" $integrator
