#!/bin/bash

pathToSuite=$(pwd)
suite="suite.mag"
integrator="integrator.sage"
t_integrator="transition-integrator.sage"


mkdir -p ${pathToSuite}"/""ode_storage/incinerator";
mkdir -p ${pathToSuite}"/""fermat_data";
touch ${pathToSuite}"/""ode_storage/incinerator/.PERIODSUITE-this-directory-is-safe-to-rm-fr";

sed -i.backup "1d" $suite
sed -i.backup "1i\\
pathToSuite:=\"${pathToSuite}"/"\";
" $suite

# set directory names in sage
cat > SAGE_CONFIG.py << EOF
pathToSuite = "$pathToSuite/"
EOF

sed -i.backup "1d" $integrator
sed -i.backup "1i\\
pathToSuite=\"${pathToSuite}"/"\";
" $integrator

sed -i.backup "1d" $t_integrator
sed -i.backup "1i\\
pathToSuite=\"${pathToSuite}"/"\";
" $t_integrator
