#
# Script that plots all the stability results.
#

#!/bin/bash

# definitions:
TEST_FOLDER=./tests_EFT

# go to the eftcamb root folder:
cd ..

# do the stability plot:
echo ; echo ; echo ;
echo 'Plotting stability.'
echo ; echo ; echo ;

cd ${TEST_FOLDER}/python/

python stab_plot.py

echo 'done.'

exit 0
