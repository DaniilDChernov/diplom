#!/bin/bash
#
# Run a simple test for MedChem_Fragment_Splitter using docker image
#

# 1. first, build the docker image. 
# NB: this step will be automatically omitted by docker once the image was built

BASE_DIR=$( dirname $(dirname $(readlink -f ${0}) ) )
echo -e "BASE_DIR=[${BASE_DIR}]"

DOCKER_DIR="${BASE_DIR}/docker"
CUR_PWD=$(pwd)
cd ${DOCKER_DIR}
echo -e "\n\tRUNNING docker build -t mcfs ."
docker build -t mcfs .
echo -e "\n\tDONE\n"
cd ${CUR_PWD}

# 2. second, run the test from files directory
echo -e "\n\tRUNNING TEST"
docker run --rm -it --name medchem_fragment_splitter -v $(pwd):/opt/mcfs \
    --user "$(id -u):$(id -g)" mcfs \
    /bin/bash -c "cd /opt/mcfs/files; ./test.sh"
echo -e "\n\tDONE\n"

