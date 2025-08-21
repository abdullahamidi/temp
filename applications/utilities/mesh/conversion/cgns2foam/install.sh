#!/bin/bash

SOURCE_DIR=$(realpath $(dirname "${BASH_SOURCE[0]}"))

function check_execute() {
    $@
    if [ $? -ne 0 ]
    then
        echo Unable to execute $@
        exit 1
    fi
}

check_execute python $SOURCE_DIR/dependency.py --unbundle
check_execute python $SOURCE_DIR/dependency.py --install
check_execute cmake $SOURCE_DIR -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$1
check_execute cmake --build .
check_execute cmake --install .

