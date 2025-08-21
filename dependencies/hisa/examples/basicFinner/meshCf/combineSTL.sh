#!/bin/sh

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

HOMEDIR=`pwd`

source ../settings.sh

# Combine stl files for cfMesh
cd constant/triSurface

echo "Combine stl files for cfMesh"

# Combine STL files
cp farfield.stl combined.stl
for i in "${parts[@]}"
do
    cat combined.stl $i\Scaled.stl > temp.stl
    mv -f temp.stl combined.stl
done

cd $HOMEDIR


# ----------------------------------------------------------------- end-of-file
