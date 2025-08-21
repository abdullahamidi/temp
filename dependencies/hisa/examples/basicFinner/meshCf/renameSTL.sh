#!/bin/sh

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

HOMEDIR=`pwd`

source ../settings.sh

# Combine stl files for cfMesh
cd constant/triSurface

echo "Rename stl files for cfMesh"

# solidRename <fileName.stl> <newSolidName>
function solidRename {
    fileName=$1
    solidName=$2
    echo Renaming solid in $fileName to $solidName 

    sed "s/\(solid *\).*/solid $solidName/"     <$fileName  >temp1
    sed "s/\(endsolid *\).*/endsolid $solidName/"   <temp1      >temp2
    mv temp2 $fileName
    rm -f temp1
} 

# Rename elements of STL
for i in "${parts[@]}"
do
    solidRename $i.stl $i
done
solidRename farfield.stl farfield

cd $HOMEDIR


# ----------------------------------------------------------------- end-of-file
