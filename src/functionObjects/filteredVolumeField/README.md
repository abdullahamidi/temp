## Configuration - controlDict
Following template can be used to use filteredVolumeField function object in OpenFOAM/HiSA/TAMS
<pre>
<code>
    filteredPressure
    {
        type            filteredVolumeField;
        libs            ("libAeroacousticsFunctionObjects.so");
        writeControl    writeTime;

        field       p;

        enableBoundingBox true;
        boundingBoxLowerLeft     (-10 -10 -10);
        boundingBoxUpperRight    (10 10 10);

        enableYFilter true;
        minY          0.00;
        maxY          0.00001;

        enableNutFilter false;
        minNut          0.0001;
        maxNut          0.0005;

        enableNuTildaFilter false;
        minNuTilda      0.0001;
        maxNuTilda      0.0005;

        enableKFilter false;
        minK            0.0001;
        maxK            0.0005;
        
        enableOmegaFilter false;
        minOmega        0.0001;
        maxOmega        0.0005;
    }
</code>
</pre>

You can collect processor data using `reconstructFilteredField` tool.

# Notes
- It may be helpful to look below link and articles to determine y values 
  https://volupe.com/simcenter-star-ccm/transition-modelling-in-simcenter-star-ccm/ 
- Don't use purge write option since it will delete whole folder.