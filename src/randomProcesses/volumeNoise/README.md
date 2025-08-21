##  Configuration - noiseDict
Add and configure below lines in system/noiseDict
<pre>
<code>
FoamFile {
    version 2.0;
    format ascii;
    class dictionary;
    location "system";
    object noiseDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

noiseModel volumeNoise;

volumeNoiseCoeffs {

    filteredFieldFunctionObjectNames ("filteredPressure");

    // Pressure column name, default = p
    p p;

    // Write interval for FFT data, default = 1
    fftWriteInterval 1;

    // Volume-weighted averaging switch, default = no (ensemble) for backwards
    // compatibility
    volumeAverage yes;


    /****** noiseModel coefficients *****/
    // Pressure reference
    pRef 0;

    // Number of samples in sampling window, default = 2^16 (=65536)
    N 64;

    // Lower frequency bounds
    fl 0;

    // Upper frequency bounds
    fu 1000;

    // Start time
    startTime 1;

    // SPL weighting; default = none
    SPLweighting dBA;

    dBRef 2e-5;

    windowModel Hanning;

    HanningCoeffs {
        // Window overlap percentage
        overlapPercent 50;
        symmetric yes;
        extended yes;

        // Optional number of windows, default = all available
        // nWindow 5;
    }

    // Collate times for ensight output - ensures geometry is only written once
    writeOptions {
        // Write Prmsf; default = yes
        writePrmsf no;
        // Write SPL; default = yes
        writeSPL yes;
        // Write PSD; default = yes
        writePSD yes;
        // Write PSDf; default = yes
        writePSDf no;
        // Write writeOctaves; default = yes
        writeOctaves no;
    }
}

</code>
</pre>

After setting up your `noiseDict` you run `noise` application that OpenFOAM provides to calculate noise levels based on your simulation data.
<pre>
<code>mpirun -n 16 noise -parallel</code>
</pre>
 Simulation data is saved by `filteredVolumeField` function object in csv format. It contains local and global cell indices, centroids, volumes, and values of the field at cells where filtering conditions are met.

 You can collect processor FFT data using `reconstructVolumeNoise` tool.
