/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Website:  www.openfoam.com                      |
|   \\  /    A nd           |                                                 |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     icoFoam;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         0.5;
deltaT          0.005;

// Deliberately less frequent than the function-object output.
writeControl    adjustableRunTime;
writeInterval   1;

purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;

functions
{
    numpyExport
    {
        type                foamToNumpy;
        libs                (numpyFunctionObjects);

        fields              (p U);
        dataType            float64;
        batchSize           2;

        writeTimes          true;
        writeCellCentres    true;
        writeCellVolumes    true;

        writeControl        adjustableRunTime;
        writeInterval       0.1;
    }
}

// ************************************************************************* //
