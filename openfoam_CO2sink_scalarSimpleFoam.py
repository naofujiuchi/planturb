#%%
import os
import zipfile
import shutil

dir_name = "tomato_case_CO2sink"
stl_name = "tomato_plant.stl"
stl_basename = os.path.splitext(stl_name)[0]

# Wind speed [m/s]
U_wind = 0.5

# Calculate gradient of CO2 concentration at the object surface
# gradCO2 [µmolCO2 molAir-1 m-1] = Pn [µmolCO2 m-2 s-1] / molar density of Air [molAir m-3] / Diffusivity of CO2 [m2 s-1]
Pn = 10 # µmolCO2 m-2 s-1
# molar density of Air [molAir m-3] 
# 20℃, 1 atm
# 41.5 mol/m3
# 25℃, 1 atm
# 40.9 mol/m3
P = 101325      # Pa
R = 8.314       # J/(mol·K)
T = 298.15      # K (25°C)
MD = P / (R * T) # mol/m³
# Diffusivity of CO2 [m2 s-1]
# 20℃, 1 atm
# 1.4e-05 m2/s
# 25℃, 1 atm
# 1.6e-05 m2/s
D_CO2 = 1.6e-05
gradCO2 = - Pn / MD / D_CO2
gradCO2

# Kinematic viscosity [m2/s]
# 20℃, 1 atm
# 1.5e-05 m2/s
# 25℃, 1 atm
# 1.8e-05 m2/s
nu = 1.8e-05

# Define the directory structure
base_dir = os.path.expanduser(f"~/Python/openfoam/{dir_name}")
dirs = [
    "0",
    "constant/triSurface",
    "constant/polyMesh",
    "system"
]

# Create the directories
for d in dirs:
    os.makedirs(os.path.join(base_dir, d), exist_ok=True)

# Header template
header_template = lambda obj: f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version: 10                                     |
|   \\  /    A nd           | Website:  https://openfoam.org                  |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/

FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""

# File contents
blockMeshDict = """
convertToMeters 1.0;

vertices (
    (0 0 0)
    (1.8 0 0)
    (1.8 1.8 0)
    (0 1.8 0)
    (0 0 2.0)
    (1.8 0 2.0)
    (1.8 1.8 2.0)
    (0 1.8 2.0)
);

blocks (
    hex (0 1 2 3 4 5 6 7) (36 36 40) simpleGrading (1 1 1)
);

edges 
(
);

boundary
(
    inlet
    {
        type patch;
        faces
        (
            (0 4 7 3)
        );
    }
    outlet
    {
        type patch;
        faces
        (
            (1 2 6 5)
        );
    }
    top
    {
        type wall;
        faces
        (
            (3 7 6 2)
        );
    }
    bottom
    {
        type wall;
        faces
        (
            (0 1 5 4)
        );
    }
    sides
    {
        type wall;
        faces
        (
            (0 3 2 1)
            (4 5 6 7)
        );
    }
);

mergePatchPairs 
(
);
"""

snappyHexMeshDict_1 = """
castellatedMesh true;
snap            true;
addLayers       false;

geometry
{{
    "{stl_name}"
    {{
        type triSurfaceMesh;
        name tomato;
        scale 1.0;
    }}
}}
""".format(stl_name=stl_name)

snappyHexMeshDict_2 = """
castellatedMeshControls
{{
    maxLocalCells 1000000;
    maxGlobalCells 2000000;
    minRefinementCells 10;
    nCellsBetweenLevels 3;

    features 
    (
        {{
            file "{stl_basename}.eMesh";
            level 2;
        }}
    );

    refinementSurfaces
    {{
        tomato
        {{
            level (2 2);
        }}
    }}

    refinementRegions {{}}

    locationInMesh (0.1 0.1 0.1);

    allowFreeStandingZoneFaces true;

    resolveFeatureAngle 30;
}}
""".format(stl_basename=stl_basename)

snappyHexMeshDict_3 = """
snapControls
{
    nSmoothPatch 3;
    tolerance 2.0;
    nSolveIter 30;
    nRelaxIter 5;
}

addLayersControls
{
    relativeSizes true;
    layers
    {
        tomato
        {
            nSurfaceLayers 1;
        }
    }
    expansionRatio 1.0;
    finalLayerThickness 0.3;
    minThickness 0.1;
    nGrow 0;
    featureAngle 60;
    slipFeatureAngle 30;
    nRelaxIter 5;
    nSmoothSurfaceNormals 1;
    nSmoothNormals 3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedianAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter 50;
}

meshQualityControls 
{
    maxNonOrtho 65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave 80;
    minFlatness 0.5;
    minVolRatio 0.01;
    minTetQuality 1e-9;
    minAreaFace 1e-13;
    minTwist 0.05;
    minDeterminant 0.001;
    minFaceWeight 0.05;
    minVolWeight 0.05;
    minTriangleTwist -1;
    minEdgeLength 1e-9;
    nSmoothScale 4;
    errorReduction 0.75;
    minVol 1e-13;
}

mergeTolerance 1e-6;
"""

snappyHexMeshDict = snappyHexMeshDict_1 + snappyHexMeshDict_2 + snappyHexMeshDict_3

surfaceFeaturesDict = """
surfaceFeaturesDict
{{
    surfaces
    (
        "{stl_name}"
    );

    includedAngle 150;
}}
""".format(stl_name=stl_name)

U_file_1 = """
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform ({U_wind} 0 0);

boundaryField
{{
    inlet
    {{
        type fixedValue;
        value uniform ({U_wind} 0 0);
    }}
""".format(U_wind=U_wind)

U_file_2 = """
    outlet
    {
        type zeroGradient;
    }
    top
    {
        type noSlip;
    }
    bottom
    {
        type noSlip;
    }
    sides
    {
        type noSlip;
    }
    tomato
    {
        type noSlip;
    }
}
"""

U_file = U_file_1 + U_file_2

p_file = """
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0;

boundaryField
{
    inlet
    {
        type zeroGradient;
    }
    outlet
    {
        type fixedValue;
        value uniform 0;
    }
    top
    {
        type zeroGradient;
    }
    bottom
    {
        type zeroGradient;
    }
    sides
    {
        type zeroGradient;
    }
    tomato
    {
        type zeroGradient;
    }
}
"""

CO2_file_1 = """
dimensions      [0 0 0 0 0 0 0];
internalField   uniform 400;

boundaryField
{
    inlet
    {
        type fixedValue;
        value uniform 400;
    }
    outlet
    {
        type zeroGradient;
    }
    top
    {
        type zeroGradient;
    }
    bottom
    {
        type zeroGradient;
    }
    sides
    {
        type zeroGradient;
    }
"""

CO2_file_2 = """
    tomato
    {{
        type fixedGradient;
        gradient uniform {gradCO2};  // CO2 sink: negative value. Adjust the value to represent the sink strength
        // type fixedValue;
        // value uniform 200;
    }}
}}
""".format(gradCO2=gradCO2)

CO2_file = CO2_file_1 + CO2_file_2

T_file = """
dimensions      [0 0 0 1 0 0 0];
internalField   uniform 300;

boundaryField
{
    inlet
    {
        type fixedValue;
        value uniform 300;
    }
    outlet
    {
        type zeroGradient;
    }
    top
    {
        type zeroGradient;
    }
    bottom
    {
        type zeroGradient;
    }
    sides
    {
        type zeroGradient;
    }
    tomato
    {
        type zeroGradient;
    }
}
"""

nut_file = """
dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    inlet
    {
        type fixedValue;
        value uniform 0;
    }
    outlet
    {
        type zeroGradient;
    }
    walls
    {
        type nutkWallFunction;
        value uniform 0;
    }
    top
    {
        type nutkWallFunction;
        value uniform 0;
    }
    bottom
    {
        type nutkWallFunction;
        value uniform 0;
    }
    sides
    {
        type nutkWallFunction;
        value uniform 0;
    }
    tomato
    {
        type nutkWallFunction;
        value uniform 0;
    }
}
"""

k_file = """
dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0.001;

boundaryField
{
    inlet
    {
        type fixedValue;
        value uniform 0.001;
    }
    outlet
    {
        type zeroGradient;
    }
    walls
    {
        type kqRWallFunction;
        value uniform 0.001;
    }
    top
    {
        type kqRWallFunction;
        value uniform 0.001;
    }
    bottom
    {
        type kqRWallFunction;
        value uniform 0.001;
    }
    sides
    {
        type kqRWallFunction;
        value uniform 0.001;
    }
    tomato
    {
        type kqRWallFunction;
        value uniform 0.001;
    }
}
"""

epsilon_file = """
dimensions      [0 2 -3 0 0 0 0];

internalField   uniform 0.001;

boundaryField
{
    inlet
    {
        type fixedValue;
        value uniform 0.001;
    }
    outlet
    {
        type zeroGradient;
    }
    walls
    {
        type epsilonWallFunction;
        value uniform 0.001;
    }
    top
    {
        type epsilonWallFunction;
        value uniform 0.001;
    }
    bottom
    {
        type epsilonWallFunction;
        value uniform 0.001;
    }
    sides
    {
        type epsilonWallFunction;
        value uniform 0.001;
    }
    tomato
    {
        type epsilonWallFunction;
        value uniform 0.001;
    }
}
"""

fvSchemes = """
ddtSchemes
{
    default         steadyState;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;
    div(phi,U)      Gauss upwind;
    div(phi,CO2)    Gauss upwind;
    div(phi,T)      Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
    div(phi,k)      Gauss upwind;
    div(phi,epsilon) Gauss upwind;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         orthogonal;
}

fluxRequired
{
    default         no;
    p               ;
}
"""

fvSolution = """
solvers
{
    p
    {
        solver          GAMG;
        tolerance       1e-06;
        relTol          0.01;
        smoother        DICGaussSeidel;
        cacheAgglomeration true;
        nCellsInCoarsestLevel 10;
        agglomerator    faceAreaPair;
        mergeLevels     1;
    }
    U
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-05;
        relTol          0;
    }
    CO2
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-05;
        relTol          0;
    }
    T
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-05;
        relTol          0;
    }
    k
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-05;
        relTol          0;
    }
    epsilon
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-05;
        relTol          0;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 0;
    residualControl
    {
        p               1e-2;
        U               1e-3;
        "(k|epsilon|omega)" 1e-3;
    }
}

relaxationFactors
{
    fields
    {
        p               0.3;
    }
    equations
    {
        U               0.7;
        k               0.7;
        epsilon         0.7;
    }
}
"""

transportProperties = """
D_CO2 D_CO2 [0 2 -1 0 0 0 0] {D_CO2}; // Diffusivity of CO2 (m^2/s)
nu nu [0 2 -1 0 0 0 0] {nu}; // Kinematic viscosity (m^2/s)
""".format(D_CO2=D_CO2, nu=nu)


physicalProperties = """

transportModel  Newtonian;
nu nu [0 2 -1 0 0 0 0] {nu}; // Kinematic viscosity (m^2/s)
rho             rho [ 1 -3 0 0 0 0 0 ] 1.225;   // Density (kg/m^3)
D_CO2 D_CO2 [0 2 -1 0 0 0 0] {D_CO2}; // Diffusivity of CO2 (m^2/s)
DT              DT [ 0 2 -1 0 0 0 0 ] 1.0e-05;  // Thermal diffusivity (m^2/s)
""".format(D_CO2=D_CO2, nu=nu)

momentumTransport = """

simulationType  RAS;
transportModel  Newtonian;

nu              nu [ 0 2 -1 0 0 0 0 ] {nu};  // Kinematic viscosity (m^2/s)
rho             rho [ 1 -3 0 0 0 0 0 ] 1.225;   // Density (kg/m^3)

RAS
{
    RASModel        kEpsilon;
    turbulence      on;
    printCoeffs     on;
}
"""

turbulenceProperties = """

simulationType  RAS;

RAS
{
    RASModel        kEpsilon;
    turbulence      on;
    printCoeffs     on;
}
"""

controlDict = """
application     scalarSimpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         1000;
deltaT          1;
writeControl    timeStep;
writeInterval   100;
purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;

functions
{
    averageCO2
    {
        type            surfaceFieldValue;
        functionObjectLibs ("libfieldFunctionObjects.so");
        enabled         true;
        writeControl    timeStep;
        writeInterval   1;
        regionType      patch;
        name            tomato;
        operation       areaAverage;
        fields
        (
            CO2
        );
        writeFields     yes;
        surfaceFormat   raw;
    }
//    totalCO2Gradient
//    {
//        type            surfaceFieldValue;
//        functionObjectLibs ("libfieldFunctionObjects.so");
//        enabled         true;
//        writeControl    timeStep;
//        writeInterval   1;
//        regionType      patch;
//        name            tomato;
//        operation       areaIntegrate;
//        fields
//        (
//            "snGrad(CO2)"
//        );
//        surfaceFormat   raw;
//    }
}
"""

# Write files to the directory

files = {
    "0/U": U_file,
    "0/p": p_file,
    "0/CO2": CO2_file,
    "0/nut": nut_file,
    "0/k": k_file,
    "0/epsilon": epsilon_file,
    "0/T": T_file,
    "system/controlDict": controlDict,
    "system/fvSchemes": fvSchemes,
    "system/fvSolution": fvSolution,
    "system/blockMeshDict": blockMeshDict,
    "system/snappyHexMeshDict": snappyHexMeshDict,
    "system/surfaceFeaturesDict": surfaceFeaturesDict,
    "constant/transportProperties": transportProperties,
    "constant/physicalProperties": physicalProperties,
    "constant/momentumTransport": momentumTransport,
    "constant/turbulenceProperties": turbulenceProperties,
}

for path, content in files.items():
    with open(os.path.join(base_dir, path), "w") as f:
        f.write(header_template(path.split("/")[-1]) + content.strip())

# Copy stl file to constant/triSurface
os.makedirs(os.path.join(base_dir, "constant/triSurface"), exist_ok=True)
shutil.copy(os.path.join(os.getcwd(), stl_name), os.path.join(base_dir, "constant/triSurface", stl_name))

# Zip the case directory
zip_path = f"{dir_name}.zip"
with zipfile.ZipFile(zip_path, 'w') as zipf:
    for root, _, files in os.walk(base_dir):
        for file in files:
            file_path = os.path.join(root, file)
            arcname = os.path.relpath(file_path, base_dir)
            zipf.write(file_path, arcname)

zip_path