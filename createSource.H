#include "zeroGradientFvPatchFields.H"
#include "uniformDimensionedFields.H"
#include "dimensionedTypes.H"
#include "Time.H"

Info << "Reading Source Properties\n" << endl;
IOdictionary sourcePropertiesDict
(
    IOobject
    (
        "sourceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info<< "Making Source\n" << endl;
heatSourceClass arcSource(sourcePropertiesDict, alphaMetalTotal, U, phi);

dimensionedScalar qSourceUnits
(
	"qSourceUnits",
    dimensionSet(1,-1,-3,0,0,0,0),
    1.0
);

volScalarField qSource
(
    IOobject
    (
        "qSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    arcSource.applyField(alphaMetalTotal)*qSourceUnits,
    zeroGradientFvPatchScalarField::typeName
);
/*
Info << "Reading Source Properties\n" << endl;
IOdictionary velPropertiesDict
(
    IOobject
    (
        "velProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info<< "Making Source\n" << endl;
velSourceClass velSource(velPropertiesDict, U, phi);

dimensionedScalar velSourceUnits
(
    "velSourceUnits",
    dimensionSet(1,-3,-2,0,0,0,0),
    1.0
);

dimensionedVector velSourceZero
(
    "velSourceZero",
    dimensionSet(1,-2,-2,0,0,0,0),
    Zero
);

volVectorField velSourceZeroField
(
    IOobject
    (
        "velSourceZeroField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("velSourceZeroDM", dimensionSet(1,-2,-2,0,0,0,0), Zero)
);

volVectorField uSource
(
    IOobject
    (
        "uSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    velSource.applyField()*velSourceUnits,
    zeroGradientFvPatchScalarField::typeName
);


Info << "Reading Source Properties\n" << endl;
IOdictionary fillerPropertiesDict
(
    IOobject
    (
        "fillerProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info<< "Making filler source\n" << endl;
alphaSourceClass alphaSource(fillerPropertiesDict, U, phi);

dimensionedScalar alphaSourceUnits
(
    "alphaSourceUnits",
    dimensionSet(0,0,-1,0,0,0,0),
    1.0
);

volScalarField sFiller
(
    IOobject
    (
        "sFiller",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alphaSource.applyField()*alphaSourceUnits,
    zeroGradientFvPatchScalarField::typeName
);
*/
// AM sources