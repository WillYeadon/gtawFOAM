/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

	Class to create up to three heat sources with a elliptic paraboloid shape
	Modified from http://dx.doi.org/10.1016/j.ijheatmasstransfer.2016.04.064

\*---------------------------------------------------------------------------*/

#include "heatSourceClass/heatSourceClass.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceFields.H"
#include "fvc.H"

Foam::heatSourceClass::heatSourceClass(
        const dictionary& sourceProperties,
        const volScalarField& alphaMetalTotal,
        const volVectorField& U,
        const surfaceScalarField& phi
):

    sourcePropertiesDict(sourceProperties),
    alphaMetalTotal_(alphaMetalTotal),
    U_(U),
    phi_(phi),
    
    sourceVel("sourceVel", dimLength, sourcePropertiesDict.lookupOrDefault<vector>("sourceVel", Zero)),
    sourcePos0("sourcePos0", dimLength, sourcePropertiesDict.lookupOrDefault<vector>("sourcePos0", Zero)),
    sourcePos1("sourcePos1", dimLength, sourcePropertiesDict.lookupOrDefault<vector>("sourcePos1", Zero)),
    sourcePos2("sourcePos2", dimLength, sourcePropertiesDict.lookupOrDefault<vector>("sourcePos2", Zero)),
    sourceOmega("sourceOmega", dimLength, sourcePropertiesDict.lookupOrDefault<scalar>("sourceOmega", 1e-3)),
    sourcePenn("sourcePenn", dimLength, sourcePropertiesDict.lookupOrDefault<scalar>("sourcePenn", 1e-3)),
    sourceModOmega("sourceModOmega", dimLength, sourcePropertiesDict.lookupOrDefault<scalar>("sourceModOmega", 1e-3)),
    sourceModPenn("sourceModPenn", dimLength, sourcePropertiesDict.lookupOrDefault<scalar>("sourceModPenn", 1e-3)),
    sourceModFraction("sourceModFraction", dimless, sourcePropertiesDict.lookupOrDefault<scalar>("sourceModFraction", 0.1)),
    sourceStartTime("sourceStartTime", dimTime, sourcePropertiesDict.lookupOrDefault<scalar>("sourceStartTime", 0.0)),
    sourcePauseTime("sourcePauseTime", dimTime, sourcePropertiesDict.lookupOrDefault<scalar>("sourcePauseTime", 0.0)),
    sourceEndTime("sourceEndTime", dimTime, sourcePropertiesDict.lookupOrDefault<scalar>("sourceEndTime", 1e3)),
    source3EndTime("source3EndTime", dimTime, sourcePropertiesDict.lookupOrDefault<scalar>("source3EndTime", 1e3)),
    sourcePower("sourcePower", dimLength, sourcePropertiesDict.lookupOrDefault<scalar>("sourcePower", 0.0)),
    
    sourceC_C("sourceC_C", dimless, sourcePropertiesDict.lookupOrDefault<scalar>("sourceC_C", 0.0)),
    sourceFieldCut("sourceFieldCut", dimless, sourcePropertiesDict.lookupOrDefault<scalar>("sourceFieldCut", 0.01)),

    sourcePause(sourcePropertiesDict.lookupOrDefault<bool>("sourcePause", true)),  
    sourceMod(sourcePropertiesDict.lookupOrDefault<bool>("sourceMod", false)),  
    sourceOn(sourcePropertiesDict.lookupOrDefault<bool>("sourceOn", true)),
    sourceTwo(sourcePropertiesDict.lookupOrDefault<bool>("sourceTwo", false)),
    sourceThree(sourcePropertiesDict.lookupOrDefault<bool>("sourceThree", false)),  

    sourceFieldSize
    (
        IOobject
        (
            "sourceFieldSize",
            //runTime.timeName(),
            U_.mesh().time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("sourceFieldSize", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    sourceFieldSizeTwo
    (
        IOobject
        (
            "sourceFieldSizeTwo",
            //runTime.timeName(),
            U_.mesh().time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("sourceFieldSize", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    sourceFieldSizeThree
    (
        IOobject
        (
            "sourceFieldSizeThree",
            //runTime.timeName(),
            U_.mesh().time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("sourceFieldSize", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    sourceFieldOutput
    (
        IOobject
        (
            "sourceFieldOutput",
            U_.mesh().time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("sourceFieldOutput", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    )
{
  
}

Foam::tmp<Foam::volScalarField>
Foam::heatSourceClass::geoField(const dimensionedVector sourcePos, const dimensionedScalar omega, const dimensionedScalar penn)
{
    const volVectorField& cellCentre = U_.mesh().C();

    scalar C_d = log(sourceFieldCut.value());

    dimensionedScalar oneMeter
    (
        "oneMeter",
        dimensionSet(0,1,0,0,0,0,0),
        1.0
    );

    scalar currentTime = U_.mesh().time().value();

    bool sourcePauseSwitch;
    (currentTime < sourcePauseTime.value()) ? sourcePauseSwitch = true : sourcePauseSwitch = false;

    volScalarField expX
    (
        exp(C_d*(sqr((cellCentre.component(vector::X)-(sourcePos.component(vector::X))
        + (sourceVel.component(vector::X)*(currentTime  - sourcePauseTime.value())))/oneMeter)/sqr(omega/oneMeter)))
    );

    volScalarField expY
    (
        exp(C_d*(sqr((cellCentre.component(vector::Y)-(sourcePos.component(vector::Y))
        + (sourceVel.component(vector::Y)*(currentTime  - sourcePauseTime.value())))/oneMeter)/sqr(omega/oneMeter)))
    ); 

    volScalarField expZ
    (
      exp(C_d*sqrt(sqr((cellCentre.component(vector::Z) - (sourcePos.component(vector::Z)))/oneMeter)) / (penn/oneMeter))
    );


    if (sourcePause && sourcePauseSwitch)      
    {
      expX = exp(C_d*(sqr((cellCentre.component(vector::X)-(sourcePos.component(vector::X)))/oneMeter)/sqr(omega/oneMeter)));
      expY = exp(C_d*(sqr((cellCentre.component(vector::Y)-(sourcePos.component(vector::Y)))/oneMeter)/sqr(omega/oneMeter)));
    }

    volScalarField expTot
    (
      expX*expY*expZ
    );

  	return tmp<volScalarField>
    (
        new volScalarField
        (
          "geoField",
            expTot
    )
  );  
}


Foam::tmp<Foam::volScalarField>
Foam::heatSourceClass::genField(const volScalarField& alphaMetalTotal, const dimensionedVector sourcePos, const dimensionedScalar omega, const dimensionedScalar penn)
{
	volScalarField genField
	(
		geoField(sourcePos, omega, penn)*alphaMetalTotal
	);

    forAll(genField, i)
    {
        ((genField[i] >= sourceFieldCut.value()) && (alphaMetalTotal[i] >= sourceFieldCut.value())) ? genField[i] = scalar(1) : genField[i] = scalar(0);
    }

  	return tmp<volScalarField>
    (
        new volScalarField
        (
          "genField",
			genField
	    )
  );  

}

Foam::tmp<Foam::volScalarField>
Foam::heatSourceClass::normField(const volScalarField& sourceField)
{
	scalar fieldMax = max(sourceField).value();	
	
	volScalarField sourceFieldNorm
	(
		sourceField/fieldMax
	);

  	return tmp<volScalarField>
    (
        new volScalarField
        (
          "sourceFieldNorm",
			sourceFieldNorm
    	)
  	);
}

Foam::tmp<Foam::volScalarField>
Foam::heatSourceClass::binField(const volScalarField& sourceField)
{
    volScalarField sourceFieldBin(sourceField);

    forAll(sourceFieldBin, i)
    {
        if (sourceFieldBin[i] >= 1.0)
            sourceFieldBin[i] = 1.0;
    }

    return tmp<volScalarField>
    (
        new volScalarField
        (
          "sourceFieldBin",
            sourceFieldBin
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::heatSourceClass::floorField(const volScalarField& sourceField)
{
    volScalarField sourceFieldFloor(sourceField);

    forAll(sourceFieldFloor, i)
    {
        if (sourceFieldFloor[i] <= 1.0)
            sourceFieldFloor[i] = 0.0;
    }

    return tmp<volScalarField>
    (
        new volScalarField
        (
          "sourceFieldFloor",
            sourceFieldFloor
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::heatSourceClass::applyField(const volScalarField& alphaMetalTotal)
{
    dimensionedScalar oneMeter
    (
        "oneMeter",
        dimensionSet(0,1,0,0,0,0,0),
        1.0
    );

    Info << " Q_arc " << endl;
    scalar Q_arc = sourcePower.value();
    Info << " Generating Heat Source Field " << endl;
    sourceFieldSize = 0.0;
    sourceFieldSizeTwo = 0.0;
    sourceFieldSizeThree = 0.0;

    volScalarField sourceModField0 = genField(alphaMetalTotal, sourcePos0, sourceModOmega, sourceModPenn);
    volScalarField sourceModField1 = genField(alphaMetalTotal, sourcePos1, sourceModOmega, sourceModPenn);
    volScalarField sourceModField2 = genField(alphaMetalTotal, sourcePos2, sourceModOmega, sourceModPenn);

    /*
    Source Mod allows you to overlay two fields to make a sombrero shape 
    */

    if (!sourceMod)
    {
        sourceModField0 *= 0.0;
        sourceModField1 *= 0.0;
        sourceModField2 *= 0.0;
    }

    if (sourceOn && U_.mesh().time() < sourceEndTime) 
    {
    	sourceFieldSize = binField(genField(alphaMetalTotal, sourcePos0, sourceOmega, sourcePenn) + sourceModField0);
    }
    if (sourceTwo && (U_.mesh().time() < sourceEndTime)) 
    {
    	sourceFieldSizeTwo = binField(genField(alphaMetalTotal, sourcePos1, sourceOmega, sourcePenn) + sourceModField1);
    }
    if (sourceThree && (U_.mesh().time() < sourceEndTime) && (U_.mesh().time() < source3EndTime)) 
    {
    	sourceFieldSizeThree = binField(genField(alphaMetalTotal, sourcePos2, sourceOmega, sourcePenn) + sourceModField2);
    }

    if (U_.mesh().time() > sourceEndTime)
    {
      sourceFieldSize = 0.0;
      sourceFieldSizeTwo = 0.0;
      sourceFieldSizeThree = 0.0;
    }

    scalar heatSourceFieldMax = max(sourceFieldSize).value();
    scalar heatSourceFieldMin = min(sourceFieldSize).value();
    Info << "heatSourceField max: " << heatSourceFieldMax 
    << " min: " << heatSourceFieldMin << endl;

    Info << " delHField " << endl;
    volScalarField delHField = ((sourceC_C * Q_arc) / (Foam::sqr(sourceOmega/oneMeter)*Foam::sqrt(sourcePenn/oneMeter)))*alphaMetalTotal;

    volScalarField source0 = normField(geoField(sourcePos0, sourceOmega, sourcePenn)); 
	volScalarField source1 = normField(geoField(sourcePos1, sourceOmega, sourcePenn)); 
    volScalarField source2 = normField(geoField(sourcePos2, sourceOmega, sourcePenn)); 
    
    if (sourceMod)
    {
        source0 = normField(geoField(sourcePos0, sourceOmega, sourcePenn) + sourceModFraction*geoField(sourcePos0, sourceModOmega, sourceModPenn));
        source1 = normField(geoField(sourcePos1, sourceOmega, sourcePenn) + sourceModFraction*geoField(sourcePos0, sourceModOmega, sourceModPenn));
        source2 = normField(geoField(sourcePos2, sourceOmega, sourcePenn) + sourceModFraction*geoField(sourcePos0, sourceModOmega, sourceModPenn));
    }

    scalar source0Max = max(source0).value();
    scalar source0Min = min(source0).value();
    Info <<  "source0 max: " << source0Max << " min: " << source0Min << endl;

    sourceFieldOutput = (sourceFieldSize*source0 + sourceFieldSizeTwo*source1 + sourceFieldSizeThree*source2)*delHField;

    return sourceFieldOutput;
}