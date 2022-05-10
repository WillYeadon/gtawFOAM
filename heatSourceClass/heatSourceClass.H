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

\*---------------------------------------------------------------------------*/

#ifndef heatSourceClass_H
#define heatSourceClass_H

#include "uniformDimensionedFields.H"
#include "dimensionedTypes.H"
#include "Time.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include "fvMesh.H"

namespace Foam
{

class heatSourceClass
{

        const dictionary& sourcePropertiesDict;
        const volScalarField& alphaMetalTotal_;
        const volVectorField& U_;
        const surfaceScalarField& phi_;

        dimensionedVector sourceVel;
        dimensionedVector sourcePos0;
        dimensionedVector sourcePos1;
        dimensionedVector sourcePos2;
        dimensionedScalar sourceOmega;
        dimensionedScalar sourcePenn;
        dimensionedScalar sourceModOmega;
        dimensionedScalar sourceModPenn;
        dimensionedScalar sourceModFraction;
        dimensionedScalar sourceStartTime;
        dimensionedScalar sourcePauseTime;
        dimensionedScalar sourceEndTime;
        dimensionedScalar source3EndTime;
        dimensionedScalar sourcePower;

        dimensionedScalar sourceC_C;
        dimensionedScalar sourceFieldCut;

        bool sourcePause;
        bool sourceMod;
        bool sourceOn;        
        bool sourceTwo;
        bool sourceThree;

        volScalarField sourceFieldSize;
//        volScalarField sourceModSize;
        volScalarField sourceFieldSizeTwo;
        volScalarField sourceFieldSizeThree;
        volScalarField sourceFieldOutput;

public:

	heatSourceClass(
	        const dictionary& sourceProperties,
            const volScalarField& alphaMetalTotal,
	        const volVectorField& U,
	        const surfaceScalarField& phi
		);

	virtual ~heatSourceClass(){}

    virtual tmp<volScalarField> geoField(const dimensionedVector sourcePos, const dimensionedScalar omega, const dimensionedScalar penn);
    virtual tmp<volScalarField> genField(const volScalarField& alphaMetalTotal, const dimensionedVector sourcePos, const dimensionedScalar omega, const dimensionedScalar penn);
    virtual tmp<volScalarField> normField(const volScalarField& sourceField);
    virtual tmp<volScalarField> binField(const volScalarField& sourceField);
    virtual tmp<volScalarField> floorField(const volScalarField& sourceField);
    virtual tmp<volScalarField> applyField(const volScalarField& alphaMetalTotal);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //