/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

Class
    Foam::incompressibleThreePhaseMixture

Description

SourceFiles
    incompressibleThreePhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleThreePhaseMixture_H
#define incompressibleThreePhaseMixture_H

#include "kinematicTransportModel.H"
#include "IOdictionary.H"
#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class incompressibleThreePhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class incompressibleThreePhaseMixture
:
    public IOdictionary,
    public kinematicTransportModel
{
    // Private Data

        word phase1Name_;
        word phase2Name_;
        word phase3Name_;

        volScalarField alpha1_;
        volScalarField alpha2_;
        volScalarField alpha3_;
        volScalarField T_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        volScalarField nu_;

        autoPtr<viscosityModel> nuModel1_;
        autoPtr<viscosityModel> nuModel2_;
        autoPtr<viscosityModel> nuModel3_;

        dimensionedScalar rho1_;
        dimensionedScalar rho2_;

        bool rho2Opt_;
        dimensionedScalar rho2a_;
        dimensionedScalar rho2b_;
        dimensionedScalar rho2c_;
        dimensionedScalar rho2d_;

        dimensionedScalar rho3_;

        dimensionedScalar Tref_;
        dimensionedScalar meltingT_;
        dimensionedScalar Lf_;
        dimensionedScalar beta_;

        dimensionedScalar CTLE_;

        bool Tdep1_;
        bool Tdep2_;
        bool Tdep3_;
        bool Shomate_;

        dimensionedScalar cp1_;

        dimensionedScalar cp2_;
        dimensionedScalar cp2a_;
        dimensionedScalar cp2b_;
        dimensionedScalar cp2c_;
        dimensionedScalar cp2d_;
        dimensionedScalar cp2e_;
        dimensionedScalar cp2f_;        

        dimensionedScalar cp3_;
        dimensionedScalar cp3a_;
        dimensionedScalar cp3b_;
        dimensionedScalar cp3c_;
        dimensionedScalar cp3d_;
        dimensionedScalar cp3e_;
        dimensionedScalar cp3f_;

        dimensionedScalar smallT_;
        dimensionedScalar smallNu_;

        dimensionedScalar k1_;

        dimensionedScalar k2_;
        bool k2Opt_;
        dimensionedScalar k2a_;
        dimensionedScalar k2b_;
        dimensionedScalar k2c_;                

        dimensionedScalar k3_;
        bool k3Opt_;
        dimensionedScalar k3a_;
        dimensionedScalar k3b_;
        dimensionedScalar k3c_;             

    // Private Member Functions

        //- Calculate and return the laminar viscosity
        void calcNu();


public:

    // Constructors

        //- Construct from components
        incompressibleThreePhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    ~incompressibleThreePhaseMixture()
    {}


    // Member Functions

        const word phase1Name() const
        {
            return phase1Name_;
        }

        const word phase2Name() const
        {
            return phase2Name_;
        }

        const word phase3Name() const
        {
            return phase3Name_;
        }

        const volScalarField& alpha1() const
        {
            return alpha1_;
        }

        volScalarField& alpha1()
        {
            return alpha1_;
        }

        const volScalarField& alpha2() const
        {
            return alpha2_;
        }

        volScalarField& alpha2()
        {
            return alpha2_;
        }

        const volScalarField& alpha3() const
        {
            return alpha3_;
        }

        volScalarField& alpha3()
        {
            return alpha3_;
        }

        const volScalarField& T() const
        {
            return T_;
        }

        volScalarField& T()
        {
            return T_;
        }

        //- Return const-access to phase1 density
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& rho2() const
        {
            return rho2_;
        };

        //- Return const-access to phase3 density
        const dimensionedScalar& rho3() const
        {
            return rho3_;
        };

        const dimensionedScalar& cp1() const
        {
            return cp1_;
        }

        const dimensionedScalar& cp2() const
        {
            return cp2_;
        };

        const dimensionedScalar& cp3() const
        {
            return cp3_;
        };

        const dimensionedScalar& meltingT() const
        {
            return meltingT_;
        };

        const dimensionedScalar& Lf() const
        {
            return Lf_;
        };

        //- Return the velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the flux
        const surfaceScalarField& phi() const
        {
            return phi_;
        }

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel1() const
        {
            return nuModel1_();
        }

        //- Return const-access to phase2 viscosityModel
        const viscosityModel& nuModel2() const
        {
            return nuModel2_();
        }

        //- Return const-access to phase3 viscosityModel
        const viscosityModel& nuModel3() const
        {
            return nuModel3_();
        }

        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;

        //- Return the kinematic laminar viscosity
        tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

        tmp<volScalarField> k() const;
        tmp<surfaceScalarField> kf() const;

        tmp<volScalarField> rhoExp() const;
        tmp<volScalarField> rhoBouss() const;      
        tmp<volScalarField> rhoBoussRho() const;
        
        tmp<volScalarField> cp2TDepVSF() const;
        tmp<volScalarField> cp3TDepVSF() const;
        tmp<volScalarField> cp() const;

        tmp<volScalarField> rho() const;
        tmp<volScalarField> rhoCp() const;
        //- Correct the laminar viscosity
        void correct()
        {
            calcNu();
        }

        //- Read base transportProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
