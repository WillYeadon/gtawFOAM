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

\*---------------------------------------------------------------------------*/

#include "incompressibleThreePhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::incompressibleThreePhaseMixture::calcNu()
{
    nuModel1_->correct();
    nuModel2_->correct();
    nuModel3_->correct();

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu()/(alpha1_*rho1_ + alpha2_*rho2_ + alpha3_*rho3_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleThreePhaseMixture::incompressibleThreePhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    phase1Name_(wordList(lookup("phases"))[0]),
    phase2Name_(wordList(lookup("phases"))[1]),
    phase3Name_(wordList(lookup("phases"))[2]),

    alpha1_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase1Name_),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    alpha2_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase2Name_),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    alpha3_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase3Name_),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    T_
    (
        IOobject
        (
            "T",
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    U_(U),
    phi_(phi),

    nu_
    (
        IOobject
        (
            "nu",
            U.time().timeName(),
            U.db()
        ),
        U.mesh(),
        dimensionedScalar(dimensionSet(0, 2, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    ),

    nuModel1_
    (
        viscosityModel::New
        (
            "nu1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    nuModel2_
    (
        viscosityModel::New
        (
            "nu2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),
    nuModel3_
    (
        viscosityModel::New
        (
            "nu3",
            subDict(phase3Name_),
            U,
            phi
        )
    ),

    rho1_("rho", dimDensity, nuModel1_->viscosityProperties()),

    rho2_("rho", dimDensity, nuModel2_->viscosityProperties()),

    rho2Opt_(nuModel2_->viscosityProperties().lookupOrDefault("rho2Opt", false)),
    rho2a_("2rho2a", dimensionSet(1,-3,0,-1,0,0,0), nuModel2_->viscosityProperties().lookupOrDefault("rho2a", 0.0)),
    rho2b_("2rho2b", dimensionSet(1,-3,0,-2,0,0,0), nuModel2_->viscosityProperties().lookupOrDefault("rho2b", 0.0)),
    rho2c_("2rho2c", dimensionSet(1,-3,0,-3,0,0,0), nuModel2_->viscosityProperties().lookupOrDefault("rho2c", 0.0)),
    rho2d_("2rho2d", dimensionSet(1,-3,0,-4,0,0,0), nuModel2_->viscosityProperties().lookupOrDefault("rho2d", 0.0)),

    rho3_("rho", dimDensity, nuModel3_->viscosityProperties()),

    Tref_("Tref", dimensionSet(0,0,0,1,0,0,0), nuModel2_->viscosityProperties()),
    meltingT_("Tref", dimensionSet(0,0,0,1,0,0,0), nuModel2_->viscosityProperties()),
    Lf_("Lf", dimensionSet(0,2,-2,0,0,0,0), nuModel2_->viscosityProperties()),
    beta_("beta", dimensionSet(0,0,0,-1,0,0,0), nuModel2_->viscosityProperties()),
    CTLE_("CTLE", dimensionSet(0,0,0,-1,0,0,0), nuModel3_->viscosityProperties().lookupOrDefault("CTLE", 0.0)),

    Tdep1_(nuModel1_->viscosityProperties().lookupOrDefault("Tdep", false)),
    Tdep2_(nuModel2_->viscosityProperties().lookupOrDefault("Tdep", false)),
    Tdep3_(nuModel3_->viscosityProperties().lookupOrDefault("Tdep", false)),
    Shomate_(nuModel2_->viscosityProperties().lookupOrDefault("Shomate", false)),

    cp1_("cp", dimensionSet(0,2,-2,-1,0,0,0), nuModel1_->viscosityProperties()),

    cp2_("cp", dimensionSet(0,2,-2,-1,0,0,0), nuModel2_->viscosityProperties()),
    cp2a_("2a", dimensionSet(0,2,-2,-1,0,0,0), nuModel2_->viscosityProperties().lookupOrDefault("a", 1.0)),
    cp2b_("2b", dimensionSet(0,2,-2,-2,0,0,0), nuModel2_->viscosityProperties().lookupOrDefault("b", 1.0)),
    cp2c_("2c", dimensionSet(0,2,-2,1,0,0,0), nuModel2_->viscosityProperties().lookupOrDefault("c", 1.0)),
    cp2d_("2d", dimensionSet(0,2,-2,-3,0,0,0), nuModel2_->viscosityProperties().lookupOrDefault("d", 0.0)),
    cp2e_("2e", dimensionSet(0,2,-2,-4,0,0,0), nuModel2_->viscosityProperties().lookupOrDefault("e", 0.0)), 
    cp2f_("2f", dimensionSet(0,2,-2,-5,0,0,0), nuModel2_->viscosityProperties().lookupOrDefault("f", 0.0)), 

    cp3_("cp", dimensionSet(0,2,-2,-1,0,0,0), nuModel3_->viscosityProperties()),
    cp3a_("3a", dimensionSet(0,2,-2,-1,0,0,0), nuModel3_->viscosityProperties().lookupOrDefault("a", 1.0)),
    cp3b_("3b", dimensionSet(0,2,-2,-2,0,0,0), nuModel3_->viscosityProperties().lookupOrDefault("b", 1.0)),
    cp3c_("3c", dimensionSet(0,2,-2,1,0,0,0), nuModel3_->viscosityProperties().lookupOrDefault("c", 1.0)),
    cp3d_("3d", dimensionSet(0,2,-2,-3,0,0,0), nuModel3_->viscosityProperties().lookupOrDefault("d", 0.0)),
    cp3e_("3e", dimensionSet(0,2,-2,-4,0,0,0), nuModel3_->viscosityProperties().lookupOrDefault("e", 0.0)), 
    cp3f_("3f", dimensionSet(0,2,-2,-5,0,0,0), nuModel2_->viscosityProperties().lookupOrDefault("f", 0.0)),

    smallT_("smallT", dimensionSet(0,0,0,1,0,0,0), 1e-6),
    smallNu_("smallNu", dimensionSet(1,-3,0,0,0,0,0), 1e-6),

    k1_("k", dimensionSet(1,1,-3,-1,0,0,0), nuModel1_->viscosityProperties()),
    k2_("k", dimensionSet(1,1,-3,-1,0,0,0), nuModel2_->viscosityProperties()),
    k2Opt_(nuModel2_->viscosityProperties().lookupOrDefault("k2Opt", false)),
    k2a_("k2a", dimensionSet(1,1,-3,-1,0,0,0), nuModel2_->viscosityProperties().lookupOrDefault("k2a", 1.0)),
    k2b_("k2b", dimensionSet(1,1,-3,-2,0,0,0), nuModel2_->viscosityProperties().lookupOrDefault("k2b", 0.0)),
    k2c_("k2c", dimensionSet(1,1,-3,-3,0,0,0), nuModel2_->viscosityProperties().lookupOrDefault("k2c", 0.0)),

    k3_("k", dimensionSet(1,1,-3,-1,0,0,0), nuModel2_->viscosityProperties()),
    k3Opt_(nuModel2_->viscosityProperties().lookupOrDefault("k3Opt", false)),
    k3a_("k3a", dimensionSet(1,1,-3,-1,0,0,0), nuModel2_->viscosityProperties().lookupOrDefault("k3a", 1.0)),
    k3b_("k3b", dimensionSet(1,1,-3,-2,0,0,0), nuModel2_->viscosityProperties().lookupOrDefault("k3b", 0.0)),
    k3c_("k3c", dimensionSet(1,1,-3,-3,0,0,0), nuModel2_->viscosityProperties().lookupOrDefault("k3c", 0.0))

{
    alpha3_ == 1.0 - alpha1_ - alpha2_;
    calcNu();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::mu() const
{
    return volScalarField::New
    (
        "mu",
        alpha1_*rho1_*nuModel1_->nu()
      + alpha2_*rho2_*nuModel2_->nu()
      + alpha3_*rho3_*nuModel3_->nu()
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleThreePhaseMixture::muf() const
{
    surfaceScalarField alpha1f(fvc::interpolate(alpha1_));
    surfaceScalarField alpha2f(fvc::interpolate(alpha2_));
    surfaceScalarField alpha3f(fvc::interpolate(alpha3_));

    return surfaceScalarField::New
    (
        "mu",
        alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
      + alpha2f*rho2_*fvc::interpolate(nuModel2_->nu())
      + alpha3f*rho3_*fvc::interpolate(nuModel3_->nu())
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleThreePhaseMixture::nuf() const
{
    surfaceScalarField alpha1f(fvc::interpolate(alpha1_));
    surfaceScalarField alpha2f(fvc::interpolate(alpha2_));
    surfaceScalarField alpha3f(fvc::interpolate(alpha3_));

    return surfaceScalarField::New
    (
        "nu",
        (
            alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
          + alpha2f*rho2_*fvc::interpolate(nuModel2_->nu())
          + alpha3f*rho3_*fvc::interpolate(nuModel3_->nu())
        )/(alpha1f*rho1_ + alpha2f*rho2_ + alpha3f*rho3_)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::k() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    const volScalarField limitedAlpha2
    (
        min(max(alpha2_, scalar(0)), scalar(1))
    );

    const volScalarField limitedAlpha3
    (
        min(max(alpha3_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "k",
            limitedAlpha1*k1_ + limitedAlpha2*k2_ + limitedAlpha3*k3_ 
//             limitedAlpha1*k1() + (k2Opt_ ? k2T() : limitedAlpha2*k2()) + (k3Opt_ ? k3T() : limitedAlpha3*k3()) 
        )
    );
}
Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleThreePhaseMixture::kf() const
{
    const surfaceScalarField alpha1f
    (
        fvc::interpolate(min(max(alpha1_, scalar(0)), scalar(1)))     
    );

    const surfaceScalarField alpha2f
    (
        fvc::interpolate(min(max(alpha2_, scalar(0)), scalar(1)))     
    );

    const surfaceScalarField alpha3f
    (
        fvc::interpolate(min(max(alpha3_, scalar(0)), scalar(1)))     
    );

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "kf",
            alpha1f*k1_
          + alpha2f*k2_
          + alpha3f*k3_
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::rho() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "rho",
            limitedAlpha1*rho1_
          + rhoBouss()
          + rhoExp()
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::rhoCp() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    if (!Tdep2_)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                "rhoCp",
               (limitedAlpha1*rho1_*cp1_
              + rhoBouss()*cp2_
              + rhoExp()*cp3_)
            )
        );       

    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                "rhoCp",
               (limitedAlpha1*cp1_*rho1_
              + rhoBouss()*cp2TDepVSF()
              + rhoExp()*cp3TDepVSF())
            )
        );       
    }
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::rhoExp() const
{
    const volScalarField limitedAlpha3
    (
        min(max(alpha3_, scalar(0)), scalar(1))
    );

    const volScalarField linearExp
    (
        max(T_, smallT_)*CTLE_
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "rhoExp",
            limitedAlpha3*rho3_*(1.0 / pow((1.0 + linearExp), 3))
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::rhoBouss() const
{
    const volScalarField secondTerm
    (
        rho2a_*pow(max(T_, smallT_) - Tref_, 1)
    );

    const volScalarField thirdTerm
    (
        rho2b_*pow(max(T_, smallT_) - Tref_, 2)
    );

    const volScalarField forthTerm
    (
        rho2c_*pow(max(T_, smallT_) - Tref_, 3)
    );

    const volScalarField fifthTerm
    (
        rho2d_*pow(max(T_, smallT_) - Tref_, 4)
    );

    const volScalarField limitedAlpha2
    (
        min(max(alpha2_, scalar(0)), scalar(1))
    );

    const volScalarField betaT
    (
        beta_*(max(T_, smallT_) - Tref_)
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "rhoBouss",
            (rho2Opt_? limitedAlpha2*(rho2_ + secondTerm + thirdTerm + forthTerm + fifthTerm) 
                : limitedAlpha2*rho2_*(scalar(1.0) - betaT))
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::rhoBoussRho() const
{
    const volScalarField secondTerm
    (
        rho2a_*pow(max(T_, smallT_) - Tref_, 1)
    );

    const volScalarField thirdTerm
    (
        rho2b_*pow(max(T_, smallT_) - Tref_, 2)
    );

    const volScalarField forthTerm
    (
        rho2c_*pow(max(T_, smallT_) - Tref_, 3)
    );

    const volScalarField fifthTerm
    (
        rho2d_*pow(max(T_, smallT_) - Tref_, 4)
    );

    const volScalarField limitedAlpha2
    (
        min(max(alpha2_, scalar(0)), scalar(1))
    );

    const volScalarField betaT
    (
        beta_*(max(T_, smallT_) - Tref_)
    );


    return tmp<volScalarField>
    (
        new volScalarField
        (
            "rhoBoussRho",
//                rho2_*betaT
            (rho2Opt_ ? (rho2_ + secondTerm + thirdTerm + forthTerm + fifthTerm)
                : rho2_*(scalar(1.0) - betaT))
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::cp2TDepVSF() const
{
    scalar shomateDiv = 1;

    if (Shomate_)
        shomateDiv = 1000;

    const volScalarField thirdTerm
    (
        cp2c_*pow(max(T_, smallT_)/shomateDiv, -2)
    );

    const volScalarField forthTerm
    (
        cp2d_*pow(max(T_, smallT_)/shomateDiv, 2)
    );

    const volScalarField fifthTerm
    (
        cp2e_*pow(max(T_, smallT_)/shomateDiv, 3)
    );

    const volScalarField sixthTerm
    (
        cp2f_*pow(max(T_, smallT_)/shomateDiv, 4)
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "cp2TDepVSF",
            cp2a_ + cp2b_*((T_)/shomateDiv) + thirdTerm + forthTerm + fifthTerm + sixthTerm
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::cp3TDepVSF() const
{
    scalar shomateDiv = 1;

    if (Shomate_)
        shomateDiv = 1000;

    const volScalarField thirdTerm
    (
        cp3c_*pow(max(T_, smallT_)/shomateDiv, -2)
    );

    const volScalarField forthTerm
    (
        cp3d_*pow(max(T_, smallT_)/shomateDiv, 2)
    );

    const volScalarField fifthTerm
    (
        cp3e_*pow(max(T_, smallT_)/shomateDiv, 3)
    );

    const volScalarField sixthTerm
    (
        cp2f_*pow(max(T_, smallT_)/shomateDiv, 4)
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "cp3TDepVSF",
            cp3a_ + cp3b_*((T_)/shomateDiv) + thirdTerm + forthTerm + fifthTerm + sixthTerm
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleThreePhaseMixture::cp() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );
    const volScalarField limitedAlpha2
    (
        min(max(alpha2_, scalar(0)), scalar(1))
    );
    const volScalarField limitedAlpha3
    (
        min(max(alpha3_, scalar(0)), scalar(1))
    );

    const volScalarField outputCp2
    (
        cp2TDepVSF()
    );

// alpha fields called twice for interface
    if (!Tdep2_)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                "cp",
               (limitedAlpha1*cp1_*rho1_
              + limitedAlpha2*rhoBouss()*cp2_
              + limitedAlpha3*rhoExp()*cp3_)
              /(rho1_
              + rhoBouss()
              + rhoExp())
            )
        );       

    }
    else
    {

        Info<< " T dep Cp Mixture triggered Phase 2:" << outputCp2().weightedAverage(U_.mesh().V()) << endl;

         return tmp<volScalarField>
        (
            new volScalarField
            (
                "cp",
               (limitedAlpha1*cp1_*rho1_
              + limitedAlpha2*rhoBouss()*cp2TDepVSF()
              + limitedAlpha3*rhoExp()*cp3TDepVSF())
              /(rho1_               
              + rhoBouss()                
              + rhoExp())
            )
        );       
    }

}


bool Foam::incompressibleThreePhaseMixture::read()
{
    if (regIOobject::read())
    {
        if
        (
            nuModel1_().read(*this)
         && nuModel2_().read(*this)
         && nuModel3_().read(*this)
        )
        {
            nuModel1_->viscosityProperties().lookup("rho") >> rho1_;
            nuModel2_->viscosityProperties().lookup("rho") >> rho2_;
            nuModel3_->viscosityProperties().lookup("rho") >> rho3_;

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
