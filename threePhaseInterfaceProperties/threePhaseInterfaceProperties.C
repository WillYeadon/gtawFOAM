/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "threePhaseInterfaceProperties.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "unitConversion.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::threePhaseInterfaceProperties::convertToRad =
    Foam::constant::mathematical::pi/180.0;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::threePhaseInterfaceProperties::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb
) const
{
    const volScalarField::Boundary& alpha1 =
        mixture_.alpha1().boundaryField();
    const volScalarField::Boundary& alpha2 =
        mixture_.alpha2().boundaryField();
    const volScalarField::Boundary& alpha3 =
        mixture_.alpha3().boundaryField();
    const volVectorField::Boundary& U =
        mixture_.U().boundaryField();

    const fvMesh& mesh = mixture_.U().mesh();
    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(alpha1[patchi]))
        {
            const alphaContactAngleFvPatchScalarField& a2cap =
                refCast<const alphaContactAngleFvPatchScalarField>
                (alpha2[patchi]);

            const alphaContactAngleFvPatchScalarField& a3cap =
                refCast<const alphaContactAngleFvPatchScalarField>
                (alpha3[patchi]);

            scalarField twoPhaseAlpha2(max(a2cap, scalar(0)));
            scalarField twoPhaseAlpha3(max(a3cap, scalar(0)));

            scalarField sumTwoPhaseAlpha
            (
                twoPhaseAlpha2 + twoPhaseAlpha3 + small
            );

            twoPhaseAlpha2 /= sumTwoPhaseAlpha;
            twoPhaseAlpha3 /= sumTwoPhaseAlpha;

            fvsPatchVectorField& nHatp = nHatb[patchi];

            scalarField theta
            (
                convertToRad
              * (
                   twoPhaseAlpha2*(180 - a2cap.theta(U[patchi], nHatp))
                 + twoPhaseAlpha3*(180 - a3cap.theta(U[patchi], nHatp))
                )
            );

            vectorField nf(boundary[patchi].nf());

            // Reset nHatPatch to correspond to the contact angle

            scalarField a12(nHatp & nf);

            scalarField b1(cos(theta));

            scalarField b2(nHatp.size());

            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;

            nHatp /= (mag(nHatp) + deltaN_.value());
        }
    }
}

void Foam::threePhaseInterfaceProperties::calculateK()
{
    const volScalarField& alpha1 = mixture_.alpha1();

    const fvMesh& mesh = alpha1.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    // Cell gradient of alpha
    volVectorField gradAlpha(fvc::grad(alpha1));

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    correctContactAngle(nHatfv.boundaryFieldRef());

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
//    volVectorField nHat(gradAlpha/(mag(gradAlpha) + deltaN_));
//    nHat.boundaryField() = nHatfv.boundaryField();
//    K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::threePhaseInterfaceProperties::threePhaseInterfaceProperties
(
    const incompressibleThreePhaseMixture& mixture
)
:
    mixture_(mixture),
    cAlpha_
    (
        readScalar
        (
            mixture.U().mesh().solverDict
            (
                mixture_.alpha1().name()
            ).lookup("cAlpha")
        )
    ),

    sigma12_("sigma12", dimensionSet(1, 0, -2, 0, 0), mixture),
    sigma13_("sigma13", dimensionSet(1, 0, -2, 0, 0), mixture),
    sigma23_("sigma23", dimensionSet(1, 0, -2, 0, 0), mixture),
    sigmaMetal_("sigmaMetal", dimensionSet(1, 0, -2, 0, 0), mixture),
    sigmaMetalFlag_("sigmaMetalFlag", dimless, mixture),

    deltaN_
    (
        "deltaN",
        1e-8/pow(average(mixture.U().mesh().V()), 1.0/3.0)
    ),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            mixture.alpha1().time().timeName(),
            mixture.alpha1().mesh()
        ),
        mixture.alpha1().mesh(),
        dimensionedScalar("nHatf", dimArea, 0.0)
    ),

    K_
    (
        IOobject
        (
            "interfaceProperties:K",
            mixture.alpha1().time().timeName(),
            mixture.alpha1().mesh()
        ),
        mixture.alpha1().mesh(),
        dimensionedScalar("K", dimless/dimLength, 0.0)
    ),
    alphaMetal_
    (
        IOobject
        (
            "interfaceProperties:alphaMetal",
            mixture.alpha1().time().timeName(),
            mixture.alpha1().mesh()
        ),
        (mixture.alpha2() + mixture.alpha3())
    )
{
    calculateK();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::threePhaseInterfaceProperties::surfaceTensionForce() const
{
    const fvMesh& mesh = mixture_.alpha1().mesh();

    if (sigmaMetalFlag_.value() != 1.0)
    {
        surfaceVectorField nHatfv12Corr = nHatfv12();
        surfaceVectorField nHatfv13Corr = nHatfv13();
        surfaceVectorField nHatfv23Corr = nHatfv23();    

        correctContactAngle(nHatfv12Corr.boundaryFieldRef());
        correctContactAngle(nHatfv13Corr.boundaryFieldRef());
        correctContactAngle(nHatfv23Corr.boundaryFieldRef());

        return 
        (
            sigma12_*fvc::interpolate(-fvc::div(nHatfv12Corr & mesh.Sf()))*
            (
                fvc::interpolate(mixture_.alpha2())*fvc::snGrad(mixture_.alpha1())
                - fvc::interpolate(mixture_.alpha1())*fvc::snGrad(mixture_.alpha2())
            )

            + sigma13_*fvc::interpolate(-fvc::div(nHatfv13Corr & mesh.Sf()))*
            (
                fvc::interpolate(mixture_.alpha3())*fvc::snGrad(mixture_.alpha1())
                - fvc::interpolate(mixture_.alpha1())*fvc::snGrad(mixture_.alpha3())
            )

            + sigma23_*fvc::interpolate(-fvc::div(nHatfv23Corr & mesh.Sf()))*
            (
                fvc::interpolate(mixture_.alpha3())*fvc::snGrad(mixture_.alpha2())
                - fvc::interpolate(mixture_.alpha2())*fvc::snGrad(mixture_.alpha3())
            )
        );
    }
    else
    {
        const volScalarField& alpha1 = mixture_.alpha1();
        const volScalarField& alpha2 = mixture_.alpha2();
        const volScalarField& alpha3 = mixture_.alpha3();
        volScalarField alphaMetal = alpha2 + alpha3;


        surfaceVectorField gradAlphaMetf
        (
            fvc::interpolate(alphaMetal)*fvc::interpolate(fvc::grad(alpha1))
          - fvc::interpolate(alpha1)*fvc::interpolate(fvc::grad(alphaMetal))
        );

        // Face unit interface normal
 
        surfaceVectorField nHatfvMetalCorr = gradAlphaMetf/(mag(gradAlphaMetf) + deltaN_);
        correctContactAngle(nHatfvMetalCorr.boundaryFieldRef());

        return 
        (
            sigmaMetal_*fvc::interpolate(-fvc::div(nHatfvMetalCorr & mesh.Sf()))*
            (
                fvc::interpolate(alphaMetal)*fvc::snGrad(alpha1)
                - fvc::interpolate(alpha1)*fvc::snGrad(alphaMetal)
            )
        );        
    }
}

Foam::tmp<Foam::surfaceScalarField> 
Foam::threePhaseInterfaceProperties::nHatfMetal() const
{
    const fvMesh& mesh = mixture_.alpha1().mesh();
    const volScalarField alphaMetal = alphaMetal_;

    surfaceVectorField gradAlphaMetf
    (
        fvc::interpolate((mixture_.alpha2() + mixture_.alpha3()))*fvc::interpolate(fvc::grad(mixture_.alpha1()))
      - fvc::interpolate(mixture_.alpha1())*fvc::interpolate(fvc::grad((mixture_.alpha2() + mixture_.alpha3())))
    );

    return (gradAlphaMetf/(mag(gradAlphaMetf) + deltaN_)) & mesh.Sf();
}

Foam::tmp<Foam::surfaceVectorField> 
Foam::threePhaseInterfaceProperties::nHatfvMetal() const
{
    const volScalarField alphaMetal = alphaMetal_;

    surfaceVectorField gradAlphaMetf
    (
        fvc::interpolate((mixture_.alpha2() + mixture_.alpha3()))*fvc::interpolate(fvc::grad(mixture_.alpha1()))
      - fvc::interpolate(mixture_.alpha1())*fvc::interpolate(fvc::grad((mixture_.alpha2() + mixture_.alpha3())))
    );

    // Face unit interface normal
    return gradAlphaMetf/(mag(gradAlphaMetf) + deltaN_);
}

Foam::tmp<Foam::surfaceScalarField> 
Foam::threePhaseInterfaceProperties::nHatMetalMG() const
{
    const volScalarField alphaMetal = alphaMetal_;

    surfaceVectorField gradAlphaMetf
    (
        fvc::interpolate((mixture_.alpha2() + mixture_.alpha3()))*fvc::interpolate(fvc::grad(mixture_.alpha1()))
      - fvc::interpolate(mixture_.alpha1())*fvc::interpolate(fvc::grad((mixture_.alpha2() + mixture_.alpha3())))
    );

    return mag(gradAlphaMetf);
}

Foam::tmp<Foam::surfaceScalarField> 
Foam::threePhaseInterfaceProperties::nHatf12() const
{
    const fvMesh& mesh = mixture_.alpha1().mesh();

    surfaceVectorField gradAlpha12f
    (
        fvc::interpolate(mixture_.alpha2())*fvc::interpolate(fvc::grad(mixture_.alpha1()))
      - fvc::interpolate(mixture_.alpha1())*fvc::interpolate(fvc::grad(mixture_.alpha2()))
    );

    return (gradAlpha12f/(mag(gradAlpha12f) + deltaN_)) & mesh.Sf();
}

Foam::tmp<Foam::surfaceVectorField> 
Foam::threePhaseInterfaceProperties::nHatfv12() const
{
    surfaceVectorField gradAlpha12f
    (
        fvc::interpolate(mixture_.alpha2())*fvc::interpolate(fvc::grad(mixture_.alpha1()))
      - fvc::interpolate(mixture_.alpha1())*fvc::interpolate(fvc::grad(mixture_.alpha2()))
    );

    // Face unit interface normal
    return gradAlpha12f/(mag(gradAlpha12f) + deltaN_);
}

Foam::tmp<Foam::surfaceScalarField> 
Foam::threePhaseInterfaceProperties::nHatf13() const
{
    const fvMesh& mesh = mixture_.alpha1().mesh();

    surfaceVectorField gradAlpha13f
    (
        fvc::interpolate(mixture_.alpha3())*fvc::interpolate(fvc::grad(mixture_.alpha1()))
      - fvc::interpolate(mixture_.alpha1())*fvc::interpolate(fvc::grad(mixture_.alpha3()))
    );

    return (gradAlpha13f/(mag(gradAlpha13f) + deltaN_)) & mesh.Sf();
}

Foam::tmp<Foam::surfaceVectorField> 
Foam::threePhaseInterfaceProperties::nHatfv13() const
{
    surfaceVectorField gradAlpha13f
    (
        fvc::interpolate(mixture_.alpha3())*fvc::interpolate(fvc::grad(mixture_.alpha1()))
      - fvc::interpolate(mixture_.alpha1())*fvc::interpolate(fvc::grad(mixture_.alpha3()))
    );

    // Face unit interface normal
    return gradAlpha13f/(mag(gradAlpha13f) + deltaN_);
}

Foam::tmp<Foam::surfaceScalarField> 
Foam::threePhaseInterfaceProperties::nHatf23() const
{
    const fvMesh& mesh = mixture_.alpha1().mesh();

    surfaceVectorField gradAlpha23f
    (
        fvc::interpolate(mixture_.alpha3())*fvc::interpolate(fvc::grad(mixture_.alpha2()))
      - fvc::interpolate(mixture_.alpha2())*fvc::interpolate(fvc::grad(mixture_.alpha3()))
    );

    return (gradAlpha23f/(mag(gradAlpha23f) + deltaN_)) & mesh.Sf();
}

Foam::tmp<Foam::surfaceVectorField> 
Foam::threePhaseInterfaceProperties::nHatfv23() const
{
    surfaceVectorField gradAlpha23f
    (
        fvc::interpolate(mixture_.alpha3())*fvc::interpolate(fvc::grad(mixture_.alpha2()))
      - fvc::interpolate(mixture_.alpha2())*fvc::interpolate(fvc::grad(mixture_.alpha3()))
    );

    // Face unit interface normal
    return gradAlpha23f/(mag(gradAlpha23f) + deltaN_);
}

Foam::tmp<Foam::volScalarField>
Foam::threePhaseInterfaceProperties::nearInterface() const
{
    return max
    (
        pos0(mixture_.alpha1() - 0.01)*pos0(0.99 - mixture_.alpha1()),
        pos0(mixture_.alpha2() - 0.01)*pos0(0.99 - mixture_.alpha2())
    );
}


// ************************************************************************* //
