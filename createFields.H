// Update second
Info << "Reading Transport Properties\n" << endl;
IOdictionary transportPropertiesDict
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

bool alphaGeo(transportPropertiesDict.lookupOrDefault<bool>("alphaGeo", true));
bool phaseChange(transportPropertiesDict.lookupOrDefault<bool>("phaseChange", true));
bool LfSmoothing(transportPropertiesDict.lookupOrDefault<bool>("LfSmoothing", false));
bool LfLimit(transportPropertiesDict.lookupOrDefault<bool>("LfLimit", false));

(void)alphaGeo;
(void)phaseChange;
(void)LfSmoothing;
(void)LfLimit;

#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

immiscibleIncompressibleThreePhaseMixture mixture(U, phi);

volScalarField& T(mixture.T());
volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());
volScalarField& alpha3(mixture.alpha3());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();
const dimensionedScalar& rho3 = mixture.rho3();

const dimensionedScalar& cp1 = mixture.cp1();
const dimensionedScalar& cp2 = mixture.cp2();
const dimensionedScalar& cp3 = mixture.cp3();

const dimensionedScalar& Lf = mixture.Lf();

(void)T;
(void)cp1;
(void)cp2;
(void)cp3;

surfaceScalarField kf
(
    IOobject
    (
        "kf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(mixture.k())
);

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2 + alpha3*rho3
);
rho.oldTime();


// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha solution before it is used in the U equation.
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho1*phi
);


// Construct incompressible turbulence model
autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

pressureReference pressureReference(p, p_rgh, pimple.dict());

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pressureReference.refValue()
      - getRefCellValue(p, pressureReference.refCell())
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha2.name());

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"

dimensionedScalar phaseChangeComputationalFactor
(
    "phaseChangeComputationalFactor",
    dimensionSet(0,0,0,0,0,0,0),
    transportPropertiesDict.lookupOrDefault("phaseChangeComputationalFactor", 10.0)
);

volScalarField Su
(
    IOobject
    (
        "Su",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Su",dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

volScalarField dmdt
(
    IOobject
    (
        "dmdt",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dmdt",dimensionSet(0,0,0,0,0,0,0), 0.0)
);

volScalarField Sp
(
    IOobject
    (
        "Sp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Sp",dimVelocity/dimLength, 0)
);

volScalarField dmdtSp
(
    IOobject
    (
        "dmdtSp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dmdtSp",dimensionSet(0,0,0,0,0,0,0), 0)
);

dimensionedScalar meltingT
(
    "meltingT",
    dimensionSet(0,0,0,1,0,0,0),
    transportPropertiesDict.lookupOrDefault("meltingT", 300.0)
);

volScalarField alphaMetalTotal
(
    IOobject
    (
        "alphaMetalTotal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (alpha2 + alpha3)
);

Info<< "Creating rho*cp\n" <<endl;
volScalarField rhoCp
(
    IOobject
    (
    "rhoCp",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),
    mixture.rhoCp()
);
rhoCp.oldTime();

surfaceScalarField rhoPhiCp
(
    IOobject
    (
    "rhoPhiCp",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),
    fvc::interpolate(rhoCp)*phi
);

volScalarField rhoLf
(
    IOobject
    (
    "rhoLf",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),
    mixture.rho()*Lf//,
);
rhoLf.oldTime();

volScalarField S_latent
(
    IOobject
    (
        "S_latent",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    rhoLf*fvc::ddt(alpha2)
);

dimensionedScalar darcyComputationalFactor
(
    "darcyComputationalFactor",
    dimensionSet(1,-3,-1,0,0,0,0),
    transportPropertiesDict.lookupOrDefault("darcyComputationalFactor", 8e8)
);

volScalarField S_d
(
    IOobject
    (
        "S_d",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    -darcyComputationalFactor*(sqr(alpha3)/(pow3(1.0 - alpha3) + 1e-3))
);

bool stOn = false;
//if (mixture.sigma12Access().value() != 0 && dsigmadT.value() != 0)
//    stOn = true;
(void)stOn;

dimensionedScalar dsigmadT
(
    "dsigmadT",
    dimensionSet(1,0,-2,-1,0,0,0),
    transportPropertiesDict.lookupOrDefault("dsigmadT", 0.0)
);

volVectorField S_Marg
(
    IOobject
    (
        "S_Marg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),   
    mesh,
    dimensionedVector("S_Marg",dimensionSet(1,-2,-2,0,0,0,0), Zero)
);
