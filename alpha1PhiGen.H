alphaPhiCorrs.set
(
    0,
    new surfaceScalarField
    (
        "alpha1PhiCorr",
        fvc::flux
        (
            phi,
            alpha1,
            alphaScheme
        )
    )
);

surfaceScalarField& alpha1PhiCorr = alphaPhiCorrs[0];

/*
surfaceScalarField alpha1PhiCorr
(
	fvc::flux
	(
        phi,
        alpha1,
    	alphaScheme
    )
);
*/    

//surfaceScalarField phir(phic*mixture.nHatf12(alpha1, alpha2));
surfaceScalarField phir(phic*mixture.nHatf12());

//phir = phic*mixture.nHatf12(alpha1, alpha2);
//phir = phic*mixture.nHatf12();

alpha1PhiCorr += fvc::flux
(
    -fvc::flux(-phir, alpha2, alpharScheme),
    alpha1,
    alpharScheme
);

//phir = phic*mixture.nHatf13(alpha1, alpha3);
phir = phic*mixture.nHatf13();

alpha1PhiCorr += fvc::flux
(
    -fvc::flux(-phir, alpha3, alpharScheme),
    alpha1,
    alpharScheme
);

MULES::limit
(
    1.0/runTime.deltaT().value(),
    geometricOneField(),
    alpha1,
    phi,
    alpha1PhiCorr,
    zeroField(),
    zeroField(),
    oneField(),
    zeroField(),
	true
);
//surfaceScalarField alpha1Phi(alpha1PhiCorr);
//alpha1Phi += upwind<scalar>(mesh, phi).flux(alpha1);