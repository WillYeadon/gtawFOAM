alphaPhiCorrs.set
(
    2,
    new surfaceScalarField
    (
        "alpha1PhiCorr",
        fvc::flux
        (
            phi,
            alpha3,
            alphaScheme
        )
    )
);

surfaceScalarField& alpha3PhiCorr = alphaPhiCorrs[2];

/*
surfaceScalarField alpha3PhiCorr
(
    fvc::flux
    (
        phi,
        alpha3,
        alphaScheme
    )
);
*/
/*
phir = phic*mixture.nHatf23();

alpha3PhiCorr += fvc::flux
(
    -fvc::flux(-phir, alpha2, alpharScheme),
    alpha3,
    alpharScheme
);

phir = phic*mixture.nHatf13();
alpha3PhiCorr += fvc::flux
(
    -fvc::flux(-phir, alpha1, alpharScheme),
    alpha3,
    alpharScheme
);
*/
MULES::limit
(
    1.0/runTime.deltaT().value(),
    geometricOneField(),
    alpha3,
    phi,
    alpha3PhiCorr,
    zeroField(),
    zeroField(),
    oneField(),
    zeroField(),
	true
);

//alphaPhiCorrs[2] = alpha3PhiCorr;
//surfaceScalarField alpha3Phi(alpha3PhiCorr);
//alpha3Phi += upwind<scalar>(mesh, phi).flux(alpha3);