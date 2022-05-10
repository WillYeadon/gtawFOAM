alphaPhiCorrs.set
(
    1,
    new surfaceScalarField
    (
        "alpha1PhiCorr",
        fvc::flux
        (
            phi,
            alpha2,
            alphaScheme
        )
    )
);

surfaceScalarField& alpha2PhiCorr = alphaPhiCorrs[1];

/*
surfaceScalarField alpha2PhiCorr
(
	fvc::flux
	(
        phi,
        alpha2,
    	alphaScheme
    )
);
*/
/*
phir = phic*mixture.nHatf12();

alpha2PhiCorr += fvc::flux
(
    -fvc::flux(-phir, alpha1, alpharScheme),
    alpha2,
    alpharScheme
);
*/
//phir = phic*mixture.nHatf23(alpha2, alpha3);
phir = phic*mixture.nHatf23();

alpha2PhiCorr += fvc::flux
(
    -fvc::flux(-phir, alpha3, alpharScheme),
    alpha2,
    alpharScheme
);

MULES::limit
(
    1.0/runTime.deltaT().value(),
    geometricOneField(),
    alpha2,
    phi,
    alpha2PhiCorr,
    zeroField(),
    zeroField(),
    oneField(),
    zeroField(),
	true
);

//alphaPhiCorrs[1] = alpha2PhiCorr;

//alpha2Phi = alpha2PhiCorr;
//surfaceScalarField alpha2Phi(alpha2PhiCorr);
//alpha2Phi += upwind<scalar>(mesh, phi).flux(alpha2);