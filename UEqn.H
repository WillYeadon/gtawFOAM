    MRF.correctBoundaryVelocity(U);

    S_d = -darcyComputationalFactor*(sqr(alpha3)/(pow3(1.0 - alpha3) + 1e-3));

    volScalarField boundedAlpha1 = min(max(alpha1, scalar(0)), scalar(1));
    volScalarField boundedAlpha2 = min(max(alpha2, scalar(0)), scalar(1));

    stOn ? Info << "Marg On" << endl : Info << "Marg Off" << endl;
    surfaceScalarField snv = mixture.nHatf12();
    S_Marg = (stOn ? 1.0 : 0.0)*dsigmadT*(fvc::grad(T) - fvc::reconstruct(snv)*(fvc::reconstruct(snv) & fvc::grad(T)))*mag(fvc::grad(alpha2))*(2*mixture.rho()/ (rho1 + rho2));

//    velSource.sourceOnCheck() ? uSource = velSource.applyField()*velSourceUnits : uSource = velSourceZeroField;
    Info << "UEqn" << endl;
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U)
      - S_Marg
      - fvm::Sp(S_d, U)
//      - lorentz*lorentzUnits
//     ==
//        fvOptions(rho, U)
//        + uSource
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
//                - (vsfRho ? (rhofT ? ghf*fvc::snGrad((boundedAlpha1*rho1 + mixture.rhoBouss() - boundedAlpha2*rhoRef + mixture.rhoExp())*(scalar(1.0) - beta_2*boundedAlpha1*(T - Tref_2))) 
//                        : ghf*fvc::snGrad(rho*(scalar(1.0) - beta_2*boundedAlpha1*(T - Tref_2)))) 
//                        : ghf*fvc::snGrad(rho*(scalar(1.0) - beta*boundedAlpha2*(T - Tref))*(scalar(1.0) - beta_2*boundedAlpha1*(T - Tref_2))))    
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );
        
        Info << "UEqn End" << endl;
        fvConstraints.constrain(U);
    }
