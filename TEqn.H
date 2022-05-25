{ 
	Info << "\n TEqn \n" << endl;
	volScalarField boundedAlphaTotal = min(max(alpha1 + alpha2 + alpha3, scalar(0)), scalar(1));

	Info << "New thermal conductivity term" << endl;
	kf = mixture.kf();

	Info << "New Latent heat term" << endl;
//	if (!LfLimit)
//		rhoLf = LfSmoothing ? mixture.rho()*Lf : mixture.rhoBoussRho()*Lf;
	
//	S_latent = LfSmoothing ? rhoLf*fvc::ddt(alpha2) + alpha2*Lf*fvc::div(rhoPhi) : rhoLf*fvc::ddt(alpha2);// + alpha2*Lf*fvc::div(rhoPhi);
//	S_latent *= LfSmoothing ? 1.0 * boundedAlphaTotal : (1.0 - alpha1.oldTime())*(1.0 - alpha1);
//	S_latent *= phaseChange ? 1.0 : 0.0;

//	S_latent = mixture.rho()*mixture.Lf()*fvc::ddt(alpha2) + alpha2*mixture.Lf()*fvc::div(rhoPhi);
	S_latent = mixture.rho()*mixture.Lf()*(alpha2 - alpha2.oldTime())/oneSecond;


    Info << "New heat source term" << endl;
    qSource = arcSource.applyField(alphaMetalTotal)*qSourceUnits;
    rhoCp = mixture.rhoCp();

	fvScalarMatrix TEqn
	(
		fvm::ddt(rhoCp, T)
		+ fvm::div(rhoPhiCp, T)
// What is this for?
	    - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoPhiCp), T)
		- fvm::laplacian(kf, T)
		==
		- S_latent
//		+ qSource
	);

    TEqn.relax();
	Info << "TEqn End" << endl;
	TEqn.solve();
	Info << "\n TEqn Solved \n" << endl;

}