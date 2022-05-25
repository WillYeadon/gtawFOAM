Info << "Phase Change \n" << endl;

dmdt = ((T - mixture.meltingT())/mixture.Lf())*mixture.cp();//*(1.0 - alpha1);

forAll(dmdt, cellI)
{
    if(alphaGeo)
    {
        if (!(alpha1[cellI] > 0.25) && ((alpha2[cellI] >= 0.5) || (alpha3[cellI] >= 0.5)))
            dmdt[cellI] *= 1;
        else
			dmdt[cellI] *= 0;
    }      
    else
    {
        if ((alpha2[cellI] >= 0.25) || (alpha3[cellI] >= 0.25))
            dmdt[cellI] *= 1;
        else
            dmdt[cellI] *= 0;
    }
    Su[cellI] = dmdt[cellI]*phaseChangeComputationalFactor.value();
    dmdt[cellI] = -dmdt[cellI]*phaseChangeComputationalFactor.value();
}
    
dmdtSp = 0.0;
zeroField Sp;
zeroField divU;
