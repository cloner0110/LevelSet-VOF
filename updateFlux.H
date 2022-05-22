{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    //surfaceScalarField phic(mixture.cAlpha()*mag(phi/mesh.magSf()));

    // Add the optional isotropic compression contribution
    // if (icAlpha > 0)
    // {
    //     phic *= (1.0 - icAlpha);
    //     phic += (mixture.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    // }

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phic.boundaryField()[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

//2.3
    //  phic = min(mixture.cAlpha()*phic, max(phic));
//  surfaceScalarField phir(phic*interface.nHatf());
    surfaceScalarField phir(phic*nVecf);

    surfaceScalarField phiH
    (
        fvc::flux
        (
            phi,
            H,//alpha1,
            alphaScheme
        )
      + fvc::flux
        (
            -fvc::flux(-phir, scalar(1) - H, alpharScheme),
            H,//alpha1, 
            alpharScheme
        )
    );

    //MULES::explicitSolve(alpha1, phi, phiAlpha, 1, 0);
    MULES::explicitSolve(H, phi, phiH, 1, 0);
    
    rhoPhiH = phiH*(rho1 - rho2) + phi*rho2;
}
