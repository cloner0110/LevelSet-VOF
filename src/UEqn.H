    MRF.correctBoundaryVelocity(U);

    // Mass imbalance of the rho and rhoPhi used below. Taking it out keeps
    // the momentum equation consistent when rho comes from the level set
    // and is not transported by rhoPhi exactly; it is negligible for VOF.
    const volScalarField massImbalance(fvc::ddt(rho) + fvc::div(rhoPhi));

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + fvm::SuSp(-massImbalance, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.capillaryForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
