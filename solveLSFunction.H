// solve Level-Set function as the re-initialization equation
Info<< "solve the reinitialization equation"     
    << nl << endl;
psi == psi0;

for (int corr=0; corr<int(epsilon.value()/deltaTau.value()); corr++)
    {
        psi = psi + psi0/mag(psi0)*(double(1)-mag(fvc::grad(psi)*dimChange))*deltaTau;
        psi.correctBoundaryConditions();
    }

// update Dirac function
forAll(mesh.cells(),celli)
{
    if(mag(psi[celli]) > epsilon.value())
        delta[celli] = double(0);
    else
        delta[celli] = double(1.0)/(double(2.0)*epsilon.value())*(double(1.0)+Foam::cos(M_PI*psi[celli]/epsilon.value()));
};
// update Heaviside function
forAll(mesh.cells(),celli)
{
    if(psi[celli] < -epsilon.value())
        H[celli] = double(0);
    else if(epsilon.value() < psi[celli])
        H [celli] = double(1);
    else
        H[celli] = double(1.0)/double(2.0)*(double(1.0)+psi[celli]/epsilon.value()+Foam::sin(M_PI*psi[celli]/epsilon.value())/M_PI);
};

const volScalarField limitedH
(
    "limitedH",
    min(max(H, scalar(0)), scalar(1))
);

rho == limitedH*rho1 + (1.0 - limitedH)*rho2;

volScalarField& nuTemp = const_cast<volScalarField&>(mixture.nu()());
nuTemp == limitedH*nu1 + (1.0 - limitedH)*nu2;

H == limitedH;
