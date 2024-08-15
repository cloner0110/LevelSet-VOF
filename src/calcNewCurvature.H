// calculate normal vector
volVectorField gradPsi(fvc::grad(psi));
surfaceVectorField gradPsif(fvc::interpolate(gradPsi));
surfaceVectorField nVecfv(gradPsif/(mag(gradPsif)+scalar(1.0e-6)/dimChange));
surfaceScalarField nVecf(nVecfv & mesh.Sf());

// calculate new curvature based on psi (LS function)
C == -fvc::div(nVecf);
