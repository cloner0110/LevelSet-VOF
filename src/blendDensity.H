// Replace the VOF density and mass flux set by the alpha equation with the
// Heaviside blend, so that ddt(rho, U) and div(rhoPhi, U) see the same
// density field
if (mixture.propertiesFromLevelSet())
{
    rho == mixture.blendedDensity();
    rhoPhi = mixture.blendedMassFlux(phi);
}
