{ 
    const scalarField& V = mesh.V() ;
    totalMass = gSum(rho*V)/totalMass0;
    // Create the residual file if not already created
    if (massFilePtr.empty())
    {
        Info<< "Creating mass file." << endl;
        // File update
        if (Pstream::master())
        {
            fileName totalMass;
            word name_ = "totalMass";
            // Open new file at start up
            massFilePtr.reset(new OFstream("totalMass.dat"));
        }
    }
    if (Pstream::master())
    {
        massFilePtr() << runTime.timeName() << tab << totalMass << endl;
    }
}
