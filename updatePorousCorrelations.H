/*---------------------------------------------------------------------------*\
    Update porous media correlations for two-phase system
\*---------------------------------------------------------------------------*/

{
    // Set up references to current region fields
    fvMesh& mesh = fluidRegions[i];
    twoPhaseSystem& fluid = phaseSystemFluid[i];
    phaseModel& phase1 = fluid.phase1();
    phaseModel& phase2 = fluid.phase2();

    volVectorField& U1 = phase1.URef();
    volVectorField& U2 = phase2.URef();
    rhoThermo& thermo1 = phase1.thermoRef();
    rhoThermo& thermo2 = phase2.thermoRef();
    volScalarField& rho1 = thermo1.rho();
    volScalarField& rho2 = thermo2.rho();

    // Use phase 1 properties for porous correlations (dominant phase)
    // This is a simplification - could be improved with volume fraction weighting

    // Update porous model correlations
    porousModels[i].updateCorrelations
    (
        rho1,
        U1,
        thermo1.mu(),
        thermo1.kappa(),
        thermo1.Cp()
    );

    // Update permeability based on correlations
    forAll(mesh.cells(), cellI)
    {
        if (porosity[i][cellI] < 1.0) // Only apply in porous regions
        {
            scalar f_local = fanningFactor[i][cellI];
            scalar Re_local = ReField[i][cellI];
            scalar d_c_local = characteristicLength[i][cellI];

            // Calculate directional permeability based on friction factor
            scalar K_parallel = d_c_local * d_c_local / max(4.0 * f_local * Re_local, 1e-12);

            // Anisotropic permeability
            tensor K_tensor = tensor::I * K_parallel;

            // Reduce permeability in perpendicular directions for turbulators
            if (betaDirection[i][cellI] > 0.1)
            {
                K_tensor.yy() *= 0.1;
                K_tensor.zz() *= 0.1;
            }

            permeability[i][cellI] = K_tensor;

            // Update Forchheimer coefficient based on turbulator type
            if (betaDirection[i][cellI] >= 0.8) // Dimples
            {
                Cf[i][cellI] = 1.75 / Foam::sqrt(150.0 * porosity[i][cellI]);
            }
            else if (betaDirection[i][cellI] >= 0.1) // Offset-strip fins
            {
                Cf[i][cellI] = 0.55 * Foam::pow(porosity[i][cellI], -1.1);
            }
        }
    }

    permeability[i].correctBoundaryConditions();
    Cf[i].correctBoundaryConditions();

    Info<< "    Updated porous correlations for region " << mesh.name() << endl;
}

// ************************************************************************* //
