// ============================================================================
//  EEqns.H  —  Modified for porous media 
// ============================================================================

for (int Ecorr = 0; Ecorr < nEnergyCorrectors; Ecorr++)
{
    fluid.correctEnergyTransport();

    autoPtr<phaseSystem::heatTransferTable>
        heatTransferPtr(fluid.heatTransfer());

    phaseSystem::heatTransferTable& heatTransfer = heatTransferPtr();

    // ---------------------------------------------------------------
    //   PHASE 1 : porous media additions
    // ---------------------------------------------------------------
    if (!phase1.isothermal())
    {
        volScalarField& he1 = thermo1.he();

        // Porous-modified thermal diffusivity
        auto alphaEffPorous1 =
            thermo1.alphahe()*porosity[i];

        // Porous heat-source term (turbulators)
        volScalarField heatSource1
        (
            IOobject
            (
                "heatSource1",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimPower/dimVolume, 0.0)
        );

        // Build porous heat source (same logic as user’s chtMultiRegionSimpleFoam mod)
        forAll(mesh.cells(), cellI)
        {
            if (betaDirection[i][cellI] > 0.1)
            {
                scalar h_local = heatTransferCoeff[i][cellI];
                scalar Tf = thermo1.T()[cellI];
                scalar Ts = Tf - 10.0;                   // simplified fluid–solid ΔT
                scalar specificSurface = 1000.0;         // assumed value

                heatSource1[cellI] =
                    h_local*specificSurface*(Ts - Tf)*(1.0 - porosity[i][cellI]);
            }
        }

        // ================================
        //     Energy equation for phase 1
        // ================================
        fvScalarMatrix E1Eqn
        (
            phase1.heEqn()
         ==
            *heatTransfer[phase1.name()]         // interfacial heat transfer
          + alpha1*rho1*(U1 & g)                 // buoyancy
          + heatSource1                          // porous turbulator source
          + fvOptions(alpha1*porosity[i], rho1, he1)
        );

        // Replace thermal diffusion with porous-corrected version
        E1Eqn -= fvm::laplacian(thermo1.alphahe(), he1);
        E1Eqn += fvm::laplacian(alphaEffPorous1, he1);

        E1Eqn.relax();
        fvOptions.constrain(E1Eqn);
        E1Eqn.solve();
        fvOptions.correct(he1);
    }

    // ---------------------------------------------------------------
    //   PHASE 2 : porous media additions
    // ---------------------------------------------------------------
    if (!phase2.isothermal())
    {
        volScalarField& he2 = thermo2.he();

        // Porous-modified thermal diffusivity
        auto alphaEffPorous2 =
            thermo2.alphahe()*porosity[i];

        // Porous heat-source term (turbulators)
        volScalarField heatSource2
        (
            IOobject
            (
                "heatSource2",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimPower/dimVolume, 0.0)
        );

        forAll(mesh.cells(), cellI)
        {
            if (betaDirection[i][cellI] > 0.1)
            {
                scalar h_local = heatTransferCoeff[i][cellI];
                scalar Tf = thermo2.T()[cellI];
                scalar Ts = Tf - 10.0;
                scalar specificSurface = 1000.0;

                heatSource2[cellI] =
                    h_local*specificSurface*(Ts - Tf)*(1.0 - porosity[i][cellI]);
            }
        }

        // ================================
        //     Energy equation for phase 2
        // ================================
        fvScalarMatrix E2Eqn
        (
            phase2.heEqn()
         ==
            *heatTransfer[phase2.name()]         // interfacial heat transfer
          + alpha2*rho2*(U2 & g)                 // buoyancy
          + heatSource2                          // porous turbulator source
          + fvOptions(alpha2*porosity[i], rho2, he2)
        );

        // Replace thermal diffusion with porous-corrected diffusion
        E2Eqn -= fvm::laplacian(thermo2.alphahe(), he2);
        E2Eqn += fvm::laplacian(alphaEffPorous2, he2);

        E2Eqn.relax();
        fvOptions.constrain(E2Eqn);
        E2Eqn.solve();
        fvOptions.correct(he2);
    }

    // After both phases
    fluid.correctThermo();
    fluid.correct();
}




// orginal code


/*for (int Ecorr=0; Ecorr<nEnergyCorrectors; Ecorr++)
{
    fluid.correctEnergyTransport();

    autoPtr<phaseSystem::heatTransferTable>
        heatTransferPtr(fluid.heatTransfer());

    phaseSystem::heatTransferTable&
        heatTransfer = heatTransferPtr();

    if (!phase1.isothermal())
    {
        fvScalarMatrix E1Eqn
        (
            phase1.heEqn()
         ==
           *heatTransfer[phase1.name()]
          + alpha1*rho1*(U1&g)
          + fvOptions(alpha1, rho1, thermo1.he())
        );

        E1Eqn.relax();
        fvOptions.constrain(E1Eqn);
        E1Eqn.solve();
        fvOptions.correct(thermo1.he());
    }

    if (!phase2.isothermal())
    {
        fvScalarMatrix E2Eqn
        (
            phase2.heEqn()
         ==
           *heatTransfer[phase2.name()]
          + alpha2*rho2*(U2&g)
          + fvOptions(alpha2, rho2, phase2.thermoRef().he())
        );

        E2Eqn.relax();
        fvOptions.constrain(E2Eqn);
        E2Eqn.solve();
        fvOptions.correct(thermo2.he());
    }

    fluid.correctThermo();
    fluid.correct();
    
}*/
