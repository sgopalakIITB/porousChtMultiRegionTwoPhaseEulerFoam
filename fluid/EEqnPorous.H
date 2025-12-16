// ============================================================================
//  EEqnPorous.H   —   Two-Phase Euler Energy Equation with Porous Media
// ============================================================================

for (int Ecorr = 0; Ecorr < nEnergyCorrectors; Ecorr++)
{
    fluid.correctEnergyTransport();

    autoPtr<phaseSystem::heatTransferTable> heatTransferPtr(fluid.heatTransfer());
    phaseSystem::heatTransferTable& heatTransfer = heatTransferPtr();


    // ============================================================================
    //   PHASE 1 : POROUS ENERGY TRANSPORT
    // ============================================================================
    if (!phase1.isothermal())
    {
        volScalarField& he1 = thermo1.he();

        // porous-enhanced thermal diffusivity
        volScalarField alphaEffPorous1 =
            turbulence1.alphaEff()*porosity[0];

        // turbulator heat source
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

        forAll(mesh.cells(), cellI)
        {
            if (porosity[0][cellI] < 1.0)
            {
                // user’s turbulator model from EqunPorous.H
                heatSource1[cellI] = 0.0;
            }
        }


        // ---- Phase-1 Energy Equation ----
        fvScalarMatrix E1Eqn
        (
            phase1.heEqn()
         ==
            *heatTransfer[phase1.name()]
          + alpha1*rho1*(U1 & g)
          + heatSource1
          + fvOptions(alpha1*porosity[0], rho1, he1)
        );

        // replace diffusion with porous-corrected version
        E1Eqn -= fvm::laplacian(turbulence1.alphaEff(), he1);
        E1Eqn += fvm::laplacian(alphaEffPorous1, he1);

        E1Eqn.relax();
        fvOptions.constrain(E1Eqn);
        E1Eqn.solve();
        fvOptions.correct(he1);
    }




    // ============================================================================
    //   PHASE 2 : POROUS ENERGY TRANSPORT
    // ============================================================================
    if (!phase2.isothermal())
    {
        volScalarField& he2 = thermo2.he();

        volScalarField alphaEffPorous2 =
            turbulence2.alphaEff()*porosity[1];

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
            if (porosity[1][cellI] < 1.0)
            {
                heatSource2[cellI] = 0.0;
            }
        }


        // ---- Phase-2 Energy Equation ----
        fvScalarMatrix E2Eqn
        (
            phase2.heEqn()
         ==
            *heatTransfer[phase2.name()]
          + alpha2*rho2*(U2 & g)
          + heatSource2
          + fvOptions(alpha2*porosity[1], rho2, he2)
        );

        E2Eqn -= fvm::laplacian(turbulence2.alphaEff(), he2);
        E2Eqn += fvm::laplacian(alphaEffPorous2, he2);

        E2Eqn.relax();
        fvOptions.constrain(E2Eqn);
        E2Eqn.solve();
        fvOptions.correct(he2);
    }


    // after both phases
    fluid.correctThermo();
    fluid.correct();
}




// modified code for single phase done by siva sir


/*{
    volScalarField& he = thermo.he();

    // Enhanced thermal diffusivity for porous media
    tmp<volScalarField> alphaEffPorous = turbulence.alphaEff() * porosity[i];
    
    // Calculate enhanced heat transfer source term due to porous structure
    tmp<volScalarField> theatSource
    (
        new volScalarField
        (
            IOobject
            (
                "heatSource",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimPower/dimVolume, 0.0)
        )
    );
    volScalarField& heatSource = theatSource.ref();

    // Apply enhanced heat transfer in porous regions
    forAll(mesh.cells(), cellI)
    {
        if (porosity[i][cellI] < 1.0) // Only apply in porous regions
        {
            // Enhanced convective heat transfer due to turbulators
            // Placeholder for wall heat exchange if needed
            heatSource[cellI] = 0.0;
        }
    }

    fvScalarMatrix EEqn
    (
//        fvm::ddt(rho*porosity[i], he) + fvm::div(phi, he)
//      + fvc::ddt(rho*porosity[i], K) + fvc::div(phi, K)
      (
            fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
        )
      - fvm::laplacian(alphaEffPorous, he)
     ==
        rho*porosity[i]*(U&g)
      + rad.Sh(thermo, he)
      + heatSource
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    if (coupled)
    {
        fvMatrixAssemblyPtr->addFvMatrix(EEqn);
    }
    else
    {
        EEqn.solve(); // EEqn.solve(he.select(finalIter));
        fvOptions.correct(he);

        thermo.correct();
        rad.correct();

        Info<< "Min/max T:" << min(thermo.T()).value() << ' '
            << max(thermo.T()).value() << endl;
    }
}*/
