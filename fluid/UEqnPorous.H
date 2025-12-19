// ============================================================================
//  Two-Phase Euler Momentum Equation with Porous Media modified by vikash
// ============================================================================

Info<< "Constructing momentum equations with porous media" << endl;

fvVectorMatrix U1Eqn
(
    U1,
    rho1.dimensions()*U1.dimensions()*dimVolume/dimTime
);

fvVectorMatrix U2Eqn
(
    U2,
    rho2.dimensions()*U2.dimensions()*dimVolume/dimTime
);

// Get interphase momentum-transfer table
autoPtr<phaseSystem::momentumTransferTable> momentumTransferPtr(fluid.momentumTransfer());
phaseSystem::momentumTransferTable& momentumTransfer = momentumTransferPtr();


// ============================================================================
//  PHASE 1 : porous Darcy–Forchheimer momentum source
// ============================================================================

{
    tmp<volVectorField> tPorousSource1
    (
        new volVectorField
        (
            IOobject
            (
                "porousSource1",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimForce/dimVolume, Zero)
        )
    );

    volVectorField& porousSource1 = tPorousSource1.ref();
    auto mu1 = thermo1.mu();

    forAll(mesh.cells(), cellI)
    {
        if (porosity[0][cellI] < 1.0)
        {
            const vector& Uc = U1[cellI];
            scalar magU = mag(Uc);

            const tensor& K = permeability[0][cellI];
            scalar Kiso = (K.xx() + K.yy() + K.zz())/3.0;

            vector darcyTerm =
                -mu1[cellI] / max(Kiso, SMALL) * Uc;

            scalar Cf_local = Cf[0][cellI];
            vector forchTerm =
                -rho1[cellI] * Cf_local * magU * Uc;

            porousSource1[cellI] = darcyTerm + forchTerm;
        }
    }

    // ---- Phase-1 Momentum Equation with Porous Sources ----
    U1Eqn =
    (
        phase1.UEqn()
     ==
        *momentumTransfer[phase1.name()]
      + fvOptions(alpha1, rho1, U1)
      + porousSource1
    );

    U1Eqn.relax();
    fvOptions.constrain(U1Eqn);
    fvOptions.correct(U1);
}



// ============================================================================
//  PHASE 2 : porous Darcy–Forchheimer momentum source
// ============================================================================

{
    tmp<volVectorField> tPorousSource2
    (
        new volVectorField
        (
            IOobject
            (
                "porousSource2",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimForce/dimVolume, Zero)
        )
    );

    volVectorField& porousSource2 = tPorousSource2.ref();
    auto mu2 = thermo2.mu();

    forAll(mesh.cells(), cellI)
    {
        if (porosity[1][cellI] < 1.0)
        {
            const vector& Uc = U2[cellI];
            scalar magU = mag(Uc);

            const tensor& K = permeability[1][cellI];
            scalar Kiso = (K.xx() + K.yy() + K.zz())/3.0;

            vector darcyTerm =
                -mu2[cellI] / max(Kiso, SMALL) * Uc;

            scalar Cf_local = Cf[1][cellI];
            vector forchTerm =
                -rho2[cellI] * Cf_local * magU * Uc;

            porousSource2[cellI] = darcyTerm + forchTerm;
        }
    }

    // ---- Phase-2 Momentum Equation with Porous Sources ----
    U2Eqn =
    (
        phase2.UEqn()
     ==
        *momentumTransfer[phase2.name()]
      + fvOptions(alpha2, rho2, U2)
      + porousSource2
    );

    U2Eqn.relax();
    fvOptions.constrain(U2Eqn);
    fvOptions.correct(U2);
}
 
 
 
 
 
 
 //siva sir code
 
 
 /*   // Solve the Momentum equation with porous media terms siva sir code

    MRF.correctBoundaryVelocity(U);

    // Calculate Darcy-Forchheimer source term: S = -μ/K * U - ρCf|U|U
    tmp<volVectorField> tporousSource
    (
        new volVectorField
        (
            IOobject
            (
                "porousSource",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimForce/dimVolume, Zero)
        )
    );
    volVectorField& porousSource = tporousSource.ref();

    // Apply Darcy-Forchheimer source term
    forAll(mesh.cells(), cellI)
    {
        if (porosity[i][cellI] < 1.0) // Only apply in porous regions
        {
            const tensor& K = permeability[i][cellI];
            const scalar& Cf_local = Cf[i][cellI];
            const vector& Ucell = U[cellI];
            const scalar magU = mag(Ucell);
            
            // Darcy term: -μ/K * U (simplified as -μ/K_iso * U for diagonal K)
            scalar K_iso = (K.xx() + K.yy() + K.zz()) / 3.0;
            tmp<volScalarField> tmu = turbulence.mu();
            const volScalarField& mu_field = tmu();
            vector darcyTerm = -mu_field[cellI] / max(K_iso, SMALL) * Ucell;
            
            // Forchheimer term: -ρCf|U|U
            vector forchheimerTerm = -rho[cellI] * Cf_local * magU * Ucell;
            
            porousSource[cellI] = darcyTerm + forchheimerTerm;
        }
    }

    UEqn =
    (
//        fvm::ddt(rho, U) + 
        fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence.divDevRhoReff(U)
     ==
        fvOptions(rho, U)
      + porousSource
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (momentumPredictor)
    {
        solve
        (
            UEqn
          ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
//            U.select(finalIter)
        );

        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }*/
