// Modified code by vikash for porous media libraries.

Info<< "Constructing momentum equations" << endl;

fvVectorMatrix U1Eqn(U1, rho1.dimensions()*U1.dimensions()*dimVolume/dimTime);
fvVectorMatrix U2Eqn(U2, rho2.dimensions()*U2.dimensions()*dimVolume/dimTime);

// === Porous-media momentum source fields (added) === //
volVectorField porousSource1
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
    dimensionedVector("zero", dimPressure/dimLength, Zero)
);

volVectorField porousSource2
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
    dimensionedVector("zero", dimPressure/dimLength, Zero)
);

{
    autoPtr<phaseSystem::momentumTransferTable> momentumTransferPtr(fluid.momentumTransfer());
    phaseSystem::momentumTransferTable& momentumTransfer(momentumTransferPtr());

    // ---------------------------------------------------------
    // === Phase 1 porous momentum source terms (added) ===
    // ---------------------------------------------------------

    // dynamic viscosity of phase 1
    auto mu1 = thermo1.mu();

    // Darcy:  -mu/K * U
    auto darcy1 = -mu1 * (inv(permeability[0]) & U1);

    // Forchheimer: -Cf * rho * |U| * U
    auto forch1 = -Cf[0] * rho1 * mag(U1) * U1;

    // Turbulator directional source
    volVectorField turbulatorSource1
    (
        IOobject
        (
            "turbulatorSource1",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimPressure/dimLength, Zero)
    );

    forAll(mesh.cells(), cellI)
    {
        if (betaDirection[0][cellI] > 0.1)
        {
            scalar f_local = fanningFactor[0][cellI];
            scalar d_c_local = characteristicLength[0][cellI];
            vector U_local = U1[cellI];
            scalar rho_local = rho1[cellI];

            vector flowDir = U_local/(mag(U_local) + 1e-10);
            scalar sourceMag =
                f_local * 4.0/d_c_local * 0.5 * rho_local * magSqr(U_local);

            turbulatorSource1[cellI] = -sourceMag * flowDir;
        }
    }

    // Total porous source applied only in porous zone
    porousSource1 =
        (darcy1 + forch1 + turbulatorSource1)
        * (1.0 - porosity[0]);


    // ---------------------------------------------------------
    // === Phase 2 porous momentum source terms (added) ===
    // ---------------------------------------------------------

    auto mu2 = thermo2.mu();

    auto darcy2 = -mu2 * (inv(permeability[1]) & U2);
    auto forch2 = -Cf[1] * rho2 * mag(U2) * U2;

    volVectorField turbulatorSource2
    (
        IOobject
        (
            "turbulatorSource2",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimPressure/dimLength, Zero)
    );

    forAll(mesh.cells(), cellI)
    {
        if (betaDirection[1][cellI] > 0.1)
        {
            scalar f_local = fanningFactor[1][cellI];
            scalar d_c_local = characteristicLength[1][cellI];
            vector U_local = U2[cellI];
            scalar rho_local = rho2[cellI];

            vector flowDir = U_local/(mag(U_local) + 1e-10);
            scalar sourceMag =
                f_local * 4.0/d_c_local * 0.5 * rho_local * magSqr(U_local);

            turbulatorSource2[cellI] = -sourceMag * flowDir;
        }
    }

    porousSource2 =
        (darcy2 + forch2 + turbulatorSource2)
        * (1.0 - porosity[1]);



    // ---------------------------------------------------------
    // === Momentum equation phase 1 (modified) ===
    // ---------------------------------------------------------

    {
        U1Eqn =
        (
            phase1.UEqn()
         ==
           *momentumTransfer[phase1.name()]
         + fvOptions(alpha1, rho1, U1)
         + porousSource1       // << ADDED
        );
        U1Eqn.relax();
        fvOptions.constrain(U1Eqn);
        fvOptions.correct(U1);
    }

    // ---------------------------------------------------------
    // === Momentum equation phase 2 (modified) ===
    // ---------------------------------------------------------

    {
        U2Eqn =
        (
            phase2.UEqn()
         ==
           *momentumTransfer[phase2.name()]
         + fvOptions(alpha2, rho2, U2)
         + porousSource2       // << ADDED
        );
        U2Eqn.relax();
        fvOptions.constrain(U2Eqn);
        fvOptions.correct(U2);
    }
}




// orginal code


/*Info<< "Constructing momentum equations" << endl;

fvVectorMatrix U1Eqn(U1, rho1.dimensions()*U1.dimensions()*dimVolume/dimTime);
fvVectorMatrix U2Eqn(U2, rho2.dimensions()*U2.dimensions()*dimVolume/dimTime);

{
    autoPtr<phaseSystem::momentumTransferTable>
        momentumTransferPtr(fluid.momentumTransfer());

    phaseSystem::momentumTransferTable&
        momentumTransfer(momentumTransferPtr());

    {
        U1Eqn =
        (
            phase1.UEqn()
         ==
           *momentumTransfer[phase1.name()]
          + fvOptions(alpha1, rho1, U1)
        );
        U1Eqn.relax();
        fvOptions.constrain(U1Eqn);
        fvOptions.correct(U1);
    }

    {
        U2Eqn =
        (
            phase2.UEqn()
         ==
           *momentumTransfer[phase2.name()]
          + fvOptions(alpha2, rho2, U2)
        );
        U2Eqn.relax();
        fvOptions.constrain(U2Eqn);
        fvOptions.correct(U2);
    }
}

*/
