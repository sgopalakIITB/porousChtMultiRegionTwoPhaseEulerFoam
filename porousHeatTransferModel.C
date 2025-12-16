/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Description
        Implementation of porousHeatTransferModel
\*---------------------------------------------------------------------------*/

#include "porousHeatTransferModel.H"
#include "cellZone.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousHeatTransferModel::porousHeatTransferModel
(
    const fvMesh& mesh,
    const Time& runTime,
    volScalarField& porosity,
    volTensorField& permeability,
    volScalarField& Cf,
    volScalarField& betaDirection,
    volScalarField& characteristicLength,
    volScalarField& ReField,
    volScalarField& fanningFactor,
    volScalarField& colburnFactor,
    volScalarField& heatTransferCoeff
)
:
    mesh_(mesh),
    porosity_(porosity),
    permeability_(permeability),
    Cf_(Cf),
    betaDirection_(betaDirection),
    characteristicLength_(characteristicLength),
    ReField_(ReField),
    fanningFactor_(fanningFactor),
    colburnFactor_(colburnFactor),
    heatTransferCoeff_(heatTransferCoeff)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::porousHeatTransferModel::momentumSource
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu
) const
{
    tmp<volVectorField> tSource
    (
        new volVectorField
        (
            IOobject
            (
                "porousSource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimForce/dimVolume, Zero)
        )
    );

    volVectorField& source = tSource.ref();

    // Calculate Darcy-Forchheimer source term: S = -μ/K * U - ρCf|U|U
    forAll(mesh_.cells(), cellI)
    {
        if (porosity_[cellI] < 1.0) // Only apply in porous regions
        {
            const tensor& K = permeability_[cellI];
            const scalar& Cf = Cf_[cellI];
            const vector& Ucell = U[cellI];
            const scalar magU = mag(Ucell);
            
            // Darcy term: -μ/K * U
            // Simplified approach using isotropic permeability
            scalar K_iso = (K.xx() + K.yy() + K.zz()) / 3.0;
            vector darcyTerm = -mu[cellI] / max(K_iso, SMALL) * Ucell;
            
            // Forchheimer term: -ρCf|U|U
            vector forchheimerTerm = -rho[cellI] * Cf * magU * Ucell;
            
            source[cellI] = darcyTerm + forchheimerTerm;
        }
    }

    return tSource;
}

Foam::tmp<Foam::volScalarField> Foam::porousHeatTransferModel::energySource
(
    const volScalarField& T
) const
{
    tmp<volScalarField> tSource
    (
        new volScalarField
        (
            IOobject
            (
                "energySource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimPower/dimVolume, 0.0)
        )
    );

    volScalarField& source = tSource.ref();

    // Enhanced heat transfer source term based on correlations
    forAll(mesh_.cells(), cellI)
    {
        if (porosity_[cellI] < 1.0) // Only apply in porous regions
        {
            // Enhanced convective heat transfer due to turbulators
            // For now, set to zero - can be implemented later for specific wall heat exchange
            source[cellI] = 0.0;
        }
    }

    return tSource;
}

void Foam::porousHeatTransferModel::updateCorrelations
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const volScalarField& kappa,
    const volScalarField& Cp
)
{
    // Calculate Reynolds number
    tmp<volScalarField> tmagU = mag(U);
    const volScalarField& magU = tmagU();
    ReField_ = rho * magU * characteristicLength_ / mu;

    // Limit Reynolds number to reasonable range
    ReField_ = max(ReField_, dimensionedScalar("minRe", dimless, 1.0));
    ReField_ = min(ReField_, dimensionedScalar("maxRe", dimless, 1e6));

    // Calculate Prandtl number
    tmp<volScalarField> tPr = mu * Cp / kappa;
    const volScalarField& Pr = tPr();

    // Update correlations based on turbulator type
    forAll(mesh_.cells(), cellI)
    {
        scalar Re_local = ReField_[cellI];
        scalar beta_local = betaDirection_[cellI];
        scalar Pr_local = Pr[cellI];
        
        // Determine turbulator type based on beta value
        if (beta_local < 0.1)
        {
            // Smooth channel - use standard correlations
            fanningFactor_[cellI] = 16.0 / max(Re_local, 1.0);
            colburnFactor_[cellI] = 0.664 / Foam::sqrt(max(Re_local, 1.0));
        }
        else if (beta_local < 0.8)
        {
            // Offset-strip fins
            fanningFactor_[cellI] = Foam::fanningFactorOffsetStrip(Re_local, beta_local);
            colburnFactor_[cellI] = Foam::colburnFactorOffsetStrip(Re_local, beta_local);
        }
        else
        {
            // Dimples
            fanningFactor_[cellI] = Foam::fanningFactorDimples(Re_local, beta_local);
            colburnFactor_[cellI] = Foam::colburnFactorDimples(Re_local, beta_local);
        }
        
        // Calculate Nusselt number and heat transfer coefficient
        scalar Nu_local = Foam::calculateNusselt
        (
            colburnFactor_[cellI],
            Re_local,
            Pr_local
        );
        
        heatTransferCoeff_[cellI] = Foam::calculateHeatTransferCoeff
        (
            Nu_local,
            kappa[cellI],
            characteristicLength_[cellI]
        );
    }

    // Update boundary conditions
    ReField_.correctBoundaryConditions();
    fanningFactor_.correctBoundaryConditions();
    colburnFactor_.correctBoundaryConditions();
    heatTransferCoeff_.correctBoundaryConditions();
}

void Foam::porousHeatTransferModel::initializePorousZones(const Time& runTime)
{
    // Read porous zone dictionary if available
    IOdictionary porousDict
    (
        IOobject
        (
            "porousZones",
            runTime.constant(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    if (porousDict.headerOk())
    {
        Info<< "Reading porous zone properties from porousZones dictionary" << endl;
        
        forAllConstIter(dictionary, porousDict, iter)
        {
            const word& zoneName = iter().keyword();
            const dictionary& zoneDict = iter().dict();
            
            Info<< "  Processing porous zone: " << zoneName << endl;
            
            // Get cell zone ID
            const label zoneID = mesh_.cellZones().findZoneID(zoneName);
            
            if (zoneID != -1)
            {
                const cellZone& zone = mesh_.cellZones()[zoneID];
                const labelList& cells = zone;
                
                // Read zone properties
                const scalar zonePorosity = zoneDict.getOrDefault<scalar>("porosity", 1.0);
                
                // Set porosity in the zone
                forAll(cells, cellI)
                {
                    porosity_[cells[cellI]] = zonePorosity;
                }
                
                // Read permeability tensor
                if (zoneDict.found("permeability"))
                {
                    tensor K = zoneDict.get<tensor>("permeability");
                    forAll(cells, cellI)
                    {
                        permeability_[cells[cellI]] = K;
                    }
                }
                
                // Read Forchheimer coefficient
                if (zoneDict.found("Cf"))
                {
                    scalar zoneCf = zoneDict.get<scalar>("Cf");
                    forAll(cells, cellI)
                    {
                        Cf_[cells[cellI]] = zoneCf;
                    }
                }
                
                // Read flow direction parameter beta
                if (zoneDict.found("betaDirection"))
                {
                    scalar beta = zoneDict.get<scalar>("betaDirection");
                    forAll(cells, cellI)
                    {
                        betaDirection_[cells[cellI]] = beta;
                    }
                }
                
                // Read characteristic length
                if (zoneDict.found("characteristicLength"))
                {
                    scalar d_c = zoneDict.get<scalar>("characteristicLength");
                    forAll(cells, cellI)
                    {
                        characteristicLength_[cells[cellI]] = d_c;
                    }
                }
                
                Info<< "    Porosity: " << zonePorosity << endl;
                Info<< "    Applied to " << cells.size() << " cells" << endl;
            }
            else
            {
                WarningInFunction
                    << "Cell zone " << zoneName << " not found in mesh "
                    << mesh_.name() << endl;
            }
        }
    }
    else
    {
        Info<< "No porous zones dictionary found. Using default properties." << endl;
    }
    
    // Update boundary conditions
    porosity_.correctBoundaryConditions();
    permeability_.correctBoundaryConditions();
    Cf_.correctBoundaryConditions();
    betaDirection_.correctBoundaryConditions();
    characteristicLength_.correctBoundaryConditions();
}

Foam::tmp<Foam::volScalarField> Foam::porousHeatTransferModel::alphaEff
(
    const volScalarField& alpha
) const
{
    tmp<volScalarField> tAlphaEff
    (
        new volScalarField
        (
            IOobject
            (
                "alphaEff",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            alpha * porosity_
        )
    );

    return tAlphaEff;
}

// ************************************************************************* //