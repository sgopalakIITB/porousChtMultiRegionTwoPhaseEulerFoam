/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "porousHeatTransferModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(porousHeatTransferModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousHeatTransferModel::porousHeatTransferModel
(
    const fvMesh& mesh,
    const Time& time
)
:
    mesh_(mesh),
    porosity_
    (
        IOobject
        (
            "porosity",
            time.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("porosity", dimless, 1.0)
    ),
    permeability_
    (
        IOobject
        (
            "permeability",
            time.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("permeability", dimLength*dimLength, tensor::I*1e12)
    ),
    Cf_
    (
        IOobject
        (
            "Cf",
            time.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Cf", dimless/dimLength, 0.0)
    ),
    betaDirection_
    (
        IOobject
        (
            "betaDirection",
            time.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("betaDirection", dimless, 0.0)
    ),
    characteristicLength_
    (
        IOobject
        (
            "characteristicLength",
            time.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("characteristicLength", dimLength, 0.001)
    ),
    ReField_
    (
        IOobject
        (
            "Re",
            time.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Re", dimless, 1.0)
    ),
    fanningFactor_
    (
        IOobject
        (
            "fanningFactor",
            time.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("fanningFactor", dimless, 0.01)
    ),
    colburnFactor_
    (
        IOobject
        (
            "colburnFactor",
            time.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("colburnFactor", dimless, 0.01)
    ),
    heatTransferCoeff_
    (
        IOobject
        (
            "h",
            time.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("h", dimPower/dimArea/dimTemperature, 100.0)
    ),
    correlations_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porousHeatTransferModel::initializePorousZones(const Time& time)
{
    Info<< "Initializing porous zones for region " << mesh_.name() << endl;

    // Read porous zone dictionary if available
    IOdictionary porousDict
    (
        IOobject
        (
            "porousZones",
            time.constant(),
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
                
                // Read other properties
                if (zoneDict.found("Cf"))
                {
                    scalar zoneCf = zoneDict.get<scalar>("Cf");
                    forAll(cells, cellI)
                    {
                        Cf_[cells[cellI]] = zoneCf;
                    }
                }
                
                if (zoneDict.found("betaDirection"))
                {
                    scalar beta = zoneDict.get<scalar>("betaDirection");
                    forAll(cells, cellI)
                    {
                        betaDirection_[cells[cellI]] = beta;
                    }
                }
                
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
    
    // Update boundary conditions
    porosity_.correctBoundaryConditions();
    permeability_.correctBoundaryConditions();
    Cf_.correctBoundaryConditions();
    betaDirection_.correctBoundaryConditions();
    characteristicLength_.correctBoundaryConditions();
}


void Foam::porousHeatTransferModel::updateCorrelations
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const volScalarField& k,
    const volScalarField& cp
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
    tmp<volScalarField> tPr = mu * cp / k;
    const volScalarField& Pr = tPr();

    // Update correlations based on turbulator type
    forAll(mesh_.cells(), cellI)
    {
        scalar Re_local = ReField_[cellI];
        scalar beta_local = betaDirection_[cellI];
        scalar Pr_local = Pr[cellI];
        
        if (beta_local < 0.1)
        {
            // Smooth channel
            fanningFactor_[cellI] = 16.0 / max(Re_local, 1.0);
            colburnFactor_[cellI] = 0.664 / sqrt(max(Re_local, 1.0));
        }
        else if (beta_local < 0.8)
        {
            // Offset-strip fins
            fanningFactor_[cellI] = correlations_.fanningFactorOffsetStrip(Re_local, beta_local);
            colburnFactor_[cellI] = correlations_.colburnFactorOffsetStrip(Re_local, beta_local);
        }
        else
        {
            // Dimples
            fanningFactor_[cellI] = correlations_.fanningFactorDimples(Re_local, beta_local);
            colburnFactor_[cellI] = correlations_.colburnFactorDimples(Re_local, beta_local);
        }
        
        // Calculate heat transfer coefficient
        scalar Nu_local = correlations_.calculateNusselt
        (
            colburnFactor_[cellI],
            Re_local,
            Pr_local
        );
        
        heatTransferCoeff_[cellI] = correlations_.calculateHeatTransferCoeff
        (
            Nu_local,
            k[cellI],
            characteristicLength_[cellI]
        );
    }

    // Update boundary conditions
    ReField_.correctBoundaryConditions();
    fanningFactor_.correctBoundaryConditions();
    colburnFactor_.correctBoundaryConditions();
    heatTransferCoeff_.correctBoundaryConditions();
}


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
            dimensionedVector("zero", dimPressure/dimLength, Zero)
        )
    );
    volVectorField& source = tSource.ref();

    // Darcy term: -mu/K * U
    tmp<volVectorField> tdarcyTerm = -mu * (inv(permeability_) & U);
    const volVectorField& darcyTerm = tdarcyTerm();
    
    // Forchheimer term: -Cf * rho * |U| * U
    tmp<volVectorField> tforchTerm = -Cf_ * rho * mag(U) * U;
    const volVectorField& forchTerm = tforchTerm();
    
    // Turbulator source term
    forAll(mesh_.cells(), cellI)
    {
        if (betaDirection_[cellI] > 0.1) // Active turbulator zone
        {
            scalar f_local = fanningFactor_[cellI];
            scalar d_c_local = characteristicLength_[cellI];
            vector U_local = U[cellI];
            scalar rho_local = rho[cellI];
            
            // Momentum source in flow direction
            vector flowDir = U_local / (mag(U_local) + 1e-10);
            scalar sourceMag = f_local * 4.0/d_c_local * 0.5 * rho_local * magSqr(U_local);
            
            source[cellI] += -sourceMag * flowDir;
        }
    }
    
    // Add Darcy and Forchheimer terms (applied only in porous zones)
    source += (darcyTerm + forchTerm) * (1.0 - porosity_);

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
                "heatSource",
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

    // Heat transfer source terms for turbulators
    forAll(mesh_.cells(), cellI)
    {
        if (betaDirection_[cellI] > 0.1) // Active turbulator zone
        {
            scalar h_local = heatTransferCoeff_[cellI];
            scalar T_f = T[cellI];
            
            // Simplified temperature difference (in practice from coupled solving)
            scalar T_s = T_f - 10.0;
            
            // Specific surface area (typical for turbulators)
            scalar specificSurface = 1000.0; // m^2/m^3
            
            // Heat source term: h * a_s * (T_s - T_f)
            source[cellI] = h_local * specificSurface * (T_s - T_f) * (1.0 - porosity_[cellI]);
        }
    }

    return tSource;
}


// ************************************************************************* //