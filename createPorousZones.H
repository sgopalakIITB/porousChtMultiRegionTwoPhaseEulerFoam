/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Description
        Initialize porous zones for turbulators based on porousZones dictionary
\*---------------------------------------------------------------------------*/

forAll(fluidRegions, i)
{
    const fvMesh& mesh = fluidRegions[i];
    
    Info<< "Initializing porous zones for region " 
        << mesh.name() << endl;

    // Read porous zone dictionary if available
    IOdictionary porousDict
    (
        IOobject
        (
            "porousZones",
            runTime.constant(),
            mesh,
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
            const label zoneID = mesh.cellZones().findZoneID(zoneName);
            
            if (zoneID != -1)
            {
                const cellZone& zone = mesh.cellZones()[zoneID];
                const labelList& cells = zone;
                
                // Read zone properties
                const scalar zonePorosity = zoneDict.getOrDefault<scalar>("porosity", 1.0);
                
                // Set porosity in the zone
                forAll(cells, cellI)
                {
                    porosity[i][cells[cellI]] = zonePorosity;
                }
                
                // Read permeability tensor
                if (zoneDict.found("permeability"))
                {
                    tensor K = zoneDict.get<tensor>("permeability");
                    forAll(cells, cellI)
                    {
                        permeability[i][cells[cellI]] = K;
                    }
                }
                
                // Read Forchheimer coefficient
                if (zoneDict.found("Cf"))
                {
                    scalar zoneCf = zoneDict.get<scalar>("Cf");
                    forAll(cells, cellI)
                    {
                        Cf[i][cells[cellI]] = zoneCf;
                    }
                }
                
                // Read flow direction parameter beta
                if (zoneDict.found("betaDirection"))
                {
                    scalar beta = zoneDict.get<scalar>("betaDirection");
                    forAll(cells, cellI)
                    {
                        betaDirection[i][cells[cellI]] = beta;
                    }
                }
                
                // Read characteristic length
                if (zoneDict.found("characteristicLength"))
                {
                    scalar d_c = zoneDict.get<scalar>("characteristicLength");
                    forAll(cells, cellI)
                    {
                        characteristicLength[i][cells[cellI]] = d_c;
                    }
                }
                
                Info<< "    Porosity: " << zonePorosity << endl;
                Info<< "    Applied to " << cells.size() << " cells" << endl;
            }
            else
            {
                WarningInFunction
                    << "Cell zone " << zoneName << " not found in mesh "
                    << mesh.name() << endl;
            }
        }
    }
    else
    {
        Info<< "No porous zones dictionary found. Using default properties." << endl;
    }
    
    // Update boundary conditions
    porosity[i].correctBoundaryConditions();
    permeability[i].correctBoundaryConditions();
    Cf[i].correctBoundaryConditions();
    betaDirection[i].correctBoundaryConditions();
    characteristicLength[i].correctBoundaryConditions();
}

// ************************************************************************* //