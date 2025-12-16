/*---------------------------------------------------------------------------*\
    Enhanced createFields.H with porous media support
\*---------------------------------------------------------------------------*/

// Note: fluid/createFluidFields.H and solid/createSolidFields.H
// are already included via createFields.H

// Add porous media models for each fluid region
PtrList<porousHeatTransferModel> porousModels(fluidRegions.size());

forAll(fluidRegions, i)
{
    Info<< "*** Adding porous heat transfer model for region "
        << fluidRegions[i].name() << nl << endl;
        
    porousModels.set
    (
        i,
        new porousHeatTransferModel
        (
            fluidRegions[i],
            runTime,
            porosity[i],
            permeability[i],
            Cf[i],
            betaDirection[i],
            characteristicLength[i],
            ReField[i],
            fanningFactor[i],
            colburnFactor[i],
            heatTransferCoeff[i]
        )
    );
    
    // Initialize porous zones
    porousModels[i].initializePorousZones(runTime);
}

// ************************************************************************* //