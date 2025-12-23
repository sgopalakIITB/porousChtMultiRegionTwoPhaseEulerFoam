/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Description
        Porous heat transfer model for plate heat exchanger turbulators
        Implements Darcy-Forchheimer momentum source terms and enhanced
        thermal diffusivity for porous media
\*---------------------------------------------------------------------------*/

#ifndef porousHeatTransferModel_H
#define porousHeatTransferModel_H

#include "fvCFD.H"
#include "rhoThermo.H"
#include "correlations/plateHeatExchangerCorrelations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class porousHeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class porousHeatTransferModel
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Porosity field
        volScalarField& porosity_;

        //- Permeability tensor field
        volTensorField& permeability_;

        //- Forchheimer coefficient field
        volScalarField& Cf_;

        //- Beta direction parameter
        volScalarField& betaDirection_;

        //- Characteristic length field
        volScalarField& characteristicLength_;

        //- Reynolds number field
        volScalarField& ReField_;

        //- Fanning friction factor field
        volScalarField& fanningFactor_;

        //- Colburn factor field
        volScalarField& colburnFactor_;

        //- Heat transfer coefficient field
        volScalarField& heatTransferCoeff_;

public:

    // Constructors

        //- Construct from mesh and fields
        porousHeatTransferModel
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
        );

    //- Destructor
    virtual ~porousHeatTransferModel() = default;

    // Member Functions

        //- Calculate momentum source term (Darcy-Forchheimer)
        tmp<volVectorField> momentumSource
        (
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu
        ) const;

        //- Calculate energy source term
        tmp<volScalarField> energySource
        (
            const volScalarField& T
        ) const;

        //- Update correlations
        void updateCorrelations
        (
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const volScalarField& kappa,
            const volScalarField& Cp
        );

        //- Initialize porous zones from dictionary
        void initializePorousZones(const Time& runTime);

        //- Return porosity field
        const volScalarField& porosity() const
        {
            return porosity_;
        }

        //- Return Reynolds number field
        const volScalarField& Re() const
        {
            return ReField_;
        }

        //- Return effective thermal diffusivity
        tmp<volScalarField> alphaEff
        (
            const volScalarField& alpha
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //