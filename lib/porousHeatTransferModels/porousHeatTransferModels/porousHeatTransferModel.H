/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Class
        porousHeatTransferModel

    Description
        Base class for porous media heat transfer models for turbulators
        in plate heat exchangers

\*---------------------------------------------------------------------------*/

#ifndef porousHeatTransferModel_H
#define porousHeatTransferModel_H

#include "fvCFD.H"
#include "plateHeatExchangerCorrelations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class porousHeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class porousHeatTransferModel
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Porosity field
        volScalarField porosity_;

        //- Permeability tensor field
        volTensorField permeability_;

        //- Forchheimer coefficient field
        volScalarField Cf_;

        //- Beta direction parameter field
        volScalarField betaDirection_;

        //- Characteristic length field
        volScalarField characteristicLength_;

        //- Reynolds number field
        volScalarField ReField_;

        //- Fanning factor field
        volScalarField fanningFactor_;

        //- Colburn factor field
        volScalarField colburnFactor_;

        //- Heat transfer coefficient field
        volScalarField heatTransferCoeff_;

        //- Correlations object
        plateHeatExchangerCorrelations correlations_;


public:

    //- Runtime type information
    TypeName("porousHeatTransferModel");

    // Constructors

        //- Construct from mesh
        porousHeatTransferModel(const fvMesh& mesh, const Time& time);

    // Destructor
        virtual ~porousHeatTransferModel() = default;


    // Member Functions

        // Access

            //- Return porosity field
            const volScalarField& porosity() const { return porosity_; }

            //- Return permeability field  
            const volTensorField& permeability() const { return permeability_; }

            //- Return Forchheimer coefficient
            const volScalarField& Cf() const { return Cf_; }

            //- Return beta direction field
            const volScalarField& betaDirection() const { return betaDirection_; }

            //- Return Reynolds number field
            const volScalarField& Re() const { return ReField_; }

            //- Return Fanning factor field
            const volScalarField& fanningFactor() const { return fanningFactor_; }

            //- Return Colburn factor field
            const volScalarField& colburnFactor() const { return colburnFactor_; }

            //- Return heat transfer coefficient field
            const volScalarField& heatTransferCoeff() const { return heatTransferCoeff_; }

        // Edit

            //- Non-const access to porosity
            volScalarField& porosity() { return porosity_; }

            //- Non-const access to permeability
            volTensorField& permeability() { return permeability_; }


        // Evolution

            //- Initialize porous zones from dictionary
            void initializePorousZones(const Time& time);

            //- Update correlations and properties
            void updateCorrelations
            (
                const volScalarField& rho,
                const volVectorField& U,
                const volScalarField& mu,
                const volScalarField& k,
                const volScalarField& cp
            );

            //- Calculate momentum source term
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //