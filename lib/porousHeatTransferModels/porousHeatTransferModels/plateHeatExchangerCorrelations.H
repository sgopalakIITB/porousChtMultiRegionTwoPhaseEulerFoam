/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Class
        plateHeatExchangerCorrelations

    Description
        Correlation functions for plate heat exchanger turbulators
        - Fanning friction factor f(Re, β)
        - Colburn factor j(Re, β)
        
        Based on multi-scale CFD modeling approach for offset-strip fins
        and dimple-type turbulators

\*---------------------------------------------------------------------------*/

#ifndef plateHeatExchangerCorrelations_H
#define plateHeatExchangerCorrelations_H

#include "scalar.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class plateHeatExchangerCorrelations Declaration
\*---------------------------------------------------------------------------*/

class plateHeatExchangerCorrelations
{

public:

    // Constructors
        //- Default constructor
        plateHeatExchangerCorrelations() = default;

    // Destructor
        virtual ~plateHeatExchangerCorrelations() = default;


    // Member Functions

        //- Calculate Fanning friction factor for offset-strip fins
        static scalar fanningFactorOffsetStrip(const scalar Re, const scalar beta);

        //- Calculate Fanning friction factor for dimples
        static scalar fanningFactorDimples(const scalar Re, const scalar beta);

        //- Calculate Colburn factor for offset-strip fins
        static scalar colburnFactorOffsetStrip(const scalar Re, const scalar beta);

        //- Calculate Colburn factor for dimples
        static scalar colburnFactorDimples(const scalar Re, const scalar beta);

        //- Calculate characteristic dimension for offset-strip fins
        static scalar characteristicDimensionOffsetStrip
        (
            const scalar volume,
            const scalar surface
        );

        //- Calculate characteristic dimension for dimples
        static scalar characteristicDimensionDimples(const scalar channelHeight);

        //- Calculate Reynolds number
        static scalar calculateReynolds
        (
            const scalar rho,
            const scalar U,
            const scalar d_c,
            const scalar mu
        );

        //- Calculate Prandtl number
        static scalar calculatePrandtl
        (
            const scalar mu,
            const scalar cp,
            const scalar k
        );

        //- Calculate Nusselt number from Colburn factor
        static scalar calculateNusselt
        (
            const scalar j,
            const scalar Re,
            const scalar Pr
        );

        //- Calculate heat transfer coefficient
        static scalar calculateHeatTransferCoeff
        (
            const scalar Nu,
            const scalar k,
            const scalar d_c
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //