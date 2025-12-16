/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "plateHeatExchangerCorrelations.H"
#include "mathematicalConstants.H"
#include <cmath>

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::plateHeatExchangerCorrelations::fanningFactorOffsetStrip
(
    const scalar Re,
    const scalar beta
)
{
    // Correlation for offset-strip fins
    // f = a * Re^b where a and b depend on geometry parameter beta
    
    scalar a, b;
    
    if (beta >= 0.1 && beta <= 0.5)
    {
        // Coefficients for typical offset-strip fin geometries
        a = 0.6522 * pow(beta, -0.5403);
        b = -0.1856 * log(beta) - 0.3053;
    }
    else
    {
        // Default correlation for beta outside typical range
        a = 0.316;
        b = -0.25;
    }
    
    scalar f = a * pow(Re, b);
    
    // Limit reasonable range
    return max(0.001, min(f, 10.0));
}


Foam::scalar Foam::plateHeatExchangerCorrelations::fanningFactorDimples
(
    const scalar Re,
    const scalar beta
)
{
    // Correlation for dimple-type turbulators
    // Different correlation due to different flow physics
    
    scalar f;
    
    if (Re < 100)
    {
        // Laminar regime
        f = 16.0 / Re;
    }
    else if (Re < 2300)
    {
        // Transition regime
        f = 0.316 / pow(Re, 0.25);
    }
    else
    {
        // Turbulent regime with dimple enhancement
        scalar f_smooth = 0.316 / pow(Re, 0.25);
        scalar enhancement = 1.0 + 2.5 * beta; // Beta represents dimple depth ratio
        f = f_smooth * enhancement;
    }
    
    return max(0.001, min(f, 10.0));
}


Foam::scalar Foam::plateHeatExchangerCorrelations::colburnFactorOffsetStrip
(
    const scalar Re,
    const scalar beta
)
{
    // Heat transfer correlation for offset-strip fins
    // j = St * Pr^(2/3) where St is Stanton number
    
    scalar a, b;
    
    if (beta >= 0.1 && beta <= 0.5)
    {
        // Coefficients for offset-strip fin heat transfer
        a = 0.0816 * pow(beta, -0.3931);
        b = -0.2029 * log(beta) - 0.3815;
    }
    else
    {
        // Default correlation
        a = 0.023;
        b = -0.2;
    }
    
    scalar j = a * pow(Re, b);
    
    return max(0.0001, min(j, 1.0));
}


Foam::scalar Foam::plateHeatExchangerCorrelations::colburnFactorDimples
(
    const scalar Re,
    const scalar beta
)
{
    // Heat transfer correlation for dimple-type turbulators
    
    scalar j;
    
    if (Re < 100)
    {
        // Laminar regime
        j = 0.664 / sqrt(Re);
    }
    else if (Re < 2300)
    {
        // Transition regime
        j = 0.037 / pow(Re, 0.2);
    }
    else
    {
        // Turbulent regime with dimple enhancement
        scalar j_smooth = 0.023 / pow(Re, 0.2);
        scalar enhancement = 1.0 + 3.2 * beta; // Enhanced heat transfer due to dimples
        j = j_smooth * enhancement;
    }
    
    return max(0.0001, min(j, 1.0));
}


Foam::scalar Foam::plateHeatExchangerCorrelations::characteristicDimensionOffsetStrip
(
    const scalar volume,
    const scalar surface
)
{
    // Hydraulic diameter: d_c = 4*V_f/S
    return 4.0 * volume / surface;
}


Foam::scalar Foam::plateHeatExchangerCorrelations::characteristicDimensionDimples
(
    const scalar channelHeight
)
{
    // For dimples, characteristic dimension is channel height
    return channelHeight;
}


Foam::scalar Foam::plateHeatExchangerCorrelations::calculateReynolds
(
    const scalar rho,
    const scalar U,
    const scalar d_c,
    const scalar mu
)
{
    return max(1e-10, rho * U * d_c / max(mu, 1e-10));
}


Foam::scalar Foam::plateHeatExchangerCorrelations::calculatePrandtl
(
    const scalar mu,
    const scalar cp,
    const scalar k
)
{
    return max(0.1, min(mu * cp / max(k, 1e-10), 100.0));
}


Foam::scalar Foam::plateHeatExchangerCorrelations::calculateNusselt
(
    const scalar j,
    const scalar Re,
    const scalar Pr
)
{
    return j * Re * pow(Pr, 1.0/3.0);
}


Foam::scalar Foam::plateHeatExchangerCorrelations::calculateHeatTransferCoeff
(
    const scalar Nu,
    const scalar k,
    const scalar d_c
)
{
    return Nu * k / d_c;
}


// ************************************************************************* //