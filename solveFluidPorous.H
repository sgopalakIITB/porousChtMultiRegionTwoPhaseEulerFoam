/*---------------------------------------------------------------------------*\
    Solve fluid equations with porous media for two-phase system
\*---------------------------------------------------------------------------*/

// Note: This file is included in the main loop for solving porous fluid regions
// For the two-phase Euler solver, porous effects are included directly in
// UEqnPorous.H and EEqnPorous.H, so this file can be minimal or empty

// The actual porous momentum and energy equations are solved via:
// - UEqnPorous.H: momentum with Darcy-Forchheimer terms
// - EEqnPorous.H: energy with modified thermal conductivity

// Any additional porous-specific operations can be added here

// ************************************************************************* //
