if (finalIter)
{
    mesh.data().setFinalIteration(true);
}

// Momentum equations must be defined before pressure equation
#include "UEqns.H"
#include "pEqn.H"
#include "EEqns.H" //energy equation with porous
//#include "YEqns.H" // this is used for species transport

fluid.correctKinematics();

// Update alpha's for new U
fluid.correctTurbulence();

if (finalIter)
{
    mesh.data().setFinalIteration(false);
}

/*if (frozenFlow)
{
    #include "EEqns.H"
}
else
{
    fluid.solve();
    fluid.correct();

    #include "YEqns.H"

    if (faceMomentum)
    {
        #include "pUf/UEqns.H"
        #include "EEqns.H"
        #include "pUf/pEqn.H"
    }
    else
    {
        #include "pU/UEqns.H"
        #include "EEqns.H"
        #include "pU/pEqn.H"
    }

    fluid.correctKinematics();

    // Update alpha's for new U
    fluid.correctTurbulence();
}

if (finalIter)
{
    mesh.data().setFinalIteration(false);
}
*/
