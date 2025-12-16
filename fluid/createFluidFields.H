// Initialise fluid field pointer lists
PtrList<twoPhaseSystem> phaseSystemFluid(fluidRegions.size());

PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());

PtrList<volScalarField> p_rghFluid(fluidRegions.size());

PtrList<multivariateSurfaceInterpolationScheme<scalar>::fieldTable>
    fieldsFluid(fluidRegions.size());

List<scalar> initialMassFluid(fluidRegions.size());
List<bool> frozenFlowFluid(fluidRegions.size(), false);

List<label> pRefCellFluid(fluidRegions.size());
List<scalar> pRefValueFluid(fluidRegions.size());
PtrList<dimensionedScalar> pMinFluid(fluidRegions.size());

// --- ADDED: porous media fields
PtrList<volScalarField> KFluid(fluidRegions.size());

PtrList<volScalarField> porosity(fluidRegions.size());
PtrList<volTensorField> permeability(fluidRegions.size());
PtrList<volScalarField> Cf(fluidRegions.size());
PtrList<volScalarField> betaDirection(fluidRegions.size());

PtrList<volScalarField> ReField(fluidRegions.size());
PtrList<volScalarField> fanningFactor(fluidRegions.size());
PtrList<volScalarField> colburnFactor(fluidRegions.size());
PtrList<volScalarField> heatTransferCoeff(fluidRegions.size());
PtrList<volScalarField> characteristicLength(fluidRegions.size());

PtrList<dimensionedScalar> rhoMaxFluid(fluidRegions.size());
PtrList<dimensionedScalar> rhoMinFluid(fluidRegions.size());
// --- END ADDED

const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);

PtrList<pimpleControl> pimpleFluid(fluidRegions.size());

forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    pimpleFluid.set
    (
        i,
        new pimpleControl(fluidRegions[i])
    );

    phaseSystemFluid.set(i, twoPhaseSystem::New(fluidRegions[i]).ptr());

    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(word::null, dimLength, Zero)
        )
    );

    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (g & fluidRegions[i].C()) - ghRef
        )
    );

    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (g & fluidRegions[i].Cf()) - ghRef
        )
    );

    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    p_rghFluid[i] =
        phaseSystemFluid[i].phase1().thermo().p()
      - phaseSystemFluid[i].phase1().thermo().rho()*ghFluid[i];

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

    initialMassFluid[i] =
        fvc::domainIntegrate(phaseSystemFluid[i].rho()).value();

    const dictionary& pimpleDict =
        fluidRegions[i].solutionDict().subDict("PIMPLE");

    pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);

    pRefCellFluid[i] = -1;
    pRefValueFluid[i] = 0.0;

    if (p_rghFluid[i].needReference())
    {
        setRefCell
        (
            phaseSystemFluid[i].phase1().thermoRef().p(),
            p_rghFluid[i],
            pimpleDict,
            pRefCellFluid[i],
            pRefValueFluid[i]
        );
    }

    pMinFluid.set
    (
        i,
        new dimensionedScalar
        (
            "pMin",
            dimPressure,
            phaseSystemFluid[i]
        )
    );

    // ----------------------------
    // CORRECTED POROUS FIELD INIT
    // ----------------------------

    KFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "K",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar("K", sqr(dimVelocity), 1e-6)
        )
    );

    porosity.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "porosity",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar("porosity", dimless, 1.0)
        )
    );

    permeability.set
    (
        i,
        new volTensorField
        (
            IOobject
            (
                "permeability",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedTensor("permeability", dimLength*dimLength, tensor(1,0,0,0,1,0,0,0,1)*1e12)
        )
    );

    Cf.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Cf",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar("Cf", dimless/dimLength, 0.0)
        )
    );

    betaDirection.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "betaDirection",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar("betaDirection", dimless, 0.0)
        )
    );

    ReField.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Re",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar("Re", dimless, 1.0)
        )
    );

    fanningFactor.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "fanningFactor",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar("fanningFactor", dimless, 0.01)
        )
    );

    colburnFactor.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "colburnFactor",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar("colburnFactor", dimless, 0.01)
        )
    );

    heatTransferCoeff.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "h",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar("h", dimPower/dimArea/dimTemperature, 100.0)
        )
    );

    characteristicLength.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "characteristicLength",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar("characteristicLength", dimLength, 0.001)
        )
    );

    rhoMaxFluid.set
    (
        i,
        new dimensionedScalar("rhoMax", dimDensity, GREAT, pimpleDict)
    );

    rhoMinFluid.set
    (
        i,
        new dimensionedScalar("rhoMin", dimDensity, Zero, pimpleDict)
    );
}

