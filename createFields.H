    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"


    Info<< "Reading transportProperties\n" << endl;
    immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);


    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());
	// tu referencja pola tworzy w twophasemixture.C




    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();


    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1 + alpha2*rho2,
        alpha1.boundaryField().types()
    );
    rho.oldTime();


    // Mass flux
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)*phi
    );


    // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, mixture)
    );

    #include "readGravitationalAcceleration.H"

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        pimple.dict(),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }


    fv::IOoptionList fvOptions(mesh);


    // MULES flux from previous time-step
    surfaceScalarField phiAlpha
    (
        IOobject
        (
            "phiAlpha",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        phi*fvc::interpolate(alpha1)
    );

    // MULES Correction
    tmp<surfaceScalarField> tphiAlphaCorr0;



  //Smoothing
label smoothIter
(
    runTime.controlDict().lookupOrDefault<label>("smoothIter", 0)
);

label smoothIterPre
(
    runTime.controlDict().lookupOrDefault<label>("smoothIterPre", 0)
);

scalar lowRange
(
    runTime.controlDict().lookupOrDefault<scalar>("lowRange", 0.0000000001)
);

scalar upRange
(
    runTime.controlDict().lookupOrDefault<scalar>("upRange", 1)
);

// pole objetosciowe(dla komorki) do przechowywania pol (m2) danych komorek
volScalarField cellTotalSurface
(
	IOobject
	(
	    "cellTotalSurface",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	mesh,
	dimless		
);

cellTotalSurface = 0;
//pobierz chronione referencje
const labelUList & owner = mesh.owner();
const surfaceScalarField & fMag = mesh.magSf();
const labelUList & neighbour = mesh.neighbour();
// policz sume pol scianek danej komorki (przeleci po wl ale zrobi dla wszystkich)
forAll(owner, faceId)
{	
	cellTotalSurface[owner[faceId]] += fMag[faceId];
	cellTotalSurface[neighbour[faceId]] += fMag[faceId];
}
// przelec po WB, scianki WB nie maja wl i sasiadow	
forAll(mesh.boundary(), patchi)
{
	// lista komorek do ktorych naleza fejsy danego [patchi] WB
    const unallocLabelList& faceCells = mesh.boundary()[patchi].faceCells();
	// ref do aktualnego WB
    const fvPatch & patch = mesh.boundary()[patchi];

	// przelec po aktualnyb paczu WB alphy numerujac fejsy jako patchFacei
    forAll(alpha1.boundaryField()[patchi], patchFacei)
    {
	cellTotalSurface[faceCells[patchFacei]] += patch.magSf()[patchFacei];
    }
}


volScalarField alpha1_tmp
(
	IOobject
	(
	    "alpha1_tmp",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	mesh,
	dimless		
);
// ============== K ================
volScalarField K
(
	IOobject
	(
		"curvature", 
		runTime.timeName(), 
		mesh, 
		IOobject::READ_IF_PRESENT, 
		IOobject::AUTO_WRITE
	), 
	mesh, 
	dimensionedScalar("K", dimless/dimLength, 0.0)
);


// wygladz pole alpha przed rozw. rown. tr. na alpha
	
/*for(label nl=0; nl<2; ++nl)
{
	surfaceScalarField alphaSmothing("smothing", fvc::interpolate(alpha1) * mesh.magSf());

	alpha1.internalField() = 0;
	forAll(owner, faceId)
	{
		alpha1[owner[faceId]] += alphaSmothing[faceId] / cellTotalSurface[owner[faceId]];
		alpha1[neighbour[faceId]] += alphaSmothing[faceId] / cellTotalSurface[neighbour[faceId]];
	}

	forAll(mesh.boundary(), patchi)
	{
	    const unallocLabelList& faceCells = mesh.boundary()[patchi].faceCells();
	    const fvPatch & patch = mesh.boundary()[patchi];

	    forAll(alphaSmothing.boundaryField()[patchi], patchFacei)
	    {
		label cellId = faceCells[patchFacei];
		alpha1[cellId] += alphaSmothing.boundaryField()[patchi][patchFacei] / cellTotalSurface[cellId];
	    }
	}
	alpha2 = 1. - alpha1;
}
*/
	
	
