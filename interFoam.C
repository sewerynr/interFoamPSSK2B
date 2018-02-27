/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    interFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

    For a two-fluid approach see twoPhaseEulerFoam.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "fixedFluxPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"
    #include "createPrghCorrTypes.H"
    #include "correctPhi.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    Info<< "\nStarting time loop\n" << endl;

	Info<< "Time = " << runTime.timeName() << nl << endl;

	for(label i=0; i<smoothIterPre; ++i)
	{
		Info<< "Initial alpha field smoothing "<< i << endl;
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
		mixture.correct();
	}

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "alphaControls.H"

            //waring: alpha1 is replaced with smoothed version of alpha
            #include "alphaEqnSubCycle.H"
		
		forAll(alpha1, cellId)
		{
			if (alpha1[cellId] < lowRange ) 
			{
			 	alpha1[cellId] = 0.0;
			}
			else if (alpha1[cellId] > upRange )
			{
				alpha1[cellId] = 1.0;
			}
		}
		alpha2 = 1. - alpha1;
            mixture.correct();

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }

	    //come back from smoothed to non-smoothed alphas	
	    alpha1.internalField() = alpha1_tmp.internalField();
            alpha2 = 1. - alpha1;
	    rho == alpha1*rho1 + alpha2*rho2;
	    mixture.correct();
	
 //- Correct the transport and interface properties immiscibleIncompressibleTwoPhaseMixture() de fakto liczy  calculateK(); z class interfaceProperties

        }
 Info<< "Time = " << runTime.timeName() << nl << endl;
	K = mixture.K();
	//mixture.K().wirte();
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
