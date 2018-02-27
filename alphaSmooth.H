
alpha1_tmp.internalField() = alpha1.internalField();

for(label nl=0; nl<smoothIter; ++nl)
{
	Info<< "!!!$$$ Internal smoother !!!$$$"<< nl << endl;
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

	//alpha1.correctBoundaryConditions();

	
	alpha2 = 1. - alpha1;
}
