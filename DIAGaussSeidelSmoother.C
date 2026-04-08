#include "DIAGaussSeidelSmoother.H"
#include "DiagTensor.H"
#include "Field.H"
#include "FieldField.H"
#include "className.H"
#include "label.H"
#include "lduInterfaceFieldPtrsList.H"
#include "lduMatrix.H"
#include <set>
#include <vector>
#include "GaussSeidelSmoother.H"
#include "scalar.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(DIAGaussSeidelSmoother, 0);

    lduMatrix::smoother::addsymMatrixConstructorToTable<DIAGaussSeidelSmoother>
        addDIAGaussSeidelSmootherSymMatrixConstructorToTable_;

    lduMatrix::smoother::addasymMatrixConstructorToTable<DIAGaussSeidelSmoother>
        addDIAGaussSeidelSmootherAsymConstructorToTable_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DIAGaussSeidelSmoother::DIAGaussSeidelSmoother(
    const word &fieldName,
    const lduMatrix &matrix,
    const FieldField<Field, scalar> &interfaceBouCoeffs,
    const FieldField<Field, scalar> &interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList &interfaces)
    : lduMatrix::smoother(
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces
    ),
    structuredMesh_(false),
    Nx_(0),
    Ny_(0),
    Nz_(0)
{
    const label nCells = matrix_.diag().size();

    const labelUList upperAddr = matrix_.lduAddr().upperAddr();
    const labelUList lowerAddr = matrix_.lduAddr().lowerAddr();

    // unique offsets
    std::set<label> offsets;
    for(label f = 0; f < lowerAddr.size(); f++)
    {
        offsets.insert(upperAddr[f] - lowerAddr[f]);
    }
    std::vector<label> sortedOffsets(offsets.begin(), offsets.end());

    // checking if the size of the set (offsets) is 2 (for 2d) or 3 (for 3d)
    if(sortedOffsets.size() == 2)
    {
        if(sortedOffsets[0] == 1 && nCells % sortedOffsets[1] == 0)
        {
            Nx_ = sortedOffsets[1];
            Ny_ = nCells / Nx_;
            Nz_ = 1;
            structuredMesh_ = true;
        }
    }
    else if(sortedOffsets.size() == 3)
    {
        if(sortedOffsets[0] == 1 && sortedOffsets[2] % sortedOffsets[1] == 0
            && nCells % sortedOffsets[2] == 0)
        {
            Nx_ = sortedOffsets[1];
            Ny_ = sortedOffsets[2] / sortedOffsets[1];
            Nz_ = nCells / sortedOffsets[2];
            structuredMesh_ = true;
        }
    }

    // Logging the outcome
    if(structuredMesh_)
    {
        Info << "DIAGaussSeidel: structured hex mesh detected, Nx=" << Nx_
             << " Ny=" << Ny_ << " Nz=" << Nz_ << " (nCells=" << nCells << ")" << endl;
    }
    else {
        Info << "DIAGaussSeidel: structured hex mesh detection failed, nCells=" << nCells << endl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// void Foam::DIAGaussSeidelSmoother::smooth(
//     const word &fieldName_,
//     scalarField &psi,
//     const lduMatrix &matrix_,
//     const scalarField &source,
//     const FieldField<Field, scalar> &interfaceBouCoeffs_,
//     const lduInterfaceFieldPtrsList &interfaces_,
//     const direction cmpt,
//     const label nSweeps)
// {
// }

// Define the virtual smooth method (the 4-arg const one). Its body should just call Foam::GaussSeidelSmoother::smooth(...) — the stock one — passing fieldName_, psi, matrix_, source, interfaceBouCoeffs_, interfaces_, cmpt, nSweeps. All the underscore-suffixed members are inherited from the protected base class.
void Foam::DIAGaussSeidelSmoother::smooth(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    if(structuredMesh_)
    {
        // Stage 2: fast path placeholder — still forwards to stock GS.
        // Stage 3 will replace this with the DIA kernel.
        static bool firstCall = true;
        if(firstCall)
        {
            Info<< "DIAGaussSeidel: fast path would run here "
                            << "(Nx=" << Nx_ << " Ny=" << Ny_ << " Nz=" << Nz_ << ")"
                            << endl;
                        firstCall = false;
        }
    }
   GaussSeidelSmoother::smooth(
       fieldName_,
       psi,
       matrix_,
       source,
       interfaceBouCoeffs_,
       interfaces_,
       cmpt,
       nSweeps
   );
}
