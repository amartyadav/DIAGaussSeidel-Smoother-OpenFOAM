#include "DIAGaussSeidelSmoother.H"
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
    Nz_(0),
    useDIA_(false),
    upperCoeffI_(0.0),
    upperCoeffJ_(0.0),
    upperCoeffK_(0.0),
    lowerCoeffI_(0.0),
    lowerCoeffJ_(0.0),
    lowerCoeffK_(0.0)
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
    // DEBUG: dump first 5 faces per direction to see actual coefficient values
    if (structuredMesh_)
    {
        const scalarField& upper = matrix_.upper();
        label countI = 0, countJ = 0, countK = 0;
        // Info<< "DEBUG: sampling upper coefficients by direction" << endl;
        for (label f = 0; f < upperAddr.size() && (countI < 5 || countJ < 5 || countK < 5); f++)
        {
            const label offset = upperAddr[f] - lowerAddr[f];
            if (offset == 1 && countI < 5)
            {
                // Info<< "  I-face " << f << ": upper=" << upper[f] << endl;
                countI++;
            }
            else if (offset == Nx_ && countJ < 5)
            {
                // Info<< "  J-face " << f << ": upper=" << upper[f] << endl;
                countJ++;
            }
            else if (offset == Nx_ * Ny_ && countK < 5)
            {
                // Info<< "  K-face " << f << ": upper=" << upper[f] << endl;
                countK++;
            }
        }
    }
    // checking coefficient uniformity
    bool allUniform = false;
    if(structuredMesh_ && !matrix_.asymmetric())
    {
        const scalarField& upper = matrix_.upper();
        bool firstI = true, firstJ = true, firstK = true;

        allUniform = true;
        label nFaces = upperAddr.size();

        for(label f = 0; f < nFaces; f++)
        {
            const label offset = upperAddr[f] - lowerAddr[f];
            const scalar coeff = upper[f];

            if(offset == 1)
            {
                if(firstI)
                {
                    upperCoeffI_ = coeff;
                    firstI = false;
                }
                else
                {
                    const scalar tol = 1e-12 * std::max(std::abs(upperCoeffI_), SMALL);
                    if(std::abs(coeff - upperCoeffI_) > tol)
                    {
                        // Info<< "  I-direction non-uniform: face " << f
                        //         << " coeff=" << coeff
                        //         << " vs canonical=" << upperCoeffI_
                        //         << " diff=" << (coeff - upperCoeffI_) << endl;
                        allUniform = false;
                        break;
                    }
                }
            }
            else if(offset == Nx_)
            {
                if(firstJ)
                {
                    upperCoeffJ_ = coeff;
                    firstJ = false;
                }
                else
                {
                    const scalar tol = 1e-12 * std::max(std::abs(upperCoeffJ_), SMALL);
                    if(std::abs(coeff - upperCoeffJ_) > tol)
                    {
                        // Info<< "  J-direction non-uniform: face " << f
                        //         << " coeff=" << coeff
                        //         << " vs canonical=" << upperCoeffJ_
                        //         << " diff=" << (coeff - upperCoeffJ_) << endl;
                        allUniform = false;
                        break;
                    }
                }
            }
            else if(offset == Nx_*Ny_)
            {
                if(firstK)
                {
                    upperCoeffK_ = coeff;
                    firstK = false;
                }
                else
                {
                    const scalar tol = 1e-12 * std::max(std::abs(upperCoeffK_), SMALL);
                    if(std::abs(coeff - upperCoeffK_) > tol)
                    {
                        // Info<< "  K-direction non-uniform: face " << f
                        //         << " coeff=" << coeff
                        //         << " vs canonical=" << upperCoeffK_
                        //         << " diff=" << (coeff - upperCoeffK_) << endl;
                        allUniform = false;
                        break;
                    }
                }
            }
        }
        if(allUniform)
        {
            lowerCoeffI_ = upperCoeffI_;

            lowerCoeffJ_ = upperCoeffJ_;
            lowerCoeffK_ = upperCoeffK_;
        }
    }

    useDIA_ = structuredMesh_ && !matrix.asymmetric() && allUniform;

    if (useDIA_)
    {
        Info<< "DIAGaussSeidel: fast path enabled, Nx=" << Nx_
            << " Ny=" << Ny_ << " Nz=" << Nz_
            << " upperCoeffs=(" << upperCoeffI_
            << ", " << upperCoeffJ_
            << ", " << upperCoeffK_ << ")"
            << " (nCells=" << nCells << ")" << endl;
    }
    else if (!structuredMesh_)
    {
        Info<< "DIAGaussSeidel: mesh not structured, falling back (nCells="
            << nCells << ")" << endl;
    }
    else if (matrix_.asymmetric())
    {
        Info<< "DIAGaussSeidel: asymmetric matrix, falling back" << endl;
    }
    else
    {
        Info<< "DIAGaussSeidel: non-uniform coefficients, falling back" << endl;
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
    if (useDIA_)
    {
        // ---- Setup ----
        // Grab raw pointers for the hot loop
        // - psiPtr: non-const because we write to psi
        // - sourcePtr: const
        // - diagPtr: const (diagonal from matrix_)
        // - nCells from psi.size() or matrix_.diag().size()

        scalar* __restrict__ psiPtr = psi.begin();
        // const scalar* __restrict__ sourcePtr = source.begin();
        const scalar* __restrict__ diagPtr = matrix_.diag().begin();
        const label nCells = psi.size();

        // Allocate bPrime as a local scalarField of size nCells
        // Grab bPrimePtr = bPrime.begin()
        scalarField bPrime(nCells);
        scalar* bPrimePtr = bPrime.begin();

        // Cache Nx, Ny, Nz into local labels (avoids repeated member access in hot loop)
        // Cache the six coefficients (upperI/J/K, lowerI/J/K) into local scalars
        label Nx = Nx_;
        label Ny = Ny_;
        label Nz = Nz_;

        scalar upperI = upperCoeffI_;
        scalar upperJ = upperCoeffJ_;
        scalar upperK = upperCoeffK_;
        scalar lowerI = lowerCoeffI_;
        scalar lowerJ = lowerCoeffJ_;
        scalar lowerK = lowerCoeffK_;

        // Compute the two strides we'll use:
        //   jStride = Nx_       (offset for idx + Nx)
        //   kStride = Nx_ * Ny_ (offset for idx + Nx*Ny)
        label jStride = Nx;
        label kStride = Nx * Ny;

        // ---- Interface boundary coefficient negation ----
        // Copy the mBouCoeffs const_cast + forAll + negate() block verbatim
        // from stock GaussSeidelSmoother.C
        // --- --- --- ---
        // Parallel boundary initialisation.  The parallel boundary is treated
            // as an effective jacobi interface in the boundary.
            // Note: there is a change of sign in the coupled
            // interface update.  The reason for this is that the
            // internal coefficients are all located at the l.h.s. of
            // the matrix whereas the "implicit" coefficients on the
            // coupled boundaries are all created as if the
            // coefficient contribution is of a source-kind (i.e. they
            // have a sign as if they are on the r.h.s. of the matrix.
            // To compensate for this, it is necessary to turn the
            // sign of the contribution.

            FieldField<Field, scalar> &mBouCoeffs =
                const_cast<FieldField<Field, scalar> &>(
                    interfaceBouCoeffs_);

            forAll(mBouCoeffs, patchi)
            {
                if (interfaces_.set(patchi))
                {
                    mBouCoeffs[patchi].negate();
                }
            }

        // ---- Sweep loop ----
        for (label sweep = 0; sweep < nSweeps; sweep++)
        {
            // Reset bPrime to source (scalarField assignment)
            bPrime = source;

            // Update interfaces (initMatrixInterfaces + updateMatrixInterfaces)
            // Copy verbatim from stock GaussSeidelSmoother.C
            matrix_.initMatrixInterfaces(
                        mBouCoeffs,
                        interfaces_,
                        psi,
                        bPrime,
                        cmpt);

                    matrix_.updateMatrixInterfaces(
                        mBouCoeffs,
                        interfaces_,
                        psi,
                        bPrime,
                        cmpt);


            // The actual DIA sweep
            for (label idx = 0; idx < nCells; idx++)
            {
                // Extract i, j, k coordinates
                //   j = idx % Nx_
                //   i = (idx / Nx_) % Ny_
                //   k = idx / kStride
                // (These are unsigned integer divisions - very cheap)
                label j = idx % Nx;
                label i = (idx / Nx) % Ny;
                label k = idx / kStride;

                scalar psii = bPrimePtr[idx];

                // Upper triangle pull (read forward neighbours)
                // Only if the forward neighbour exists (not at boundary):
                //   if (j < Nx_ - 1):  psii -= upperI * psiPtr[idx + 1]
                //   if (i < Ny_ - 1):  psii -= upperJ * psiPtr[idx + jStride]
                //   if (k < Nz_ - 1):  psii -= upperK * psiPtr[idx + kStride]
                if (j < Nx - 1) {
                    psii -= upperI * psiPtr[idx + 1];
                }
                if (i < Ny - 1) {
                    psii -= upperJ * psiPtr[idx + jStride];
                }
                if (k < Nz - 1) {
                    psii -= upperK * psiPtr[idx + kStride];
                }

                // Divide by diagonal
                psii /= diagPtr[idx];

                // Lower triangle push (update forward neighbours' bPrime)
                // Same conditionals:
                //   if (j < Nx_ - 1):  bPrimePtr[idx + 1] -= lowerI * psii
                //   if (i < Ny_ - 1):  bPrimePtr[idx + jStride] -= lowerJ * psii
                //   if (k < Nz_ - 1):  bPrimePtr[idx + kStride] -= lowerK * psii
                if (j < Nx - 1) {
                    bPrimePtr[idx + 1] -= lowerI * psii;
                }
                if (i < Ny - 1) {
                    bPrimePtr[idx + jStride] -= lowerJ * psii;
                }
                if (k < Nz - 1) {
                    bPrimePtr[idx + kStride] -= lowerK * psii;
                }

                // Commit
                psiPtr[idx] = psii;
            }
        }

        // ---- Restore interface coefficients ----
        // Same negate loop again, to undo the earlier negation
        // Copy verbatim from stock GS
        forAll(mBouCoeffs, patchi)
            {
                if (interfaces_.set(patchi))
                {
                    mBouCoeffs[patchi].negate();
                }
            }
    }
    else
    {
        // fallback - existing forwarding call
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
}
