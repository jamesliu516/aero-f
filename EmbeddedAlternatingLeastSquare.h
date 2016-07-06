//
// Created by lei on 2/3/16.
//

#ifndef _EMBEDDED_ALTERNATING_LEAST_SQUARE_H
#define _EMBEDDED_ALTERNATING_LEAST_SQUARE_H

#include <NonlinearRom.h>
#include <vector>
#include <DistInfo.h>

template<int dim>
class EmbeddedAlternatingLeastSquare : public NonlinearRom<dim> {

public: // TODO: change back to protected when done testing

    VecSet<DistSVec<char, dim> > *mask; //<! create by C++ new when needed
    int nSnapShotFiles;
    int numSnapshots; //<! total number of snapshots across all files
    std::vector<int> numMasters; //<! numMasters[i] = #{master nodes} in subdomain i
    int reducedDimension;
    /**
     * Wrapper to NonlinearRom::readSnapshotsFiles
     * @param keyword determines what snapshots is used.
     * @param preprocess defaults to false, not sure what it is.
     * @return the number of files being read.
     * The result is stored in private variable NonlinearRom::snapshots
     */
    int readSnapshotsFilesHelper(char *keyword, bool preprocess = false);

    /**
     * read maskSnapshot from file if existed
     * note that maskSnapshot is of type VecSet<DistVec<bool> >.
     * The result is stored in private variable mask.
     * TODO: read/write vector<bool> seems to produce problem
     */
    int readStateMaskFile();

    /**
     * strips distSVec of slave nodes, put it into a continuous region of memory, pointed by mem
     */
    void freeSlaves(void *&mem, const VecSet<DistSVec<double, dim> > &X, const int M, const int N);
    void freeSlaves(void *&mem, const VecSet<DistSVec<char, dim> > &X, const int M, const int N);

    /**
     * Refills distSVec with mem, a continuous region of memory, put it in X
     */
    void summonSlaves(void *&mem, VecSet<DistSVec<double, dim> > &X, const int M, const int N);
    void summonZombies(void *&mem, VecSet<DistSVec<double, dim> > &X, const int M, const int N);

    /**
     * convert a matrix stored in a continuous region of memory from row-major to column-major
     */
    void transpose(void* &buff1, void* &buff2, int nrow, int ncol);

    /**
     * returns the number of master nodes in a vector
     * res[i] = # of master nodes in subdomain i
     */
    std::vector<int> countMasters(DistInfo &distinfo);

    /**
     * test if two VecSet<DistSVec<double, dim> >& X, Y are equal
     */
    bool isEqualStateMatrices(const VecSet<DistSVec<double, dim> > &X,
                         const VecSet<DistSVec<double, dim> > &Y);
    bool isEqualMaskMatrices(const VecSet<DistSVec<char, dim> > &X,
                         const VecSet<DistSVec<char, dim> > &Y);
    void outputBasis(const VecSet<DistSVec<double, dim> > &U);
    void readBasisFiles(VecSet<DistSVec<double, dim> > &U);

public:

    EmbeddedAlternatingLeastSquare(Communicator *, IoData &, Domain &/*, DistGeoState * */);

    ~EmbeddedAlternatingLeastSquare();
    using NonlinearRom<dim>::determinePath;
    /**
     * Compute the reduced order basis given snapshot and mask.
     * Input: mask, snapshots
     */
    void ReducedOrderBasisConstruction(int _dim);

    /**
     * testing if snapshot I/O is correct; testing if armadillo compile correctly.
     */
    void testingSnapshotIO();
    void testingALS();
};


#include "EmbeddedAlternatingLeastSquare.C"

#endif // _EMBEDDED_ALTERNATING_LEAST_SQUARE_H
