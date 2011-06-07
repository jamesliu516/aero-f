#ifndef VECTOR_SET_OP_H_
#define VECTOR_SET_OP_H_

#include <DistVector.h>
#include <VectorSet.h>
#include <DistInfo.h>
#include <cassert>

//------------------------------------------------------------------------------

template <int dim>
void transMatVecProd (const VecSet< DistSVec<double, dim> > &matrix, const
		DistSVec<double, dim> &vector, double *targetBuffer){

	// see DistVec<double>::operator*(const DistVec<double> &vector)

	assert(&vector.info() == &matrix.size());

	int numVec = matrix.numVectors();

	const DistInfo &distInfo = vector.info();

	int iSub;

	double *totalAllRes = reinterpret_cast<double *>(alloca(sizeof(double) *
				distInfo.numGlobSub * numVec));

	for (int iVec = 0; iVec < numVec; ++iVec) {

		const DistSVec<double, dim> &matVector = matrix[iVec];
		double &res = targetBuffer[iVec];
		res = 0;

#ifndef MPI_OMP_REDUCTION
		double *allres = &totalAllRes[iVec * distInfo.numGlobSub];

		for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) allres[iSub] = 0;
#endif

		if (distInfo.masterFlag) {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
			for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

				int locOffset = distInfo.subOffset[iSub];
				int locLen = distInfo.subLen[iSub];

				double locres = 0;

				for (int i = 0; i < locLen; ++i)
					if (distInfo.masterFlag[locOffset+i])
						for (int j = 0; j < dim; ++j)
							locres += matVector.data()[locOffset+i][j] * vector.data()[locOffset+i][j];

#ifdef MPI_OMP_REDUCTION
				res += locres;
#else
				allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif

			}

		} 
		else {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
			for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

				int locOffset = distInfo.subOffset[iSub];
				int locLen = distInfo.subLen[iSub];

				double locres = 0;

				for (int i = 0; i < locLen; ++i)
					for (int j = 0; j < dim; ++j)
						locres += matVector.data()[locOffset+i][j] * vector.data()[locOffset+i][j];

#ifdef MPI_OMP_REDUCTION
				res += locres;
#else
				allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif

			}
		}
	}

#ifdef MPI_OMP_REDUCTION
	distInfo.com->globalSum(numVec, targetBuffer);
#else
	distInfo.com->globalSum(distInfo.numGlobSub * numVec, totalAllRes);

	for (int iVec = 0; iVec < numVec; ++iVec) {
		targetBuffer[iVec] = 0;
		for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) {
			targetBuffer[iVec] += totalAllRes[iVec * distInfo.numGlobSub + iSub];
		}
	}
#endif

}

template <int dim>
void transMatVecProdRestrict(const VecSet< DistSVec<double, dim> > &matrix, const
		DistSVec<double, dim> &vector, double *targetBuffer,
		const std::vector<std::vector<int> > & sampleNodes) {

	// see DistVec<double>::operator*(const DistVec<double> &vector)

	assert(&vector.info() == &matrix.size());

	int numVec = matrix.numVectors();

	const DistInfo &distInfo = vector.info();

	int iSub, i;

	double *totalAllRes = reinterpret_cast<double *>(alloca(sizeof(double) *
				distInfo.numGlobSub * numVec));

	for (int iVec = 0; iVec < numVec; ++iVec) {

		const DistSVec<double, dim> &matVector = matrix[iVec];
		double &res = targetBuffer[iVec];
		res = 0;

#ifndef MPI_OMP_REDUCTION
		double *allres = &totalAllRes[iVec * distInfo.numGlobSub];

		for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) allres[iSub] = 0;
#endif

		if (distInfo.masterFlag) {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
			for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

				int locOffset = distInfo.subOffset[iSub];
				int locLen = distInfo.subLen[iSub];

				double locres = 0;

				for (int iNode = 0; iNode < sampleNodes[iSub].size(); ++iNode) {
					i = sampleNodes[iSub][iNode];
					if (distInfo.masterFlag[locOffset+i])
						for (int j = 0; j < dim; ++j)
							locres += matVector.data()[locOffset+i][j] * vector.data()[locOffset+i][j];
				}

#ifdef MPI_OMP_REDUCTION
				res += locres;
#else
				allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif

			}

		} 
		else {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
			for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

				int locOffset = distInfo.subOffset[iSub];
				int locLen = distInfo.subLen[iSub];

				double locres = 0;

				for (int iNode = 0; iNode < sampleNodes[iSub].size(); ++iNode) {
					i = sampleNodes[iSub][iNode];
					for (int j = 0; j < dim; ++j)
						locres += matVector.data()[locOffset+i][j] * vector.data()[locOffset+i][j];
				}

#ifdef MPI_OMP_REDUCTION
				res += locres;
#else
				allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif

			}
		}
	}

#ifdef MPI_OMP_REDUCTION
	distInfo.com->globalSum(numVec, targetBuffer);
#else
	distInfo.com->globalSum(distInfo.numGlobSub * numVec, totalAllRes);

	for (int iVec = 0; iVec < numVec; ++iVec) {
		targetBuffer[iVec] = 0;
		for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) {
			targetBuffer[iVec] += totalAllRes[iVec * distInfo.numGlobSub + iSub];
		}
	}
#endif

}
template <int dim>
void transMatMatProd (const VecSet< DistSVec<double, dim> > &matrix1, 
		const VecSet< DistSVec<double, dim> > &matrix2, double *targetBuffer){

	// see DistVec<double>::operator*(const DistVec<double> &matrix2)

	assert(&matrix1.size() == &matrix2.size());

	int numVec1 = matrix1.numVectors();
	int numVec2 = matrix2.numVectors();

	const DistInfo &distInfo = matrix1.size();

	int iSub;

	double *totalAllRes = reinterpret_cast<double *>(alloca(sizeof(double) *
				distInfo.numGlobSub * numVec1 * numVec2));

	for (int iVec1 = 0; iVec1 < numVec1; ++iVec1) {
		for (int iVec2 = 0; iVec2 < numVec2; ++iVec2) {

			const DistSVec<double, dim> &matVector1 = matrix1[iVec1];
			const DistSVec<double, dim> &matVector2 = matrix2[iVec2];
			double &res = targetBuffer[iVec1 + iVec2 * numVec1];
			res = 0;

#ifndef MPI_OMP_REDUCTION
			double *allres = &totalAllRes[iVec1 * distInfo.numGlobSub + iVec2 *
				(numVec1 *distInfo.numGlobSub)];

			for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) allres[iSub] = 0;
#endif

			if (distInfo.masterFlag) {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
				for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

					int locOffset = distInfo.subOffset[iSub];
					int locLen = distInfo.subLen[iSub];

					double locres = 0;

					for (int i = 0; i < locLen; ++i)
						if (distInfo.masterFlag[locOffset+i])
							for (int j = 0; j < dim; ++j)
								locres += matVector1.data()[locOffset+i][j] * matVector2.data()[locOffset+i][j];

#ifdef MPI_OMP_REDUCTION
					res += locres;
#else
					allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif

				}

			} 
			else {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
				for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

					int locOffset = distInfo.subOffset[iSub];
					int locLen = distInfo.subLen[iSub];

					double locres = 0;

					for (int i = 0; i < locLen; ++i)
						for (int j = 0; j < dim; ++j)
							locres += matVector1.data()[locOffset+i][j] * matVector2.data()[locOffset+i][j];

#ifdef MPI_OMP_REDUCTION
					res += locres;
#else
					allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif

				}
			}
		}
	}

#ifdef MPI_OMP_REDUCTION
	distInfo.com->globalSum(numVec1 * numVec2, targetBuffer);
#else
	distInfo.com->globalSum(distInfo.numGlobSub * numVec1 * numVec2, totalAllRes);

	for (int iVec1 = 0; iVec1 < numVec1 ; ++iVec1) {
		for (int iVec2 = 0; iVec2 < numVec2; ++iVec2) {
			targetBuffer[iVec1 + iVec2 * numVec1] = 0;
			for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) {
				targetBuffer[iVec1 + iVec2 * numVec1] += totalAllRes[iVec1 *
					distInfo.numGlobSub + iVec2 * (numVec1 * distInfo.numGlobSub) + iSub];
			}
		}
	}
#endif

}

template <int dim>
void transMatMatProdRestrict(const VecSet< DistSVec<double, dim> > &matrix1,
		const VecSet< DistSVec<double, dim> > &matrix2, double *targetBuffer,
		const std::vector<std::vector<int> > &sampleNodes){

	// see DistVec<double>::operator*(const DistVec<double> &matrix2)

	assert(&matrix1.size() == &matrix2.size());

	int numVec1 = matrix1.numVectors();
	int numVec2 = matrix2.numVectors();

	const DistInfo &distInfo = matrix1.size();

	int iSub, i;

	double *totalAllRes = reinterpret_cast<double *>(alloca(sizeof(double) *
				distInfo.numGlobSub * numVec1 * numVec2));

	for (int iVec1 = 0; iVec1 < numVec1; ++iVec1) {
		for (int iVec2 = 0; iVec2 < numVec2; ++iVec2) {

			const DistSVec<double, dim> &matVector1 = matrix1[iVec1];
			const DistSVec<double, dim> &matVector2 = matrix2[iVec2];
			double &res = targetBuffer[iVec1 + iVec2 * numVec1];
			res = 0;

#ifndef MPI_OMP_REDUCTION
			double *allres = &totalAllRes[iVec1 * distInfo.numGlobSub + iVec2 *
				(numVec1 *distInfo.numGlobSub)];

			for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) allres[iSub] = 0;
#endif

			if (distInfo.masterFlag) {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
				for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

					int locOffset = distInfo.subOffset[iSub];
					int locLen = distInfo.subLen[iSub];

					double locres = 0;
					for (int iNode = 0; iNode < sampleNodes[iSub].size(); ++iNode) {
						i = sampleNodes[iSub][iNode];

						if (distInfo.masterFlag[locOffset+i])
							for (int j = 0; j < dim; ++j)
								locres += matVector1.data()[locOffset+i][j] * matVector2.data()[locOffset+i][j];
					}

#ifdef MPI_OMP_REDUCTION
					res += locres;
#else
					allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif

				}

			} 
			else {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
				for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

					int locOffset = distInfo.subOffset[iSub];
					int locLen = distInfo.subLen[iSub];

					double locres = 0;

					for (int iNode = 0; iNode < sampleNodes[iSub].size(); ++iNode) {
						i = sampleNodes[iSub][iNode];
						for (int j = 0; j < dim; ++j)
							locres += matVector1.data()[locOffset+i][j] * matVector2.data()[locOffset+i][j];
					}

#ifdef MPI_OMP_REDUCTION
					res += locres;
#else
					allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif

				}
			}
		}
	}

#ifdef MPI_OMP_REDUCTION
	distInfo.com->globalSum(numVec1 * numVec2, targetBuffer);
#else
	distInfo.com->globalSum(distInfo.numGlobSub * numVec1 * numVec2, totalAllRes);

	for (int iVec1 = 0; iVec1 < numVec1 ; ++iVec1) {
		for (int iVec2 = 0; iVec2 < numVec2; ++iVec2) {
			targetBuffer[iVec1 + iVec2 * numVec1] = 0;
			for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) {
				targetBuffer[iVec1 + iVec2 * numVec1] += totalAllRes[iVec1 *
					distInfo.numGlobSub + iVec2 * (numVec1 * distInfo.numGlobSub) + iSub];
			}
		}
	}
#endif

}
#endif
