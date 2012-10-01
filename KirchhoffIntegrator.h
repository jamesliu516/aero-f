#ifndef _KIRCHHOFF_INTEGRATOR_H
#define _KIRCHHOFF_INTEGRATOR_H


//------------------------------------------------------------------------------
#include <complex>
#include <vector>


class Domain;
class IoData;


template<class Scalar, int dim> class DistSVec;


class KirchhoffIntegrator {

  private:

    // Do not define these functions
    KirchhoffIntegrator(const KirchhoffIntegrator &ref);
    KirchhoffIntegrator & operator=(const KirchhoffIntegrator &ref);

  protected:

    IoData &d_iod;
    Domain *d_domain_p;
    DistSVec<double,3> *d_X_p;

    std::vector<int> d_globalToNodeID;

    enum TypeGamma {SPHERE = 0, CYLINDER = 1} d_SurfType;
    double d_R;

    ///////////////////////
    // Utility functions
    ///////////////////////


    // Function to evaluate the spherical bessel h_n at a point x.
    std::complex<double> besselh
      (
        const int n, const double x
      ) const;


    // Function to evaluate the derivative of a spherical bessel h_n.
    std::complex<double> besselh_prime
      (
        const int n, const double x
      ) const;


    // This function computes the normalized spherical harmonics
    // at a point (theta, phi).
    // For a given value n, all values between -n and n are computed.
    // They are stored in the output array Yn.
    // Note: Yn may be resized.
    void sphericalHarmonic
      (
        const int n,
        const double theta,
        const double phi,
        std::vector< std::complex<double> > &Yn
      ) const;


    // This function compute the L2-norm squared of the pressure 
    // on the surface.
    void getL2NormSquare
    (
      std::complex<double> *nodal,
      const int numFreq,
      const int *vecLen,
      double *l2Norm
    );
  
  
    // This function converts nodal values on a spherical surface
    // to the coefficients in a spherical harmonic series expansion.
    void convertToSHSeries
      (
        std::complex<double> *nodal,
        const int numVec,
        const int *vecLen,
        std::vector< std::complex<double> > &series,
        int &nmax
      );


    // This function computes the coefficients in a spherical harmonic series
    // expansion for the normal derivative of the pressure.
    // The pressure is also specified by the coefficients in a spherical
    // harmonic series.
    void getdpdnSHseries
     (
       std::complex<double> *coeff,
       std::complex<double> *coeff_dudn,
       const double *freqSamples,
       const int numFrequencies,
       const int nMax
      ) const;

  
    // This function computes the Kirchhoff integral on a spherical
    // surface.
    void integrateOnSphere
      (
        std::complex<double> *pvalues,
        int *vecLen,
        std::complex<double> *dpdn,
        int nmax,
        const double *freqSamples,
        const int numFreq
      );

  
    // This function computes the far-field pattern
    // when the surface is a sphere.
    void ffpDataOnSphere
    (
      std::complex<double> *pvalues,
      int *vecLen,
      std::complex<double> *dpdn,
      const int nmax,
      const double *freqSamples,
      const int numFreq
     );
  
  
    // This function evaluates a series in spherical harmonics
    // at a point (theta, phi) with the input coefficients.
    std::complex<double> evaluateSHS
      (
        std::complex<double> *coeff,
        int nmax,
        double theta,
        double phi
      ) const;
  
  
    // This function returns quadrature points on a face.
    //
    // \note UH (08/2012) It is implemented only for triangles
    // with a 7-point quadrature rule.
    void getQuadrature
      (
        std::vector<double> &_xigauss,
        std::vector<double> &w
      ) const;


    // This function converts cartesian coordinates to spherical coordinates.
    // theta is between 0 and pi, phi between 0 and 2*pi
    //
    // \note The cut is for {x > 0, y = 0}
    void cart2sph_x
      (
        double xo, double yo, double zo,
        double &ro, double &to, double &po
      ) const;


    // This function converts cartesian coordinates to spherical coordinates.
    // theta is between 0 and pi, phi between -pi/2 and 3*pi/2
    //
    // \note The cut is for {x > 0, y = 0}
    void cart2sph_y
      (
        double xo, double yo, double zo,
        double &ro, double &to, double &po
      ) const;


  public:


    //! Constructor
    ///
    ///
    KirchhoffIntegrator
      (
        IoData &iod, 
        Domain *domain
      );


    //! Destructor 
    ///
    ///
    ~KirchhoffIntegrator();


    //! Function to compute the Kirchhoff integral at different points.
    ///
    /// This function evaluates the Kirchhoff integral from snapshots
    /// of pressure nodal values.
    ///
    void computeIntegral();


};


//------------------------------------------------------------------------------


#endif

