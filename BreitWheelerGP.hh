#ifndef BreitWheelerGP_HH
#define BreitWheelerGP_HH

#include "G4VDiscreteProcess.hh"
#include "PhotonField.hh"

#include "gp.h"
#include "gp_utils.h"
#include "rprop.h"
#include <Eigen/Dense>


class BreitWheelerGP: public G4VDiscreteProcess
{
public:
    explicit BreitWheelerGP(PhotonField* field,
        const G4String& name = "BreitWheelerGP",
        G4ProcessType type = fUserDefined);

    ~BreitWheelerGP() override;

    G4bool IsApplicable(const G4ParticleDefinition& particle);

    /* Method to set the paramters used by the gaussian process */
    void setParamsGP(bool save, int trainSize, double errorMax);

    /* Method to caculate the mean free path for interacting gamma. */
    G4double GetMeanFreePath(const G4Track& track, G4double, 
            G4ForceCondition*) override;

    G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
            const G4Step& aStep) override;


public:

    /* Returns the total cross-section for the breit wheeler process */
    double crossSection(double comEnergy);

    /* Returns the differential cross sectionof the breit wheeler process 
       used for sampling output energies */
    double diffCrossSection(double comEnergy, double theta);

    /* Returns a photon energy sampled from the d_tau / d_energy. 
       uses a basic random dart approach. */
    double samplePhotonEnergy();

    /* Returns centre of mass energy for the interaction */ 
    double sampleComEnergy(double photonEnergy, double gammaEnergy);

    /* Returns an angle of the scattered pair in the COM frame with the 
       collision axis along */
    double samplePairAngle(double comEnergy);

    /* The following method takes in a theta a phi value in the roated frame
       and gives back the index of the theta and phi in the old frame. The
       return value is an 2 element array of ints */
    void RotateThetaPhi(double theta, double phi, int thetaIndex,
            int phiIndex);
    /* The following functions are basic numerical methods used in calculating 
       tec cross-section and properties of the electron / positron producsts */

    /* Basic trapezium method for integration */
    double trapezium(double* variable, double* integrand, int resolusion);

    /* Basic linear intepolation method in 1D */
    double interpolate1D(double* samplePoints, double* sampleValues,
            int sampleSize, double queryPoint);

    /* Method foc cacluating the array index. Returns a double by interpolating
       between the two closest indices. */
    double arrayIndex(double* samplePoints, double queryPoint, int sampleSize);

private:
    PhotonField* m_field;          // Photon field which gamma interacts with
    G4RotationMatrix m_rotateForward;  // Axis to rotate gamma to z axis
    G4RotationMatrix m_rotateBackward;        // Angle to rate gamma to z axis 

    // Meber data for GP which has to be defined here annoyingly
    libgp::GaussianProcess m_gp_gausProc =
        libgp::GaussianProcess(3, "CovSum ( CovSEiso, CovNoise)");
    libgp::RProp m_gp_optimiser; // Class to optimise GP
    bool m_gp_on;   // bool setting if GP is to be used
    bool m_gp_save;  // bool setting if GP will be saved
    double** m_gp_input;  // Training energy for the GP
    double* m_gp_mfp;       // Training mfp for the GP
    int m_gp_trainSize;     // The number of data points which is trained
    double m_gp_errorMax;   // Max error before full method is used.
    int trainCount;         // count giving points in training set

};
#endif