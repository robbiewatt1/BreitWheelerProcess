#ifndef BREITWHEELER_HH
#define BREITWHEELER_HH

#include "G4VDiscreteProcess.hh"
#include "PhotonField.hh"

class BreitWheeler: public G4VDiscreteProcess
{
public:
    explicit BreitWheeler(PhotonField* field,
        const G4String& name = "BreitWheeler",
        G4ProcessType type = fUserDefined);

    ~BreitWheeler() override;

    G4bool IsApplicable(const G4ParticleDefinition& particle);

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
    PhotonField* m_field;    // Photon field which gamma interacts with
    int* m_rotatedIndex;     // index array for rotated photon distrobution
    double** m_comEnergy;    // center of mass energy squared array
    double** m_comEnergyInt; // centre of mass integrand    
    double* m_energyInt;     // Integrand to be integrated over energy
};
#endif