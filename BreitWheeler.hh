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

private:

    /* Returns the total cross-section for the breit wheeler process */
    double crossSection(double comEnergy);

    /* Returns the differential cross sectionof the breit wheeler process 
       used for sampling output energies */
    double diffCrossSection(double comEnergy, double theta);

    /* Basic trapezium method for integration */
    double trapezium(double* variable, double* integrand, int resolusion);

private:
    PhotonField* m_field;
};
#endif