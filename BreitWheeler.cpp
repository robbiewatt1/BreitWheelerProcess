
#include "BreitWheeler.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

BreitWheeler::BreitWheeler(PhotonField* field, const G4String& name,
        G4ProcessType type):
G4VDiscreteProcess(name, type)
{
    m_field = field;
}

BreitWheeler::~BreitWheeler()
{
}

G4bool BreitWheeler::IsApplicable(const G4ParticleDefinition& particle)
{
    return ( &particle == G4Gamma::Gamma() );
}

G4double BreitWheeler::GetMeanFreePath(const G4Track& track, G4double,
         G4ForceCondition*)
{
    // The interacting particle properties
    const G4DynamicParticle *aDynamicGamma = aTrack.GetDynamicParticle();
    G4double gammaEnergy = aDynamicGamma->GetKineticEnergy();
    
    // The photon field properties
    int resolusion = getResolution();
    double* photonEnergy   = m_field->getEnergy();
    double* photonAngle    = m_field->getAngle();
    double** photonDensity = m_field->getDensity();

    // rotate the photon distrobution function into gamma axis

    double* energyInt = new double [resolusion];
    double* s = new double [resolusion];
    double* sInt = new double [resolusion];
    // Integral over photon energy epsilon
    for (int i = 0; i < resolusion; i++)
    {
        // integrate over s, such that s > 1
        for (int j = 0; j < resolusion; j++)
        {
           s[j] = photonEnergy[i] * gammaEnergy * (1.0 - std::cos(photonAngle[j]))
                / (electron_mass_c2 * electron_mass_c2);
            if (s[j] > 1.0)
            {
                sInt[j] = crossSection(s[j]) * photonDensity[i][j] * s[j];
            } else
            {
                sInt[j] = 0;
            }
        }
        energyInt[i] = trapezium(s, sInt, resolusion) / (photonEnergy[i] 
                    * photonEnergy[i]);
    }

    delete [] s;
    delete [] sInt;
}

G4VParticleChange* BreitWheeler::PostStepDoIt(const G4Track& aTrack,
        const G4Step& aStep)
{
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

double BreitWheeler::crossSection(double s)
{
    double beta = std::sqrt(1.0 - 1.0 / s);
    return (1.0 - beta * beta) * ((3.0 - beta * beta * beta * beta)
            * std::log((1.0 + beta) / (1.0 - beta)) - 2.0 * beta 
            * (2.0 - beta * beta));
}

 double BreitWheeler::trapezium(double* variable, double* integrand,
        int resolusion)
 {
    double result(0);
    for (int i = 0; i < resolusion-1; i++)
    {
        result += (variable[i+1] - variable[i]) * (integrand[i] 
                + integrand[i+1]) / 2.0;
    }
    return result;
 }