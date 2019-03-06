#include "BreitWheeler.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Gamma.hh"

#include <cmath>

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
    return (&particle == G4Gamma::Gamma());
}

G4double BreitWheeler::GetMeanFreePath(const G4Track& track, G4double,
         G4ForceCondition*)
{
    // The interacting particle properties
    const G4DynamicParticle *aDynamicGamma = track.GetDynamicParticle();
    double gammaEnergy = aDynamicGamma->GetKineticEnergy();
    double gammaAngle  = aDynamicGamma->GetMomentumDirection().
            angle(G4ThreeVector(0, 0, 1));

    // The photon field properties
    int resolusion         = m_field->getResolution();
    double* photonEnergy   = m_field->getEnergy();
    double* photonAngle    = m_field->getAngle();
    double** photonDensity = m_field->getDensity();

    /* Here a rotation is found that gives the photon distrobution
       function in the fame with the gamma direction along the z
       axis. This gives a new index for the photon density 
       angle dimension */
    int angleIndex(0);
    for (int i = 0; i < resolusion-1; ++i)
    {
        if (photonAngle[i+1] > gammaAngle)
        {
            if (photonAngle[i] > 2.0 * gammaAngle - photonAngle[i+1])
            {
                angleIndex = i;
            } else
            {
                angleIndex = i + 1;
            }
            break;
        }
    }
    int* rotatedIndex = new int [resolusion];
    for (int i = 0; i < resolusion; ++i)
    {
        if (i >= angleIndex)
        {
            rotatedIndex[i] = i - angleIndex;
        } else
        {
            rotatedIndex[i] = resolusion - (angleIndex - i);
        }
    }

    double* energyInt = new double [resolusion];
    double* comEnergy = new double [resolusion];
    double* comEnergyInt = new double [resolusion];
    // Integral over photon energy epsilon
    for (int i = 0; i < resolusion; i++)
    {
        // integrate over s
        for (int j = 0; j < resolusion; j++)
        {
           comEnergy[j] = photonEnergy[i] * gammaEnergy * (1.0 - 
                    std::cos(photonAngle[j])) / (electron_mass_c2 
                    * electron_mass_c2);
            if (comEnergy[j] > 1.0)
            {
                comEnergyInt[j] = crossSection(comEnergy[j])
                    * photonDensity[i][rotatedIndex[j]] * comEnergy[j];
            } else
            {
                comEnergyInt[j] = 0;
            }
        }
        if (photonEnergy[i] > electron_mass_c2 * electron_mass_c2 
                / gammaEnergy)
        {
            energyInt[i] = trapezium(comEnergy, comEnergyInt, resolusion)
                    / (photonEnergy[i] * photonEnergy[i]);
        } else
        {
            energyInt[i] = 0;
        }
    }
    double meanPath = gammaEnergy * gammaEnergy / 
            (trapezium(photonEnergy, energyInt, resolusion) * pi 
                * classic_electr_radius * classic_electr_radius
                * electron_mass_c2 * electron_mass_c2);
    
    delete [] rotatedIndex;
    delete [] energyInt;
    delete [] comEnergy;
    delete [] comEnergyInt;
   
    return meanPath;
}

G4VParticleChange* BreitWheeler::PostStepDoIt(const G4Track& aTrack,
        const G4Step& aStep)
{
    aParticleChange.Initialize(aTrack);

    
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

double BreitWheeler::crossSection(double comEnergy)
{
    double beta = std::sqrt(1.0 - 1.0 / comEnergy);
    return (1.0 - beta * beta) * ((3.0 - beta * beta * beta * beta)
            * std::log((1.0 + beta) / (1.0 - beta)) - 2.0 * beta 
            * (2.0 - beta * beta));
}

double BreitWheeler::diffCrossSection(double comEnergy, double theta)
{
    double beta = std::sqrt(1.0 - 1.0 / comEnergy);
    double sinT = std::sin(theta);
    double cosT = std::cos(theta);
    return (beta / comEnergy) * (1.0 + 2.0 * beta * beta * sinT * sinT
            - beta * beta * beta * beta - beta * beta * beta * beta
            * sinT * sinT * sinT * sinT) / ((1.0 - beta * beta * cosT * cosT)
            * (1.0 - beta * beta * cosT * cosT));

}

 double BreitWheeler::trapezium(double* variable, double* integrand,
        int resolusion)
 {
    double result(0);
    for (int i = 0; i < resolusion-1; i++)
    {
        result += std::abs(variable[i+1] - variable[i]) * (integrand[i] 
                + integrand[i+1]) / 2.0;
    }
    return result;
 }