#include "BreitWheeler.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"

#include <cmath>
#include <fstream>

BreitWheeler::BreitWheeler(PhotonField* field, const G4String& name,
        G4ProcessType type):
G4VDiscreteProcess(name, type)
{
    m_field = field;
    m_rotatedIndex = new int [m_field->getResolution()];
    m_energyInt    = new double [m_field->getResolution()];
    m_comEnergy    = new double* [m_field->getResolution()];
    m_comEnergyInt = new double* [m_field->getResolution()];
    for (int i = 0; i < m_field->getResolution(); i++)
    {
        m_comEnergy[i] = new double [m_field->getResolution()];
        m_comEnergyInt[i] = new double [m_field->getResolution()];
    }
}

BreitWheeler::~BreitWheeler()
{
    delete [] m_rotatedIndex;
    delete [] m_energyInt;
    for (int i = 0; i < m_field->getResolution(); i++)
    {
        delete [] m_comEnergy[i];
        delete [] m_comEnergyInt[i];
    }
    delete [] m_comEnergy;
    delete [] m_comEnergyInt;
}

G4bool BreitWheeler::IsApplicable(const G4ParticleDefinition& particle)
{
    return (&particle == G4Gamma::Gamma());
}

G4double BreitWheeler::GetMeanFreePath(const G4Track& track, G4double,
         G4ForceCondition*)
{
    G4Material* aMaterial = track.GetMaterial();
    if(aMaterial->GetMaterialPropertiesTable()->GetConstProperty("Radiation")
            < 0.0) return 1e99;

    /* Get the interacting particle properties */
    const G4DynamicParticle *aDynamicGamma = track.GetDynamicParticle();
    double gammaEnergy = aDynamicGamma->GetKineticEnergy();
    G4ThreeVector gammaDirection = aDynamicGamma->GetMomentumDirection();
 //   double gammaTheat = std::acos(gammaDirection[2]);
 //   double gammaPhi   = std::atan(gammaDirection[1] / gammaDirection[0]);
  
    /* Find the roation matricies that rotate the gamma ray onto
       the z axis and return the gamma ray back */
    G4ThreeVector rotationAxis = gammaDirection.
            cross(G4ThreeVector(0, 0, 1)).unit();
    double rotationAngle = gammaDirection.angle(G4ThreeVector(0, 0, 1)); 
    G4RotationMatrix m_rotateForward = G4RotationMatrix(m_rotationAxis,
            m_rotationAngle);
    G4RotationMatrix m_rotateBackward = G4RotationMatrix(m_rotationAxis,
            m_rotationAngle);

    /* Get the photon field properties */
    int resolusion          = m_field->getResolution();
    double* photonEnergy    = m_field->getEnergy();
    double* photonTheta     = m_field->getTheta();
    double* photonPhi       = m_field->getTheta();
    double*** photonDensity = m_field->getDensity();

    double* energyInt = new double [resolusion];
    double* comInt    = new double [resolusion];
    double* comAxis   = new double [resolusion];
    double* phiInt    = new double [resolusion];
    // Integral over photon energy epsilon
    for (int i = 0; i < resolusion; i++)
    {
        // Integrate over s
        for (int j = 0; j < resolusion; j++)
        {
            // Integrate over phi
            for (int k = 0; k < resolusion; k++)
            {
                int[2] indexPrimne = RotateThetaPhi(photonTheta[j], photonPhi[k],
                        photonTheta, photonPhi, resolusion, resolusion);
                phiInt[k] = photonDensity[i][indexPrimne[0]][indexPrimne[1]];
            }

            comEnergy[j] = photonEnergy[i] * gammaEnergy * (1.0 - 
                    std::cos(photonAngle[j])) / (electron_mass_c2 
                    * electron_mass_c2);
            if (comEnergy[j] > 1.0)
            {
                comEnergyInt[j] = crossSection(comEnergy[j]) 
                        * comEnergy[j] * trapezium(photonPhi, phiInt,
                            resolusion);
            } else
            {
                comEnergyInt[j] = 0;
            }
        }
        if (photonEnergy[i] > electron_mass_c2 * electron_mass_c2 
                / gammaEnergy)
        {
            energyInt[i] = trapezium(comEnergy, comEnergyInt, 
                    resolusion) / (photonEnergy[i] * photonEnergy[i]);
        } else
        {
            energyInt[i] = 0;
        }
    }
    double meanPath = gammaEnergy * gammaEnergy / 
            (trapezium(photonEnergy, energyInt, resolusion) * pi 
                * classic_electr_radius * classic_electr_radius
                * electron_mass_c2 * electron_mass_c2);
    return meanPath;
}


G4VParticleChange* BreitWheeler::PostStepDoIt(const G4Track& aTrack,
        const G4Step& aStep)
{
    /*
    aParticleChange.Initialize(aTrack);
    const G4DynamicParticle *aDynamicGamma = aTrack.GetDynamicParticle();
    G4ThreeVector gammaMomentum = aDynamicGamma->GetMomentum();
    double gammaEnergy     = aDynamicGamma->GetKineticEnergy();
    double photonEnergy    = samplePhotonEnergy();
    double comEnergy       = sampleComEnergy(photonEnergy, gammaEnergy);
    double photonAnglePhi  = G4UniformRand() * 2.0 * pi;
    double photonAngleThe  = std::acos(1.0 - 2.0 * electron_mass_c2
            * electron_mass_c2 * comEnergy / (gammaEnergy * photonEnergy));
    */
    /* Find particles properties in the COM frame*/
 /*
    double pairEnergy   = 0.5 * std::sqrt(comEnergy);
    double pairMomentum = std::sqrt(pairEnergy * pairEnergy
                - electron_mass_c2 * electron_mass_c2);    
    double pairAngleThe = samplePairAngle(comEnergy);
    double pairAnglePhi = std::acos(2.0 * G4UniformRand() - 1.0);
    double pairPx = std::sin(pairAnglePhi) * std::cos(pairAngleThe)
            * pairMomentum;
    double pairPy = std::sin(pairAnglePhi) * std::sin(pairAngleThe)
            *pairMomentum;
    double pairPz = std::cos(pairAnglePhi) * pairMomentum;
    G4ThreeVector electronMomentum = G4ThreeVector(pairPx, pairPy, pairPz);
    G4ThreeVector positronMomentum = G4ThreeVector(-pairPx, -pairPy, -pairPz);
*/
    /* Apply lorentz transform into rorated lab frame */
/*
    double lorentzFact = (photonEnergy + gammaEnergy) / std::sqrt(comEnergy);
    double electronEnergy = lorentzFact * pairEnergy + std::sqrt(lorentzFact
            * lorentzFact - 1.0) * electronMomentum[2];
    double positronEnergy = lorentzFact * pairEnergy + std::sqrt(lorentzFact
            * lorentzFact - 1.0) * positronMomentum[2];
    electronMomentum[2] = std::sqrt(lorentzFact * lorentzFact - 1.0)
            * pairEnergy + lorentzFact * electronMomentum[2];
    positronMomentum[2] = std::sqrt(lorentzFact * lorentzFact - 1.0)
            * pairEnergy + lorentzFact * positronMomentum[2];
*/
    /* Apply rotation into frame with gamma along z axis */
 /*
    double thetaGamma = std::atan(photonEnergy * std::sin(photonAngleThe)
                / (gammaEnergy + photonEnergy * std::cos(photonAngleThe)));
    electronMomentum[0] = std::cos(thetaGamma) * electronMomentum[0]
            - std::sin(thetaGamma) * electronMomentum[2];
    positronMomentum[0] = std::cos(thetaGamma) * positronMomentum[0]
            - std::sin(thetaGamma) * positronMomentum[2];
    electronMomentum[2] = std::sin(thetaGamma) * electronMomentum[0]
            + std::cos(thetaGamma) * electronMomentum[2];
    positronMomentum[2] = std::sin(thetaGamma) * positronMomentum[0]
            + std::cos(thetaGamma) * positronMomentum[2];
*/
    /* Apply rotation around the gamma axis by -photonAnglePhi */
/*
    electronMomentum[0] = std::cos(photonAnglePhi) * electronMomentum[0]
            + std::sin(photonAnglePhi) * electronMomentum[1];
    positronMomentum[0] = std::cos(photonAnglePhi) * positronMomentum[0]
            + std::sin(photonAnglePhi) * positronMomentum[1];
    electronMomentum[1] = - std::sin(photonAnglePhi) * electronMomentum[0]
            + std::cos(photonAnglePhi) * electronMomentum[1];            
    positronMomentum[1] = - std::sin(photonAnglePhi) * positronMomentum[0]
            + std::cos(photonAnglePhi) * positronMomentum[1];

    /* Apply final rotation into simulation frame */
/*
    double thetaLabXZ = - std::atan(gammaMomentum[0] / gammaMomentum[2]);
    double thetaLabYZ = std::atan(gammaMomentum[1] / gammaMomentum[2] *
            std::sqrt(1 + gammaMomentum[0] * gammaMomentum[0]
                / (gammaMomentum[2] * gammaMomentum[2])));
    electronMomentum[0] = std::cos(thetaLabXZ) * electronMomentum[0]
            - std::sin(thetaLabXZ) * (- std::sin(thetaLabYZ)
                * electronMomentum[1] + std::cos(thetaLabYZ)
                * electronMomentum[2]);
    positronMomentum[0] = std::cos(thetaLabXZ) * positronMomentum[0]
            - std::sin(thetaLabXZ) * (- std::sin(thetaLabYZ)
                * positronMomentum[1] + std::cos(thetaLabYZ)
                * positronMomentum[2]);
    electronMomentum[1] = std::cos(thetaLabYZ) * electronMomentum[1]
            + std::sin(thetaLabYZ) * electronMomentum[2];
    positronMomentum[1] = std::cos(thetaLabYZ) * positronMomentum[1]
            + std::sin(thetaLabYZ) * positronMomentum[2];
    electronMomentum[2] = std::sin(thetaLabXZ) * electronMomentum[0]
            + std::cos(thetaLabXZ) * (- std::sin(thetaLabYZ)
                * electronMomentum[1] + std::cos(thetaLabYZ)
                * electronMomentum[2]);
    positronMomentum[2] = std::sin(thetaLabXZ) * positronMomentum[0]
            + std::cos(thetaLabXZ) * (- std::sin(thetaLabYZ)
                * positronMomentum[1] + std::cos(thetaLabYZ)
                * positronMomentum[2]);
*/
    /* Add new electron and positron and kill the gamma ray */
/*
    G4DynamicParticle* electron = new G4DynamicParticle(
            G4Electron::Electron(), electronMomentum);            
    G4DynamicParticle* positron = new G4DynamicParticle(
            G4Positron::Positron(), positronMomentum);
    aParticleChange.SetNumberOfSecondaries(2);
    aParticleChange.AddSecondary(electron);
    aParticleChange.AddSecondary(positron);
    aParticleChange.ProposeMomentumDirection(G4ThreeVector(0,0,0));
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    */
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

 double BreitWheeler::samplePhotonEnergy()
 {
    double maxDensity = *(std::max_element(m_energyInt, 
            &m_energyInt[m_field->getResolution()-1]));
    double randEnergy;
    double randDensity;
    double photonDensity;
    do 
    {
        randEnergy = m_field->getMinEnergy() + G4UniformRand() *
                (m_field->getMaxEnergy() - m_field->getMinEnergy());
        randDensity = G4UniformRand() * maxDensity;
        photonDensity = interpolate1D(m_field->getEnergy(), m_energyInt,
                m_field->getResolution(), randEnergy);

    } while (randDensity > photonDensity);
    return randEnergy;
}

double BreitWheeler::sampleComEnergy(double photonEnergy, double gammaEnergy)
{
    double fracIndex = arrayIndex(m_field->getEnergy(), photonEnergy,
            m_field->getResolution());
    int index  = fracIndex;
    double frac = fracIndex - index;

    double maxCon = photonEnergy * gammaEnergy / (electron_mass_c2
                * electron_mass_c2);
    double maxDensity = *(std::max_element(m_comEnergyInt[index],
            &m_comEnergyInt[index][m_field->getResolution()]));
    double randCom;
    double randDensity;
    double comDensity;
    do
    {
        randCom = 1.0 + (maxCon - 1.0) * G4UniformRand();
        randDensity = maxDensity * G4UniformRand();
        comDensity = (1.0 - frac) * interpolate1D(m_comEnergy[index],
                m_comEnergyInt[index], m_field->getResolution(), randCom)
            + (frac) * interpolate1D(m_comEnergy[index+1],
                m_comEnergyInt[index+1], m_field->getResolution(), randCom);
    } while (randDensity > comDensity);
    return randCom;
}

double BreitWheeler::samplePairAngle(double comEnergy)
{
    double maxDensity = diffCrossSection(comEnergy, 2.0 * pi);
    double randAngle;
    double randDensity;
    double angleDensity;
    do
    {
        randAngle    = G4UniformRand() * 2.0 * pi;
        randDensity  = G4UniformRand() * maxDensity;
        angleDensity = diffCrossSection(comEnergy, randAngle);
    } while (randDensity > angleDensity);
    return randAngle;
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

double BreitWheeler::interpolate1D(double* samplePoints, double* sampleValues,
        int sampleSize, double queryPoint)
{
    double queryValue;
    if (queryPoint < samplePoints[0])
    {
        queryValue = sampleValues[0] + (sampleValues[1] - sampleValues[0]) 
                * ((queryPoint - samplePoints[0]) / (samplePoints[1] 
                            - samplePoints[0]));

    } else if (queryPoint > samplePoints[sampleSize-1])
    {
        int end = sampleSize - 1;
        queryValue = sampleValues[end] + (sampleValues[end] - sampleValues[end-1])
                * (queryPoint - samplePoints[end-1]) / (samplePoints[end] 
                            - samplePoints[end-1]);
    } else
    {
        int lowIndex(0), highIndex(0);
        for (int i = 0; i < sampleSize; i++)
        {
            if (queryPoint < samplePoints[i])
            {
                lowIndex = i - 1;
                highIndex = i;
                break;
            }
        }
        queryValue = sampleValues[lowIndex] + (queryPoint - samplePoints[lowIndex])
                * (sampleValues[highIndex] - sampleValues[lowIndex])
                / (samplePoints[highIndex] - samplePoints[lowIndex]);
    }
    return queryValue;
}

double BreitWheeler::arrayIndex(double* samplePoints, double queryPoint,
        int sampleSize)
{
    if (queryPoint < samplePoints[0] || queryPoint >
            samplePoints[sampleSize-1])
    {
        std::cerr << "Error: Sample point out of bounds!" << std::endl;
        exit(1);
    }
    double index(0);
    for (int i = 0; i < sampleSize; i++)
    {
        if (samplePoints[i] > queryPoint)
        {
            index = (i - 1) + (queryPoint - samplePoints[i-1]) / 
                    (samplePoints[i] - samplePoints[i-1]);
            break;
        }
    }
    return index;
}

int* RotateThetaPhi(double theta, double phi, double* thetaAxis,
        double* phiAxis, int thetaResolusion, phiResolusion);
{
    // find components in cartesian
    G4ThreeVector vector = G4ThreeVector(std::sin(theta) * std::cos(phi),
            std::sin(theta) * std::sin(phi), std::cos(theta));
    G4ThreeVector vectorPrime = m_rotationMatrix(vector1);
    double thetaPrime = std::acos(vector2[2]);
    double phiPrime   = std::atan(vector2[1] / vector2[0]);

    int[2] index = {arrayIndex(thetaAxis, thetaPrime, thetaResolusion),
            arrayIndex(phiaAxis, phiPrime, phiResolusion)}

    return index;
}