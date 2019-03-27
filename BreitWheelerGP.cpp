#include "BreitWheelerGP.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"

#include <cmath>
#include <fstream>

BreitWheelerGP::BreitWheelerGP(PhotonField* field, const G4String& name,
        G4ProcessType type):
G4VDiscreteProcess(name, type), trainCount(0)
{
    m_field = field;
    m_gp_optimiser.init();
    Eigen::VectorXd params(m_gp_gausProc.covf().get_param_dim());
    params << -1, -1, -1;
    m_gp_gausProc.covf().set_loghyper(params);
}

BreitWheelerGP::~BreitWheelerGP()
{
    if (m_gp_save == true)
    {
        m_gp_gausProc.write("./load.gp");
    }
    delete [] m_gp_input;
    delete [] m_gp_mfp;
}

G4bool BreitWheelerGP::IsApplicable(const G4ParticleDefinition& particle)
{
    return (&particle == G4Gamma::Gamma());
}

void BreitWheelerGP::setParamsGP(bool save, int trainSize, double errorMax)
{
    m_gp_save = save;
    m_gp_trainSize = trainSize;
    m_gp_errorMax = errorMax;
    m_gp_input = new double* [trainSize];
    for (int i = 0; i < trainSize; ++i)
    {
        m_gp_input[i] = new double [3];
    }
    m_gp_mfp = new double [trainSize];
}

G4double BreitWheelerGP::GetMeanFreePath(const G4Track& track, G4double,
         G4ForceCondition*)
{
    G4Material* aMaterial = track.GetMaterial();
    if(aMaterial->GetMaterialPropertiesTable()->GetConstProperty("Radiation")
            < 0.0) return 1e99;

    /* Get the interacting particle properties */
    const G4DynamicParticle *aDynamicGamma = track.GetDynamicParticle();
    double gammaEnergy = aDynamicGamma->GetKineticEnergy();
    G4ThreeVector gammaDirection = aDynamicGamma->GetMomentumDirection();
    double thetaPrime = std::acos(gammaDirection[2]);
    double phiPrime   = std::atan(gammaDirection[1] / gammaDirection[0]);


    /* Generate input for GP with close to normalised values */
    double input[] = {gammaEnergy / 1e3, thetaPrime / (2.0 * pi),
            phiPrime / pi};
    double variance = m_gp_gausProc.var(input);
    double meanPath;

    if (variance < m_gp_errorMax)
    {  
        /* We use the regressed value from the GP. Need to convert 
           the mean free path back to proper units */
        meanPath = std::pow(10.0, (100 * m_gp_gausProc.f(input)));
    } else
    {
        /* We use the slow full method of integration */

        /* Find the roation matricies that rotate the gamma ray onto
           the z axis and return the gamma ray back */
        G4ThreeVector rotationAxis = gammaDirection.
                cross(G4ThreeVector(0, 0, 1)).unit();
        double rotationAngle = gammaDirection.angle(G4ThreeVector(0, 0, 1)); 
        m_rotateForward = G4RotationMatrix(rotationAxis,
                rotationAngle);

        /* Get the photon field properties */
        int resolusion          = m_field->getResolution();
        double* photonEnergy    = m_field->getEnergy();
        double* photonTheta     = m_field->getTheta();
        double* photonPhi       = m_field->getPhi();
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
                    int thetaIndex, phiIndex;
                    RotateThetaPhi(photonTheta[j], photonPhi[k], thetaIndex,
                            phiIndex);
                    phiInt[k] = photonDensity[i][thetaIndex][phiIndex];
                }

                comAxis[j] = photonEnergy[i] * gammaEnergy * (1.0 - 
                        std::cos(photonTheta[j])) / (electron_mass_c2 
                        * electron_mass_c2);
                if (comAxis[j] > 1.0)
                {
                    comInt[j] = crossSection(comAxis[j]) 
                            * comAxis[j] * trapezium(photonPhi, phiInt,
                                resolusion);
                } else
                {
                    comInt[j] = 0;
                }
            }
            if (photonEnergy[i] > electron_mass_c2 * electron_mass_c2 
                    / gammaEnergy)
            {
                energyInt[i] = trapezium(comAxis, comInt, 
                        resolusion) / (photonEnergy[i] * photonEnergy[i]);
            } else
            {
                energyInt[i] = 0;
            }
        }
        meanPath = gammaEnergy * gammaEnergy / (trapezium(photonEnergy,
                energyInt, resolusion) * pi * classic_electr_radius
                * classic_electr_radius * electron_mass_c2 * electron_mass_c2);

        /* add the new point to the training set, These are roughly normilised
           to the maximum values expected */
        m_gp_input [trainCount][0] = gammaEnergy / 1000;
        m_gp_input [trainCount][1] = thetaPrime / (2.0 * pi);
        m_gp_input [trainCount][2] = phiPrime / pi;
        m_gp_mfp[trainCount] = std::log10(meanPath) / 100;
        trainCount++;

        delete [] energyInt;
        delete [] comInt;
        delete [] comAxis;
        delete [] phiInt;
    }

    /* Check if we need to carry out the training again */
    if(trainCount == m_gp_trainSize)
    {
        for (int i = 0; i < m_gp_trainSize; i++)
        {
            double input[] = {m_gp_input[i][0], m_gp_input[i][1],
                    m_gp_input[i][2]};
            m_gp_gausProc.add_pattern(input, m_gp_mfp[trainCount]);
        }
        m_gp_optimiser.maximize(&m_gp_gausProc, 50, 0);
        trainCount = 0;
    }
    return meanPath;
}


G4VParticleChange* BreitWheelerGP::PostStepDoIt(const G4Track& aTrack,
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
    */
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

double BreitWheelerGP::crossSection(double comEnergy)
{
    double beta = std::sqrt(1.0 - 1.0 / comEnergy);
    return (1.0 - beta * beta) * ((3.0 - beta * beta * beta * beta)
            * std::log((1.0 + beta) / (1.0 - beta)) - 2.0 * beta 
            * (2.0 - beta * beta));
}

double BreitWheelerGP::diffCrossSection(double comEnergy, double theta)
{
    double beta = std::sqrt(1.0 - 1.0 / comEnergy);
    double sinT = std::sin(theta);
    double cosT = std::cos(theta);
    return (beta / comEnergy) * (1.0 + 2.0 * beta * beta * sinT * sinT
            - beta * beta * beta * beta - beta * beta * beta * beta
            * sinT * sinT * sinT * sinT) / ((1.0 - beta * beta * cosT * cosT)
            * (1.0 - beta * beta * cosT * cosT));
}

 double BreitWheelerGP::samplePhotonEnergy()
 {
    /*
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
*/
}

double BreitWheelerGP::sampleComEnergy(double photonEnergy, double gammaEnergy)
{
    /*
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
*/
}

double BreitWheelerGP::samplePairAngle(double comEnergy)
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

 double BreitWheelerGP::trapezium(double* variable, double* integrand,
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

double BreitWheelerGP::interpolate1D(double* samplePoints, double* sampleValues,
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

double BreitWheelerGP::arrayIndex(double* samplePoints, double queryPoint,
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

void BreitWheelerGP::RotateThetaPhi(double theta, double phi, int thetaIndex,
        int phiIndex)
{
    // find components in cartesian
    G4ThreeVector vector = G4ThreeVector(std::sin(theta) * std::cos(phi),
            std::sin(theta) * std::sin(phi), std::cos(theta));
    G4ThreeVector vectorPrime = m_rotateForward(vector);
    double thetaPrime = std::acos(vectorPrime[2]);
    double phiPrime   = std::atan(vectorPrime[1] / vectorPrime[0]);

    thetaIndex = arrayIndex(m_field->getTheta(), thetaPrime,
            m_field->getResolution());
    phiIndex   = arrayIndex(m_field->getPhi(), phiPrime,
            m_field->getResolution());
}