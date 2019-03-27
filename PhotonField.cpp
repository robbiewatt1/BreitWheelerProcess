#include "PhotonField.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <fstream>
#include <iostream>
#include <cmath>

PhotonField::PhotonField(double temp, double energyMin, double energyMax,
        int resolution):
m_resolution(resolution)
{
    m_energy  = new double [m_resolution];
    m_theta   = new double [m_resolution];
    m_phi     = new double [m_resolution];
    m_density = new double** [m_resolution];

    double energyDelta = (energyMax - energyMin) / m_resolution;
    double thetaDelta  = 2.0 * pi / m_resolution;
    double phiDelta  = pi / m_resolution;
    for (int i = 0; i < m_resolution; i++)
    {
        m_theta[i]    = i * thetaDelta;
        m_phi[i]      = i * phiDelta;
        m_energy[i] = i * energyDelta + energyMin;
        m_density[i] = new double* [m_resolution];
        for (int j = 0; j < m_resolution; j++)
        {
            m_density[i][j] = new double [m_resolution];
        }
    }
    thermalField(temp);
}

PhotonField::PhotonField(const std::string fileName, int resolution):
m_resolution(resolution)
{
    m_energy  = new double [m_resolution];
    m_theta   = new double [m_resolution];
    m_phi     = new double [m_resolution];
    m_density = new double** [m_resolution];
    for (int i = 0; i < m_resolution; i++)
    {
        m_density[i] = new double* [m_resolution];
        for (int j = 0; j < m_resolution; j++)
        {
            m_density[i][j] = new double [m_resolution];
        }
    }
    fileField(fileName);
}

PhotonField::~PhotonField()
{
    delete [] m_energy;
    delete [] m_theta;
    delete [] m_phi;
    for (int i = 0; i < m_resolution; i++)
    {
        for (int j = 0; j < m_resolution; j++)
        {
            delete [] m_density[i][j];
        }
        delete [] m_density[i];
    }
    delete [] m_density;
}

void PhotonField::fileField(const std::string &fileName)
{
    /*
    std::ifstream file(fileName);
    if (!file)
    {
        std::cerr << "File: " << fileName << " not found." << std::endl;
        std::cerr << "Aborting program." << std::endl;
        exit(1);
    }

    int fileResolution;
    double energyMin, energyMax, angleMin, angleMax;
    file >> fileResolution >> energyMin >> energyMax >> angleMin >> angleMax;

    double energyDelta = (energyMax - energyMin) / m_resolution;
    double angleDelta  = (angleMax  - angleMin)  / m_resolution;
    for (int i = 0; i < m_resolution; i++)
    {
        m_energy[i] = i * energyDelta + energyMin;
        m_angle[i]  = i * angleDelta + angleMin;
        for (int j = 0; j < m_resolution; ++j)
        {
            file >> m_density[i][j];
        }
    }
    */
}

void PhotonField::thermalField(double temp)
{
    for (int i = 0; i < m_resolution; i++)
    {
        for (int j = 0; j < m_resolution; j++)
        {
            for (int k = 0; j < m_resolution; k++)
            {            
                m_density[i][j][k] = 1.0 / (hbar_Planck * hbar_Planck
                    * hbar_Planck * c_light * c_light * c_light * pi * pi)
                    * m_energy[i] * m_energy[i] / std::exp(m_energy[i]
                        / temp - 1.0) / (m_resolution * m_resolution *
                        2.0 * pi * pi);
            }
        }
    }
}