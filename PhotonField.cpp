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
    m_angle   = new double [m_resolution];
    m_density = new double* [m_resolution];

    double energyDelta = (energyMax - energyMin) / m_resolution;
    double angleDelta  = (2.0 * pi) / m_resolution;
    for (int i = 0; i < m_resolution; i++)
    {
        m_angle[i]   = i * angleDelta;
        m_energy[i]  = i * energyDelta + energyMin;
        m_density[i] = new double [m_resolution];
    }
    thermalField(temp);
}

PhotonField::PhotonField(const std::string fileName, int resolution):
m_resolution(resolution)
{
    m_energy  = new double [m_resolution];
    m_angle   = new double [m_resolution];
    m_density = new double* [m_resolution];
    for (int i = 0; i < m_resolution; i++)
    {
        m_density[i] = new double [m_resolution];
    }
    fileField(fileName);
}

PhotonField::~PhotonField()
{
    delete [] m_energy;
    delete [] m_angle;
    for (int i = 0; i < m_resolution; ++i)
    {
        delete [] m_density[i];
    }
    delete [] m_density;
}

void PhotonField::fileField(const std::string &fileName)
{
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
}

void PhotonField::thermalField(double temp)
{
    for (int i = 0; i < m_resolution; i++)
    {
        for (int j = 0; j < m_resolution; j++)
        {
            m_density[i][j] = 1.0 / (hbar_Planck * hbar_Planck * hbar_Planck
                    * c_light * c_light * c_light * pi * pi) * m_energy[i]
                    * m_energy[i] / std::exp(m_energy[i] / temp - 1.0)
                    / m_resolution;
        }
    }
}

/*
double PhotonField::interpolate1D(double* samplePoints, double* sampleValues,
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
*/