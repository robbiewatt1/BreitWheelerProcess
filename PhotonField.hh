#ifndef PHOTON_FIELD
#define PHOTON_FIELD

#include <string>

class PhotonField
{
public:
    /* Constructor to generate a thermal isotropic field field. All input 
       paramters should be given in MeV mm ns. */
    PhotonField(double temp, double energyMin, double energyMax,
            int resolution = 1000);

    /* Constructor to generate a photon field from a file. The energy axis
       should be in MeV and the density axis in mm^-3 MeV^-1 */
    PhotonField(const std::string fileName, int resolution);

    /* Destructor dealocates m_energyAxis and m_photonDensity arrays. */ 
    ~PhotonField();

    /* Following getter methods allow acces to the photon field 
       energy, angle and density variables. */
    int getResolution() const {return m_resolution;}

    double* getEnergy() const {return m_energy;}

    double* getAngle() const {return m_angle;}

    double** getDensity() const {return m_density;}

private:

    /* Reads the photon field from a file. Then interpolates onto the required
       grid. Energy values should be in MeV and density m^-3 MeV^-1. */
    void fileField(const std::string &fileName);

    /* Generates the photon field assuming a black body spectrum. The 
        temperature should be given in MeV. */
    void thermalField(double temp);

    /* First order 1D intepolation method, samplePoints and sampleValues are 
       the know points and values and queryPoint is where data is wanted. This
       method is used to get fields are on the required grid. */
    double interpolate1D(double* samplePoints, double* sampleValues,
            int sampleSize, double queryPoint);

private:
    double* m_energy;
    double* m_angle;
    double** m_density;
    int     m_resolution;
};
#endif
