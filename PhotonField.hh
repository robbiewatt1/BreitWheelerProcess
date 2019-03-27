#ifndef PHOTON_FIELD
#define PHOTON_FIELD

#include <string>

class PhotonField
{
public:
    /* Constructor to generate a thermal isotropic field field. All input 
       paramters should be given in MeV mm ns. */
    PhotonField(double temp, double energyMin, double energyMax,
            int resolution = 200);

    /* Constructor to generate a photon field from a file. The energy axis
       should be in MeV and the density axis in mm^-3 MeV^-1 */
    PhotonField(const std::string fileName, int resolution);

    /* Destructor dealocates m_energyAxis and m_photonDensity arrays. */ 
    ~PhotonField();

    /* Following getter methods allow acces to the photon field 
       energy, angle and density variables. */
    int getResolution() const {return m_resolution;}

    double* getEnergy() const {return m_energy;}

    double* getTheta() const {return m_theta;}

    double* getPhi() const {return m_phi;}

    double*** getDensity() const {return m_density;}

    double getMaxEnergy() const {return m_energy[m_resolution-1];}

    double getMinEnergy() const {return m_energy[0];}

private:

    /* Reads the photon field from a file. Then interpolates onto the required
       grid. Energy values should be in MeV and density m^-3 MeV^-1. */
    void fileField(const std::string &fileName);

    /* Generates the photon field assuming a black body spectrum. The 
        temperature should be given in MeV. */
    void thermalField(double temp);

private:
    double* m_energy;
    double* m_theta;
    double* m_phi;
    double*** m_density;
    int     m_resolution;
};
#endif
