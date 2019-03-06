#include "PhotonField.hh"

#include <iostream>

int main(int argc, char* argv[])
{
	PhotonField field = PhotonField(0.1, 0.001, 2);

	double* energy = field.getEnergy();
	double** density = field.getDensity();
	int resolution = field.getResolution();

	for (int i = 0; i < resolution; ++i)
	{
		std::cout << energy[i] << "\t" << density[i] << std::endl;
	}

	return 0;
}