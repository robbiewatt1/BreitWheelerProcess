#include "PhotonField.hh"
#include "BreitWheeler.hh"
#include "G4Track.hh"
#include <iostream>
#include "G4ForceCondition.hh"

int main(int argc, char* argv[])
{
	PhotonField* field = new PhotonField(0.1, 0.001, 2);
	BreitWheeler* process = new BreitWheeler(field);
	G4Track a;
	G4double b;
	G4ForceCondition* c;

	std::cout << process->GetMeanFreePath(a, b, c) << std::endl;

	return 0;
}