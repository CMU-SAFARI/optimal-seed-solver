#include <iostream>
#include <string>
#include "RefDB.h"
#include <cstdlib>

using namespace std;

int main(int argc, const char* argv[]) {
	RefDB refDB;
	if (argc != 2) {
		cerr << "Please also put the reference file name" << endl;
		exit(1);
	}

	refDB.loadRefFile(argv[1]);

	int chromo_number = -1;
	int location;
	int length;

	cout << "Num of chromosomes: " << refDB.getNumOfChromo() << endl;

	while (1) {
		cout << "please input chromo number: " << endl;
		int tempChromoNum;
		cin >> tempChromoNum;
		if (tempChromoNum != chromo_number) {
			chromo_number = tempChromoNum;
			refDB.loadChromo(chromo_number);
		}

		cout << "please input location: " << endl;
		cin >> location;
		cout << "please input length: " << endl;
		cin >> length;

		cout << refDB.getRefSeq(location, length) << endl;
	}
	return 0;
}
