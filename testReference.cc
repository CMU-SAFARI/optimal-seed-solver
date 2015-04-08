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

	for (int i = 0; i < refDB.getNumOfChromo(); i++) {
		refDB.loadChromo(i);
		for (int j = 0; j < refDB.getChromoLength(); j++) {
			char base = refDB.getRefBase(j, 0);
			if (base != 'A' && base != 'C' && base != 'G' && base != 'T' && base != 'N') {
				cout << i << "-" << j << ": " << base << endl;
				cout << refDB.getRefSeq(j - 5, 11) << endl;
			}

			if (base == 'N') {
				if (refDB.getRefBase(j+1, 0) != 'N') {
					cout << i << "-" << j << ": N" << endl;
					cout << refDB.getRefSeq(j - 5, 11) << endl;
					j++;
				}
				else if (refDB.getRefBase(j+2, 0) != 'N') {
					cout << i << "-" << j << ": NN" << endl;
					cout << refDB.getRefSeq(j - 5, 12) << endl;
					j = j + 2;
				}
				else if (refDB.getRefBase(j+3, 0) != 'N') {
					cout << i << "-" << j << ": NNN" << endl;
					cout << refDB.getRefSeq(j - 5, 13) << endl;
					j = j + 3;
				}
				else if (refDB.getRefBase(j+4, 0) != 'N') {
					cout << i << "-" << j << ": NNNN" << endl;
					cout << refDB.getRefSeq(j - 5, 14) << endl;
					j = j + 4;
				}
				else if (refDB.getRefBase(j+5, 0) != 'N') {
					cout << i << "-" << j << ": NNNN" << endl;
					cout << refDB.getRefSeq(j - 5, 15) << endl;
					j = j + 5;
				}
				else {
					while (j < refDB.getChromoLength() && refDB.getRefBase(j+1, 0) == 'N' )
						j++;
				}
			}
		}
	
	}

	return 0;
}
