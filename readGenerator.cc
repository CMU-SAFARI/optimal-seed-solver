/*
 * readGenerator.cc
 *
 *  Created on: Aug 12, 2013
 *      Author: hxin
 */

#include "RefDB.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>

using namespace std;

int main(int argc, char* argv[]) {\
	if (argc != 6) {
		cout << "usage: >readGenerator refFile readFile chromoNum readLen readNum" << endl;
		return 1;
	}
	ofstream readFile(argv[2]);

	RefDB refDB;
	refDB.loadRefFile(argv[1]);
	refDB.loadChromo(atoi(argv[3]) );

	int readLen = atoi(argv[4]);
	int readNum = atoi(argv[5]);

	srand(1);
	string read;
	int location;

	for (int i = 0; i < readNum; i++) {
		do {
			location = rand() % (refDB.getChromoLength() - readLen);
			read = refDB.getRefSeq(location, readLen);
		}
		while (read.find('N') != string::npos);

//		cout << location << " "<< read << endl;
		readFile << location << " "<< read << endl;
	}

	readFile.close();
}


