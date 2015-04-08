#include "spacedSeedSolver.h"
#include <fstream>

SpacedSeedSolver::SpacedSeedSolver() {}

SpacedSeedSolver::~SpacedSeedSolver() {}

void SpacedSeedSolver::generateHash(string refName, string mask) {
	refDB.loadRefFile(refName);

	this->mask = mask;

	unsigned int oneCount = 0;
	for (int i = 0; i < mask.length(); i++) {
		if (mask[i] == '1')
			oneCount++;
	}

	freqTable.resize( (1 << (2 * oneCount) ), 0);

	seed = string(oneCount, 'N');

	for (unsigned int chromoCount = 0; chromoCount < refDB.getNumOfChromo(); chromoCount++) {
		refDB.loadChromo(chromoCount);
		for (unsigned int refPos = 0; refPos < refDB.getChromoLength() - mask.length(); refPos++) {
			string spacedSeed = refDB.getRefSeq(refPos, mask.length() );
			if (spacedSeed.find('N') == string::npos && spacedSeed.find('M') == string::npos && spacedSeed.find('R') == string::npos) {
				oneCount = 0;
				for (int i = 0; i < spacedSeed.length(); i++) {
				   if (mask[i] == '1') {
						seed[oneCount] = spacedSeed[i];
						oneCount++;
				   }
				}
				unsigned int hash = hashGen.calculateHash(seed, seed.length() );
				freqTable[hash]++;
			}
		}
	}	
}

void SpacedSeedSolver::storeHash(string tableName) {
	ofstream storeFile;
	storeFile.open(tableName.c_str() );

	storeFile<< mask << endl;

	for (int i = 0; i < freqTable.size(); i++)
		storeFile << freqTable[i] << endl;

	storeFile.close();
}

void SpacedSeedSolver::loadHash(string tableName) {
	ifstream storeFile;
	storeFile.open(tableName.c_str() );

	storeFile >> mask;

	unsigned int oneCount = 0;
	for (int i = 0; i < mask.length(); i++)
		if (mask[i] == '1')
			oneCount++;

	freqTable.resize( (1 << (2 * oneCount) ), 0);
	
	seed = string(oneCount, 'N');

	for (int i = 0; i < freqTable.size(); i++)
		storeFile >> freqTable[i];

	storeFile.close();
}

unsigned int SpacedSeedSolver::solveDNA(string DNA, unsigned int seedNum) {
	unsigned int totalFreq = 0;
	
	for (int i = 0; i < seedNum; i++) {
		string spacedSeed = DNA.substr(i * mask.length(), mask.length() );
		unsigned int oneCount = 0;
		for (int i = 0; i < spacedSeed.length(); i++) {
		   if (mask[i] == '1') {
				seed[oneCount] = spacedSeed[i];
				oneCount++;
		   }
		}
		unsigned int hash = hashGen.calculateHash(seed, seed.length() );
		totalFreq += freqTable[hash];
	}
	return totalFreq;
}

