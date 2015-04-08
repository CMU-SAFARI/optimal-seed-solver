#include "analyzeBasic.h"
#include <cassert>
#include <fstream>

AnalyzeBasic::AnalyzeBasic() {}

AnalyzeBasic::~AnalyzeBasic() {}

void AnalyzeBasic::generateHash(string refName, unsigned int seedLength) {
	refDB.loadRefFile(refName);

	this->seedLength = seedLength;

	string seed = string(seedLength, 'N');
	
	freqTable.resize( (1 << (2 * seedLength) ), 0);

	for (unsigned int chromoCount = 0; chromoCount < refDB.getNumOfChromo(); chromoCount++) {
		refDB.loadChromo(chromoCount);
		for (unsigned int refPos = 0; refPos < refDB.getChromoLength() - seedLength; refPos++) {
			seed = refDB.getRefSeq(refPos, seedLength);
			if (seed.find('N') == string::npos && seed.find('M') == string::npos && seed.find('R') == string::npos) {
				unsigned int hash = hashGen.calculateHash(seed, seed.length() );
				freqTable[hash]++;
			}
		}
	}	
}

void AnalyzeBasic::storeHash(string tableName) {
	ofstream storeFile;
	storeFile.open(tableName.c_str() );

	storeFile << seedLength << endl;

	for (int i = 0; i < freqTable.size(); i++)
		storeFile << freqTable[i] << endl;

	storeFile.close();
}

void AnalyzeBasic::loadHash(string tableName) {
	ifstream storeFile;
	storeFile.open(tableName.c_str() );

	storeFile >> seedLength;

	freqTable.resize( (1 << (2 * seedLength) ), 0);
	
	for (int i = 0; i < freqTable.size(); i++)
		storeFile >> freqTable[i];

	storeFile.close();
}

unsigned int AnalyzeBasic::checkFreq(string seed) {
	assert (seed.length() >= seedLength);
	unsigned int hash = hashGen.calculateHash(seed, seed.length() );
	return freqTable[hash];
}

unsigned int AnalyzeBasic::getLength() {
	return seedLength;
}
