	#include "predictorSolver.h"
	#include <fstream>
	#include <cassert>

	PredictorSolver::PredictorSolver() {
		hashLength = 0;
		numThresholds = 0;
		numTotalEntries = 0;
		
		table = NULL;
		thresholds = NULL;
	}

	PredictorSolver::~PredictorSolver() {
		if (table != NULL) {
			delete [] table[0];
			delete [] table;
			table = NULL;
		}

		if (thresholds != NULL) {
			delete [] thresholds;
			thresholds = NULL;
		}
	}

	void PredictorSolver::init(unsigned int hashLength, unsigned int numThresholds, unsigned int* thresholds) {
		this->hashLength = hashLength;
		this->numThresholds = numThresholds;
		numTotalEntries = 1 << (2 * hashLength);

		this->thresholds = new unsigned int [numThresholds];
		table = new unsigned int* [numTotalEntries];
		table[0] = new unsigned int [numTotalEntries * numThresholds];

		for (int i = 1; i < numTotalEntries; i++)
			table[i] = table[i-1] + numThresholds;

		for (int i = 0; i < numThresholds; i++)
			this->thresholds[i] = thresholds[i];
	}

	void PredictorSolver::generate(string treeFileName, unsigned int numThresholds, unsigned int * thresholds) {
		HashTree treeGen;
		treeGen.loadTree(treeFileName);

		init(treeGen.getHashLength(), numThresholds, thresholds);

		for (unsigned int i = 0; i < numTotalEntries; i++) {
			for (unsigned int j = 0; j < numThresholds; j++)
				table[i][j] = treeGen.estimateFreqLength(i, thresholds[j]);
		}
	}

	void PredictorSolver::storeTable(string tableFileName) {
		ofstream output;
		output.open(tableFileName.c_str(), ios::out | ios::binary);

		output.write( (char*) &hashLength, sizeof(unsigned int) );
		output.write( (char*) &numThresholds, sizeof(unsigned int) );

		for (int i = 0; i < numThresholds; i++) {
			output.write( (char*) &thresholds[i], sizeof(unsigned int) );
		}

		for (int i = 0; i < numTotalEntries; i++) {
			for (int j = 0; j < numThresholds; j++) {
				output.write( (char*) &table[i][j], sizeof(unsigned int) );
			}
		}
		output.close();
	}

	void PredictorSolver::loadTable(string tableFileName) {
		ifstream input;
		input.open(tableFileName.c_str(), ios::in | ios::binary);

		unsigned int tempHashLength;
		unsigned int tempNumThresholds;

		input.read( (char*) &tempHashLength, sizeof(unsigned int) );
		input.read( (char*) &tempNumThresholds, sizeof(unsigned int) );

		unsigned int * tempThresholds = new unsigned int [tempNumThresholds];

		for (int i = 0; i < tempNumThresholds; i++) {
			input.read( (char*) &tempThresholds[i], sizeof(unsigned int) );
		}

		init(tempHashLength, tempNumThresholds, tempThresholds);

		delete [] tempThresholds;
		tempThresholds = NULL;

		for (int i = 0; i < numTotalEntries; i++) {
			for (int j = 0; j < numThresholds; j++) {
				input.read( (char*) &table[i][j], sizeof(unsigned int) );
			}
		}

		input.close();
	}

	unsigned int PredictorSolver::getLength(unsigned int hashVal, int thresholdIdx) {
		assert(thresholdIdx < numThresholds);
		assert(hashVal < numTotalEntries);

		return table[hashVal][thresholdIdx];
	}

	unsigned int PredictorSolver::getThreshold(unsigned int thresholdIdx) {
		return thresholds[thresholdIdx];
	}

	int PredictorSolver::solveIdx(string& DNA, unsigned int numOfSeeds) {
		unsigned int beginIdx = 0;
		unsigned int endIdx = numThresholds - 1;
		unsigned int midIdx = (beginIdx + endIdx) / 2;

		while (beginIdx < midIdx) {
			if (idxIsValid(DNA, numOfSeeds, midIdx) ) {
				beginIdx = midIdx;
			}
			else {
				endIdx = midIdx - 1;
			}
			midIdx = (beginIdx + endIdx) / 2;
		}

		if (idxIsValid(DNA, numOfSeeds, midIdx) == false) {
			assert(midIdx == 0);
			return -1;
		}
		else
			return midIdx;
	}

	int PredictorSolver::solveLength(string& DNA, unsigned int numOfSeeds, int idx, unsigned int* lengthVec) {
		if (idx < 0)
			return -1;

		unsigned int seedIter = 0;
		unsigned int requiredLength = 0;
		unsigned int hashVal;

		while (seedIter < numOfSeeds && requiredLength + hashLength <= DNA.length() ) {
			hashVal = hashGen.calculateHash(DNA.substr(requiredLength, hashLength), hashLength);
			lengthVec[seedIter] = table[hashVal][idx];
			requiredLength += lengthVec[seedIter];
			seedIter++;
		}

		return numOfSeeds * thresholds[idx];
	}

	unsigned int PredictorSolver::getNumThresholds() {
		return numThresholds;
	}

	unsigned int PredictorSolver::getHashLength() {
	return hashLength;
}

unsigned int PredictorSolver::getNumTotalEntries() {
	return numTotalEntries;
}

bool PredictorSolver::idxIsValid(string& DNA, unsigned int numOfSeeds, unsigned int idx) {
	unsigned int seedIter = 0;
	unsigned int requiredLength = 0;
	unsigned int hashVal;

	while (seedIter < numOfSeeds && requiredLength + hashLength <= DNA.length() ) {
		hashVal = hashGen.calculateHash(DNA.substr(requiredLength, hashLength), hashLength);
		requiredLength += table[hashVal][idx];
		seedIter++;
	}

	if (seedIter == numOfSeeds && requiredLength <= DNA.length() )
		return true;
	else
		return false;
}

