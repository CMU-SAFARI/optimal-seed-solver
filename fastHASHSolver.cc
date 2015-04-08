#include "fastHASHSolver.h"
#include <cassert>
#include <algorithm>

FastHASHSolver::FastHASHSolver() {
}

FastHASHSolver::~FastHASHSolver() {
}

void FastHASHSolver::loadTree(string treeFileName) {
	tree.loadTree(treeFileName);
}

void FastHASHSolver::init(unsigned int readLength, unsigned int seedNum, unsigned int seedLength) {
	this->seedNum = seedNum;
	this->seedLength = seedLength;
	this->readLength = readLength;

	seedFreq.resize(readLength / seedLength);
}

unsigned int FastHASHSolver::solveDNA(string DNA) {
	assert(DNA.length() == readLength);

	for (int i = 0; i < readLength / seedLength; i++) {
		string seed = DNA.substr(i * seedLength, seedLength);

		tree.query(seed);
		seedFreq[i] = tree.frequency();
	}

	sort(seedFreq.begin(), seedFreq.end() );

#ifdef DEBUG	
	for (int i = 0; i < readLength / seedLength; i++)
		cout << seedFreq[i] << " ";

	cout << endl;
#endif

	unsigned int totalFreq = 0;
	for (int i = 0; i < seedNum; i++)
		totalFreq += seedFreq[i];

	return totalFreq;
}

