#include "basicSolver.h"
#include <cassert>
#include <algorithm>

BasicSolver::BasicSolver() {
}

BasicSolver::~BasicSolver() {
}

void BasicSolver::loadTree(string treeFileName) {
	tree.loadTree(treeFileName);
}

void BasicSolver::init(unsigned int readLength, unsigned int seedNum, unsigned int seedLength) {
	this->seedNum = seedNum;
	this->seedLength = seedLength;
	this->readLength = readLength;
}

unsigned int BasicSolver::solveDNA(string DNA) {
	assert(DNA.length() == readLength);

	unsigned int interval = readLength / seedNum;

	unsigned int totalFreq = 0;

	for (int i = 0; i < seedNum; i++) {
		string seed = DNA.substr(i * interval, seedLength);

		tree.query(seed);
		totalFreq += tree.frequency();
	}

	return totalFreq;
}

