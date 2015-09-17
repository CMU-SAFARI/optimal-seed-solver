#include "optimalSolverLN.h"
#include <cstring>
#include <cassert>
#include <iostream>
#include <climits>
//#define DEBUG

OptimalSolverLN::OptimalSolverLN() {
	matrix = NULL;
	base = NULL;
	defaultMatrix = NULL;
	defaultBase = NULL;
	readLength = 0;
	seedNum = 0;
	matrixSize = 0;
	baseSize = 0;
	L0Loaded = false;
	minLength = 1;
}

OptimalSolverLN::~OptimalSolverLN() {
	if  (matrix != NULL) {
		delete [] matrix;
		matrix = NULL;
	}
	if  (base != NULL) {
		delete [] base;
		base = NULL;
	}
	if  (defaultMatrix != NULL) {
		delete [] defaultMatrix;
		defaultMatrix = NULL;
	}
	if  (defaultBase != NULL) {
		delete [] defaultBase;
		defaultBase = NULL;
	}

	if (level.size() != 0)
		level.clear();
	if (level0.size() != 0)
		level0.clear();
	if (seeds.size() != 0)
		seeds.clear();
}

void OptimalSolverLN::loadTree(string treeFileName) {
	tree.loadTree(treeFileName);
	minLength = tree.getHashLength();
}

void OptimalSolverLN::generateTree(string refFileName) {
	tree.generateTree(refFileName);
}

void OptimalSolverLN::init(int readLength, int seedNum) {
	//Initialize the matrix
	if  (matrix != NULL) {
		delete [] matrix;
		matrix = NULL;
	}
	if  (base != NULL) {
		delete [] base;
		base = NULL;
	}
	if  (defaultMatrix != NULL) {
		delete [] defaultMatrix;
		defaultMatrix = NULL;
	}
	if  (defaultBase != NULL) {
		delete [] defaultBase;
		defaultBase = NULL;
	}

	this->readLength = readLength;

	if (readLength >= seedNum * minLength)
		this->seedNum = seedNum;
	else
		this->seedNum = readLength / minLength;

	//Position count for the level
	int lvPosCount;

	//Matrix and base overall size
	matrixSize = 0;
	baseSize = 0;

	//For the base case, we need to consider all possible seeds
	lvPosCount = readLength + 1 - minLength;
	baseSize = (lvPosCount + 1) * lvPosCount / 2;

	//For the matrix, we need to consider all cases, while leaving
	//the last minLength out (because we don't need to)
	if (seedNum > 1)
		matrixSize = (readLength + 1 - minLength * seedNum) * (seedNum - 1);
	else
		matrixSize = readLength + 1 - minLength;
	/*
	for (int i = 0; i < this->seedNum - 1; i++) {
		//Notice that for the rest levels, we do not need to consider the last
		//hashLength positions
		lvPosCount = readLength + 1 - minLength * (i + 2);
		matrixSize += lvPosCount;
	}
	*/

#ifdef DEBUG
	cout << "matrixSize: " << matrixSize << endl;
#endif
	//Allocate space for the base
	base = new Cell[baseSize];
	defaultBase = new Cell[baseSize];

	//Allocate space for the matrix
	matrix = new Cell[matrixSize];
	defaultMatrix = new Cell[matrixSize];

	//Fill the level0 pointers
	lvPosCount = readLength + 1 - minLength;
	level0.resize(lvPosCount);
	int levelIdx = 0;
	int baseProgress = 0;
	//Initialize the first one
	for (int i = lvPosCount; i > 0; i--) {
		level0[levelIdx] = base + baseProgress;
		baseProgress += i;
		levelIdx++;
	}

	//Just to check the baseSize
	assert (baseProgress == baseSize);

	if (seedNum > 1) {
		//Fill the level pointers. No need to have the last level
		level.resize(seedNum - 1);
		levelIdx = 0;
		int matrixProgress = 0;
		for (int i = 0; i < seedNum - 1; i++) {
			lvPosCount = readLength + 1 - minLength * seedNum;
			level[i] = matrix + matrixProgress;
			matrixProgress += lvPosCount;
		}

		//Just to check the matrixSize
		assert (matrixProgress == matrixSize);
	}
	else {
		level.resize(1);
		level[0] = matrix;
	}


	//Initialize all the seeds
	seeds.resize(seedNum);

	//Initialize statistics
	processedReads = 0;
	divTravel.resize(readLength, 0);
}

void OptimalSolverLN::reset() {
	memcpy(matrix, defaultMatrix, matrixSize);
	memcpy(base, defaultBase, baseSize);
}

void OptimalSolverLN::feedL0() {
	reset();
	tree.setHashLength(10);
	int lvPosCount = readLength + 1 - minLength;

	for (int i = 0; i < lvPosCount; i++) {
		for (int j = 0; j < lvPosCount - i; j++) {
			cin >> level0[i][j].start;
			cin >> level0[i][j].end;
			cin >> level0[i][j].frequency;
			cin >> level0[i][j].isleaf;
#ifdef DEBUG
			cout << "level0[" << i << "][" << j << "]:" << level0[i][j].start << "-" 
				<<level0[i][j].end << ":" << level0[i][j].frequency << " ";
			cout.flush();
		}
		cout << endl;
#else
		}
#endif

	}
	L0Loaded = true;
}

void OptimalSolverLN::setMinLength(int minLength) {
	this->minLength = minLength;
}

void OptimalSolverLN::loadL0(string& DNA) {
	bool edge_left;
	bool edge_right;
	int lvPosCount = readLength + 1 - minLength;

	if (!L0Loaded) {
		assert(DNA.length() == static_cast<unsigned int>(readLength) );

		reset();

		//Calculate the first level

#ifdef DEBUG
		cout << "*L0*: " << endl << "*subl 0*: ";
#endif
		for (int i = 0; i < lvPosCount; i++) {
			string seed = DNA.substr(i, minLength);
			tree.query(seed);
			level0[0][i].start = i;
			level0[0][i].end = i + minLength - 1;
			level0[0][i].frequency = tree.frequency();
			level0[0][i].isleaf = tree.isLeaf();

#ifdef DEBUG
			cout << level0[0][i].start << "-" << level0[0][i].end << ":" << level0[0][i].frequency << " ";
#endif
		}
#ifdef DEBUG
		cout << endl;
#endif

		for (int i = 1; i < lvPosCount; i++) {
#ifdef DEBUG
			cout << "*subl" << i << "*: ";
#endif
			for (int j = 0; j < lvPosCount - i; j++) {
				edge_left = level0[i-1][j].start == static_cast<unsigned int>(j);
				edge_right = level0[i-1][j+1].end == static_cast<unsigned int>(j + i + minLength - 1);

				// If left edge but not right, copy left seed
				if (edge_left && !edge_right) {
					level0[i][j].isleaf = level0[i-1][j].isleaf;
					level0[i][j].frequency = level0[i-1][j].frequency;
					level0[i][j].start = level0[i-1][j].start;
					level0[i][j].end = level0[i-1][j].end;
				}
				// If right edge but not left, copy right seed
				else if (!edge_left && edge_right) {
					level0[i][j].isleaf = level0[i-1][j+1].isleaf;
					level0[i][j].frequency = level0[i-1][j+1].frequency;
					level0[i][j].start = level0[i-1][j+1].start;
					level0[i][j].end = level0[i-1][j+1].end;
				}
				// If both are edging, then we have to test potentially three cases
				else if (edge_left && edge_right) {
					// If left is not leaf, always merge
					if (!level0[i-1][j].isleaf) {
						string seed = DNA.substr(j, minLength + i);

						//TODO: Optimize query
						tree.query(seed);
						level0[i][j].isleaf = tree.isLeaf();
						level0[i][j].frequency = tree.frequency();
						level0[i][j].start = j;
						level0[i][j].end = j + minLength + i - 1;
					}
					// Check left seed vs right seed if not merging.
					else {
						// Prioritize left position over right position
						if (level0[i-1][j].frequency <= level0[i-1][j+1].frequency) {
							level0[i][j].isleaf = level0[i-1][j].isleaf;
							level0[i][j].frequency = level0[i-1][j].frequency;
							level0[i][j].start = level0[i-1][j].start;
							level0[i][j].end = level0[i-1][j].end;
						}
						else {
							level0[i][j].isleaf = level0[i-1][j+1].isleaf;
							level0[i][j].frequency = level0[i-1][j+1].frequency;
							level0[i][j].start = level0[i-1][j+1].start;
							level0[i][j].end = level0[i-1][j+1].end;
						}
					}

				}
				// None is edging. The left and right seed should be equal.
				else {
					assert(level0[i-1][j].start == level0[i-1][j+1].start);
					assert(level0[i-1][j].end == level0[i-1][j+1].end);
					assert(level0[i-1][j].frequency == level0[i-1][j+1].frequency);
					assert(level0[i-1][j].isleaf == level0[i-1][j+1].isleaf);
					level0[i][j].start = level0[i-1][j].start;
					level0[i][j].end = level0[i-1][j].end;
					level0[i][j].frequency = level0[i-1][j].frequency;
					level0[i][j].isleaf = level0[i-1][j].isleaf;
				}

#ifdef DEBUG
				cout << level0[i][j].start << "-" << level0[i][j].end << ":" << level0[i][j].frequency << "|" << level0[i][j].isleaf << " ";
#endif
			}
#ifdef DEBUG
			cout << endl;
#endif
		}
	}
}


void OptimalSolverLN::fillMatrix(string& DNA) {
	//Load all substrings
	loadL0(DNA);

#ifdef DEBUG
	cout << "l-" << 0 << "  ";
#endif

	//Fill level 0
	int lvPosCount = readLength - minLength * seedNum;
	for (int pos = 0; pos <= lvPosCount; pos++) {
		level[0][pos].start = level0[pos][0].start;
		level[0][pos].end = level0[pos][0].end;
		level[0][pos].frequency = level0[pos][0].frequency;
		//Fill up again for tracing back
		level[0][pos].rstart = level0[pos][0].start;
		level[0][pos].rend = level0[pos][0].end;
		level[0][pos].rfreq = level0[pos][0].frequency;
#ifdef DEBUG
		cout << 0 << "-" << pos + minLength - 1 << ":" << level[0][pos].start << "-"
			<< level[0][pos].end << ":" << level[0][pos].frequency << " ";
#endif
	}
#ifdef DEBUG
	cout << endl;
#endif

	int opt_div;

	//Now calculate the rest levels
	for (int l = 1; l < seedNum - 1; l++) {

		lvPosCount = readLength - minLength * seedNum;
		opt_div = solveFirstOptimal(lvPosCount, lvPosCount, l);
#ifdef DEBUG
		cout << "lvPosCount: " << lvPosCount << " opt_div: " << opt_div << endl;
#endif
		level[l][lvPosCount].lstart = level[l-1][opt_div].start;
		level[l][lvPosCount].lend = level[l-1][opt_div].end;
		level[l][lvPosCount].lfreq = level[l-1][opt_div].frequency;
		level[l][lvPosCount].rstart = level0[lvPosCount - opt_div][opt_div + l * minLength].start;
		level[l][lvPosCount].rend = level0[lvPosCount - opt_div][opt_div + l * minLength].end;
		level[l][lvPosCount].rfreq = level0[lvPosCount - opt_div][opt_div + l * minLength].frequency;
		level[l][lvPosCount].start = level[l][lvPosCount].lstart;
		level[l][lvPosCount].end = level[l][lvPosCount].rend;
		level[l][lvPosCount].frequency = level[l][lvPosCount].lfreq + level[l][lvPosCount].rfreq;

		for (int pos = lvPosCount - 1; pos >= 0; pos--) {

			int prev_opt_div = level[l][pos+1].rstart - l * minLength;
#ifdef DEBUG
			cout << "prev_opt_div: " << prev_opt_div << endl;
#endif
			//Testing knob
			if (true && prev_opt_div <= pos && level0[pos - prev_opt_div][prev_opt_div + l * minLength].frequency == level[l][pos+1].rfreq) {
#ifdef TEST
				cout << "same as previous; pos: " << pos << endl;
				cout.flush();
#endif
				level[l][pos].lstart = level[l][pos+1].lstart;
				level[l][pos].lend = level[l][pos+1].lend;
				level[l][pos].lfreq = level[l][pos+1].lfreq;
				level[l][pos].rstart = level0[pos - prev_opt_div][opt_div + l * minLength].start;
				level[l][pos].rend = level0[pos - prev_opt_div][prev_opt_div + l * minLength].end;
				level[l][pos].rfreq = level0[pos - prev_opt_div][prev_opt_div + l * minLength].frequency;
				level[l][pos].start = level[l][pos].lstart;
				level[l][pos].end = level[l][pos].rend;
				level[l][pos].frequency = level[l][pos].lfreq + level[l][pos].rfreq;
				divTravel[0]++;
			}
			else {
				if (opt_div > pos)
					opt_div = pos;

				opt_div = solveFirstOptimal(opt_div, pos, l);

				level[l][pos].lstart = level[l-1][opt_div].start;
				level[l][pos].lend = level[l-1][opt_div].end;
				level[l][pos].lfreq = level[l-1][opt_div].frequency;
				level[l][pos].rstart = level0[pos - opt_div][opt_div + l * minLength].start;
				level[l][pos].rend = level0[pos - opt_div][opt_div + l * minLength].end;
				level[l][pos].rfreq = level0[pos - opt_div][opt_div + l * minLength].frequency;
				level[l][pos].start = level[l][pos].lstart;
				level[l][pos].end = level[l][pos].rend;
				level[l][pos].frequency = level[l][pos].lfreq + level[l][pos].rfreq;
			}
		}

#ifdef TEST
		cout << "freq: ";
		for (int pos = 0; pos <= lvPosCount; pos++) {
			cout << level[l][pos].frequency << " ";
		}
		cout << endl;
		cout << "rstart: ";
		for (int pos = 0; pos <= lvPosCount; pos++) {
			cout << level[l][pos].rstart << " ";
		}
		cout << endl;
		cout << "rfreq: ";
		for (int pos = 0; pos <= lvPosCount; pos++) {
			cout << level[l][pos].rfreq << " ";
		}
		cout << endl;
#endif
#ifdef DEBUG
		cout << "l-" << l << "  ";
		for (int pos = 0; pos <= lvPosCount; pos++) {
			cout << 0 << "-" << pos + (l + 1) * minLength - 1 << ":" << level[l][pos].lstart << "-"
				<< level[l][pos].lend << "(" << level[l][pos].lfreq << ")"
				<< level[l][pos].rstart << "-" << level[l][pos].rend << "("
				<< level[l][pos].rfreq << "):" << level[l][pos].frequency << " ";

		}
		cout << endl;
#endif
	}

	L0Loaded = false;
}

unsigned int OptimalSolverLN::calculateLastDiv() {
	// The last level

	if (seedNum > 1) {
		int lvPosCount = readLength - seedNum * minLength;
		finalDiv = solveFirstOptimal(lvPosCount, lvPosCount, seedNum - 1);
	}
	else {
		finalDiv = readLength;
	}

#ifdef DEBUG
	cout << "opt_div: " << finalDiv << endl;
#endif

	return finalDiv;
}
	   
unsigned int OptimalSolverLN::calcualteFreq() {
	int bestFreq;

	if (seedNum > 1) {
		int lvPosCount = readLength - seedNum * minLength;
		bestFreq = level[seedNum - 2][finalDiv].frequency + level0[lvPosCount - finalDiv][finalDiv + (seedNum - 1) * minLength].frequency;
	} else
		bestFreq = level0[readLength - minLength][0].frequency;
	
#ifdef DEBUG
	cout << "bestFreq: " << bestFreq << endl;
#endif

	return bestFreq;
}

unsigned int OptimalSolverLN::solveDNA(string& DNA) {
	fillMatrix(DNA);
	calculateLastDiv();
	processedReads++;
	return calcualteFreq();
}

int OptimalSolverLN::solveFirstOptimal(int opt_div, int pos, int l) {
	bool print = false;
	if (pos == 27 || pos == 55 || pos == 70)
		print = true;

	int lend = level[l-1][opt_div].end;
	int lfreq = level[l-1][opt_div].frequency;
	int rfreq = level0[pos - opt_div][opt_div + l * minLength].frequency;
	int minFreq = lfreq + rfreq;
	int prev_lfreq = lfreq;
	int prev_rfreq = rfreq;
	int next_lfreq = lfreq;

#ifdef TEST
	if (print) {
		cout << "div: " << opt_div << " lfreq: " << lfreq << " rfreq: " << rfreq
			<< " prev_lfreq: " << prev_lfreq << " prev_rfreq: " << prev_rfreq
			<< " total: " << lfreq + rfreq << endl;
	}
#endif
	int div;
	int totalTravel = 0;
	for (div = opt_div - 1; div >= 0; div--) {
		lfreq = level[l-1][div].frequency;
		
		//divider sprinting left
		if (lfreq == prev_lfreq && div > 0 && lfreq == level[l-1][div-1].frequency) {
			continue;
		}

		rfreq = level0[pos - div][div + l * minLength].frequency;

#ifdef TEST
		if (print) {
			cout << "div: " << div << " lfreq: " << lfreq << " rfreq: " << rfreq
				<< " prev_lfreq: " << prev_lfreq << " prev_rfreq: " << prev_rfreq
				<< " total: " << lfreq + rfreq << endl;
		}
#endif

		if (lfreq + rfreq <= minFreq) {
			minFreq = lfreq + rfreq;
			opt_div = div;
		}

		//early divider termination
		if (lfreq - prev_lfreq > prev_rfreq) {
			div--;
			break;
		}

#ifdef DEBUG
		cout << " div_after: " << div << endl;
#endif

		prev_lfreq = lfreq;
		prev_rfreq = rfreq;
		totalTravel++;
	}
	
	divTravel[totalTravel]++;

#ifdef TEST
	if (print)
		cout << " opt_div: " << opt_div << endl;
#endif

	return opt_div;
}

void OptimalSolverLN::backtrack() {
	if (seedNum > 1) {
		int seedIdx = seedNum - 1;
		int level0lv = readLength - seedNum * minLength - finalDiv;
		int level0pos = finalDiv + (seedNum - 1) * minLength;
		seeds[seedIdx].start = level0[level0lv][level0pos].start;
		seeds[seedIdx].end = level0[level0lv][level0pos].end;
		seeds[seedIdx].frequency = level0[level0lv][level0pos].frequency;
	
		int opt_div = finalDiv;
	
		for (seedIdx = seedNum - 2; seedIdx >= 0; seedIdx--) {
			seeds[seedIdx].start = level[seedIdx][opt_div].rstart;
			seeds[seedIdx].end = level[seedIdx][opt_div].rend;
			seeds[seedIdx].frequency = level[seedIdx][opt_div].rfreq;
			opt_div = level[seedIdx][opt_div].lend + 1 - seedIdx * minLength;
		}
	}
	else {
		seeds[0].start = level0[readLength - minLength][0].start;
		seeds[0].end = level0[readLength - minLength][0].end;
		seeds[0].frequency = level0[readLength - minLength][0].frequency;
	}
}

bool OptimalSolverLN::compFreq(Cell left, Cell right) {
	return (left.frequency < right.frequency);
}

bool OptimalSolverLN::compLength(Cell left, Cell right) {
	return ( (left.end + 1 - left.start) < (right.end + 1 - right.start) );
}

void OptimalSolverLN::sortOfFreq() {
	sortSeeds(this->compFreq);
}

void OptimalSolverLN::sortOfLength() {
	sortSeeds(this->compLength);
}

template<class T>
void OptimalSolverLN::sortSeeds(T relation) {
	sort(seeds.begin(), seeds.end(), relation);
}

void OptimalSolverLN::printSeeds(ostream& stream) {
	for (int i = 0; i < seedNum; i++) {
		stream << "Seed[" << i << "]: start: " << seeds[i].start;
		stream << " end: " << seeds[i].end;
		stream << " frequency: " << seeds[i].frequency;
		stream << endl;
	}
}

void OptimalSolverLN::printFreqs(ostream& stream) {
	for (int i = 0; i < seedNum; i++)
		stream << seeds[i].frequency << " ";
	stream << endl;
}

void OptimalSolverLN::printLength(ostream& stream) {
	for (int i = 0; i < seedNum; i++)
		stream << seeds[i].end + 1 - seeds[i].start << " ";
	stream << endl;
}

void OptimalSolverLN::printStats(ostream& stream) {
	unsigned long long sumOfTravel = 0;
	stream << "divTravel distribution: ";
	for (int i = 0; i < readLength; i++) {
		stream << divTravel[i] << " ";
		sumOfTravel += divTravel[i] * i;
	}
	stream << endl;

	unsigned long long prefixesPerRead = (readLength + 1 - seedNum * minLength) * (seedNum - 2) + 1;
	unsigned long long totalPrefixes = processedReads * prefixesPerRead;
	stream << "divTravel per prefix: " << (double) sumOfTravel / totalPrefixes << endl;
	stream << endl;
}
