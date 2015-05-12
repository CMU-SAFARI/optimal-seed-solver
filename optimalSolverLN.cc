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
}

void OptimalSolverLN::loadTree(string treeFileName) {
	tree.loadTree(treeFileName);
	minLength = tree.getHashLength();
	maxLength = tree.getFullLength();
}

void OptimalSolverLN::generateTree(string refFileName) {
	tree.generateTree(refFileName);
}

void OptimalSolverLN::init(int readLength, int seedNum) {
	//Initialize the matrix
	if (matrix != NULL) {
		delete [] matrix;
		matrix = NULL;
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
	for (int i = 0; i < this->seedNum - 1; i++) {
		//Notice that for the rest levels, we do not need to consider the last
		//hashLength positions
		lvPosCount = readLength + 1 - minLength * (i + 2);
		matrixSize += lvPosCount;
	}

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

	//Fill the level pointers. No need to have the last level
	level.resize(seedNum - 1);
	levelIdx = 0;
	int matrixProgress = 0;
	for (int i = 0; i < seedNum - 1; i++) {
		lvPosCount = readLength + 1 - minLength * (i + 2);
		level[i] = matrix + matrixProgress + lvPosCount;
		matrixProgress += lvPosCount;
	}

	//Just to check the matrixSize
	assert (matrixProgress == matrixSize);
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
			cout << "level0[" << i << "][" << j << "].start=" << level0[i][j].start << endl;
			cout.flush();
#endif
		}
	}
			L0Loaded = true;
}

void OptimalSolverLN::loadL0(string& DNA) {
	bool edge_left;
	bool edge_right;
	int lvPosCount = readLength + 1 - minLength;

	if (!L0Loaded) {
		assert(DNA.length() == readLength);

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
				edge_left = level0[i-1][j].start == j;
				edge_right = level0[i-1][j+1].end == j + i + minLength - 1;

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
				cout << level[0][i][j].start << "-" << level[0][i][j].end << ":" << level[0][i][j].frequency << "|" << level[0][i][j].isleaf << " ";
#endif
			}
#ifdef DEBUG
			cout << endl;
#endif
		}
	}
}


unsigned int OptimalSolverLN::solveDNA(string DNA) {
	//Load all substrings
	loadL0(DNA);

	//Fill level 0
	int lvPosCount = readLength + 1 - minLength * 2;
	for (int pos = 0; pos < lvPosCount; pos++) {
		level[0][pos].start = level0[pos][0].start;
		level[0][pos].end = level0[pos][0].end;
		level[0][pos].frequency = level0[pos][0].frequency;
	}

	int opt_div;

	//Now calculate the rest levels
	for (int l = 1; l < seedNum - 1; l++) {

		lvPosCount = readLength - minLength * (l + 2);
		opt_div = solve_first_optimal(lvPosCount, lvPosCount, l);
		level[l][lvPosCount].lstart = level[l-1][opt_div].start;
		level[l][lvPosCount].lend = level[l-1][opt_div].end;
		level[l][lvPosCount].lfreq = level[l-1][opt_div].frequency;
		level[l][lvPosCount].rstart = level0[lvPosCount - opt_div][opt_div + l * minLength].start;
		level[l][lvPosCount].rend = level0[lvPosCount - opt_div][opt_div + l * minLength].end;
		level[l][lvPosCount].rfreq = level0[lvPosCount - opt_div][opt_div + l * minLength].frequency;
		level[l][lvPosCount].frequency = level[l][lvPosCount].lfreq + level[l][lvPosCount].rfreq;
	
		for (int pos = lvPosCount - 1; pos >= 0; pos--) {

			if (opt_div > pos)
				opt_div = pos;
		
			opt_div = solve_first_optimal(opt_div, pos, l);
		
			level[l][pos].lstart = level[l-1][opt_div].start;
			level[l][pos].lend = level[l-1][opt_div].end;
			level[l][pos].lfreq = level[l-1][opt_div].frequency;
			level[l][pos].rstart = level0[pos - opt_div][opt_div + l * minLength].start;
			level[l][pos].rend = level0[pos - opt_div][opt_div + l * minLength].end;
			level[l][pos].rfreq = level0[pos - opt_div][opt_div + l * minLength].frequency;
			level[l][pos].frequency = level[l][pos].lfreq + level[l][pos].rfreq;
		}
#ifdef DEBUG
		for (int pos = lvPosCount; pos >= 0; pos--) {
			cout << pos << "-" << pos + l * minLength - 1 << ":" << level[l][pos].lstart << "-"
					<< level[l][pos].lend << "(" << level[l][pos].lfreq << ")"
					<< level[l][pos].rstart << "-" << level[l][pos].rend << "("
					<< level[l][pos].rfreq << "):" << level[l][pos].frequency << " ";
		
		}
		cout << endl;
#endif
	}

	// The last level... In fact we don't need to feel it up.
	// The same as previous divider, we choose the divider as the starting position of the right seed
	unsigned int bestFreq = UINT_MAX;

	if (seedNum > 1) {
		lvPosCount = readLength + 1 - seedNum * minLength;
		opt_div = solve_first_optimal(lvPosCount + 1, lvPosCount, seedNum - 1);
		bestFreq = level[seedNum - 2][opt_div - (seedNum - 2) * minLength - 1].frequency + level0[lvPosCount - opt_div][opt_div].frequency;
	}
	else {
		opt_div = readLength;
		bestFreq = level0[readLength - minLength][0].frequency;
	}

#ifdef DEBUG
	cout << "bestDiv: " << bestDiv << endl;
	cout << "bestFreq: " << bestFreq << endl;
#endif
	L0Loaded = false;

	return bestFreq;
}

int OptimalSolverLN::solve_first_optimal(int opt_div, int pos, int l) {
	int lend = level[l-1][opt_div].end;
	int lfreq = level[l-1][opt_div].frequency;
	int rfreq = level0[pos - opt_div][opt_div + l * minLength].frequency;
	int minFreq = lfreq + rfreq;
	int prev_lfreq = lfreq;
	int pref_rfreq = rfreq;

	for (int div = pos - 1; div >= 0; div--) {
		lfreq = level[l-1][opt_div].frequency;
		rfreq = level0[pos - opt_div][opt_div + l * minLength].frequency;

		if (lfreq + rfreq <= minFreq)
			opt_div = div;

		//early divider termination
		if (lfreq - prev_lfreq < pref_rfreq)
			break;

		//divider sprinting left
		if (div + l * minLength > lend + 1) {
			div = lend + 1 - l * minLength;
		}
	}

	return opt_div;
}

