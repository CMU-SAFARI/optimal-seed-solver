#include "optimalSolver.h"
#include <cstring>
#include <cassert>
#include <iostream>
#include <climits>
//#define DEBUG

OptimalSolver::OptimalSolver() {
	matrix = NULL;
	defaultCells = NULL;
	readLength = 0;
	seedNum = 0;
	matrixSize = 0;
	L0Loaded = false;
}

OptimalSolver::~OptimalSolver() {
	if  (matrix != NULL) {
		delete [] matrix;
		matrix = NULL;
	}
	if  (defaultCells != NULL) {
		delete [] defaultCells;
		defaultCells = NULL;
	}

	if (level.size() != 0)
		level.clear();
}

void OptimalSolver::loadTree(string treeFileName) {
	tree.loadTree(treeFileName);
}

void OptimalSolver::generateTree(string refFileName) {
	tree.generateTree(refFileName);
}

void OptimalSolver::init(int readLength, int seedNum) {
	//Initialize the matrix
	if (matrix != NULL) {
		delete [] matrix;
		matrix = NULL;
	}

	this->readLength = readLength;

	if (readLength >= seedNum * tree.getHashLength() )
		this->seedNum = seedNum;
	else
		this->seedNum = readLength / tree.getHashLength();

	//Position count for the level
	int lvPosCount;

	//Matrix and base overall size
	matrixSize = 0;
	baseSize = 0;

	//For the base case, we need to consider all possible seeds
	lvPosCount = readLength + 1 - tree.getHashLength();
	baseSize = (lvPosCount + 1) * lvPosCount / 2;

	//For the matrix, we need to consider all cases, while leaving
	//the last minLength out (because we don't need to)
	for (int i = 0; i < this->seedNum - 1; i++) {
		//Notice that for the rest levels, we do not need to consider the last
		//hashLength positions
		lvPosCount = readLength + 1 - tree.getHashLength() * (i + 2);
		matrixSize += lvPosCount;
	}

	//Allocate space for the base
	base = new Cell[base];
	defaultBase = new Cell[base];
	
	//Allocate space for the matrix
	matrix = new Cell[matrixSize];
	defaultMatrix = new Cell[matrixSize];

	

	//Progress counter
	int matrixProgress = 0;
	//sublevel counter
	int levelIdx = 0;
	//Fill the level pointers
	level.resize(this->seedNum);
	//Find the first row
	lvPosCount = readLength + 1 - tree.getHashLength();
	level[0].resize(lvPosCount);
	for (int i = lvPosCount; i > 0; i--) {
		level[0][levelIdx] = matrix + matrixProgress;
		matrixProgress += i;
		levelIdx++;
	}

	//Find the reminding levels
	for (int idx = 1; idx < this->seedNum - 1; idx++) {
		lvPosCount = readLength + 1 - tree.getHashLength() * (idx + 2);
		level[idx].resize(lvPosCount);
		//reset the sublevel counter
		levelIdx = 0;
		for (int i = lvPosCount; i > 0; i--) {
			level[idx][levelIdx] = matrix + matrixProgress;
			matrixProgress += i;
			levelIdx++;
		}
	}

	//Set the last level
	level[this->seedNum - 1].resize(1);
	level[this->seedNum - 1][0] = matrix + matrixProgress;

	//Just to check the matrixSize
	assert (matrixProgress + readLength + 1 - this->seedNum * tree.getHashLength() == matrixSize);
}

void OptimalSolver::reset() {
	memcpy(matrix, defaultCells, matrixSize);
}

void OptimalSolver::feedL0() {
	reset();
	tree.setHashLength(10);
	int lvPosCount = readLength + 1 - tree.getHashLength();

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

unsigned int OptimalSolver::solveDNA(string DNA) {
	int lvPosCount = readLength + 1 - tree.getHashLength();
	bool edge_left;
	bool edge_right;

	if (!L0Loaded) {
		assert(DNA.length() == readLength);

		reset();

		//Calculate the first level

#ifdef DEBUG
		cout << "*L0*: " << endl << "*subl 0*: ";
#endif
		for (int i = 0; i < lvPosCount; i++) {
			string seed = DNA.substr(i, tree.getHashLength() );
			tree.query(seed);
			level0[0][i].start = i;
			level0[0][i].end = i + tree.getHashLength() - 1;
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
				edge_left = level[0][i-1][j].start == j;
				edge_right = level[0][i-1][j+1].end == j + i + tree.getHashLength() - 1;

				// If left edge but not right, copy left seed
				if (edge_left && !edge_right) {
					level[0][i][j].isleaf = level[0][i-1][j].isleaf;
					level[0][i][j].frequency = level[0][i-1][j].frequency;
					level[0][i][j].start = level[0][i-1][j].start;
					level[0][i][j].end = level[0][i-1][j].end;
				}
				// If right edge but not left, copy right seed
				else if (!edge_left && edge_right) {
					level[0][i][j].isleaf = level[0][i-1][j+1].isleaf;
					level[0][i][j].frequency = level[0][i-1][j+1].frequency;
					level[0][i][j].start = level[0][i-1][j+1].start;
					level[0][i][j].end = level[0][i-1][j+1].end;
				}
				// If both are edging, then we have to test potentially three cases
				else if (edge_left && edge_right) {
					// If left is not leaf, always merge
					if (!level[0][i-1][j].isleaf) {
						string seed = DNA.substr(j, tree.getHashLength() + i);

						//TODO: Optimize query
						tree.query(seed);
						level[0][i][j].isleaf = tree.isLeaf();
						level[0][i][j].frequency = tree.frequency();
						level[0][i][j].start = j;
						level[0][i][j].end = j + tree.getHashLength() + i - 1;
					}
					// Check left seed vs right seed if not merging.
					else {
						// Prioritize left position over right position
						if (level[0][i-1][j].frequency <= level[0][i-1][j+1].frequency) {
							level[0][i][j].isleaf = level[0][i-1][j].isleaf;
							level[0][i][j].frequency = level[0][i-1][j].frequency;
							level[0][i][j].start = level[0][i-1][j].start;
							level[0][i][j].end = level[0][i-1][j].end;
						}
						else {
							level[0][i][j].isleaf = level[0][i-1][j+1].isleaf;
							level[0][i][j].frequency = level[0][i-1][j+1].frequency;
							level[0][i][j].start = level[0][i-1][j+1].start;
							level[0][i][j].end = level[0][i-1][j+1].end;
						}
					}

				}
				// None is edging. The left and right seed should be equal.
				else {
					assert(level[0][i-1][j].start == level[0][i-1][j+1].start);
					assert(level[0][i-1][j].end == level[0][i-1][j+1].end);
					assert(level[0][i-1][j].frequency == level[0][i-1][j+1].frequency);
					assert(level[0][i-1][j].isleaf == level[0][i-1][j+1].isleaf);
					level[0][i][j].start = level[0][i-1][j].start;
					level[0][i][j].end = level[0][i-1][j].end;
					level[0][i][j].frequency = level[0][i-1][j].frequency;
					level[0][i][j].isleaf = level[0][i-1][j].isleaf;
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

	//Now calculate the middle levels
	for (int l = 1; l < seedNum - 1; l++) {
#ifdef DEBUG
		cout << "*L" << l << "*: " << endl << "*subl 0*: ";
#endif
		//Execpt the first level, all the rest levels do not have the last hashLength entries
		lvPosCount = readLength + 1 - tree.getHashLength() * (l + 2);
		//Fill up the first layer
		for (int j = 0; j <lvPosCount; j++) {
			level[l][0][j].start = j;
			level[l][0][j].end = j + (l + 1) * tree.getHashLength() - 1;

			// Get left seed. Left seed from previous level
			level[l][0][j].lstart = j;
			level[l][0][j].lend = j + l * tree.getHashLength() - 1;
			level[l][0][j].lfreq = level[l-1][0][j].frequency;

			// Get right seed. Right seed from 0 level
			level[l][0][j].rstart = j + l * tree.getHashLength();
			level[l][0][j].rend = j + (l + 1) * tree.getHashLength() - 1;
			level[l][0][j].rfreq = level[0][0][j+l*tree.getHashLength()].frequency;

			//Get freq from level 0 for the right part
			level[l][0][j].frequency = level[l][0][j].lfreq + level[l][0][j].rfreq;
#ifdef DEBUG
			cout << level[l][0][j].start << "|"
				<< level[l][0][j].lstart << "-" << level[l][0][j].lend << ":" << level[l][0][j].lfreq << "|"
				<< level[l][0][j].rstart << "-" << level[l][0][j].rend << ":" << level[l][0][j].rfreq
				<< "|" << level[l][0][j].end << ":" << level[l][0][j].frequency << "+";
#endif
		}
#ifdef DEBUG
		cout << endl;
#endif

		//Fill up the rest layers, similar to before. Check left_edge and right_edge
		for (int i = 1; i < lvPosCount; i++) {

#ifdef DEBUG
			cout << "*subl" << i << "*: ";
#endif

			for (int j = 0; j < lvPosCount - i; j++) {
				// Test if we touch edges
				edge_left = level[l][i-1][j].start == j;
				edge_right = level[l][i-1][j+1].end == j + i + (l + 1) * tree.getHashLength() - 1;

				// If left edge but not right, copy left seed
				if (edge_left && !edge_right) {
					level[l][i][j].frequency = level[l][i-1][j].frequency;
					level[l][i][j].start = level[l][i-1][j].start;
					level[l][i][j].end = level[l][i-1][j].end;
					level[l][i][j].lstart = level[l][i-1][j].lstart;
					level[l][i][j].lend = level[l][i-1][j].lend;
					level[l][i][j].lfreq = level[l][i-1][j].lfreq;
					level[l][i][j].rstart = level[l][i-1][j].rstart;
					level[l][i][j].rend = level[l][i-1][j].rend;
					level[l][i][j].rfreq = level[l][i-1][j].rfreq;
				}
				// If right edge but not left, copy right seed
				else if (!edge_left && edge_right) {
					level[l][i][j].frequency = level[l][i-1][j+1].frequency;
					level[l][i][j].start = level[l][i-1][j+1].start;
					level[l][i][j].end = level[l][i-1][j+1].end;
					level[l][i][j].lstart = level[l][i-1][j+1].lstart;
					level[l][i][j].lend = level[l][i-1][j+1].lend;
					level[l][i][j].lfreq = level[l][i-1][j+1].lfreq;
					level[l][i][j].rstart = level[l][i-1][j+1].rstart;
					level[l][i][j].rend = level[l][i-1][j+1].rend;
					level[l][i][j].rfreq = level[l][i-1][j+1].rfreq;
				}
				// If both are edging, then we have test potentially three cases
				else if (edge_left && edge_right) {
					// First get the best config from the two seeds.
					// Prioritize left position over right position
					if (level[l][i-1][j].frequency <= level[l][i-1][j+1].frequency) {
						level[l][i][j].frequency = level[l][i-1][j].frequency;
						level[l][i][j].start = level[l][i-1][j].start;
						level[l][i][j].end = level[l][i-1][j].end;
						level[l][i][j].lstart = level[l][i-1][j].lstart;
						level[l][i][j].lend = level[l][i-1][j].lend;
						level[l][i][j].lfreq = level[l][i-1][j].lfreq;
						level[l][i][j].rstart = level[l][i-1][j].rstart;
						level[l][i][j].rend = level[l][i-1][j].rend;
						level[l][i][j].rfreq = level[l][i-1][j].rfreq;
					}
					else {
						level[l][i][j].frequency = level[l][i-1][j+1].frequency;
						level[l][i][j].start = level[l][i-1][j+1].start;
						level[l][i][j].end = level[l][i-1][j+1].end;
						level[l][i][j].lstart = level[l][i-1][j+1].lstart;
						level[l][i][j].lend = level[l][i-1][j+1].lend;
						level[l][i][j].lfreq = level[l][i-1][j+1].lfreq;
						level[l][i][j].rstart = level[l][i-1][j+1].rstart;
						level[l][i][j].rend = level[l][i-1][j+1].rend;
						level[l][i][j].rfreq = level[l][i-1][j+1].rfreq;
					}

					// We set the divider as the beginning for the second seed.
					for (int d = level[l][i-1][j].lend + 1; d <= level[l][i-1][j+1].rstart; d++) {
						//						cout << "&d:" << d << "&";
						//						cout << "left:" << "level[" << l-1 <<"][" << d - j - l * tree.getHashLength() 
						//							<< "][" << j <<"]|" << level[l-1][d - j - l * tree.getHashLength()][j].start << "-"
						//
						//							<< level[l-1][d - j - l * tree.getHashLength()][j].end << ":"
						//							<< level[l-1][d - j - l * tree.getHashLength()][j].frequency << "|";
						//
						//						cout << "right:" << "level[0][" << j + i + tree.getHashLength() - d 
						//							<< "][" << d <<"]|" << level[0][j + i + tree.getHashLength() - d][d].start << "-"
						//
						//							<< level[0][j + i +  tree.getHashLength() - d][d].end << ":"
						//							<< level[0][j + i + tree.getHashLength() - d][d].frequency;
						//
						//						cout << "	";

						// Stop when left stoped touching edge
						if (level[l-1][d - j - l * tree.getHashLength()][j].start != j)
							break;
				
						// Continue until the right reaches edge
						if (level[0][j + i + l * tree.getHashLength() - d][d].end !=
								j + i + (l + 1) * tree.getHashLength() - 1)
							continue;

						// If smaller or equal && prioritize left
						if (level[l-1][d - j - l * tree.getHashLength()][j].frequency +
								level[0][j + i + l * tree.getHashLength() - d][d].frequency
								< level[l][i][j].frequency
								||
								(level[l-1][d - j - l * tree.getHashLength()][j].frequency +
								 level[0][j + i + l * tree.getHashLength() - d][d].frequency
								 == level[l][i][j].frequency && level[l][i][j].start != j)
						   ) {
							level[l][i][j].start = j;
							level[l][i][j].end = j + i + (l + 1) * tree.getHashLength() - 1;
							level[l][i][j].lstart = j;
							level[l][i][j].lend = level[l-1][d - j - l * tree.getHashLength()][j].end;
							level[l][i][j].lfreq = level[l-1][d - j - l * tree.getHashLength()][j].frequency;
							level[l][i][j].rstart = level[0][j + i + l * tree.getHashLength() - d][d].start;
							level[l][i][j].rend = j + i + (l + 1) * tree.getHashLength() - 1;
							level[l][i][j].rfreq = level[0][j + i + l * tree.getHashLength() - d][d].frequency;
							level[l][i][j].frequency = level[l][i][j].lfreq + level[l][i][j].rfreq;
						}
					}
				}
				// None is edging. The left and right seed should be equal.
				else {
					//cout << "i: " << i << " j: " << j << endl;
					assert(level[l][i-1][j].start == level[l][i-1][j+1].start);
					assert(level[l][i-1][j].end == level[l][i-1][j+1].end);
					assert(level[l][i-1][j].frequency == level[l][i-1][j+1].frequency);
					assert(level[l][i-1][j].lstart == level[l][i-1][j+1].lstart);
					assert(level[l][i-1][j].lend == level[l][i-1][j+1].lend);
					assert(level[l][i-1][j].lfreq == level[l][i-1][j+1].lfreq);
					assert(level[l][i-1][j].rstart == level[l][i-1][j+1].rstart);
					assert(level[l][i-1][j].rend == level[l][i-1][j+1].rend);
					assert(level[l][i-1][j].rfreq == level[l][i-1][j+1].rfreq);
					level[l][i][j].start = level[l][i-1][j].start;
					level[l][i][j].end = level[l][i-1][j].end;
					level[l][i][j].frequency = level[l][i-1][j].frequency;
					level[l][i][j].lstart = level[l][i-1][j].lstart;
					level[l][i][j].lend = level[l][i-1][j].lend;
					level[l][i][j].lfreq = level[l][i-1][j].lfreq;
					level[l][i][j].rstart = level[l][i-1][j].rstart;
					level[l][i][j].rend = level[l][i-1][j].rend;
					level[l][i][j].rfreq = level[l][i-1][j].rfreq;
				}
#ifdef DEBUG
				cout << level[l][i][j].start << "|"
					<< level[l][i][j].lstart << "-" << level[l][i][j].lend << ":" << level[l][i][j].lfreq << "|"
					<< level[l][i][j].rstart << "-" << level[l][i][j].rend << ":" << level[l][i][j].rfreq
					<< "|" << level[l][i][j].end << ":" << level[l][i][j].frequency << "+ ";
#endif
			}
#ifdef DEBUG
			cout << endl;
#endif
		}

	}

	// The last level... In fact we don't need to feel it up.
	// The same as previous divider, we choose the divider as the starting position of the right seed
	unsigned int leftDiv = (seedNum - 1) * tree.getHashLength();
	unsigned int rightDiv = readLength - tree.getHashLength();

	unsigned int bestFreq = UINT_MAX;
	unsigned int bestDiv;
	unsigned int frequency;

	if (seedNum > 1) {
		// Check leftDiv first
		frequency = 0;
		frequency += level[seedNum-2][leftDiv-(seedNum-1)*tree.getHashLength()][0].frequency;
		frequency += level[0][readLength-leftDiv-tree.getHashLength()][leftDiv].frequency;

#ifdef DEBUG
		cout << "leftDiv: " << leftDiv << " f: " << frequency << "	";
#endif

		if (frequency < bestFreq) {
			bestFreq = frequency;
			bestDiv = leftDiv;
		}

		// Check rightDiv next
		frequency = 0;
		frequency += level[seedNum-2][rightDiv-(seedNum-1)*tree.getHashLength()][0].frequency;
		frequency += level[0][readLength-rightDiv-tree.getHashLength()][rightDiv].frequency;

#ifdef DEBUG
		cout << "rightDiv: " << rightDiv << " f: " << frequency << "	" << endl;
#endif

		if (frequency < bestFreq) {
			bestFreq = frequency;
			bestDiv = rightDiv;
		}

		while (leftDiv <= rightDiv) {
			// Update the 2 divs
			if (level[0][readLength-leftDiv-tree.getHashLength()][leftDiv].start != leftDiv)
				leftDiv = level[0][readLength-leftDiv-tree.getHashLength()][leftDiv].start;
			else
				leftDiv++;

			if (level[seedNum-2][rightDiv-(seedNum-1)*tree.getHashLength()][0].end != rightDiv - 1)
				rightDiv = level[seedNum-2][rightDiv-(seedNum-1)*tree.getHashLength()][0].end + 1;
			else
				rightDiv--;

			// Check leftDiv first
			frequency = 0;
			frequency += level[seedNum-2][leftDiv-(seedNum-1)*tree.getHashLength()][0].frequency;
			frequency += level[0][readLength-leftDiv-tree.getHashLength()][leftDiv].frequency;

#ifdef DEBUG
			cout << "leftDiv: " << leftDiv << " f: " << frequency << "	";
#endif

			if (frequency < bestFreq) {
				bestFreq = frequency;
				bestDiv = leftDiv;
			}

			// Check rightDiv next
			frequency = 0;
			frequency += level[seedNum-2][rightDiv-(seedNum-1)*tree.getHashLength()][0].frequency;
			frequency += level[0][readLength-rightDiv-tree.getHashLength()][rightDiv].frequency;

#ifdef DEBUG
			cout << "rightDiv: " << rightDiv << " f: " << frequency << "	" << endl;
#endif

			if (frequency < bestFreq) {
				bestFreq = frequency;
				bestDiv = rightDiv;
			}

		}
	}
	else {
		bestDiv = readLength;
		bestFreq = level[0][readLength - tree.getHashLength()][0].frequency;
	}

#ifdef DEBUG
	cout << "bestDiv: " << bestDiv << endl;
	cout << "bestFreq: " << bestFreq << endl;
#endif
	//	
	//	// Print seeds from right to left
	//	string seed;
	//	int levelLength = readLength;
	//	int levelDiv = bestDiv;
	//	int seedPos;
	//	int seedLength;
	//
	//	for (int l = seedNum - 2; l >= 0; l--) {
	//		seedPos = level[0][levelLength - levelDiv - tree.getHashLength()][levelDiv].start;
	//		seedLength = level[0][levelLength - levelDiv - tree.getHashLength()][levelDiv].end - seedPos + 1;
	//		seed = DNA.substr(seedPos, seedLength);
	//
	//		cout << "seed: " << seed << " at " << seedPos << " with frequency of "
	//			<< level[0][levelLength - levelDiv - tree.getHashLength()][levelDiv].frequency << endl;
	//		
	//		cout << "Checking seed[0][" << levelLength - levelDiv - tree.getHashLength() << "][" << levelDiv
	//			<< "]|" << level[0][levelLength - levelDiv - tree.getHashLength()][levelDiv].start
	//			<< "-" << level[0][levelLength - levelDiv - tree.getHashLength()][levelDiv].end
	//			<< ":" << level[0][levelLength - levelDiv - tree.getHashLength()][levelDiv].frequency << endl;
	//
	//		if (l != 0) {
	//			levelLength = level[l][levelDiv - (l+1) * tree.getHashLength()][0].end + 1;	
	//			levelDiv = level[l][levelDiv - (l+1) * tree.getHashLength()][0].lend + 1;
	//		}
	//	}
	//	// Print the last, the left most seed
	//	seedPos = level[0][levelDiv - tree.getHashLength()][0].start;
	//	seedLength = level[0][levelDiv - tree.getHashLength()][0].end + 1
	//		- level[0][levelDiv - tree.getHashLength()][0].start;
	//	seed = DNA.substr(seedPos, seedLength);
	//	cout << "seed: " << seed << " at " << 0 << " with frequency of "
	//		<< level[0][levelDiv - tree.getHashLength()][0].frequency << endl;
	//		
	//	cout << "Checking seed[0][" << levelDiv - tree.getHashLength() 
	//		<< "[0]|" << level[0][levelDiv - tree.getHashLength()][0].start
	//		<< "-" << level[0][levelDiv - tree.getHashLength()][0].end
	//		<< ":" << level[0][levelDiv - tree.getHashLength()][0].frequency << endl;

	L0Loaded = false;

	return bestFreq;
}

