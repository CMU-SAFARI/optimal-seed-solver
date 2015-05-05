#include "HashTree.h"
#include <vector>

using namespace std;

struct Cell {
	Cell() {
		start = -1;
		end = -1;
		isleaf = false;
		frequency = -1;
		lstart = -1;
		lend = -1;
		lfreq = -1;
		rstart = -1;
		rend = -1;
		rfreq = -1;
	};

	unsigned int start;
	unsigned int end;
	bool isleaf;
	unsigned int frequency;

	unsigned int lstart;
	unsigned int lend;
	unsigned int lfreq;
	unsigned int rstart;
	unsigned int rend;
	unsigned int rfreq;
};

class OptimalSolverLN {
public:
	OptimalSolverLN();
	~OptimalSolverLN();
	void loadTree(string treeFileName);
	void generateTree(string treeFileName);
	void init(int readLength, int seedNum);
	void reset();
	unsigned solveDNA(string DNA);
	void print();

private:
	//Internal functions
	void loadL0(string DNA);

	//This is for debugging
	void feedL0();
	bool L0Loaded;
	
	HashTree tree;

	int readLength;
	int seedNum;

	Cell* matrix;
	Cell* base;
	vector<Cell*> level;
	vector<Cell*> level0;
	Cell* defaultMatrix;
	Cell* defaultBase;
	int matrixSize;
	int baseSize;
};
