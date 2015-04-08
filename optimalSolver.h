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

class OptimalSolver {
public:
	OptimalSolver();
	~OptimalSolver();
	void loadTree(string treeFileName);
	void generateTree(string treeFileName);
	void init(int readLength, int seedNum);
	void reset();
	unsigned solveDNA(string DNA);
	void print();
	void feedL0();

private:
	HashTree tree;

	bool L0Loaded;
	int readLength;
	int seedNum;
	Cell* matrix;
	vector<vector<Cell*> > level;
	Cell* defaultCells;
	int matrixSize;
};
