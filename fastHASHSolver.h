#include "HashTree.h"
#include <vector>
#include <string>

using namespace std;

class FastHASHSolver {
public:
	FastHASHSolver();
	~FastHASHSolver();

	void loadTree(string treeFileName);
	void init(unsigned int readLength, unsigned int seedNum, unsigned int seedLength);
	unsigned int solveDNA(string DNA);

private:
	vector<unsigned int> seedFreq;
	HashTree tree;
	unsigned int readLength;
	unsigned int seedNum;
	unsigned int seedLength;

};

