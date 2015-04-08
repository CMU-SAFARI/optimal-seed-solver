#include "HashTree.h"
#include <vector>
#include <string>

using namespace std;

class BasicSolver {
public:
	BasicSolver();
	~BasicSolver();

	void loadTree(string treeFileName);
	void init(unsigned int readLength, unsigned int seedNum, unsigned int seedLength);
	unsigned int solveDNA(string DNA);

private:
	HashTree tree;
	unsigned int readLength;
	unsigned int seedNum;
	unsigned int seedLength;

};

