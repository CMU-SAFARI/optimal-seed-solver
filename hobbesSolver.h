#include "HashTree.h"
#include <vector>

using namespace std;


class HobbesSolver {
public:
	HobbesSolver();
	~HobbesSolver();
	void loadTree(string treeFileName);
	void generateTree(string treeFileName);
	void init(int readLength, int seedNum, int seedLength);
	void reset();
	unsigned int solveDNA(string DNA);
	void print();

private:
	HashTree tree;
	int *invertedList;
	int **dynamicMatrix;
	int readLength;
	int seedNum;
	int seedLength;
	int *defaultInvertedList;
	int *defaultDynamicMatrix;
	int dynamicMatrixWidth;
};
