#include "HashTree.h"
#include <vector>
#include <string>

using namespace std;

class templateClass{
public:
	templateClass();
	void generateTree(string treeFileName);
	void loadTree(string treeFileName);
	void setThresholdVec(const vector<int>& thresholdVec);
	void analyzeRead(string Read);
	void printInfo();

private:
	string read;
	HashTree tree;
	vector<int> thresholdVec;
	vector<int> lengthVec;
	vector<int> frequencyVec;
};
