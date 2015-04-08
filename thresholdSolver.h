#include "HashTree.h"
#include <vector>
#include <string>

using namespace std;

class ThresholdSolver{
public:
	ThresholdSolver();
	void generateTree(string treeFileName);
	void loadTree(string treeFileName);

	void setThresholdVec(const vector<int>& thresholdVec);

	void init(int readLength, int seedNum);

	void analyzeRead(string Read);
	void printInfo();
	void printFreqTable();
	void printTotalFreqTable();

    vector<vector <int> > missTable;
  	vector<vector <int> > penaltyTable;

  	vector<vector <int> > freqTable;
	vector<vector <int> > seedLenTable;

	vector<int> avgFreqTable;
	vector<int> avgSeedLenTable;
	vector<int> lettersLeft;
	vector<int> totalFrequency;

  	int threshold;
  	int seedLength;
  	int seedNum;

  	bool verbose;

	void load(string file);
	int solve(string read);

private:
	HashTree tree;
	int index;
	string read;
	vector<int> thresholdVec;
	vector<int> lengthVec;
	vector<int> frequencyVec;
};
