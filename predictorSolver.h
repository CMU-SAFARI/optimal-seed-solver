#include "HashTree.h"
#include "KmerHash.h"
#include <string>

using namespace std;

class PredictorSolver {
public:
	PredictorSolver ();
	~PredictorSolver ();
	
	void init(unsigned int hashLength, unsigned int numThresholds, unsigned int * thresholds);

	unsigned int getLength(unsigned int hashVal, int thresholdIdx);
	unsigned int getNumThresholds();
	unsigned int getHashLength();
	unsigned int getNumTotalEntries();
	unsigned int getThreshold(unsigned int thresholdIdx);

	int solveIdx(string& DNA, unsigned int numOfSeeds);
	int solveLength(string& DNA, unsigned int numOfSeeds, int idx, unsigned int* lengthVec);

	void generate(string treeFileName, unsigned int numThresholds, unsigned int * thresholds);
	void storeTable(string tableFileName);
	void loadTable(string tableFileName);

private:
	bool idxIsValid(string& DNA, unsigned int numOfSeeds, unsigned int idx);

	//hash generator
	KmerHash hashGen;

	//data
	unsigned int ** table;
	//Important: thresholds have to follow a decreasing order
	unsigned int * thresholds;
	unsigned int numThresholds;
	unsigned int hashLength;
	unsigned int numTotalEntries;
};

