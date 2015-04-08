#include "RefDB.h"
#include "KmerHash.h"

using namespace std;

class SpacedSeedSolver {
public:
	SpacedSeedSolver();
	~SpacedSeedSolver();
	void loadHash(string tableName);
	void storeHash(string tableName);
	void generateHash(string refName, string mask);
	unsigned int solveDNA(string DNA, unsigned int seedNum);
	void print();

private:
	RefDB refDB;
	KmerHash hashGen;
	string seed;

	vector <unsigned int> freqTable;
	string mask;
};
