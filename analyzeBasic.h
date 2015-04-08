#include "RefDB.h"
#include "KmerHash.h"

using namespace std;

class AnalyzeBasic {
public:
	AnalyzeBasic();
	~AnalyzeBasic();
	void loadHash(string tableName);
	void storeHash(string tableName);
	void generateHash(string refName, unsigned int seedLength);
	unsigned int checkFreq(string seed);
	unsigned int getLength();
	void print();

private:
	RefDB refDB;
	KmerHash hashGen;

	unsigned int seedLength;
	vector <unsigned int> freqTable;
};
