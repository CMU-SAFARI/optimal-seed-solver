/*
 * refDB.h
 *
 *  Created on: May 30, 2013
 *      Author: hxin
 */

#ifndef REFDB_H_
#define REFDB_H_

#include <string>
#include <fstream>
#include <vector>

using namespace std;

class RefDB {
public:
	RefDB();
	~RefDB();
	void loadRefFile(string refName);
	bool loadChromo(int chromoNum);
	void unloadAll();

	unsigned int getChromoLength();
	unsigned int getNumOfChromo();
	string getChromoName();
	string getRefSeq(unsigned int refLoc, unsigned int length);
	char getRefBase(unsigned int refLoc, unsigned int length);

private:
	ifstream refFile;
	vector<streampos> refPos;
	string chromoName;
	string chromosome;
	int curChromoNum;
};

#endif /* REFDB_H_ */

