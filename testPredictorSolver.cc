#include "predictorSolver.h"
#include "HashTree.h"
#include "KmerHash.h"
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <cstring>
#include <map>

using namespace std;

int main(int argc, const char* argv[]) {
	//Solve reads
	PredictorSolver solver;
	HashTree tree;

	map <int, unsigned int> predictedFreqCounter;
	map <int, unsigned int> realFreqCounter;
	map <int, unsigned int> differenceCounter;
	map <double, unsigned int> precentDifferenceCounter;

	if (argc != 4) {
		cerr << "USAGE: ./testPredictorSolver <Table File> <Tree File> <Read File>" << endl;
		exit(1);
	}

	cout << "Loading table " << argv[1] << " into memory.\n";
	solver.loadTable(argv[1]);
	cout << "Finished loading. \n";

	cout << "Loading tree " << argv[2] << " into memory.\n";
	tree.loadTree(argv[2]);	
	cout << "Finished loading. \n";

	char benchFileName[80];

	strcpy(benchFileName, argv[3]);

	char* benchName = strtok(benchFileName, ".");
	
	//unsigned int hashVal;
	int thresholdIdx;

	//unsigned int length;
	//unsigned int threshold;

	unsigned int curLength;
	int seedNum;
	int predictedFreq;
	int realFreq;
	unsigned int *lengthVec = NULL;
	cin >> seedNum;

	ifstream readFile;
	readFile.open(argv[3]);
	string read;

	ofstream predictOutput;
	ofstream realOutput;
	ofstream diffOutput;
	ofstream percentDiffOutput;

	predictOutput.open( (string(benchName) + string(".predict_predict") ).c_str() );
	realOutput.open( (string(benchName) + string(".predict_real") ).c_str() );
	diffOutput.open( (string(benchName) + string(".predict_diff") ).c_str() );
	percentDiffOutput.open( (string(benchName) + string(".predict_\%diff") ).c_str() );

	cout << "benchName: " << string(benchName) << endl;
	
	while (seedNum > 0) {
		//clean up last iteration data
		predictedFreqCounter.clear();
		realFreqCounter.clear();
		differenceCounter.clear();
		precentDifferenceCounter.clear();

		if (lengthVec != NULL) {
			delete [] lengthVec;
			lengthVec = NULL;
		}
		lengthVec = new unsigned int [seedNum];

		readFile.clear();
		readFile.seekg(0, ios::beg);

		//Name
		getline(readFile, read);
		//Read
		getline(readFile, read);
//		cout << "DNA: " << read << endl;

		//TODO: Solve
		thresholdIdx = solver.solveIdx(read, seedNum);
		predictedFreq = solver.solveLength(read, seedNum, thresholdIdx, lengthVec);

		if (thresholdIdx < 0)
			realFreq = -1;
		else {
			curLength = 0;
			realFreq = 0;
			for (int i = 0; i < seedNum; i++) {
				//cout << "seed: " << read.substr(curLength, lengthVec[i]) << endl;
				tree.query(read.substr(curLength, lengthVec[i]), solver.getThreshold(thresholdIdx) );
				realFreq += tree.frequency();
				curLength += lengthVec[i];
			}
		}

		if (realFreqCounter.find(realFreq) == realFreqCounter.end() )
			realFreqCounter[realFreq] = 1;
		else
			realFreqCounter[realFreq]++;

		if (predictedFreqCounter.find(predictedFreq) == predictedFreqCounter.end() )
			predictedFreqCounter[predictedFreq] = 1;
		else
			predictedFreqCounter[predictedFreq]++;

		int freqDiff = realFreq - predictedFreq;
		if (differenceCounter.find(freqDiff) == differenceCounter.end() )
			differenceCounter[freqDiff] = 1;
		else
			differenceCounter[freqDiff]++;

		double percentFreqDiff = (double) freqDiff / predictedFreq;
		if (precentDifferenceCounter.find(percentFreqDiff) == precentDifferenceCounter.end() )
			precentDifferenceCounter[percentFreqDiff] = 1;
		else
			precentDifferenceCounter[percentFreqDiff]++;

		//Trash
		getline(readFile, read);
		getline(readFile, read);
			
		//Name
		getline(readFile, read);

		//int debugCount = 1;

		while (readFile.good()) {
			//cout << "debugCount: " << debugCount << endl;

			//Read
			getline(readFile, read);
//			cout << "DNA: " << read << endl;

			//TODO: Solve
			if (read.find('N') == string::npos) {
				thresholdIdx = solver.solveIdx(read, seedNum);
				predictedFreq = solver.solveLength(read, seedNum, thresholdIdx, lengthVec);

				curLength = 0;
				
				if (thresholdIdx < 0)
					realFreq = -1;
				else {
					realFreq = 0;
					for (int i = 0; i < seedNum; i++) {
						//cout << "seed: " << read.substr(curLength, lengthVec[i]) << endl;
						tree.query(read.substr(curLength, lengthVec[i]), solver.getThreshold(thresholdIdx) );
						realFreq += tree.frequency();
						curLength += lengthVec[i];
					}
				}

				//Gathering data
				if (realFreqCounter.find(realFreq) == realFreqCounter.end() )
					realFreqCounter[realFreq] = 1;
				else
					realFreqCounter[realFreq]++;

				if (predictedFreqCounter.find(predictedFreq) == predictedFreqCounter.end() )
					predictedFreqCounter[predictedFreq] = 1;
				else
					predictedFreqCounter[predictedFreq]++;

				unsigned int freqDiff = realFreq - predictedFreq;
				if (differenceCounter.find(freqDiff) == differenceCounter.end() )
					differenceCounter[freqDiff] = 1;
				else
					differenceCounter[freqDiff]++;

				double percentFreqDiff = (double) freqDiff / predictedFreq;
				if (precentDifferenceCounter.find(percentFreqDiff) == precentDifferenceCounter.end() )
					precentDifferenceCounter[percentFreqDiff] = 1;
				else
					precentDifferenceCounter[percentFreqDiff]++;

			}

			//Trash
			getline(readFile, read);
			getline(readFile, read);
			
			//Name
			getline(readFile, read);

			//debugCount++;
		}

		cout << "seedNum: " << seedNum << endl;
		predictOutput << "seedNum: " << seedNum << endl;
		realOutput << "seedNum: " << seedNum << endl;
		diffOutput << "seedNum: " << seedNum << endl;
		percentDiffOutput << "seedNum: " << seedNum << endl;

		for (map<int, unsigned int>::iterator iter = predictedFreqCounter.begin(); iter != predictedFreqCounter.end(); iter++)
			predictOutput << iter->first << ": " << iter->second << endl;
		
		for (map<int, unsigned int>::iterator iter = realFreqCounter.begin(); iter != realFreqCounter.end(); iter++)
			realOutput << iter->first << ": " << iter->second << endl;

		for (map<int, unsigned int>::iterator iter = differenceCounter.begin(); iter != differenceCounter.end(); iter++)
			diffOutput << iter->first << ": " << iter->second << endl;

		for (map<double, unsigned int>::iterator iter = precentDifferenceCounter.begin(); iter != precentDifferenceCounter.end(); iter++)
			percentDiffOutput << iter->first << ": " << iter->second << endl;

		cin >> seedNum;
	}
	
	predictOutput.close();
	realOutput.close();
	diffOutput.close();
	percentDiffOutput.close();
/*
	//load table
	PredictorSolver solver;

	if (argc != 2) {
		cerr << "Usage: $> tableFileName" << endl;
		exit(1);
	}

	solver.loadTable(argv[1]);

	unsigned int hashVal;
	unsigned int thresholdIdx;

	unsigned int length;
	unsigned int threshold;

	for (hashVal = 0; hashVal < solver.getNumTotalEntries(); hashVal++) {
		for (thresholdIdx = 0; thresholdIdx < solver.getNumThresholds(); thresholdIdx++) {
			length = solver.getLength(hashVal, thresholdIdx);
			assert(length >= solver.getHashLength() );
		}
	}
	
	while (1) {
		cout << "Please input hashVal: " << endl;
		cin >> hashVal;
		cout << "Please input thresholdIdx: " << endl;
		cin >> thresholdIdx;

		threshold = solver.getThreshold(thresholdIdx);
		length = solver.getLength(hashVal, thresholdIdx);

		cout << "threshold: " << threshold << " length: " << length << endl;
	}
*/	
/*
	//Generate table
	PredictorSolver solver;

	if (argc != 3) {
		cerr << "Usage: $> treeFileName tableFileName" << endl;
		exit(1);
	}
	cout << "please input thresholdNum: " << endl;
	unsigned int thresholdNum;
	cin >> thresholdNum;
	unsigned int* thresholds = new unsigned int [thresholdNum];
	for (int i = 0; i < thresholdNum; i++) {
		cout << "please input thresholdNum: " << endl;
		cin >> thresholds[i];
	}
	solver.generate(argv[1], thresholdNum, thresholds);
	solver.storeTable(argv[2]);
*/
	return 0;
}
