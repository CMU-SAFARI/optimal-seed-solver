#include "optimalSolverLN.h"
#include <vector>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cstring>
#include <iostream>
#include <map>

using namespace std;

int main(int argc, const char* argv[]) {
	map <unsigned int, unsigned long long> frequencyCounter;

	OptimalSolverLN solver;

	if (argc == 3) {
		cout << "Loading tree " << argv[1] << " into memory.\n";
		solver.loadTree(argv[1]);
		cout << "Finished loading. \n";

		int seedNum;
		unsigned int frequency;
		//cout << "Please input number of seeds: ";
		cin >> seedNum;

		ifstream readFile;
		readFile.open(argv[2]);
		string read;

		char benchFileName[80];
		strcpy(benchFileName, argv[2]);
		char* benchName = strtok(benchFileName, ".");

		ofstream output;
		ofstream freqOutput;
		output.open( (string(benchName) + string(".optimalLN") ).c_str() );
		freqOutput.open( (string(benchName) + string(".optimalLN_freq") ).c_str() );

		while (seedNum > 0) {
			cout << "seedNum: " << seedNum << endl;
			output << "seedNum: " << seedNum << endl;
			freqOutput << "seedNum: " << seedNum << endl;

			frequencyCounter.clear();

			readFile.clear();
			readFile.seekg(0, ios::beg);

			//Name
			getline(readFile, read);
			//Read
			getline(readFile, read);
			solver.init(read.length(), seedNum);
			//cout << "DNA: " << read << endl;

			//Solve
			frequency = solver.solveDNA(read);
			if (frequencyCounter.find(frequency) == frequencyCounter.end() )
				frequencyCounter[frequency] = 1;
			else
				frequencyCounter[frequency]++;

			//Trash
			getline(readFile, read);
			getline(readFile, read);

			//Name
			getline(readFile, read);

			while (!readFile.eof()) {
				//Read
				getline(readFile, read);
				//			cout << "DNA: " << read << endl;

				//Solve
				if (read.find('N') == string::npos) {
					frequency = solver.solveDNA(read);
					solver.backtrack();
					solver.sortOfFreq();
					solver.printFreqs(freqOutput);
					solver.printLength(freqOutput);
					if (frequencyCounter.find(frequency) == frequencyCounter.end() )
						frequencyCounter[frequency] = 1;
					else
						frequencyCounter[frequency]++;
				}

				//Trash
				getline(readFile, read);
				getline(readFile, read);

				//Name
				getline(readFile, read);
			}

			for (map<unsigned int, unsigned long long>::iterator iter = frequencyCounter.begin(); iter != frequencyCounter.end(); iter++)
				output << iter->first << ": " << iter->second << endl;

			cin >> seedNum;
		}
		output.close();
		freqOutput.close();

	}
	else if (argc == 1) {
		cout << "input seedNum and readLength: " << endl;
		int seedNum;
		int readLength;
		cin >> seedNum;
		string dummy = "";
		while (seedNum > 0) {
			cin >> readLength;
			solver.setMinLength(10);
			solver.init(readLength, seedNum);
			solver.feedL0();
			solver.solveDNA(dummy);
			solver.backtrack();
			solver.sortOfFreq();
			solver.printFreqs(cout);
			solver.printStats(cout);
			cin >> seedNum;
		}
	}
	else
	{
		cerr << "USAGE: ./testOptimalSolver <Tree File> <Read File>" << endl;
		exit(1);
	}

	//	cout << "Loading tree " << argv[1] << " into memory.\n";
	//	solver.loadTree(argv[1]);
	//	cout << "Finished loading. \n";
	//
	//	int seedNum;
	//	unsigned int frequency;
	//	//cout << "Please input number of seeds: ";
	//	cin >> seedNum;
	//
	//	ifstream readFile;
	//	readFile.open(argv[2]);
	//	string read;
	//
	//	while (seedNum > 0) {
	//		frequencyCounter.clear();
	//
	//		readFile.clear();
	//		readFile.seekg(0, ios::beg);
	//
	//		//Name
	//		getline(readFile, read);
	//		//Read
	//		getline(readFile, read);
	//		solver.init(read.length(), seedNum);
	//		//		cout << "DNA: " << read << endl;
	//
	//		//Solve
	//		frequency = solver.solveDNA(read);
	//		if (frequencyCounter.find(frequency) == frequencyCounter.end() )
	//			frequencyCounter[frequency] = 1;
	//		else
	//			frequencyCounter[frequency]++;
	//
	//		//Trash
	//		getline(readFile, read);
	//		getline(readFile, read);
	//
	//		//Name
	//		getline(readFile, read);
	//
	//		while (!readFile.eof()) {
	//			//Read
	//			getline(readFile, read);
	//			//			cout << "DNA: " << read << endl;
	//
	//			//Solve
	//			if (read.find('N') == string::npos) {
	//				frequency = solver.solveDNA(read);
	//				if (frequencyCounter.find(frequency) == frequencyCounter.end() )
	//					frequencyCounter[frequency] = 1;
	//				else
	//					frequencyCounter[frequency]++;
	//			}
	//
	//			//Trash
	//			getline(readFile, read);
	//			getline(readFile, read);
	//
	//			//Name
	//			getline(readFile, read);
	//		}
	//
	//		cout << "seedNum: " << seedNum << endl;
	//		for (map<unsigned int, unsigned long long>::iterator iter = frequencyCounter.begin(); iter != frequencyCounter.end(); iter++)
	//			cout << iter->first << ": " << iter->second << endl;
	//
	//		cin >> seedNum;
	//	}

	return 0;
}
