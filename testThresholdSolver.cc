#include "thresholdSolver.h"
#include <vector>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cstring>
#include <map>

using namespace std;

int main(int argc, const char* argv[]) {
	map <int, unsigned long long> frequencyCounter;

	ThresholdSolver solver;

	if (argc != 3) {
		cerr << "USAGE: ./testOptimalSolver <Tree File> <Read File>" << endl;
		exit(1);
	}

	cout << "Loading tree " << argv[1] << " into memory.\n";
	solver.loadTree(argv[1]);
	cout << "Finished loading. \n";

	int seedNum;
	int frequency;
	//cout << "Please input number of seeds: ";
	cin >> seedNum;

	ifstream readFile;
	readFile.open(argv[2]);
	string read;
	
	char benchFileName[80];
	strcpy(benchFileName, argv[2]);
	char* benchName = strtok(benchFileName, ".");
	
	ofstream output;
	output.open( (string(benchName) + string(".threshold") ).c_str() );

	int threshold = 0;
	vector<int> thresholds;

	while (seedNum > 0) {
		cin >> threshold;
		thresholds.clear();
		thresholds.resize(seedNum, threshold);
		solver.setThresholdVec(thresholds);

		frequencyCounter.clear();

		readFile.clear();
		readFile.seekg(0, ios::beg);

		//Name
		getline(readFile, read);
		//Read
		getline(readFile, read);
		//cout << "DNA: " << read << endl;

		//Solve
		frequency = solver.solve(read);
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
			//cout << "DNA: " << read << endl;

			//Solve
			if (read.find('N') == string::npos) {
				frequency = solver.solve(read);
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

		cout << "seedNum: " << seedNum << " | threshold: " << threshold << endl;
		output << "seedNum: " << seedNum << " | threshold: " << threshold << endl;
		for (map<int, unsigned long long>::iterator iter = frequencyCounter.begin(); iter != frequencyCounter.end(); iter++)
			output << iter->first << ": " << iter->second << endl;
		
		cin >> seedNum;
	}

	return 0;
}
