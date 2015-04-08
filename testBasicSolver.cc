#include "basicSolver.h"
#include <vector>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cstring>
#include <map>

using namespace std;

int main(int argc, const char* argv[]) {
	int frequency;
	map <int, unsigned long long> frequencyCounter;

	BasicSolver solver;

	if (argc != 3) {
		cerr << "USAGE: ./testOptimalSolver <Tree File> <Read File>" << endl;
		exit(1);
	}

	cout << "Loading tree " << argv[1] << " into memory.\n";
	solver.loadTree(argv[1]);
	cout << "Finished loading. \n";

	int seedNum;
	//cout << "Please input number of seeds: ";
	cin >> seedNum;

	int seedLength;

	ifstream readFile;
	readFile.open(argv[2]);
	string read;
	
	char benchFileName[80];
	strcpy(benchFileName, argv[2]);
	char* benchName = strtok(benchFileName, ".");
	
	ofstream output;
	output.open( (string(benchName) + string(".analyzeBasic") ).c_str() );

	while (seedNum > 0) {	
		//cout << "Please input seed length: ";
		cin >> seedLength;
		
		frequencyCounter.clear();

		//Start reading the file
		readFile.clear();
		readFile.seekg(0, ios::beg);

		//Name
		getline(readFile, read);
		//Read
		getline(readFile, read);
		//cout << "DNA: " << read << endl;
		solver.init(read.length(), seedNum, seedLength);

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
			//cout << "DNA: " << read << endl;

			//Solve
			if (read.find('N') == string::npos) {
				frequency = solver.solveDNA(read);
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

		cout << "seedNum: " << seedNum << " | seedLength: " << seedLength << endl;
		output << "seedNum: " << seedNum << " | seedLength: " << seedLength << endl;
		for (map<int, unsigned long long>::iterator iter = frequencyCounter.begin(); iter != frequencyCounter.end(); iter++)
			output << iter->first << ": " << iter->second << endl;

		cin >> seedNum;
	}
	
	return 0;
}
