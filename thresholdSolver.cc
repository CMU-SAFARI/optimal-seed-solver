#include "thresholdSolver.h"

ThresholdSolver::ThresholdSolver() {
	this->verbose = false;
	this->seedNum = -1;
}

void ThresholdSolver::generateTree(string treeFileName) {
	tree.generateTree(treeFileName);
}

void ThresholdSolver::loadTree(string treeFileName) {
	tree.loadTree(treeFileName);
}

void ThresholdSolver::init(int readLength, int seedNum)
{
	this->seedNum = seedNum;
	vector< int > thresholdVec (seedNum);
	for(int i = 0; i < thresholdVec.size(); i++)
		thresholdVec[i] = readLength;
	this->thresholdVec = thresholdVec;
	this->frequencyVec.resize(thresholdVec.size(), -1);
	this->lengthVec.resize(thresholdVec.size(), 0);
}


void ThresholdSolver::setThresholdVec(const vector<int>& thresholdVec) {
	this->thresholdVec = thresholdVec;
	this->frequencyVec.resize(thresholdVec.size(), -1);
	this->lengthVec.resize(thresholdVec.size(), 0);
}

void ThresholdSolver::analyzeRead(string read) {
	int minLength = tree.getHashLength();
	int maxLength = tree.getFullLength();

	int pos = 0;

	for (int i = 0; i < thresholdVec.size(); i++) {
		string seed = read.substr(pos, maxLength);

		tree.query(seed, thresholdVec[i]);
		frequencyVec[i] = tree.frequency();
		lengthVec[i] = tree.effectiveLength();
		pos += tree.effectiveLength();

		if (pos > read.length() - minLength)
			break;
	}
}
void ThresholdSolver::load(string file)
{
	std::ifstream test;
	test.open(file.c_str(), ios::in);
	string line;
	int num = 0;
	if(test.is_open())
	{
		while(std::getline(test,line))
			num++;
	}

	int SIZE = thresholdVec.size();
	std::ifstream f;
	f.open(file.c_str(), ios::in);

	if(f.is_open())
	{
		int i = 0;

		vector< vector <int> > freqTable (num/4, vector<int>(SIZE));
		vector< vector <int> > seedLenTable (num/4, vector<int>(SIZE));

		vector<int> avgFreqTable (SIZE);
		vector<int> avgSeedLenTable (SIZE);
		vector<int> lettersLeft (num/4);

		vector<int> totalFrequency (100000);
		for(int j=0; j<100000; j++) totalFrequency[j] = 0;

		int avgLettersLeft = 0;

		int pass = 0;
		int total = 0;

		this->freqTable = freqTable;
		this->seedLenTable = seedLenTable;
		this->lettersLeft = lettersLeft;
		this->avgFreqTable = avgFreqTable;
		this->avgSeedLenTable = avgSeedLenTable;
		this->totalFrequency = totalFrequency;

		while(!f.eof())
		{
			string id, id2, read, plus, quality;
			f >> id; f >> id2; f >> read; f >> plus; f >> quality;

			this->index = i;
			if(read.find('N') != string::npos)
			{
				if(verbose)
				{
					cout << "Read #" << i << "failed." << endl;
					cout << "Found N ... fail" << endl;
				}
			}
			else if(read.length() > 0){
				int ret = solve(read);
				if(ret != -1)
				{
					for(int j = 0; j < SIZE; j++)
					{
						this->avgFreqTable[j] += this->freqTable[i][j];
						this->avgSeedLenTable[j] += this->seedLenTable[i][j];
					}
					pass++;
				}
				else
					avgLettersLeft += this->lettersLeft[i];
				int s = this->totalFrequency.size();
				if(ret + 1 >= s)
					this->totalFrequency.resize(2*s, 0);
				total++;
				this->totalFrequency[ret + 1]++;
			}

			i++;
		}
		if(verbose)
		{
			cout << "# of Valid reads: " << pass << "/" << total << endl << endl;

			if(pass > 0)
			{
				cout << "For valid reads ... " << endl;
				cout << "Average Frequency Table:" << endl;
				for(int j = 0; j < SIZE; j++)
					cout << this->avgFreqTable[j]/pass << " ";
				cout << endl;
				cout << "Average SeedLength Table:" << endl;
				for(int j = 0; j < SIZE; j++)
					cout << this->avgSeedLenTable[j]/pass << " ";
				cout << endl;
				cout << endl;
			}
			if(total - pass > 0)
			{
				cout << "For invalid reads ... " << endl;
				cout << "Average lettersLeft: " << avgLettersLeft/(total - pass) << endl;
				cout << endl << endl;
			}
		}
	}
	if(verbose)
		printFreqTable();
}
int ThresholdSolver::solve(string read)
{
	int index = this->index;
	if(verbose)
	{
		cout << "Reading read #" << index << ": " << read << endl;
		cout.flush();
		cout << "Original length: " << read.length() << endl;
	}

	int minLength = tree.getHashLength();
	int maxLength = tree.getFullLength();

	int pos = 0;
	int count = 0;

	bool pass = true;

	int ll = 0;
	int totalfreq = 0;

	this->read = read;

	for (int i = 0; i < thresholdVec.size(); i++) {

		int currLength = maxLength;

		if(read.length() - pos < maxLength)
		{
			currLength = read.length() - pos;
			if(currLength <= minLength)
			{
				// Read too small
				if(verbose)
					cout << "Read too small" << endl;
				pass = false;
				break;
			}
		}

		string seed = read.substr(pos, currLength);
		if(verbose)
			cout << "Position " << pos << "." << endl;

		tree.query(seed, thresholdVec[i]);
		int len = tree.effectiveLength();
		int freq = tree.frequency();

		if(verbose)
			cout << "Effective length: " << len << " Frequency: " << freq << endl;

		if(freq <= thresholdVec[i])
		{
			//Case 1
			//this works
			if(verbose)
				cout << "Frequency is <= threshold." << endl;
			pos += len;
			//freqTable[index][count] = freq;
			//seedLenTable[index][count] = len;
		}
		else if(tree.isLeaf())
		{
			//Case 2
			//hits leaf before threshold - ok
			if(verbose)
				cout << "Hit leaf before threshold." << endl;
			pos += len;
			//freqTable[index][count] = freq;
			//seedLenTable[index][count] = len;
		}
		else
		{
			//Case 3
			//does not meet threhold, not leaf - fail
			if(verbose)
				cout << "Leaf does not meet threshold = fail." << endl;
			pass = false;
			break;
		}

		count++;
		if (pos == read.length())
		{
			if(verbose)
				cout << "Position at read length. Success." << endl;
			break;
		}
		if (pos > read.length() - minLength)
		{
			if(verbose)
				cout << "Next read size too small." << endl;
			pass = false;
			break;
		}
		totalfreq += freq;
	}

	if(!pass)
	{
		//Did not pass - calculate letters left
		/*
		ll = read.length() - pos;
		if(verbose)
			cout << ll << " letters left." << endl;
		lettersLeft[index] = ll;
		*/
		return -1;
	}
	if(verbose)
		cout << "Finished read at position " << pos << " from " << read.length() << " with totalFrequency " << totalfreq << ".\n" << endl;

	return totalfreq;
}

void ThresholdSolver::printFreqTable() {
	cout << "Printing Frequency Table." << endl;
	for (int i = 0; i < freqTable.size(); i++)
	{
		for (int j = 0; j<freqTable[i].size(); j++)
			cout << freqTable[i][j] << " ";
		cout << "" << endl;
	}
}

void ThresholdSolver::printTotalFreqTable() {
	cout << "seedNum: " << seedNum << endl;
	for(int i = 0; i < totalFrequency.size(); i++)
	{
		if(totalFrequency[i] != 0)
			cout << i-1 << ": " << totalFrequency[i] << endl;
	}
}

void ThresholdSolver::printInfo() {
	cout << "Thresholds: ";
	for (int i = 0; i < thresholdVec.size(); i++)
		cout << "\t|" << thresholdVec[i];

	cout << endl << "Lengths:";
	for (int i = 0; i < thresholdVec.size(); i++)
		cout << "\t|" << lengthVec[i];

	cout << endl << "Frequencys:";
	for (int i = 0; i < thresholdVec.size(); i++)
		cout << "\t|" << frequencyVec[i];

	int pos = 0;
	cout << endl << "Seeds:" << endl;
	for (int i = 0; i < thresholdVec.size(); i++) {
		cout << read.substr(pos, lengthVec[i]) << endl;
		pos += lengthVec[i];
		if (lengthVec[i] == 0)
			break;
	}
}
