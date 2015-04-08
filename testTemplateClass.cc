#include "templateClass.h"
#include <vector>
#include <cstdlib>

using namespace std;

int main(int argc, const char* argv[]) {
	templateClass test;
	if (argc != 2) {
		cerr << "Please also put the tree file name" << endl;
		exit(1);
	}

	test.loadTree(argv[1]);

	while (1) {
		cout << "please input DNA string: " << endl;
		string DNA;
		cin >> DNA;

		cout << "please input number of thresholds: ";
		int input;
		cin >> input;
		
		cout << "please input thresholds, one by one:" << endl;

		vector<int> thresholds(input, 1);
		for (int i = 0; i < thresholds.size(); i++)
			cin >> thresholds[i];

		test.setThresholdVec(thresholds);
		test.analyzeRead(DNA);
		test.printInfo();
	}

	return 0;
}
