#include "templateClass.h"

templateClass::templateClass() {}
void templateClass::generateTree(string treeFileName) {
	tree.generateTree(treeFileName);
}

void templateClass::loadTree(string treeFileName) {
	tree.loadTree(treeFileName);
}

void templateClass::setThresholdVec(const vector<int>& thresholdVec) {
	this->thresholdVec = thresholdVec;
	this->frequencyVec.resize(thresholdVec.size(), -1);
	this->lengthVec.resize(thresholdVec.size(), 0);
}

void templateClass::analyzeRead(string read) {
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

	this->read = read;
}

void templateClass::printInfo() {
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
