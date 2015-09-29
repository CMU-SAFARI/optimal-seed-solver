#include "HashTree.h"
#include "KmerHash.h"
#include <cstdlib>

using namespace std;

int main(int argc, const char* argv[]) {
	HashTree tree;
	KmerHash hash;
	unsigned int DNAHash;
	unsigned int estLength;


	if (argc != 2 && argc != 3) {
		cerr << "Help:" << endl;
		cerr << "To generate the tree, please run: $ ./testHashTree referenceName.fasta treeName.tree" << endl;
		cerr << "To query the tree, please run: $ ./testHashTree treeName.tree" << endl;
		exit(1);
	}

	if (argc == 3) {
		tree.generateTree(argv[1]);
		tree.storeTree(argv[2]);
		exit(0);
	}
	tree.loadTree(argv[1]);
	hash.setKmerSize(tree.getHashLength() );

	while (1) {
		cout << "Please input mode: (f) frequency test; (s) seed test (a) average test (t) table test: " << endl;
		char option;
		cin >> option;
		string DNA;
		int threshold;
		int depth;
		int avgfreq;
		int index;
		switch (option) {
			case 's':

				cout << "please input DNA string: " << endl;
				cin >> DNA;

				cout << "please input threshold: " << endl;
				cin >> threshold;

				tree.test(DNA);
				tree.query(DNA, threshold);
				cout << "eff length: " << tree.effectiveLength() << " frequency: " << tree.frequency() << " isLeaf: " << tree.isLeaf() << endl;
				break;

			case 'f':

				cout << "please input DNA string: " << endl;
				cin >> DNA;

				cout << "please input threshold: " << endl;
				cin >> threshold;

				DNAHash = hash.calculateHash(DNA, tree.getHashLength() );
				cout << "DNA Hash: " << DNAHash << endl;
				estLength = tree.estimateFreqLength(DNAHash, threshold/*, DNA*/);
				cout << "Estimated extension length to achieve " << threshold << " locations of DNA " << DNA << " is " << estLength << endl;
				break;

			case 'a':
				cout << "please tree hash: " << endl;
				cin >> index;

				cout << "please input depth: " << endl;
				cin >> depth;

				avgfreq = tree.averageFreqDepth(index, depth);
				cout << "Average weighted frequency at depth " << depth << " is " << avgfreq << endl;
				break;

			case 't':
				cout << "Testing table generation." << endl;

				cout << "Calculating table." << endl;
				tree.calculateLenToFreqTable();

				cout << "Printing Table." << endl;
				tree.printHashLenTable("table.txt");

				cout << "Storing table." << endl;
				tree.storeHashLenTable("table.tb");

				cout << "Loading table." << endl;
				tree.readHashLenTable("table.tb");

				cout << "Printing table." << endl;
				tree.printHashLenTable("table2.txt");
				break;

			default:
				cout << "Unrecognized option" << endl;
				break;
		}
	}

	return 0;
}
