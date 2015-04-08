/*
 * KmerHashMain.cc
 *
 *  Created on: Aug 10, 2013
 *      Author: hxin
 */

#include "KmerHash.h"
#include <iostream>
#include <string>

using namespace std;

int main() {
	KmerHash kmer;

	while (1) {
		string Kmer;
		int hash;
		char input;
		cout << "1: get kmer, 2 get hash, q quit: ";
		cin >> input;
		switch (input) {
		case '1':
			cout << "set KmerLengh" << endl;
			cin >> hash;
			kmer.setKmerSize(hash);
			cout << "Input hash: " << endl;
			cin >> hash;
			kmer.setHashVal(hash);

			cout << "Kmer: " << kmer.getKmer() << endl;
			break;

		case '2':
			cout << "Input Kmer" << endl;
			cin >> Kmer;
			kmer.calculateHash(Kmer, Kmer.size() );

			cout << "Hash: " << kmer.getHash() << endl;
			break;

		case 'q':
			return 0;

		default:
			break;
		}
	}
}
