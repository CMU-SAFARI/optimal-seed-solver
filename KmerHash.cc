#include "KmerHash.h"
#include <cassert>
#include <iostream>
#include <cstdlib>

KmerHash::KmerHash() {
	KmerSize = 10;
	hashVal = 0;
}

KmerHash::KmerHash(const string& Kmer, unsigned int KmerSize) {
	//cout << "KmerSize: " <<  KmerSize << endl;
	hashVal = 0;
	calculateHash(Kmer, KmerSize);
}

KmerHash::~KmerHash() {
}

KmerHash::KmerHash(const KmerHash& rhs) {
	KmerSize = rhs.KmerSize;
	hashVal = rhs.hashVal;
}

void KmerHash::setKmerSize(unsigned int kmerSize) {
	this->KmerSize = kmerSize;
}

void KmerHash::setHashVal(unsigned int hashVal) {
	this->hashVal = hashVal;
}

string KmerHash::getKmer() const {
	string Kmer(KmerSize, 'N');

	cout << "kmer size: " << KmerSize << endl;

	for (int i = 0; i < KmerSize; i++) {

		unsigned int bpVal = (hashVal >> (2 * (KmerSize - i) ) ) & 3;
		switch (bpVal) {
		case 0:
			Kmer[i] = 'A';
			break;
		case 1:
			Kmer[i] = 'C';
			break;
		case 2:
			Kmer[i] = 'G';
			break;
		case 3:
			Kmer[i] = 'T';
			break;
		default:
			cerr << "unrecognized character: " << Kmer[i] << endl;
			break;
		}

	}
	return Kmer;
}

unsigned int KmerHash::calculateHash(const string& Kmer, unsigned int KmerSize) {
	assert (KmerSize <= Kmer.size() );
	
	//Clear hashVal
	hashVal = 0;
	//Fill different Lv of the hash
	unsigned long bpVal = 0;
	this->KmerSize = KmerSize;

	unsigned int tempHashVal = 0;

	for (int i = 0; i < this->KmerSize; i++) {

		//cout << Kmer[i];

		switch (Kmer[i]) {
		case 'A':
			bpVal = 0;
			break;
		case 'C':
			bpVal = 1;
			break;
		case 'G':
			bpVal = 2;
			break;
		case 'T':
			bpVal = 3;
			break;
		case 'N':
		default:
			bpVal = 0;
		}

		hashVal = (hashVal << 2) | bpVal;
	}
	return hashVal;
}

//bool KmerHash::compare(const KmerHash& lhs, const KmerHash& rhs) {
//	return (lhs.hashVal < rhs.hashVal);
//}

unsigned int KmerHash::getHash() const {
	return hashVal;
}
