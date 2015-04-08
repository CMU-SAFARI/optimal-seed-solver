/*
 * KmerHash.h
 *
 *  Created on: Dec 2, 2012
 *      Author: hxin
 */

#ifndef KMERHASH_H_
#define KMERHASH_H_

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class KmerHash {
public:
	KmerHash();
	KmerHash(const string& DNA, unsigned KmerSize);
	~KmerHash();
	KmerHash(const KmerHash& rhs);

	void setKmerSize(unsigned int kmerSize);
	void setHashVal(unsigned int hashVal);
	unsigned int getHash() const;
	
	//Get kmer
	string getKmer() const;


	//Calculate Hash
	unsigned int calculateHash(const string& DNA, unsigned int KmerSize);

private:
	unsigned int hashVal;
	unsigned int KmerSize;
};

#endif
