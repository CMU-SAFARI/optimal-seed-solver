#ifndef HASHTREE_H_
#define HASHTREE_H_

#define MAX_LENGTH 30

#include "KmerHash.h"
#include "RefDB.h"
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;

class TreeNode {
public:
	TreeNode();
	TreeNode(unsigned int locNum);
	//Copy constructor. Copies the entire children branch
	TreeNode(const TreeNode& rhs);
	virtual ~TreeNode();

	//this is a virtual function
	virtual void setContent(unsigned int chrNum, unsigned int location){};

	//set member variables
	void setLocNum(unsigned int locNum);
	void addLoc();
	void setA(TreeNode* A);
	void setC(TreeNode* C);
	void setT(TreeNode* T);
	void setG(TreeNode* G);
	void deleteChildren();

	//get member variables
	unsigned int getLocNum();
	TreeNode * getA();
	TreeNode * getC();
	TreeNode * getG();
	TreeNode * getT();

	//these are virtual functions
	virtual unsigned int getChrNum() {return (unsigned int) -1;};
	virtual unsigned int getLocation() {return (unsigned int) -1;};

	//check if node is leaf
	virtual bool isLeaf() {return false;};

private:
	unsigned int locNum;

	//Ignore, testing purpose only
	unsigned int chrNum;
	unsigned int location;

	TreeNode* A;
	TreeNode* C;
	TreeNode* G;
	TreeNode* T;
};

class LeafNode : public TreeNode {
public:
	LeafNode();
	//*** Need to be implemented !!!!//
	LeafNode(unsigned int locNum);
	LeafNode(unsigned int chrNum, unsigned int location);
	//Copy constructor. Copies the entire children branch
	LeafNode(const LeafNode& rhs);

	~LeafNode();

	virtual void setContent(unsigned int chrNum, unsigned int location);

	unsigned int getChrNum();
	unsigned int getLocation();

	bool isLeaf() {return true;}

private:
	unsigned int chrNum;
	unsigned int location;

};

class HashTree {
public:
	HashTree(unsigned int hashLength = 10, unsigned int fullLength = 30);
	~HashTree();

	void generateTree(string refName);

	//Storing and loading the tree from fs; Richard and Sunny
	void storeTree(string dbName);
	void loadTree(string dbName);

	//Get data base property
	unsigned int getHashLength();
	unsigned int getFullLength();

	//Access the database
	void query(string DNA, unsigned int threshold = 1);

	//Functions to be called after query
	unsigned int frequency();
	unsigned int effectiveLength();
	bool isLeaf();

	//Estimate length with target frequency
	unsigned int estimateFreqLength(unsigned int hashVal, unsigned int freqThreshold/*, string DNA*/);
	unsigned int averageFreqDepth(unsigned int hash, unsigned int depth);

	// Calculate the hashes vs. lengths table
	void calculateLenToFreqTable();

	void printHashLenTable(string filename);

	void storeHashLenTable(string filename);

	void readHashLenTable(string filename);

	void calculateFreqMap();

	void calculatePrefixSumFreqMap();

	int initTreeGetRoot();

	int getLeft();

	int getRight();

	int searchCDFFreqMap(double val);

	//Testing function
	void test(string DNA);

	//Setting variable
	void setHashLength(unsigned int hashLength) {this->hashLength = hashLength;};

private:
	//Adding an DNA piece to the tree. Creating and incrementing the location counting along the path
	void addEntry(const string& DNA, unsigned int chrNum, unsigned int location);
	void checkDouble(string& DNA, unsigned int chrNum, unsigned int location);

	void storeTreeHelp(TreeNode * toOut, ofstream &f);
	void readBinaryTree(TreeNode *&p, ifstream &fin);

	//Estimate frequency
	unsigned long long estimateFreqLengthHelper(TreeNode* node, unsigned int freqThreshold, unsigned long long freqSum, unsigned int curLength/*, string DNA*/);
	unsigned long long averageFreqDepthHelp(TreeNode* node, unsigned int depth, unsigned long long weightedSum, unsigned int currDepth/*, string DNA*/);

	//Data
	unsigned int hashLength;
	unsigned int fullLength;

	map<unsigned int, unsigned long long> freqMap;
	map<unsigned int, unsigned long long> prefixSumFreqMap;
	map<unsigned int, double> cdfFreqMap;

	double start_pos;
	double end_pos;

	RefDB refdb;
	RefDB* checkdb;
	KmerHash hashgen;
	TreeNode** root;

	//Table
	vector< vector<unsigned int> > hashLenTable;

	//For query
	unsigned int hashVal;
	unsigned int length;
	TreeNode* queryCursor;
	bool hitLeaf;
};

#endif
