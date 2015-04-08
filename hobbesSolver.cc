#include "hobbesSolver.h"
#include <cstring>
#include <cassert>
#include <iostream>
#include <climits>

HobbesSolver::HobbesSolver() {
	invertedList = NULL;
	defaultInvertedList = NULL;
	dynamicMatrix = NULL;
	defaultDynamicMatrix = NULL;
	readLength = 0;
	seedNum = 0;
	seedLength = 0;
}

HobbesSolver::~HobbesSolver() {
	if  (invertedList != NULL) {
		delete [] invertedList;
		invertedList = NULL;
	}
	if  (defaultInvertedList != NULL) {
		delete [] defaultInvertedList;
		defaultInvertedList = NULL;
	}
		if  (dynamicMatrix != NULL) {
		delete [] dynamicMatrix;
		dynamicMatrix = NULL;
	}
	if  (defaultDynamicMatrix != NULL) {
		delete [] defaultDynamicMatrix;
		defaultDynamicMatrix = NULL;
	}
}

void HobbesSolver::loadTree(string treeFileName) {
	tree.loadTree(treeFileName);
}

void HobbesSolver::generateTree(string refFileName) {
	tree.generateTree(refFileName);
}

void HobbesSolver::init(int readLength, int seedNum, int seedLength) {
	//Initialize the matrix
	if  (invertedList != NULL) {
		delete [] invertedList;
		invertedList = NULL;
	}
	if  (defaultInvertedList != NULL) {
		delete [] defaultInvertedList;
		defaultInvertedList = NULL;
	}
		if  (dynamicMatrix != NULL) {
		delete [] dynamicMatrix;
		dynamicMatrix = NULL;
	}
	if  (defaultDynamicMatrix != NULL) {
		delete [] defaultDynamicMatrix;
		defaultDynamicMatrix = NULL;
	}

	this->readLength = readLength;
	this->seedNum = seedNum;
	this->seedLength = seedLength;
	dynamicMatrixWidth = readLength + 1 - seedNum * seedLength;

	invertedList = new int [readLength - seedLength + 1];
	defaultInvertedList = new int [readLength - seedLength + 1];

	dynamicMatrix = new int* [seedNum + 1];
	dynamicMatrix[0] = new int [dynamicMatrixWidth * (seedNum + 1)];
	defaultDynamicMatrix = new int [dynamicMatrixWidth * (seedNum + 1)];

	for (int i = 1; i <= seedNum; i++)
		dynamicMatrix[i] = dynamicMatrix[i-1] + dynamicMatrixWidth;

	for (int i = 0; i < readLength - seedLength + 1; i++)
		defaultInvertedList[i] = 0;
	
	for (int i = 0; i < dynamicMatrixWidth * (seedNum + 1); i++)
		defaultDynamicMatrix[i] = 0;

}

void HobbesSolver::reset() {
	memcpy(dynamicMatrix[0], defaultDynamicMatrix, dynamicMatrixWidth * ( this->seedNum + 1));
	memcpy(invertedList, defaultInvertedList, readLength - this->seedLength + 1);
}

unsigned int HobbesSolver::solveDNA(string DNA) {
	assert(DNA.length() == readLength);
	reset();

	for (int i = 0; i < readLength - seedLength + 1; i++)
	{
		string seed = DNA.substr(i, seedLength);
		tree.query(seed);
		invertedList[i] = tree.frequency();
		cout << invertedList[i];
		cout << " ";
	}

	cout << "\n";

	for (int counter = 1; counter < seedNum + 1; counter++)
	{
		dynamicMatrix[counter][0] = invertedList[(counter-1) * this->seedLength] + dynamicMatrix[counter - 1][0];
		//cout << dynamicMatrix[counter][0];
		//cout << " ";
	}
	cout << "\n";
	cout << "\n";
	for (int heightCount = 1; heightCount < seedNum + 1; heightCount++)
	{
		for (int widthCount = 1; widthCount < dynamicMatrixWidth; widthCount++)
		{
			int indexNum = (heightCount - 1) * seedLength + widthCount;
			dynamicMatrix[heightCount][widthCount] = dynamicMatrix[heightCount-1][widthCount] + invertedList[indexNum];
			if (dynamicMatrix[heightCount][widthCount-1] < dynamicMatrix[heightCount][widthCount])
				dynamicMatrix[heightCount][widthCount] = dynamicMatrix[heightCount][widthCount-1];
			cout << dynamicMatrix[heightCount][widthCount-1];
			cout << " ";			
		}
		cout << dynamicMatrix[heightCount][dynamicMatrixWidth-1];
		cout << "\n";
	}
	cout << "\n";
	cout << dynamicMatrix[seedNum][dynamicMatrixWidth-1];
	cout << "\n";
	return dynamicMatrix[seedNum][dynamicMatrixWidth-1];
}

