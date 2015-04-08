#include "HashTree.h"

TreeNode::TreeNode() {
	locNum = 0;
	A = NULL;
	C = NULL;
	G = NULL;
	T = NULL;
}

TreeNode::TreeNode(unsigned int locNum) {
	this->locNum = locNum;
	A = NULL;
	C = NULL;
	G = NULL;
	T = NULL;
}

TreeNode::TreeNode(const TreeNode& rhs) {
	if (A != NULL) {
		delete A;
		A = NULL;
	}
	if (C != NULL) {
		delete C;
		C = NULL;
	}
	if (G != NULL) {
		delete G;
		G = NULL;
	}
	if (T != NULL) {
		delete T;
		T = NULL;
	}

	if (rhs.A != NULL)
		A = new TreeNode(*rhs.A);
	if (rhs.C != NULL)
		C = new TreeNode(*rhs.C);
	if (rhs.C != NULL)
		G = new TreeNode(*rhs.G);
	if (rhs.T != NULL)
		T = new TreeNode(*rhs.T);

	this->locNum = rhs.locNum;
}

TreeNode::~TreeNode() {
	locNum = 0;
	deleteChildren();
}

void TreeNode::setLocNum(unsigned int locNum) {
	this->locNum = locNum;
}

void TreeNode::addLoc() {
	this->locNum++;
}

void TreeNode::setA(TreeNode* A) {
	this->A = A;
}

void TreeNode::setC(TreeNode* C) {
	this->C = C;
}

void TreeNode::setG(TreeNode* G) {
	this->G = G;
}

void TreeNode::setT(TreeNode* T) {
	this->T = T;
}

void TreeNode::deleteChildren() {

	if (A != NULL) {
		delete A;
		A = NULL;
	}

	if (C != NULL) {
		delete C;
		C = NULL;
	}

	if (G != NULL) {
		delete G;
		G = NULL;
	}

	if (T != NULL) {
		delete T;
		T = NULL;
	}
}

void LeafNode::setContent(unsigned int chrNum, unsigned int location) {
	this->chrNum = chrNum;
	this->location = location;
}

unsigned int TreeNode::getLocNum() {
	return locNum;
}

TreeNode * TreeNode::getA() {
	return A;
}

TreeNode * TreeNode::getC() {
	return C;
}

TreeNode * TreeNode::getG() {
	return G;
}

TreeNode * TreeNode::getT() {
	return T;
}

LeafNode::LeafNode() : TreeNode() {
	chrNum = 0;
	location = 0;
}

LeafNode::LeafNode(unsigned int locNum) : TreeNode(locNum) {
}

LeafNode::LeafNode(unsigned int chrNum, unsigned int location) : TreeNode() {
	this->chrNum = chrNum;
	this->location = location;
}

LeafNode::LeafNode(const LeafNode &rhs) : TreeNode(rhs) {
	this->chrNum = rhs.chrNum;
	this->location = rhs.location;
}

LeafNode::~LeafNode() {
	chrNum = 0;
	location = 0;
}

unsigned int LeafNode::getChrNum() {
	return chrNum;
}

unsigned int LeafNode::getLocation() {
	return location;
}

HashTree::HashTree(unsigned int hashLength, unsigned int fullLength) {
	this->hashLength = hashLength;
	this->fullLength = fullLength;

	root = new TreeNode*[1 << (2 * hashLength)];

	for (int i = 0; i < 1 << (2 * hashLength); i++)
		root[i] = NULL;

	checkdb = NULL;
}

HashTree::~HashTree() {
	for (int i = 0; i < 1 << (2 * hashLength); i++) {
		if (root[i] != NULL) {
			delete root[i];
			root[i] = NULL;
		}
	}

	if (checkdb != NULL) {
		delete[] checkdb;
		checkdb = NULL;
	}

	if (root != NULL) {
		delete[] root;
		root = NULL;
	}
}

void HashTree::addEntry(const string& DNA, unsigned int chrNum, unsigned int location) {
	unsigned int hashVal = hashgen.calculateHash(DNA, hashLength);

	//cout << "DNA: " << DNA << endl;
	//cout << "hashVal: " << hashVal << endl;

	TreeNode* cursor = root[hashVal];
	TreeNode* parentCursor = NULL;
	TreeNode* nextCursor = NULL;

	//cout << "***location: " << location << endl;

	/*
	   if (DNA == "CCTGGTGCTCCCACAAAGGAGAAGGGCTGA") {
	   cout << "found something" << endl;
	   cout << "hashVal: " << hashVal << endl;
	   }
	 */

	//The depth depends on the length of the string: DNA
	for (int i = hashLength; i < DNA.length(); i++) {
		//cout << i << ' ' << DNA[i] << ' ';
		//Create a root node if there isn't one
		if (cursor == NULL) {
			cursor = new LeafNode();
			root[hashVal] = cursor;
			cursor->setContent(chrNum, location);
			//Initialize the counter
			cursor->addLoc();
			return;
		}

		//cout << "i: " << i << " of DNA.length() - 1 = " << DNA.length() - 1 << endl;

		//This is a leaf node but it should no longer be one
		if (cursor->isLeaf() ) {
			//cout << "i: " << i << " Replacing!" << endl;
			// Creating a trunk
			TreeNode* newTrunk = new TreeNode(1);
			//cout << "cursor->getLocNum: " << cursor->getLocNum() << endl;
			//cout.flush();
			assert(cursor->getLocNum() == 1);

			assert(cursor->getChrNum() < refdb.getNumOfChromo() );

			//Push the LeafNode 1 level lower
			switch (checkdb[cursor->getChrNum()].getRefBase(cursor->getLocation(), i) ) {
				case 'A':
					newTrunk->setA(new LeafNode(cursor->getChrNum(), cursor->getLocation() ) );
					newTrunk->getA()->addLoc();
					break;
				case 'C':
					newTrunk->setC(new LeafNode(cursor->getChrNum(), cursor->getLocation() ) );
					newTrunk->getC()->addLoc();
					break;
				case 'G':
					newTrunk->setG(new LeafNode(cursor->getChrNum(), cursor->getLocation() ) );
					newTrunk->getG()->addLoc();
					break;
				case 'T':
					newTrunk->setT(new LeafNode(cursor->getChrNum(), cursor->getLocation() ) );
					newTrunk->getT()->addLoc();
					break;
				case 'M':
					newTrunk->setA(new LeafNode(cursor->getChrNum(), cursor->getLocation() ) );
					newTrunk->getA()->addLoc();
					newTrunk->setC(new LeafNode(cursor->getChrNum(), cursor->getLocation() ) );
					newTrunk->getC()->addLoc();
					break;
				case 'R':
					newTrunk->setA(new LeafNode(cursor->getChrNum(), cursor->getLocation() ) );
					newTrunk->getA()->addLoc();
					newTrunk->setG(new LeafNode(cursor->getChrNum(), cursor->getLocation() ) );
					newTrunk->getG()->addLoc();
					break;
				default:
					cerr << "unrecognized character: " << DNA[i] << endl;
					break;
			}

			//Check if we need to update root or a trunk
			if (parentCursor == NULL)
				root[hashVal] = newTrunk;
			else {
				switch (DNA[i-1]) {
					case 'A':
						parentCursor->setA(newTrunk);
						break;
					case 'C':
						parentCursor->setC(newTrunk);
						break;
					case 'G':
						parentCursor->setG(newTrunk);
						break;
					case 'T':
						parentCursor->setT(newTrunk);
						break;
				}
			}

			delete cursor;
			cursor = newTrunk;
			newTrunk = NULL;
		}

		//Incrementing the counter along the path
		cursor->addLoc();

		switch (DNA[i]) {
			case 'A':
				nextCursor = cursor->getA();
				break;
			case 'C':
				nextCursor = cursor->getC();
				break;
			case 'G':
				nextCursor = cursor->getG();
				break;
			case 'T':
				nextCursor = cursor->getT();
				break;
			default:
				cerr << "unrecognized character: " << DNA[i] << endl;
				break;
		}

		if (nextCursor == NULL) {
			//cout << "WTF!!!!!" << "DNA[" << i << "]: " << DNA[i] << endl;

			nextCursor = new LeafNode(chrNum, location);
			//Incrementing the counter of the new node
			nextCursor->addLoc();

			switch (DNA[i]) {
				case 'A':
					cursor->setA(nextCursor);
					break;
				case 'C':
					cursor->setC(nextCursor);
					break;
				case 'G':
					cursor->setG(nextCursor);
					break;
				case 'T':
					cursor->setT(nextCursor);
					break;
				default:
					cerr << "unrecognized character: " << DNA[i] << endl;
					break;
			}

			return;
		}

		//if this is the last node, we have to increment the counter
		if (i == DNA.length() - 1)
			nextCursor->addLoc();

		parentCursor = cursor;
		cursor = nextCursor;
	}
}

void HashTree::checkDouble(string& DNA, unsigned int chrNum, unsigned int location) {
	bool leaf = true;

	if (DNA.find('M') != string::npos) {
		int pos = DNA.find('M');
		DNA[pos] = 'A';
		checkDouble(DNA, chrNum, location);
		DNA[pos] = 'C';
		checkDouble(DNA, chrNum, location);
		DNA[pos] = 'M';

		leaf = false;
	}

	if (DNA.find('R') != string::npos) {
		int pos = DNA.find('R');
		DNA[pos] = 'A';
		checkDouble(DNA, chrNum, location);
		DNA[pos] = 'G';
		checkDouble(DNA, chrNum, location);
		DNA[pos] = 'R';

		leaf = false;
	}

	//Add DNA to database unless we have resolved all cases
	if (leaf)
		addEntry(DNA, chrNum, location);

}

void HashTree::generateTree(string refName) {
	refdb.loadRefFile(refName);

	//Load all databases. This only works for 1TB machine.
	checkdb = new RefDB[refdb.getNumOfChromo()];
	for (int i = 0; i < refdb.getNumOfChromo(); i++) {
		checkdb[i].loadRefFile(refName);
		checkdb[i].loadChromo(i);
	}

	string DNA;

	for (unsigned int chromoIter = 0; chromoIter < refdb.getNumOfChromo(); chromoIter++) {
		refdb.loadChromo(chromoIter);

		for (unsigned int location = 0; location < refdb.getChromoLength() - hashLength; location++){
			DNA = refdb.getRefSeq(location, fullLength);
			if (DNA.find('N') == string::npos)
				checkDouble(DNA, chromoIter, location);
		}
	}
}

unsigned int HashTree::getHashLength() {
	return hashLength;
}

unsigned int HashTree::getFullLength() {
	return fullLength;
}

void HashTree::query(string DNA, unsigned int threshold) {
	hashVal = hashgen.calculateHash(DNA, hashLength);
	queryCursor = root[hashVal];
	length = hashLength;

	TreeNode* nextCursor;

	hitLeaf = false;

	if (queryCursor == NULL || queryCursor->isLeaf() || queryCursor->getLocNum() <= threshold ) {
		hitLeaf = true;
		return;
	}

	while (length < DNA.length() ) {


		switch (DNA[length]) {
			case 'A':
				nextCursor = queryCursor->getA();
				break;
			case 'C':
				nextCursor = queryCursor->getC();
				break;
			case 'G':
				nextCursor = queryCursor->getG();
				break;
			case 'T':
				nextCursor = queryCursor->getT();
				break;
			default:
				cerr << "unrecognized character: " << DNA[length] << endl;
				break;
		}

		queryCursor = nextCursor;
		length++;

		if (queryCursor == NULL || queryCursor->isLeaf() || queryCursor->getLocNum() <= threshold ) {
			hitLeaf = true;
			return;
		}
	}
}

unsigned int HashTree::effectiveLength() {
	//Always return length. In query, if no root seed is found it will be
	//initialized to hashLength
	return length;
}

unsigned int HashTree::frequency() {
	if (queryCursor == NULL)
		return 0;
	else
		return queryCursor->getLocNum();
}

bool HashTree::isLeaf() {
	return hitLeaf;
}

unsigned int HashTree::estimateFreqLength(unsigned int hash, unsigned int freqThreshold/*, string DNA*/) {
	assert(hash < (1 << (2 * hashLength) ) );

	TreeNode* base = root[hash];

	if (base == NULL || base->getLocNum() <= freqThreshold)
		return hashLength;

	unsigned long long freqSum = 0;
	freqSum = estimateFreqLengthHelper(base, freqThreshold, freqSum, hashLength/*, DNA*/);

	return (freqSum - 1) / base->getLocNum() + 1;
}

unsigned long long HashTree::estimateFreqLengthHelper(TreeNode* node, unsigned int freqThreshold, unsigned long long freqSum, unsigned int curLength/*, string DNA*/) {
	if (node == NULL) {
		//cout << DNA << " curLength: " << curLength << " NULL." << endl;
		return freqSum;
	}

	if (node->getLocNum() <= freqThreshold || node->isLeaf() ) {
		//cout << DNA << " curLength: " << curLength << " frequency: " << node->getLocNum() << endl;
		return freqSum + (unsigned long long) curLength * node->getLocNum();
	}

	freqSum = estimateFreqLengthHelper(node->getA(), freqThreshold, freqSum, curLength + 1/*, DNA + 'A'*/);
	freqSum = estimateFreqLengthHelper(node->getC(), freqThreshold, freqSum, curLength + 1/*, DNA + 'C'*/);
	freqSum = estimateFreqLengthHelper(node->getG(), freqThreshold, freqSum, curLength + 1/*, DNA + 'G'*/);
	freqSum = estimateFreqLengthHelper(node->getT(), freqThreshold, freqSum, curLength + 1/*, DNA + 'T'*/);

	return freqSum;
}

unsigned int HashTree::averageFreqDepth(unsigned int hash, unsigned int depth)
{
	assert(hash < (1 << (2 * hashLength) ) );

	TreeNode* base = root[hash];

	if (base == NULL)
		return 0;

	unsigned long long weightedSum = 0;
	weightedSum = averageFreqDepthHelp(base, depth, weightedSum, 0);

	return weightedSum / base->getLocNum();
}

unsigned long long HashTree::averageFreqDepthHelp(TreeNode* node, unsigned int depth, unsigned long long weightedSum, unsigned int currDepth/*, string DNA*/) {
	if (node == NULL) {
		return weightedSum;
	}

	if(currDepth == depth || (node->isLeaf() && currDepth < depth))
	{
		unsigned long long freq =  node->getLocNum();
		return freq * freq + weightedSum;
	}
	else
	{
		weightedSum = averageFreqDepthHelp(node->getA(), depth, weightedSum, currDepth + 1/*, DNA + 'A'*/);
		weightedSum = averageFreqDepthHelp(node->getC(), depth, weightedSum, currDepth + 1/*, DNA + 'C'*/);
		weightedSum = averageFreqDepthHelp(node->getG(), depth, weightedSum, currDepth + 1/*, DNA + 'G'*/);
		weightedSum = averageFreqDepthHelp(node->getT(), depth, weightedSum, currDepth + 1/*, DNA + 'T'*/);
	}

	return weightedSum;
}

void HashTree::calculateLenToFreqTable()
{
	unsigned int tableWidth = fullLength - hashLength + 1;
	unsigned int tableHeight = 1 << (2 * hashLength);

	vector< vector<unsigned int> > hashLenTable (tableHeight, vector<unsigned>(tableWidth));

	for(int i=0; i<tableHeight; i++)
		for(int j=0; j<tableWidth; j++)
			hashLenTable[i][j] = averageFreqDepth(i, j);

	this->hashLenTable = hashLenTable;
}

void HashTree::printHashLenTable(string filename)
{
	std::ofstream f;
	f.open(filename.c_str(), ios::out | ios::binary);

	unsigned int tableWidth = fullLength - hashLength + 1;
	unsigned int tableHeight = 1 << (2 * hashLength);
	if(&hashLenTable != NULL)
	{
		for(int i=0; i<tableHeight; i++)
		{
			f << "Index: " << i << " ";
			for(int j=0; j<tableWidth; j++)
			{
				f << hashLenTable[i][j] << "\t";
			}
			f << "\n";
		}
	}
	else
		cout << "Empty table. " << endl;

	f.close();
}

void HashTree::storeHashLenTable(string filename)
{
	unsigned int tableWidth = fullLength - hashLength + 1;
	unsigned int tableHeight = 1 << (2 * hashLength);

	if(&hashLenTable == NULL)
		cout << "Empty table. " << endl;

	std::ofstream f;
	f.open(filename.c_str(), ios::out | ios::binary);

	int hl[] = {(int)(fullLength - hashLength + 1)};
	int fl[] = {(int)(1 << (2 * hashLength))};

	f.write((char*)&hl, sizeof(int));
	f.write((char*)&fl, sizeof(int));

	for(int i=0; i<tableHeight; i++)
		for(int j=0; j<tableWidth; j++)
		{
			int h[] = {(int)(hashLenTable[i][j])};
			f.write((char*)&h, sizeof(int));
		}
	f.close();
}

void HashTree::readHashLenTable(string filename)
{
	std::ifstream f;
	f.open(filename.c_str(),  ios::in | ios::binary);
	if(f.is_open())
	{
		f.seekg(0, ios::beg);

		int curr[1];
		f.read((char*)(&curr), sizeof(int));
		unsigned int tableWidth = (unsigned int)(curr[0]);

		int curr2[1];
		f.read((char*)(&curr2), sizeof(int));
		unsigned int tableHeight = (unsigned int)(curr2[0]);

		vector< vector <unsigned int> > hashLenTable (tableHeight,
			vector<unsigned int>(tableWidth));

		for(int i=0; i<tableHeight; i++)
			for(int j=0; j<tableWidth; j++)
			{
				int curr[1];
				f.read((char*)(&curr), sizeof(int));
				hashLenTable[i][j] = (unsigned int)(curr[0]);
			}

		this->hashLenTable = hashLenTable;
		f.close();
		return;
	}
}


void HashTree::calculateFreqMap()
{
	/*
	 * Map which stores <freq, numfreq > pairs
	 */
	map <unsigned int, unsigned long long> freqMap;
	unsigned int cfreq;
	freqMap.clear();
	if(&hashLenTable != NULL)
	{
		unsigned int counter = 0;
		/* Iterate through the table */
		for(int i = 0; i < 1 << (2 * hashLength); i++)
		{
			for(int j = 0; j < (fullLength - hashLength + 1); j++)
			{
				counter++;

				cfreq = hashLenTable[i][j];
				/* Increase counter */
				if (freqMap.find(cfreq) == freqMap.end())
					freqMap[cfreq] = 1;
				else
					freqMap[cfreq]++;

				/* We only want the first element with frequency 1
				   since we are guaranteed the rest are 1. */
				if(cfreq == 1)
					break;
			}
		}
		this->freqMap = freqMap;
	}
}

void HashTree::calculatePrefixSumFreqMap()
{
	map <unsigned int, unsigned long long> prefixSumFreqMap;
	prefixSumFreqMap.clear();
	int total = 0;
	if(&freqMap != NULL)
	{
		/* Add first element */
		map<unsigned int, unsigned long long>::iterator i = freqMap.begin();
		prefixSumFreqMap[i->first] = freqMap[i->first];
		i++;
		map<unsigned int, unsigned long long>::iterator prev = freqMap.begin();
		/* Loop for prefix sum */
		for (; i!=freqMap.end(); i++, prev++)
		{
			prefixSumFreqMap[i->first] = prefixSumFreqMap[prev->first] +
				freqMap[i->first];
			total = prefixSumFreqMap[i->first];
		}

		this->prefixSumFreqMap = prefixSumFreqMap;

		/* Calculate cumulative distribution of frequencies */
		map <unsigned int, double> cdfFreqMap;

		map<unsigned int, unsigned long long>::iterator iter = prefixSumFreqMap.begin();
		while(iter != prefixSumFreqMap.end())
		{
			cdfFreqMap[iter->first] = prefixSumFreqMap[iter->first]/total;
			iter++;
		}
		this->cdfFreqMap = cdfFreqMap;
	}
}

int HashTree::initTreeGetRoot()
{
	if(&prefixSumFreqMap != NULL)
	{
		start_pos = 0;
		end_pos = 1;
		return searchCDFFreqMap((start_pos + end_pos)/2);
	}
}

int HashTree::getLeft()
{
	end_pos = (start_pos + end_pos)/2;
	return searchCDFFreqMap((start_pos + end_pos)/2);
}

int HashTree::getRight()
{
	start_pos = (start_pos + end_pos)/2;
	return searchCDFFreqMap((start_pos + end_pos)/2);
}

int HashTree::searchCDFFreqMap(double val)
{
	int ret_freq = -1;
	for(map<unsigned int, double>::iterator i = cdfFreqMap.begin(); i != cdfFreqMap.end(); i++)
	{
		if(val <= cdfFreqMap[i->first])
		{
			ret_freq = i->second;
		}
	}
	return ret_freq;
}

void HashTree::test(string DNA) {
	hashVal = hashgen.calculateHash(DNA, hashLength);

	queryCursor = root[hashVal];

	if (queryCursor == NULL) {
		cout << "NULL pointer detected" << endl;
		cout << "DNA: " << DNA << endl;
		cout << "hashLength: " << hashLength << endl;

	}

	TreeNode* nextCursor = queryCursor;

	int i;
	//The depth depends on the length of the string: DNA
	for (i = hashLength; i < DNA.length(); i++) {
		cout << "TreeNode: " << DNA.substr(0, i)
			<< " of frequency: " << queryCursor->getLocNum()
			<< " from chromosome: " << queryCursor->getChrNum()
			<< " at location: " << queryCursor->getLocation()
			<< " isLeaf: " << queryCursor->isLeaf() << endl;

		if (!queryCursor->isLeaf() ) {
			switch (DNA[i]) {
				case 'A':
					nextCursor = queryCursor->getA();
					break;
				case 'C':
					nextCursor = queryCursor->getC();
					break;
				case 'G':
					nextCursor = queryCursor->getG();
					break;
				case 'T':
					nextCursor = queryCursor->getT();
					break;
				default:
					cerr << "unrecognized character: " << DNA[i] << endl;
					break;
			}
		}

		queryCursor = nextCursor;

		if (queryCursor == NULL || queryCursor->isLeaf() )
			break;
	}

	cout << "LastNode: " << DNA.substr(0, i);
		if (queryCursor != NULL)
			cout << " of frequency: " << queryCursor->getLocNum()
				<< " from chromosome: " << queryCursor->getChrNum()
				<< " at location: " << queryCursor->getLocation()
				<< " isLeaf: " << queryCursor->isLeaf() << endl;
		else
			cout << " is NULL." << endl;
}

void HashTree::storeTreeHelp(TreeNode * toOut, ofstream &f)
{
	if (toOut == NULL || toOut -> getLocNum() == 0)
	{
		int c[] = {0};
		f.write((char*)(&c), sizeof(int));
	}
	else
	{
		int c[] = {(int) toOut->getLocNum()};
		f.write((char*)&c, sizeof(int));
		storeTreeHelp(toOut->getA(), f);
		storeTreeHelp(toOut->getC(), f);
		storeTreeHelp(toOut->getG(), f);
		storeTreeHelp(toOut->getT(), f);
	}
}

void HashTree::storeTree(string dbName)
{
	std::ofstream f;
	f.open(dbName.c_str(),ios::out | ios::binary);

	int hl[] = {(int) hashLength};
	int fl[] = {(int) fullLength};

	f.write((char*)&hl, sizeof(int));
	f.write((char*)&fl, sizeof(int));

	for (int toPrint=0; toPrint< 1 << (2 * hashLength); toPrint++)
	{
		storeTreeHelp(root[toPrint], f);
	}
	f.close();
}

void HashTree::readBinaryTree(TreeNode *&p, ifstream &fin) {
	//TODO: handle exceptions
	int curr[1];
	fin.read((char*)curr, sizeof(int));
	unsigned int currNum = (unsigned int)(curr[0]);
	if(currNum == 0)
		return;
	else
	{
		TreeNode *a = NULL;
		TreeNode *c = NULL;
		TreeNode *g = NULL;
		TreeNode *t = NULL;

		readBinaryTree(a, fin);
		readBinaryTree(c, fin);
		readBinaryTree(g, fin);
		readBinaryTree(t, fin);
		if (a == NULL && c == NULL && g == NULL && t == NULL) {
			p = new LeafNode(currNum);
		}
		else {
			p = new TreeNode(currNum);
			p->setA(a);
			p->setC(c);
			p->setG(g);
			p->setT(t);
		}
	}
}

void HashTree::loadTree(string dbName)
{
	std::ifstream f;
	f.open(dbName.c_str(),ios::in | ios::binary);
	if(f.is_open())
	{
		f.seekg(0, ios::beg);

		int curr[1];
		f.read((char*)(&curr), sizeof(int));
		unsigned int hashLength = (unsigned int)(curr[0]);

		int curr2[1];
		f.read((char*)(&curr2), sizeof(int));
		unsigned int fullLength = (unsigned int)(curr2[0]);

		this->hashLength = hashLength;
		this->fullLength = fullLength;

		unsigned int tbsize = 1 << (2 * hashLength);
		this->root = new TreeNode*[tbsize];

		for (int i = 0; i < tbsize; i++)
			this->root[i] = NULL;

		for(int i = 0; i < tbsize; i++)
		{
			readBinaryTree(this->root[i], f);
		}

		f.close();
	}
}
