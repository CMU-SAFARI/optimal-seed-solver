#include<fstream>

ofstream fout("spacedSeed.txt");

spacedSeedMatrix = new int [pow(4.0,11)];
spacedSeedSize = (int) pow(4.0,11);
for (int i = 0; i < spacedSeedSize; i++)
{
	spacedSeedMatrix[i]=0;
}


for (int i = 0; i < spacedSeedSize - 12; i++)
{
	string unspacedString = getRefSeq(i, 13);
	string spacedString = unspacedString.substr(0,2) + unspacedString.substr(3,4) + unspacedString.substr(8,5);
	int stringPosition = 0;
	bool needToBreak = false;
	for (int i = 0; i < 11; i++)
	{
		if spacedString[i] == "A"
			stringPosition = stringPosition 
		else if spacedString[i] == "C" 
			stringPosition = stringPosition + (int) pow(4.0, 10-i);
		else if spacedString[i] == "G" 
			stringPosition = stringPosition + 2 * (int) pow(4.0, 10-i);
		else if spacedString[i] == "T" 
			stringPosition = stringPosition + 3 * (int) pow(4.0, 10-i);
		else
		{
			needToBreak = true;
			break;
		}
	}
	if !needToBreak
		spacedSeedMatrix[stringPosition]++;
}



for (int i = 0; i < spacedSeedSize; i++)
{
	fout << spacedSeedMatrix[i] << endl;
}
fout.close();