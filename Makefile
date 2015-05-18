XX = g++
CPPFLAGS=-g -DDEBUG
EXE=getHash getReference readGenerator testHashTree testTemplateClass testThresholdSolver \
    testOptimalSolver testPredictorSolver testHobbesSolver testReference testFastHASHSolver \
    testGemsSolver testSpacedSeedSolver testBasicSolver testAnalyzeBasic testOptimalSolverLN

all: $(EXE)

KmerHash.o: KmerHash.cc KmerHash.h
	$(CXX) $(CPPFLAGS) -c $<
	
RefDB.o: RefDB.cc RefDB.h
	$(CXX) $(CPPFLAGS) -c $<

HashTree.o: HashTree.cc HashTree.h
	$(CXX) $(CPPFLAGS) -c $<

templateClass.o: templateClass.cc templateClass.h
	$(CXX) $(CPPFLAGS) -c $<

thresholdSolver.o: thresholdSolver.cc thresholdSolver.h 
	$(CXX) $(CPPFLAGS) -c $<
	
optimalSolver.o: optimalSolver.cc optimalSolver.h 
	$(CXX) $(CPPFLAGS) -c $<

optimalSolverLN.o: optimalSolverLN.cc optimalSolverLN.h 
	$(CXX) $(CPPFLAGS) -c $<

predictorSolver.o: predictorSolver.cc predictorSolver.h 
	$(CXX) $(CPPFLAGS) -c $<

hobbesSolver.o: hobbesSolver.cc hobbesSolver.h 
	$(CXX) $(CPPFLAGS) -c $<

fastHASHSolver.o: fastHASHSolver.cc fastHASHSolver.h
	$(CXX) $(CPPFLAGS) -c $<

spacedSeedSolver.o: spacedSeedSolver.cc spacedSeedSolver.h
	$(CXX) $(CPPFLAGS) -c $<

analyzeBasic.o: analyzeBasic.cc analyzeBasic.h
	$(CXX) $(CPPFLAGS) -c $<

basicSolver.o: basicSolver.cc basicSolver.h
	$(CXX) $(CPPFLAGS) -c $<

getHash: KmerHashMain.cc KmerHash.o
	$(CXX) $(CPPFLAGS) -o $@ $^
	
getReference: getReference.cc RefDB.o
	$(CXX) $(CPPFLAGS) -o $@ $^
	
readGenerator: readGenerator.cc RefDB.o
	$(CXX) $(CPPFLAGS) -o $@ $^

testReference: testReference.cc RefDB.o
	$(CXX) $(CPPFLAGS) -o $@ $^

testHashTree: testHashTree.cc HashTree.o RefDB.o KmerHash.o
	$(CXX) $(CPPFLAGS) -o $@ $^
	
testTemplateClass: testTemplateClass.cc templateClass.o HashTree.o RefDB.o KmerHash.o
	$(CXX) $(CPPFLAGS) -o $@ $^

testThresholdSolver: testThresholdSolver.cc thresholdSolver.o HashTree.o RefDB.o KmerHash.o
	$(CXX) $(CPPFLAGS) -o $@ $^

testOptimalSolver: testOptimalSolver.cc optimalSolver.o HashTree.o RefDB.o KmerHash.o
	$(CXX) $(CPPFLAGS) -o $@ $^

testOptimalSolverLN: testOptimalSolverLN.cc optimalSolverLN.o HashTree.o RefDB.o KmerHash.o
	$(CXX) $(CPPFLAGS) -o $@ $^

testPredictorSolver: testPredictorSolver.cc predictorSolver.o HashTree.o RefDB.o KmerHash.o
	$(CXX) $(CPPFLAGS) -o $@ $^

testHobbesSolver: testHobbesSolver.cc hobbesSolver.o HashTree.o RefDB.o KmerHash.o
	$(CXX) $(CPPFLAGS) -o $@ $^

testFastHASHSolver: testFastHASHSolver.cc fastHASHSolver.o HashTree.o RefDB.o KmerHash.o
	$(CXX) $(CPPFLAGS) -o $@ $^

testSpacedSeedSolver: testSpacedSeedSolver.cc spacedSeedSolver.o HashTree.o RefDB.o KmerHash.o

	$(CXX) $(CPPFLAGS) -o $@ $^

testGemsSolver: testGemsSolver.cc fastHASHSolver.o thresholdSolver.o HashTree.o RefDB.o KmerHash.o
	$(CXX) $(CPPFLAGS) -o $@ $^

testBasicSolver: testBasicSolver.cc basicSolver.o HashTree.o RefDB.o KmerHash.o
	$(CXX) $(CPPFLAGS) -o $@ $^

testAnalyzeBasic: testAnalyzeBasic.cc analyzeBasic.o RefDB.o KmerHash.o
	$(CXX) $(CPPFLAGS) -o $@ $^

.PHONY: clean copy

clean:
	rm *.o $(EXE)

copy: run*.py HashTable RepeatAnalyser SeedTree getHash getReference readGenerator mapSimulator LongHashTable
	cp $^ ~/data/hashtable/
testPredictorSolver: testPredictorSolver.cc predictorSolver.o HashTree.o RefDB.o KmerHash.o
	$(CXX) $(CPPFLAGS) -o $@ $^

testHobbesSolver: testHobbesSolver.cc hobbesSolver.o HashTree.o RefDB.o KmerHash.o
	$(CXX) $(CPPFLAGS) -o $@ $^

testFastHASHSolver: testFastHASHSolver.cc fastHASHSolver.o HashTree.o RefDB.o KmerHash.o
	$(CXX) $(CPPFLAGS) -o $@ $^

testSpacedSeedSolver: testSpacedSeedSolver.cc spacedSeedSolver.o HashTree.o RefDB.o KmerHash.o

	$(CXX) $(CPPFLAGS) -o $@ $^

testGemsSolver: testGemsSolver.cc fastHASHSolver.o thresholdSolver.o HashTree.o RefDB.o KmerHash.o
	$(CXX) $(CPPFLAGS) -o $@ $^

testBasicSolver: testBasicSolver.cc basicSolver.o HashTree.o RefDB.o KmerHash.o
	$(CXX) $(CPPFLAGS) -o $@ $^

testAnalyzeBasic: testAnalyzeBasic.cc analyzeBasic.o RefDB.o KmerHash.o
	$(CXX) $(CPPFLAGS) -o $@ $^

.PHONY: clean copy

clean:
	rm *.o $(EXE)

copy: run*.py HashTable RepeatAnalyser SeedTree getHash getReference readGenerator mapSimulator LongHashTable
	cp $^ ~/data/hashtable/
