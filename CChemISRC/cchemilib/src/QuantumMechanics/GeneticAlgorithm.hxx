#ifndef GENETICALGHORITH
#define GENETICALGHORITH

#include "Rand.h"
#include "SEField.h"

class GeneticAlgorithm
{
/*
 *		maxGens: Maximum number of generations to evaluate.
 *		popSize: Size of the population
 *		pCross: Probability of cross between individuals.
 *		pMutat: Probability of mutation in an idividual.
 *		low: Vector that containts the lower bounds of the solution space.
 *		up: Vector that containts the upper bounds of the solution space.
 */
private:
	SEField* seField;
	double objFunction(vector<double> vPoints, int proc);
	int maxGens;
	int popSize;
	double pCross;
	double pMutat;
	vector<double> low;
	vector<double> up;
	vector<double> lowHistogram;
	vector<double> upHistogram;
	vector< vector<double> > histogram;
	vector< vector< vector< double> > > histogram2D;
	void resetLowUp();
	void initLowUpHistogram();
	void initHistogram();
	void updateHistogram();
	void saveHistogram();
	void initHistogram2D();
	void updateHistogram2D();
	void saveHistogram2D();
	void saveParameters();
	void initGAParameters();
	void evaluatePopulation();
	void initStat();
	void updateStat();
	void finalizeStat();
	void printStat();
	void elitist();
	void applyMutation(int itr);
	void makeSelection();

	vector< vector<double> > x;
	vector<double> fun;
	vector<double> FitnessFun;
	vector<double> xmin;
	vector<double> xmean;
	vector<double> sigmax;
	double fmin;
	int sizeHistogram;
public:
	GeneticAlgorithm(int sizeHistogram, int maxG, int popS, double pC, double pM,  SEField* seField);
	virtual ~GeneticAlgorithm();
	double run();
};

#endif
