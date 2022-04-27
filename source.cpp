#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <vector>
#include <time.h>
#include <fstream>

/*for the alimentation of code, please contatct me by mail: nostalgiadecent@gmail.com*/

//using srand, it's can make the numbers more randomly

using namespace std;

#define ChooseNumOfIncidentAngle 9 //To chose an incidence angle
#define PI 3.1415926

/*formula choice*/
const string Formula[4] = { "Ursenbach", "precise Shuey", "Shuey", "Zoeppritz" };
const int ForwardNum = 0; //chose the model of forward 
const int InverseNum = 0; //chose the model of inversion
const int DataNum = 0;	  //chose the model
/*formula choice*/

/*used for stretch fitness formula*/
const double c = 1.5;


/*incidence angle*/
const double IncidentAngle[ChooseNumOfIncidentAngle] = { 5, 10, 15, 20, 25, 30, 35, 40, 45 };
/*incidence angle*/

const int InversionTimes = 2; //time of Inversion
const int PopSize = 100;	  
const int Generation = 100;	  
const int ChromoNumber = 4;	  
const double PcFixed = 0.7;	  //Crossover probability
const double PmFixed = 0.1;	  //Mutation probability

/*use for dynamically change the range of crossover and mutation probabilities*/
const double Pcmax = 0.9;  
const double Pcmin = 0.6;  
const double Pmmax = 0.2;  
const double Pmmin = 0.05; 

/*Determine how many unconditionally will join the next population with the highest fitness each time*/
const int AmountOfBiggestPops = 1; 

/*used for Cataclysm Algorithm*/
double AverageFitnessInEveryGeneration[PopSize];
double MaxFitnessInEveryGeneration[PopSize];
int TimeOfCatastrophe = 0;
/*used for Cataclysm Algorithm*/

/*Typical values ​​of lithological parameters of 8 groups of mid-deep and deep strata in China and the Gulf of Mexico*/
const double Xmax[4] = { 0.5, 0.2, 0.5, 2.2 };
const double Xmin[4] = { 0, 0, 0, 1.8 };
const double Vp1[8] = { 2750, 2750, 2550, 2550, 4250, 4250, 3950, 3950 };
const double Vp2[8] = { 3000, 2350, 3000, 2350, 4500, 4200, 4500, 4225 };
const double Vs1[8] = { 1200, 1200, 1025, 1025, 2491, 2491, 2233, 2233 };
const double Vs2[8] = { 1550, 1547, 1550, 1545, 2650, 2633, 2650, 2645 };
const double Rou1[8] = { 2.240, 2.240, 2.203, 2.203, 2.502, 2.502, 2.458, 2.458 };
const double Rou2[8] = { 2.290, 2.160, 2.290, 2.160, 2.540, 2.496, 2.540, 2.499 };
/*Typical values ​​of lithological parameters of 8 groups of mid-deep and deep strata in China and the Gulf of Mexico*/

/*The truth value of the four data*/
const double TruthValuea = (Vp2[DataNum] - Vp1[DataNum]) / (0.5 * (Vp2[DataNum] + Vp1[DataNum]));
const double TruthValueb = (Rou2[DataNum] - Rou1[DataNum]) / (0.5 * (Rou2[DataNum] + Rou1[DataNum]));
const double TruthValuec = (Vs2[DataNum] - Vs1[DataNum]) / (0.5 * (Vs2[DataNum] + Vs1[DataNum]));
const double TruthValued = (Vp2[DataNum] + Vp1[DataNum]) / ((Vs2[DataNum] + Vs1[DataNum]));
/*The truth value of the four data*/

/* Four-parameter variation range*/
double getRanChromoSudubianhualv()
{
	double ch = 0;
	ch = 0 + (1.0 * rand() / RAND_MAX) * 0.5;
	return ch;
}

double getRanChromoMidubianhualv()
{
	double ch = 0;
	ch = 0 + (1.0 * rand() / RAND_MAX) * 0.2;
	return ch;
}

double getRanChromoSudubibianhualv1() //1.8~2.2(middle-deep)
{
	double ch = 0;
	ch = 1.8 + (1.0 * rand() / RAND_MAX) * 0.4;
	return ch;
}

double getRanChromoSudubibianhualv2() //1.5~1.9(deep)
{
	double ch = 0;
	ch = 1.5 + (1.0 * rand() / RAND_MAX) * 0.4;
	return ch;
}
/* Four-parameter variation range*/

/* Trigonometric function angle*/
double sind(double x)
{
	double a;
	a = sin(x * PI / 180);
	return a;
}

double cosd(double x)
{
	double b;
	b = cos(x * PI / 180);
	return b;
}

double tand(double x)
{
	double b;
	b = tan(x * PI / 180);
	return b;
}
/* Trigonometric function angle package*/

/* Calculate standard deviation*/
double SDCalculator(double AdditionSD) 
{
	double SD;
	SD = sqrt(AdditionSD / (InversionTimes - 1));
	return SD;
}
/* Calculate standard deviation*/

const int CrossAmount = PcFixed * PopSize;// Population size for storing crossover individuals

/*individuals*/
typedef struct individual
{
	double chromo[4];
} individual;
struct individual pop[PopSize];
struct individual nextPop[PopSize];
struct individual CopyPop[PopSize]; 
struct individual CrossPop[CrossAmount];
int GenPop = 1;			 //the first generation
int InversionRecord = 1; //the first inversion
double BiggestIndividualInEveryGeneration[InversionTimes][Generation];
double GapBetweenEveryInversion[InversionTimes][Generation];
int NextPopNum = 0;
double LineProb[PopSize];
double Fit[PopSize];
double Prob[PopSize];
/*individuals*/

double ForwardCalculator(double i1, double a, double b, double c, double d) 
//A collection of forward formulas, which formula is used for inversion according to the number entered at the beginning
{
	double Result;
	double Lin = 2 * c + b;
	double cosPhi = sqrt(1 - pow(1 / d, 2) * pow((sind(i1) / (1 - a / 2)), 2));
	double i2 = asin((Vp2[DataNum] / Vp1[DataNum]) * sind(i1)) / PI * 180; 
	double sita = (i1 + i2) / 2;										   
	double Q = (1 + a / 2) * cosd(i1) + (1 - a / 2) * cosd(i2);
	if (ForwardNum == 0)
	{
		Result = 4 * cosd(i1) * cosd(i2) / (Q * Q) * (a / (2 * cosd(i1) * cosd(i2)) - 2 * sind(i1) * sind(i2) * (1 / (d * d)) * Lin * (1 + (-1 / d) * (cosPhi * cosd(i2)) / (2 + a) + 4 / (d * d) * sind(i1) * sind(i1) / ((1 - 0.5 * a) * (1 - 0.5 * a) * Q * Q) * a + pow((1 / d) * sind(i1), 3) * sind(i2) / (2 * pow((1 - 0.5 * a), 3) * cosPhi * cosd(i2)) * Lin - (2 * a / (Q * Q) + (1 / d) * sind(i1) * sind(i2) / (2 * (1 - 0.5 * a) * cosPhi * cosd(i2)) * b)) + 0.5 * b * (1 - 0.25 * a * a) * (1 - 1 / (Q * Q) * a * b - (1 / d) * sind(i1) * sind(i2) / (2 * (1 - 0.5 * a) * cosPhi * cosd(i2)) * b)); // Ursenbach公式算出的结果
	}
	if (ForwardNum == 1)
	{
		Result = 0.5 * b - 4 / (d * d) * (c + 0.5 * b) * sind(sita) * sind(sita) + (0.5 + 0.5 * tand(sita) * tand(sita)) * a; //precise Shuey
	}

	if (ForwardNum == 2)
	{
		Result = 0.5 * b - 4 / (d * d) * (c + 0.5 * b) * sind(sita) * sind(sita) + (0.5 + 0.5 * sind(sita) * sind(sita)) * a; //simole Shuey
	}

	return Result;
}

double InverseCalculator(double i1, double a, double b, double c, double d) //Ursenbach for inversion
{

	double Result;
	double Lin = 2 * c + b;
	double cosPhi = sqrt(1 - pow(1 / d, 2) * pow((sind(i1) / (1 - a / 2)), 2));
	double i2 = asin((Vp2[DataNum] / Vp1[DataNum]) * sind(i1)) / PI * 180; 
	double sita = (i1 + i2) / 2;										   
	double Q = (1 + a / 2) * cosd(i1) + (1 - a / 2) * cosd(i2);
	if (InverseNum == 0)
	{
		Result = 4 * cosd(i1) * cosd(i2) / (Q * Q) * (a / (2 * cosd(i1) * cosd(i2)) - 2 * sind(i1) * sind(i2) * (1 / (d * d)) * Lin * (1 + (-1 / d) * (cosPhi * cosd(i2)) / (2 + a) + 4 / (d * d) * sind(i1) * sind(i1) / ((1 - 0.5 * a) * (1 - 0.5 * a) * Q * Q) * a + pow((1 / d) * sind(i1), 3) * sind(i2) / (2 * pow((1 - 0.5 * a), 3) * cosPhi * cosd(i2)) * Lin - (2 * a / (Q * Q) + (1 / d) * sind(i1) * sind(i2) / (2 * (1 - 0.5 * a) * cosPhi * cosd(i2)) * b)) + 0.5 * b * (1 - 0.25 * a * a) * (1 - 1 / (Q * Q) * a * b - (1 / d) * sind(i1) * sind(i2) / (2 * (1 - 0.5 * a) * cosPhi * cosd(i2)) * b)); // Ursenbach公式算出的结果
	}
	if (InverseNum == 1)
	{
		Result = 0.5 * b - 4 / (d * d) * (c + 0.5 * b) * sind(sita) * sind(sita) + (0.5 + 0.5 * tand(sita) * tand(sita)) * a; //precise Shuey
	}

	if (InverseNum == 2)
	{
		Result = 0.5 * b - 4 / (d * d) * (c + 0.5 * b) * sind(sita) * sind(sita) + (0.5 + 0.5 * sind(sita) * sind(sita)) * a; //Simple Shuey
	}

	return Result;
}

double ReturnMax(double a, double b)
{
	double temp;
	if (a >= b)
		temp = a;
	else
		temp = b;
	return temp;
}

double GetFitOfIndividual(individual Frog) //get the fitness of each individuals
{
	double opt = 0; 
	double cosPhi;
	double a, b, c, d;
	double Robs[ChooseNumOfIncidentAngle];
	double Rmodel[ChooseNumOfIncidentAngle];
	double NormalIndividualFit;
	a = Frog.chromo[0];
	b = Frog.chromo[1];
	c = Frog.chromo[2];
	d = Frog.chromo[3];

	for (int q = 0; q < ChooseNumOfIncidentAngle; q++)
	{
		Robs[q] = ForwardCalculator(IncidentAngle[q], TruthValuea, TruthValueb, TruthValuec, TruthValued);
		//forward calculation
	}
	for (int w = 0; w < ChooseNumOfIncidentAngle; w++)
	{
		Rmodel[w] = InverseCalculator(IncidentAngle[w], a, b, c, d); //inversion calculation
	}


	for (int e = 0; e < ChooseNumOfIncidentAngle; e++)
	{
		opt = opt + pow((Rmodel[e] - Robs[e]), 2); 
		//Calculate the difference between the reflection coefficient of each angle and the true value, the smaller the gap, the greater the fitness
	}

	
	opt = sqrt(opt);
	double IndividualFit = 1 / opt;
	//The inverse of the objective function is the fitness (the general fitness of the individual is calculated here)
	return IndividualFit;
}

double GetTotalFitOfPop(struct individual pop[])
//Get the total fitness of pop (used to calculate individual probability)
{
	double TotalFit = 0;
	for (int i = 0; i < PopSize; i++)
		TotalFit = TotalFit + GetFitOfIndividual(pop[i]);
	return TotalFit;
}

double GetProbOfIndividual(individual m) 
//get individual probability of each population
{
	double IndividualProb;
	double TotalFit;
	double IndividualFit;
	TotalFit = GetTotalFitOfPop(pop);
	IndividualFit = GetFitOfIndividual(m);
	IndividualProb = IndividualFit / TotalFit;
	return IndividualProb;
}

void PrintPop(struct individual pop[]) // print out the individuals of population
{
	double TotalFit;
	for (int i = 0; i < PopSize; i++)
	{
		TotalFit = GetTotalFitOfPop(pop);
		cout << "The " << GenPop << "th Generation and the " << i + 1 << "th individual is:";
		for (int j = 0; j < ChromoNumber; j++)
		{
			cout << pop[i].chromo[j] << "  ";
		}
		cout << "Its fitness is:" << GetFitOfIndividual(pop[i]) << "  ";
		cout << "The probability of chosen is: " << GetProbOfIndividual(pop[i]) * 100 << "%" << endl;
	}
}
//print pop

void GenerateRandomIndividuals(struct individual pop[]) //Get random individuals and generate initial population
{
	int i = 0;
	double TotalFit;
	double IndividualProb;
	for (i = 0; i < PopSize; i++)
	{
		pop[i].chromo[0] = getRanChromoSudubianhualv();
		pop[i].chromo[1] = getRanChromoMidubianhualv();
		pop[i].chromo[2] = getRanChromoSudubianhualv();
		pop[i].chromo[3] = getRanChromoSudubibianhualv1();
	}
}

void GetInitialPop(struct individual pop[]) //print the first generation
{
	GenerateRandomIndividuals(pop);
	PrintPop(pop);
}

bool InLegality(individual in) //used for judging legality
{
	double x1, x2, x3, x4;
	x1 = in.chromo[0];
	x2 = in.chromo[1];
	x3 = in.chromo[2];
	x4 = in.chromo[3];
	if ((x1 > 0.5 || x1 < 0) || (x2 > 0.2 || x2 < 0) || (x3 > 0.5 || x3 < 0) || (x4 > 2.2 || x4 < 1.8))
		return 0;
	else
		return 1;
}



double AverageFitness(struct individual pop[]) //calculate the average fitness
{
	double a = GetTotalFitOfPop(pop);
	double AF = a / PopSize;
	return AF;
}
//used to calculate the average fitness of the population
double ReturnBigger(double a, double b) 
{
	double c;
	if (a > b)
		c = a;
	else
		c = b;
	return c;
}

double FitnessMax(struct individual pop[]) //return the biggest fitness of population
{
	double Fit[PopSize];
	double Max = 0;
	for (int i = 0; i < PopSize; i++)
	{
		Fit[i] = GetFitOfIndividual(pop[i]);
	}
	for (int j = 0; j < PopSize; j++)
	{
		if (Max < Fit[j])
		{
			Max = Fit[j];
		}
	}
	return Max;
}
//return maximum fitness



double GetFitOfFinalValue(double a, double b, double c, double d)
{
	double opt = 0; 

	double Robs[ChooseNumOfIncidentAngle];
	double Rmodel[ChooseNumOfIncidentAngle];

	for (int q = 0; q < ChooseNumOfIncidentAngle; q++)
	{
		Robs[q] = ForwardCalculator(IncidentAngle[q], TruthValuea, TruthValueb, TruthValuec, TruthValued); //forward
	}
	for (int w = 0; w < ChooseNumOfIncidentAngle; w++)
	{
		Rmodel[w] = InverseCalculator(IncidentAngle[w], a, b, c, d); //inversion
	}
	
	for (int e = 0; e < ChooseNumOfIncidentAngle; e++)
	{
		opt = opt + pow((Rmodel[e] - Robs[e]), 2);
		//calculate the difference between the reflection coefficient of each angle and the true value, the smaller the gap, the greater the fitness
	}

	
	opt = sqrt(opt);
	double IndividualFit = 1 / opt;
	//The inverse of the objective function is the fitness (the general fitness of the individual is calculated here)
	return IndividualFit;
}
//the fitness used to calculate the last value



/*Catastrophic function, if the adaptation is too small, it will be catastrophic*/
void Catastrophe()
{
	if (GenPop == 20)
	{
		if (MaxFitnessInEveryGeneration[GenPop - 1] < 1000)
		{
			for (int i = 0; i < PopSize; i++)
			{
				pop[i].chromo[0] = getRanChromoSudubianhualv();
				pop[i].chromo[1] = getRanChromoMidubianhualv();
				pop[i].chromo[2] = getRanChromoSudubianhualv();
				pop[i].chromo[3] = getRanChromoSudubibianhualv1();
			}
			TimeOfCatastrophe++;
			GenPop = GenPop - 19;
		}
	}
}
/*Catastrophic function, if the adaptation is too small, it will be catastrophic*/


double StepSize()
{
	double d = 0;
	int a;
	double Criterion;
	for (int i = 0; i < 100; i++)
	{
		Criterion = rand() / double(RAND_MAX);
		double Pd = 0.01;
		if (Criterion < Pd)
		{
			a = 1;
		}
		else
		{
			a = 0;
		}
		d = d + (a / pow(2, i));
	}
	return d;
}


/*change the fitness */
void ChangeFit()
{
	double AF;
	double MaxFit;
	double a;
	double b;
	if (GenPop > 0)
	{
		AF = AverageFitness(pop);
		MaxFit = FitnessMax(pop);
		a = AF * (c - 1) / (MaxFit - AF);
		b = AF * (MaxFit - c * AF) / (MaxFit - AF);

		for (int NF = 0; NF < PopSize; NF++)
		{
			Fit[NF] = a * Fit[NF] + b;
		}
	}
}
/*change the fitness */


void Select(struct individual pop[]) //selection function (in this module, the value in pop is selected into nextpop)
{
	double Criterion; //判据
	Criterion = rand() / double(RAND_MAX);
	for (int j = 0; j < PopSize; j++)
	{
		if (Criterion <= LineProb[j])
		{
			cout << "The" << NextPopNum + 1 << "th Select, we've chosen the individual" << j + 1 << ", its fitness is" << Fit[j] << "the probability of chosen is:" << Prob[j] * 100 << "%" << endl;
			nextPop[NextPopNum] = pop[j];
			NextPopNum = NextPopNum + 1;
			break;
		}
	} //The roulette is over, the next population is full, at this time nextpopnum=100
}

void Cross(struct individual pop[])
{
	cout << "Begin to cross!" << endl;
	int AmountOfCross = 0; //used to record the number of crosses
	double Criterion;	   
	double r1, r2, r3, r4;

	int m, n;
	int z;
	int i = 0; //criterion for individuals

	while (i < PopSize)
	{

		Criterion = rand() / double(RAND_MAX);
		i++;

		if (Criterion < PcFixed)
		{
			m = (rand() % (PopSize)); //Randomly generate numbers from 0 to 100
			n = (rand() % (PopSize)); //Randomly generate numbers from 0 to 100

			r1 = rand() / double(RAND_MAX);
			r2 = rand() / double(RAND_MAX);
			r3 = rand() / double(RAND_MAX);
			r4 = rand() / double(RAND_MAX);

			nextPop[i - 1].chromo[0] = r1 * pop[n].chromo[0] + (1 - r1) * pop[m].chromo[0];
			nextPop[i].chromo[0] = r1 * pop[m].chromo[0] + (1 - r1) * pop[n].chromo[0];

			nextPop[i - 1].chromo[1] = r2 * pop[n].chromo[1] + (1 - r2) * pop[m].chromo[1];
			nextPop[i].chromo[1] = r2 * pop[m].chromo[1] + (1 - r2) * pop[n].chromo[1];

			nextPop[i - 1].chromo[2] = r3 * pop[n].chromo[2] + (1 - r3) * pop[m].chromo[2];
			nextPop[i].chromo[2] = r3 * pop[m].chromo[2] + (1 - r3) * pop[n].chromo[2];

			nextPop[i - 1].chromo[3] = r4 * pop[n].chromo[3] + (1 - r4) * pop[m].chromo[3];
			nextPop[i].chromo[3] = r4 * pop[m].chromo[3] + (1 - r4) * pop[n].chromo[3];

			if (InLegality(nextPop[i - 1]) && InLegality(nextPop[i])) //&&InLegalityFit(nextPop[i-1])&&InLegalityFit(nextPop[i])
				//The late crossover(more than 5 populations) is difficult to generate individuals with greater fitness, so we should delete it
			{
				AmountOfCross = AmountOfCross + 2;
				cout << "This crossover produced two new individuals, numbered as follows:" << i << "，" << i + 1 << ", and their fitness is:" << GetFitOfIndividual(nextPop[i - 1]) << "，" << GetFitOfIndividual(nextPop[i]) << endl;
				//At this time i=1, the next one should generate nextpop[2] and nextpop[3], so i should be incremented by 1
				i = i + 1;
			}
			else //One or more of two individuals is illegal
			{
				i--;
			}
		}
		else //The random number is greater than the probability of crossover
		{
			i--;
		}

		if (AmountOfCross >= (PopSize * PcFixed))
			break;
	}

	cout << "We have produced " << AmountOfCross << " new individuals" << endl;
}

void Mutation(struct individual pop[]) //mutation function (at this time nextpop still retains the results of the previous crossover)
{

	cout << "开始变异" << endl;
	int AmountOfMutation = 0; 
	double Criterion;		  
	double r;
	double k;

	int a, b;
	int z;
	double d;

	int i = 0; //criteria for individuals

    /* Adaptive algorithm constant */
	double A = 9.903438; 
	double AverageFit = AverageFitness(pop);
	double Max = FitnessMax(pop);
	double Pm;
	double IndividualFit;
	/* Adaptive algorithm constant */

	while (i < PopSize)
	{

		Criterion = rand() / double(RAND_MAX);
		a = (rand() % (PopSize)); //Select individuals first, then decide whether to mutate

		if (GenPop < 30)
		{
			k = 0.5;
		}
		if ((GenPop >= 30) && (GenPop < 50))
		{
			k = 0.05;
		}
		if ((GenPop >= 50) && (GenPop < 80))
		{
			k = 0.01;
		}
		if (GenPop >= 80)
		{
			k = 0.005;
		}

		/*Adaptive algorithm*/

		IndividualFit = GetFitOfIndividual(pop[a]);
		if (IndividualFit >= AverageFit)
		{
			Pm = (Pmmax - Pmmin) / (1 + exp(A * (2 * (IndividualFit - AverageFit) / (Max - AverageFit) - 1))) + Pmmin;
		}
		else
		{
			Pm = Pmmin;
		}
		/*Adaptive algorithm*/


		if (Criterion < Pm) //Determining whether a single individual mutates
		{
			z = rand() % 2; //randomly generate 0 or 1
			d = StepSize(); 

			if (z == 0)
			{
				nextPop[a].chromo[0] = pop[a].chromo[0] + k * (Xmax[0] - Xmin[0]) * d;
			}
			else
			{
				nextPop[a].chromo[0] = pop[a].chromo[0] - k * (Xmax[0] - Xmin[0]) * d;
			}

			z = rand() % 2; 
			d = StepSize();

			if (z == 0)
			{
				nextPop[a].chromo[1] = pop[a].chromo[1] + k * (Xmax[1] - Xmin[1]) * d;
			}
			else
			{
				nextPop[a].chromo[1] = pop[a].chromo[1] - k * (Xmax[1] - Xmin[1]) * d;
			}

			z = rand() % 2; 
			d = StepSize(); 

			if (z == 0)
			{
				nextPop[a].chromo[2] = pop[a].chromo[2] + k * (Xmax[2] - Xmin[2]) * d;
			}
			else
			{
				nextPop[a].chromo[2] = pop[a].chromo[2] - k * (Xmax[2] - Xmin[2]) * d;
			}

			z = rand() % 2; 
			d = StepSize(); 

			if (z == 0)
			{
				nextPop[a].chromo[3] = pop[a].chromo[3] + k * (Xmax[3] - Xmin[3]) * d;
			}
			else
			{
				nextPop[a].chromo[3] = pop[a].chromo[3] - k * (Xmax[3] - Xmin[3]) * d;
			}

			if (InLegality(nextPop[a])) 
			{
				AmountOfMutation++;
				i++;
				cout << "This mutation have produced a new individual, numbered as:" << a + 1 << ", and its fitness is: " << GetFitOfIndividual(nextPop[a]) << endl;
			}
			else 
			//If the newly generated individual does not meet the requirements, it will be replaced with the value of the original individual
			{
				nextPop[a] = pop[a];
			}
		}
		else //The random number is greater than the probability of crossover
		{
			i++; 
		}

		if (i == PopSize)
			break;
	}

	cout << "We have mutation for:" << AmountOfMutation << " times" << endl;
} //At this time, the latest population is nextpop

void Evolution() //Single evolution function (including print population function, selection, crossover, mutation function)
{
	/*The average fitness record of each population, used for disaster judgment*/
	AverageFitnessInEveryGeneration[GenPop - 1] = AverageFitness(pop); //Get the average fitness of each generation of the population
	MaxFitnessInEveryGeneration[GenPop - 1] = FitnessMax(pop);		   //Get the average fitness of each generation of the population
	

	Catastrophe(); //Catastrophe if the conditions are met

	/* Used to record the gap between each generation and the maximum value */
	double R = 0;
	double R1;
	double R2;
	for (int gap = 0; gap < ChooseNumOfIncidentAngle; gap++)
	{
		R1 = InverseCalculator(IncidentAngle[gap], pop[0].chromo[0], pop[0].chromo[1], pop[0].chromo[2], pop[0].chromo[3]);
		R2 = ForwardCalculator(IncidentAngle[gap], TruthValuea, TruthValueb, TruthValuec, TruthValued);
		R = R + fabs((R1 - R2) / R2);
	}
	GapBetweenEveryInversion[InversionRecord - 1][GenPop - 1] = R / 9;
	/* Used to record the gap between each generation and the maximum value */

	for (int f = 0; f < PopSize; f++)
	{
		Fit[f] = GetFitOfIndividual(pop[f]);
	}
	double TotalFit = GetTotalFitOfPop(pop);

/* Change the fitness of each individual here*/
/* Change the fitness of each individual here*/
/* Change the fitness of each individual here*/
// ChangeFit();
/* Change the fitness of each individual here*/
/* Change the fitness of each individual here*/
/* Change the fitness of each individual here*/

	for (int p = 0; p < PopSize; p++)
	{
		Prob[p] = Fit[p] / TotalFit;
	}

	LineProb[0] = Prob[0];
	for (int i = 1; i < PopSize; i++)
	{
		LineProb[i] = LineProb[i - 1] + Prob[i];
	}
	/*The linear probability of the selection step, etc.*/

	while (NextPopNum < PopSize)
	{
		Select(pop);
		if (NextPopNum == PopSize)
			break;
	}
	cout << endl;
	NextPopNum = 0;

	for (int m = 0; m < PopSize; m++) //redefine the population
	{
		pop[m] = nextPop[m];
	}

	// cout << "The selected population is:" << endl;
	// PrintPop(pop);
	cout << endl;

	Cross(pop); 

	for (int n = 0; n < PopSize; n++) 
	{
		pop[n] = nextPop[n];
	}

	// cout << "The crossed population is" << endl;
	// PrintPop(pop);
	cout << endl;

	Mutation(pop); 

	for (int j = 0; j < PopSize; j++)
	{
		pop[j] = nextPop[j]; //Redefine the next population as a new population
	}

	// cout << "The mutated population is" << endl;
	// PrintPop(pop);
	cout << endl;

	cout << "At the end of evolution, the population is:" << endl;
	PrintPop(pop); //print the new population
}

double Inversion() 
//Single inversion program (including obtaining initial population, multiple evolutions, and outputting results)
{
	double TotalProb = 0;
	cout << "The initial population" << endl;
	GetInitialPop(pop); //Get the initial randomly generated population

	cout << endl;
	cout << endl;

	/*evolution*/
	/*evolution*/
	while (GenPop - 1 < Generation) 
	{
		cout << endl;
		cout << "The" << GenPop << " th evolution begins:" << endl;
		Evolution();
		GenPop++; 
	}
	/*evolution*/
	/*evolution*/

	GenPop = 1; //After all evolutionary algebras in the above single inversion are over, the criterion is reset to 1

	double cosPhi;
	double JudgeFit = 0;
	double a;
	double b;
	double c;
	double d;
	double Q;
	double Lin;
	double Rmodel;

	for (int Fin = 0; Fin < PopSize; Fin++) 
	//Used to compare and assign the four chromosomes of the individual with the greatest fitness to abcd
	{
		if (GetFitOfIndividual(pop[Fin]) > JudgeFit)
		{
			JudgeFit = GetFitOfIndividual(pop[Fin]);
			a = pop[Fin].chromo[0];
			b = pop[Fin].chromo[1];
			c = pop[Fin].chromo[2];
			d = pop[Fin].chromo[3];
		}
	}
	cout << endl;
	cout << "The best solution is: (" << a << "," << b << "," << c << "," << d << ")"
		<< ", and the corresponding fitness is:" << GetFitOfFinalValue(a, b, c, d) << endl;
	cout << endl;
	int j = 0;
	for (j = 0; j < ChooseNumOfIncidentAngle; j++)
	{
		cout << "its corresponding angle is " << IncidentAngle[j] << ",the calculated reflection coefficient is" << InverseCalculator(IncidentAngle[j], a, b, c, d) << endl; //inversion
	}
	cout << endl;

	cout << "******************************" << endl;
	cout << endl;
	cout << endl;
	return 0;
}

int main()
{
	srand(time(0)); //the random number is different each time
	clock_t startTime, endTime;
	startTime = clock(); //record initial time

	/*Mean error, standard deviation, relative error calculation parameters*/
	double MeanErrora, MeanErrorb, MeanErrorc, MeanErrord;								  
	double StandardDeviationa, StandardDeviationb, StandardDeviationc, StandardDeviationd; 
	double RelativeErrora, RelativeErrorb, RelativeErrorc, RelativeErrord;				   
	double AdditionREa = 0;
	double AdditionREb = 0;
	double AdditionREc = 0;
	double AdditionREd = 0;
	double AdditionSDa = 0;
	double AdditionSDb = 0;
	double AdditionSDc = 0;
	double AdditionSDd = 0;
	/*Mean error, standard deviation, relative error calculation parameters*/

	/*Inversion result record array*/ /*Used to record the average value of the four chromosomes in the final inversion result population*/
	double a[InversionTimes] = { 0 };
	double b[InversionTimes] = { 0 };
	double c[InversionTimes] = { 0 };
	double d[InversionTimes] = { 0 };
	/*Inversion result record array*/ /*Used to record the average value of the four chromosomes in the final inversion result population*/

	/*Calculation parameters of the average value of the inversion result*/
	double Totala = 0;
	double Totalb = 0;
	double Totalc = 0;
	double Totald = 0;
	double Averagea, Averageb, Averagec, Averaged;
	/*Calculation parameters of the average value of the inversion result*/

	int i = 0;
	int m = 0;
	int n = 0;
	int k = 0;

	while (i < InversionTimes)
	{
		double JudgeFit = 0;
		cout << "The " << i + 1 << " th inversion:" << endl;
		Inversion();
		InversionRecord++;

		for (int Fin = 0; Fin < PopSize; Fin++) //used to compareand assign the four chromosomes of the individual with the greatest fitness to abcd
		{

			if (GetFitOfIndividual(pop[Fin]) > JudgeFit) //Choose the best fit
			{
				JudgeFit = GetFitOfIndividual(pop[Fin]);
				a[i] = pop[Fin].chromo[0];
				b[i] = pop[Fin].chromo[1];
				c[i] = pop[Fin].chromo[2];
				d[i] = pop[Fin].chromo[3];
			}
		}
		i = i + 1;
		cout << endl;
		cout << endl;
	}

	cout << endl;
	cout << endl;
	cout << endl;
	cout << endl;
	cout << endl;
	cout << "**********Multiple inversion results**********" << endl;
	for (int j = 0; j < InversionTimes; j++)
	{
		cout << "The " << j + 1 << "th inversion result is:"
			<< "(" << a[j] << "," << b[j] << "," << c[j] << "," << d[j] << ")"
			<< ", and the corresponding fitness is:" << GetFitOfFinalValue(a[j], b[j], c[j], d[j]) << endl;
	}
	cout << "*******Multiple inversion results are displayed over********" << endl;
	for (k = 0; k < InversionTimes; k++)
	{
		Totala = Totala + a[k];
		Totalb = Totalb + b[k];
		Totalc = Totalc + c[k];
		Totald = Totald + d[k];
	}


	Averagea = Totala / InversionTimes;
	Averageb = Totalb / InversionTimes;
	Averagec = Totalc / InversionTimes;
	Averaged = Totald / InversionTimes;



	MeanErrora = fabs(TruthValuea - Averagea) / TruthValuea;
	MeanErrorb = fabs(TruthValueb - Averageb) / TruthValueb;
	MeanErrorc = fabs(TruthValuec - Averagec) / TruthValuec;
	MeanErrord = fabs(TruthValued - Averaged) / TruthValued;
	
	for (m = 0; m < InversionTimes; m++)
	{
		AdditionSDa = AdditionSDa + pow((a[m] - Averagea), 2);
		AdditionSDb = AdditionSDb + pow((b[m] - Averageb), 2);
		AdditionSDc = AdditionSDc + pow((c[m] - Averagec), 2);
		AdditionSDd = AdditionSDd + pow((d[m] - Averaged), 2);
	}

	cout << endl;
	cout << endl;
	cout << endl;

	cout << "******Multiple inversion results displaying********" << endl;
	cout << "This is the " << DataNum + 1 << "th calculation results of a model, this model realizes" << InversionTimes << "times of inversion,"
		<< "In every inversion we have " << PopSize << " individuals and " << Generation << " evolution" << endl;
	cout << endl;
	cout << "The forward we have used the formula of " << Formula[ForwardNum] << ","
		<< "and the inversion we have used the formula of " << Formula[InverseNum] << endl;
	cout << endl;
	cout << "In this model:" << endl;
	cout << endl;
	cout << "The fixed crossover probability is: " << PcFixed << ", the fixed mutation probability is:" << PmFixed << endl;
	cout << endl;
	cout << "The true value calculated from the original data is: "
		<< "(" << TruthValuea << "," << TruthValueb << "," << TruthValuec << "," << TruthValued << ")" << endl;
	cout << endl;
	for (int k = 0; k < ChooseNumOfIncidentAngle; k++)
	{
		cout << "The angle corresponding is: " << IncidentAngle[k] << " ,and the corresponding reflection coefficient is" << ForwardCalculator(IncidentAngle[k], TruthValuea, TruthValueb, TruthValuec, TruthValued) << endl;
	}
	cout << endl;
	cout << "The average value of multiple inversions is:"
		<< "(" << Averagea << "," << Averageb << "," << Averagec << "," << Averaged << ")" << endl;
	cout << endl;
	for (int k = 0; k < ChooseNumOfIncidentAngle; k++)
	{
		double Rmodel = InverseCalculator(IncidentAngle[k], Averagea, Averageb, Averagec, Averaged);		   
		double Robs = ForwardCalculator(IncidentAngle[k], TruthValuea, TruthValueb, TruthValuec, TruthValued); 
		cout << "The angle is: " << IncidentAngle[k] << " and the reflection coefficient corresponding is:" << Rmodel << ", the percentage difference between the true value is: " << (fabs(Rmodel - Robs) / Robs) * 100 << "%" << endl;
	}
	cout << endl;
	cout << "The mean error is:"
		<< "(" << MeanErrora * 100 << "%，" << MeanErrorb * 100 << "%，" << MeanErrorc * 100 << "%," << MeanErrord * 100 << "%)" << endl;
	cout << endl;
	cout << "The standard deviation is:"
		<< "(" << SDCalculator(AdditionSDa) << "," << SDCalculator(AdditionSDb)
		<< "," << SDCalculator(AdditionSDc) << "," << SDCalculator(AdditionSDd) << ")" << endl;
	cout << endl;
	cout << "The relative deviation is:"
		<< "(" << SDCalculator(AdditionSDa) / Averagea * 100 << "%"
		<< "," << SDCalculator(AdditionSDb) / Averageb * 100 << "%"
		<< "," << SDCalculator(AdditionSDc) / Averagec * 100 << "%"
		<< "," << SDCalculator(AdditionSDd) / Averaged * 100 << "%"
		<< ")" << endl;
	cout << endl;
	cout << "In eight groups of models:" << endl;
	cout << "The rate of change of longitudinal wave velocity are:";
	for (int i = 0; i < 8; i++)
	{
		cout << (Vp2[i] - Vp1[i]) / (0.5 * (Vp2[i] + Vp1[i])) << " ";
	}
	cout << endl;
	cout << "The density change rates are:";
	for (int i = 0; i < 8; i++)
	{
		cout << (Rou2[i] - Rou1[i]) / (0.5 * (Rou2[i] + Rou1[i])) << " ";
	}
	cout << endl;
	cout << "The shear wave velocity change rates are:";
	for (int i = 0; i < 8; i++)
	{
		cout << (Vs2[i] - Vs1[i]) / (0.5 * (Vs2[i] + Vs1[i])) << " ";
	}
	cout << endl;
	cout << "The velocity ratio of longitudinal and shear waves are:";
	for (int i = 0; i < 8; i++)
	{
		cout << (Vp2[i] + Vp1[i]) / (Vs2[i] + Vs1[i]) << " ";
	}
	cout << endl;

	cout << "The total number of disasters: " << TimeOfCatastrophe << endl;
	cout << endl;
	endTime = clock();
	cout << "The total time for this experiment is: " << (int)((endTime - startTime) / CLOCKS_PER_SEC) / 60 << "minute" << (int)(endTime - startTime) / CLOCKS_PER_SEC - ((int)((endTime - startTime) / CLOCKS_PER_SEC) / 60) * 60 << "seconde" << endl;

	cout << "**********This is the end of the results**********" << endl;

	/*
	cout << "The difference between the average of each group and the reflection coefficient is:" << endl;
	for (int z = 0; z < InversionTimes; z++)
	{
		cout << "The" << z + 1 << "th value of secondary inversion and the true value difference is:" << endl;
		for (int x = 0; x < Generation; x++)
		{
			cout << "The " << x + 1 << "th inversion and the real value difference is" << GapBetweenEveryInversion[z][x] << " " << endl;
		}

	}
	*/
	//used to view the trend of convergence in each generation

	return 0;
}