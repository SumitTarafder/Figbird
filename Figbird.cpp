#include <stdio.h>
#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <assert.h>
#include <bits/stdc++.h>
using namespace std;

#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>

#define MAX_REC_LEN 1024
#define MAX_READLENGTH 200

int noContigs=0;
int noReads=0;
long int contigLength=0;
int maxReadLength=0;

long int totalContigLength=0;
vector<char*> contigs;
vector<char*> contigNames;
vector<long int> contigLengths;
double *contigReadCounts;
const int MAX_GAP=100000;


FILE *contigFile;
FILE *mapFile;
FILE *summaryFile;
FILE *countFile;
FILE * minmax1;
FILE * minmax2;
FILE* draw;
FILE* cases;

FILE *filledContigFile;

char *contigFileName;
char const *mapFileName;

long int *insertCounts;
int maxInsertSize=0;
int MAX_INSERT_SIZE;
double insertSizeMean;
double insertSizeVar;
double insertSizeSD;
int insertSizeMode;
int insertCutoffMax=0;
int insertCutoffMin=0;
int insertThresholdMax=0;
int insertThresholdMin=0;
int tolerance;
int insertCountMax;
int total_read=0,discarded_read=0;
double left_coeff=0,right_coeff=0;
double leftSD, rightSD;
int partial_flag=-1;
int priorSize=0;
long int *insertCountsMapped;
long int *insertCountsUnmapped;

long int errorTypes[5][5];
long int baseCounts[5];
long int *errorPos;
long int *inPos;
long int *inLengths;
long int *delPos;
long int *delLengths;
long int *readLengths;

long int *effectiveLengths;

double errorTypeProbs[5][5];
double baseErrorRates[5];
double *errorPosDist;
double *inPosDist;
double *inLengthDist;
double *delPosDist;
double *delLengthDist;
double *insertLengthDist;
double *insertLengthDistSmoothed;

int windowSize=12;

double *noErrorProbs;
int tmpCount=0;
int toAlign;

long int erroredReads=0;
long int uniqueMappedReads=0;
long int discardedReads=0;
long int totalCount,unCount;
int factor;
int covergae_threshold1;
int covergae_threshold2;
int partial_threshold;
int overlap_threshold;
int match_count_discont;
int clip_thresh;

int num_itr;
int flag=0;
int unmapped=-1;
int *gaptofill,*finalUsedReads;
double *timereq;
int Ncount1,Ncount2;
int step2match;
int partial_limit=3000;
int unmapped_limit=3000;
int upper_limit_insert;
int iszflag=1;//cfu

//zx
char tempCigar[500], tempMD[500];
char noErrorCigar[500], noErrorMD[500];

int left_maxDistance,right_maxDistance,script_itr=0,loadmodel;
int tempmaxDistance;
int it_count=0;
int charCodes[256];
pthread_mutex_t overlaps_mutex = PTHREAD_MUTEX_INITIALIZER;

int isJump=0;
int isConservative=0;
int isGaussian=0;
double gaussianMean;
double gaussianSD;
double inputMean;

int MIN_GAP=-500;

int MIN_READ=0;
int MIN_READ_JOIN=1;
int NO_INTERVALS=0;


double *insertProbsSums;

long int gapProbs[1000];
int gapProbCutOff;
int alloc_arg;

FILE *gapFile;
FILE *gapInfoFile;
FILE *gapOutFile;
FILE * gap_prob_dist;
FILE* stat2;
//=======================================================================

struct MAP
{
	double errorProb;
	int insertSize;
    long int pos1;
    long int pos2;
    long int contigNo1;
    long int contigNo2;
    int readLength1;
    int readLength2;
    bool isSameStrand;
};

struct InsertTable
{
	int insertSize;
    int count;
};


void initInsertCounts(int max)
{
	maxInsertSize=max;
	insertCounts=new long int[maxInsertSize];
	for(int i=0;i<maxInsertSize;i++)
	{
		insertCounts[i]=1;
	}
}

void updateInsertCounts(int index)
{
    
	if(index<=0)
		return;
	if(index<maxInsertSize)
	{
		insertCounts[index]++;
	}
	else
	{
		
		if(index>MAX_INSERT_SIZE)
		{
			discardedReads++;
			return;
		}
		int tempInsertSize=max(maxInsertSize*2,index);
		long int *tempCounts=new long int[maxInsertSize];
		for(int i=0;i<maxInsertSize;i++)
		{
			tempCounts[i]=insertCounts[i];
		}
		insertCounts=new long int[tempInsertSize];
		for(int i=0;i<maxInsertSize;i++)
		{
			insertCounts[i]=tempCounts[i];
		}
		for(int i=maxInsertSize;i<tempInsertSize;i++)
		{
			insertCounts[i]=1;
		}
        
		insertCounts[index]++;
		maxInsertSize=tempInsertSize;
		delete []tempCounts;
		
	}
    
}

void initErrorTypes(int readLength)
{
	for(int i=0;i<5;i++)
		for(int j=0;j<5;j++)
			errorTypes[i][j]=1;
    
	for(int i=0;i<5;i++)
		baseCounts[i]=1;
    
	errorPos=new long int[readLength];
	inPos=new long int[readLength];
	inLengths=new long int[readLength];
	delPos=new long int[readLength];
	delLengths=new long int[readLength];
	readLengths=new long int[readLength];
	
	for(int i=0;i<readLength;i++)
	{
		errorPos[i]=1;
		inPos[i]=1;
		inLengths[i]=1;
		delPos[i]=1;
		delLengths[i]=1;
		readLengths[i]=0;
	}
}


int getLength(char *read)
{
	int i=0;
	while(read[i])
	{
		if(read[i]=='A')
			baseCounts[0]++;
		else if(read[i]=='C')
			baseCounts[1]++;
		else if(read[i]=='G')
			baseCounts[2]++;
		else if(read[i]=='T')
			baseCounts[3]++;
		else
			baseCounts[4]++;
        
		i++;
	}
	
	return i;
}

long int getContigNo(char *contigName)
{
    
	return atol(contigName);
    
    for(long int i=0;i<contigNames.size();i++)
    {
		if(strcmp(contigNames[i],contigName)==0)
			return i;
	}
	return -1;
    
}

void processErrorTypes(char *cigar, char *md, char *read, int strandNo)
{
    
	int readLength=getLength(read);
	readLengths[readLength-1]++;
    
	if(strcmp(md,noErrorCigar)!=0)
		erroredReads++;
	else
		return;
    
    
	unsigned long mdLength=strlen(md)-5;
	unsigned long tempLength=0;
    
	char *temp;
	int index=0,totalLength=0;
	
    
	int curIndex=0;
	int *inserts=new int[readLength];
	
	for(int i=0;i<readLength;i++)
	{
		inserts[i]=0;
	}
    
	unsigned long cigarLength=strlen(cigar);
    //	char *tempCigar=new char[cigarLength];
	char cigarChar;
    
	strcpy(tempCigar,cigar);
    
    
	temp=strtok(tempCigar,"IDMS^\t\n ");
    
	while(temp!=NULL)
	{
        
		tempLength=atoi(temp);
		totalLength+=strlen(temp);
		cigarChar=cigar[totalLength];
        
		if(cigarChar=='M')
		{
			index+=tempLength;
			curIndex+=tempLength;
		}
		else if(cigarChar=='I' || cigarChar=='S')
		{
			if(strandNo==0)
			{
				inPos[index]++;
				inLengths[tempLength-1]++;
                
			}
			else
			{
                
				inPos[readLength-index-1]++;
				inLengths[tempLength-1]++;
			}
            
			inserts[curIndex]=tempLength;
            
			index+=tempLength;
		}
		else if(cigarChar=='D' )
		{
			if(strandNo==0)
			{
				delPos[index]++;
				delLengths[tempLength-1]++;
                
			}
			else
			{
                
				delPos[readLength-index-1]++;
				delLengths[tempLength-1]++;
			}
		}
		totalLength++;
		temp=strtok(NULL,"IDMS^\t\n ");
	}
    
	
	strcpy(tempMD,md);
    
	strtok(tempMD,":");
	strtok(NULL,":");
    
	
	index=0,totalLength=0,tempLength=0;
    
	int f, t;
    
	while((temp=strtok(NULL,"ACGTN^\t\n "))!=NULL)
	{
		tempLength=strlen(temp);
        
        
		totalLength+=tempLength;
		
        
		if(totalLength<mdLength)
		{
			char from=md[5+totalLength];
			
			
			if(from=='^')
			{
				totalLength++;
				index+=atoi(temp);
				for(int i=totalLength;i<mdLength;i++)
				{
					from=md[5+totalLength];
					if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
						totalLength++;
					else
						break;
                    
				}
			}
			else if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
			{
				totalLength++;
				index+=atoi(temp)+1;
				
				
				
				curIndex=0;
				for(int i=0;i<index;i++)
				{
					curIndex+=inserts[i];
				}
				char to=read[index-1+curIndex];
                
				if(strandNo==0)
					errorPos[index-1+curIndex]++;
				else
					errorPos[readLength-index-curIndex]++;
                
                
				switch(from)
				{
					case 'A':
						f=0;
						break;
					case 'C':
						f=1;
						break;
					case 'G':
						f=2;
						break;
					case 'T':
						f=3;
						break;
					default:
						f=4;
				}
                
				switch(to)
				{
					case 'A':
						t=0;
						break;
					case 'C':
						t=1;
						break;
					case 'G':
						t=2;
						break;
					case 'T':
						t=3;
						break;
					default:
						t=4;
				}
                
				if(f==t)
				{
                    
				}
				else
                    
					errorTypes[f][t]++;
                
			}
			else
				break;
		}
		
	}
	delete []inserts;
    
}

double dnorm(double x,double mean, double variance)
{
	double val=1/sqrt(M_PI*2*variance);
	val*=exp(-((x-mean)*(x-mean))/(2*variance));
	return val;
}


void computeProbabilites()
{

	int errorCount=0;
    
	for(int i=0;i<5;i++)
	{
		errorCount=0;
		for(int j=0;j<5;j++)
		{
			errorCount+=errorTypes[i][j];
		}
		for(int j=0;j<5;j++)
		{
			errorTypeProbs[i][j]=(double)errorTypes[i][j]/errorCount;
		}
		
		baseErrorRates[i]=errorCount/(double)baseCounts[i];
	}
    
	double sum=0;
	for(int i=0;i<4;i++)
		sum+=baseErrorRates[i];
    
	for(int i=0;i<4;i++)
	{
		baseErrorRates[i]=4*baseErrorRates[i]/sum;
	}
    
	baseErrorRates[4]=1;
    
	for(int i=maxReadLength-1;i>0;i--)
	{
		readLengths[i-1]=readLengths[i]+readLengths[i-1];
	}
    
	errorPosDist=new double[maxReadLength];
    
	for(int i=0;i<maxReadLength;i++)
	{
		errorPosDist[i]=(double)errorPos[i]/readLengths[i];
	}
    
	inPosDist=new double[maxReadLength];
    
	for(int i=0;i<maxReadLength;i++)
	{
		inPosDist[i]=(double)inPos[i]/readLengths[i];
	}
    
	inLengthDist=new double[maxReadLength];
    
	int inCount=0;
    
	for(int i=0;i<maxReadLength;i++)
	{
		inCount+=inLengths[i];
	}
	
	for(int i=0;i<maxReadLength;i++)
	{
		inLengthDist[i]=(double)inLengths[i]/inCount;
	}
    
	delPosDist=new double[maxReadLength];
    
	for(int i=0;i<maxReadLength;i++)
	{
		delPosDist[i]=(double)delPos[i]/readLengths[i];
	}
    
	delLengthDist=new double[maxReadLength];
    
	int delCount=0;
    
	for(int i=0;i<maxReadLength;i++)
	{
		delCount+=delLengths[i];
	}
	
	for(int i=0;i<maxReadLength;i++)
	{
		delLengthDist[i]=(double)delLengths[i]/delCount;
	}
	
    
	insertLengthDist=new double[maxInsertSize];
    
	long int insCount=discardedReads;
    
	sum=0;
    
	for(int i=0;i<maxInsertSize;i++)
	{
		insCount+=(insertCounts[i]-1);
		sum+=i*(insertCounts[i]-1);
        
	}
	insertSizeMean=sum/insCount;

	sum=0;
    
	for(int i=0;i<maxInsertSize;i++)
	{
		insertLengthDist[i]=(double)insertCounts[i]/insCount;
        
		sum+=(insertCounts[i]-1)*(insertSizeMean-i)*(insertSizeMean-i);
	}
    
	insertSizeVar=sum/insCount;
    
    insertSizeSD=sqrt(insertSizeVar);
    
	noErrorProbs=new double[maxReadLength];
    
	double noErrorProb=1.0;
    
    
	for(int i=0;i<maxReadLength;i++)
	{
		noErrorProb*=(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);
		noErrorProbs[i]=noErrorProb;
  	}
    
	effectiveLengths=new long int[maxInsertSize];
    
	for(int i=0;i<maxInsertSize;i++)
	{
		effectiveLengths[i]=-1;
	}
    
	long int totalContigLength=0;
	for(int i=0;i<contigLengths.size();i++)
	{
		totalContigLength+=contigLengths[i];
	}
	effectiveLengths[0]=totalContigLength;
	
    
    insertCountsMapped=new long int[maxInsertSize];
    insertCountsUnmapped=new long int[maxInsertSize];
    
    
    for(long int i=0;i<maxInsertSize;i++)
    {
        insertCountsMapped[i]=0;
        insertCountsUnmapped[i]=0;
    }
    
    insertLengthDistSmoothed=new double[maxInsertSize];
    double windowSum=0;
    
    for(int i=0;i<windowSize;i++)
    {
        insertLengthDistSmoothed[i]=insertLengthDist[i];
        
    }
    
    for(int i=0;i<2*windowSize+1;i++)
    {
        windowSum+=insertLengthDist[i];
        
    }
    insertLengthDistSmoothed[windowSize]=windowSum/(2*windowSize+1);
    
    for(int i=windowSize+1;i<maxInsertSize-windowSize;i++)
    {
        windowSum-=insertLengthDist[i-windowSize-1];
        windowSum+=insertLengthDist[i+windowSize];
        insertLengthDistSmoothed[i]=windowSum/(2*windowSize+1);
    }
    for(int i=maxInsertSize-windowSize;i<maxInsertSize;i++)
    {
        insertLengthDistSmoothed[i]=insertLengthDist[i];
    }
    
	for(int i=0;i<maxInsertSize;i++)
	{
		insertLengthDistSmoothed[i]=insertLengthDistSmoothed[i]-1/(double)(insCount)+(1/(double)maxInsertSize)/(double)(insCount+1);
        
	}
    
    int count=0;
    
    for(int i=insertSizeMean;i<maxInsertSize;i++)
    {
        if(insertCounts[i]<=1)
        {
            count++;
            if(count==10)
            {
                insertCutoffMax=i;
                break;
            }
        }
        else
        {
            count=0;
        }
        
    }
    
    count=0;
    
    for(int i=insertSizeMean;i>=0;i--)
    {
        if(insertCounts[i]<=1)
        {
            count++;
            if(count==10)
            {
                insertCutoffMin=i;
                break;
            }
        }
        else
        {
            count=0;
        }
        
    }
    
    insertCountMax=0;
    for(int i=0;i<maxInsertSize;i++)
    {
        if(insertCounts[i]>insertCountMax)
        {
            insertCountMax=insertCounts[i];
            insertSizeMode=i;
            
        }
    }
    
    
    count=0;
    
    for(int i=insertSizeMean;i<maxInsertSize;i++)
    {
        if(insertCounts[i]<=max(insertCountMax/1000,2))
        {
            count++;
            if(count==2)
            {
                insertThresholdMax=i;
                break;
            }
        }
        else
        {
            count=0;
        }
        
    }
    
    count=0;
    
    for(int i=insertSizeMean;i>=0;i--)
    {
        if(insertCounts[i]<=max(insertCountMax/1000,2))
        {
            count++;
            if(count==2)
            {
                insertThresholdMin=i;
                break;
            }
        }
        else
        {
            count=0;
        }
        
    }

    //mean used instead of mode
    
    double insertSum=0;
    double insertCount=0;
	
    for(int i=insertCutoffMin;i<insertCutoffMax;i++)
    {
        insertCount+=insertCounts[i]-1;
        insertSum+=(insertCounts[i]-1)*i;
    }
    insertSizeMode=insertSum/insertCount;
    

    
	insertSum=0;
    insertCount=0;
    
	for(int i=insertSizeMean+1;i<maxInsertSize;i++)
	{
		insertSum=insertSum+(insertCounts[i]-1)*(i-insertSizeMean)*(i-insertSizeMean);
		insertCount+=(insertCounts[i]-1);
	}
	rightSD=sqrt(insertSum/insertCount);

	insertSum=0;
    insertCount=0;
	for(int i=max((int)(insertSizeMean-10*rightSD),0);i<insertSizeMean;i++)
	{
		insertSum=insertSum+(insertCounts[i]-1)*(insertSizeMean-i)*(insertSizeMean-i);
		insertCount+=(insertCounts[i]-1);
	}
	leftSD=sqrt(insertSum/insertCount);
    
	if((rightSD>1000 || leftSD>1000) && isGaussian==0)
	{
		cerr<<"Switching to conservative mode"<<endl;
		isConservative=1;
	}

	//cout<<"Insertsizemean = "<<insertSizeMean<<" , LeftSD = "<<leftSD<<" , rightSD = "<<rightSD<<endl;

	if(isGaussian==1)
	{
		insertThresholdMin=max((int)(gaussianMean-2.5*gaussianSD),1);
		insertThresholdMax=min((int)(gaussianMean+2.5*gaussianSD),maxInsertSize);
		insertSizeMode=gaussianMean;
        /*
         for(int i=0;i<insertThresholdMin;i++)
         {
         insertLengthDistSmoothed[i]=(1/(double)maxInsertSize)/(double)(insCount+1);
         }
         for(int i=insertThresholdMin;i<insertThresholdMax;i++)
         {
         insertLengthDistSmoothed[i]=dnorm(i,gaussianMean,gaussianSD*gaussianSD)+(1/(double)maxInsertSize)/(double)(insCount+1);
         }
         for(int i=insertThresholdMax;i<maxInsertSize;i++)
         {
         insertLengthDistSmoothed[i]=(1/(double)maxInsertSize)/(double)(insCount+1);
         }
         */
		for(int i=0;i<maxInsertSize;i++)
		{
			insertLengthDistSmoothed[i]=dnorm(i,gaussianMean,gaussianSD*gaussianSD)+(1/(double)maxInsertSize)/(double)(insCount+1);
		}
        
	}

    insertCutoffMax=insertThresholdMax;
    insertCutoffMin=insertThresholdMin;
    
//	MAX_GAP=insertCutoffMax;
    
	
}

void processMapping(char *line)
{
	
	char * temp;
	char *qname, *rname, *mapq;
	int	pos,flag;
	char * cigar, * readString; // * md, *nhstring;
	long int contigNo;
    
    
	char md[500];
	char nhstring[500];
    
	int nh;
    
	int strandNo=0;
    
    
	qname=strtok(line,"\t");
	
    
	temp=strtok(NULL,"\t");
	flag=atoi(temp);
	
    
	strandNo=(flag&16)>>4;
    
    rname=strtok(NULL,"\t");
	
    
	temp=strtok(NULL,"\t");
	pos=atoi(temp);
	
	
	cigar=strtok(NULL,"\t");
	
	
	temp=strtok(NULL,"\t");
	
	readString=strtok(NULL,"\t");
    
	int insertSize=atoi(temp);
    
	while((temp=strtok(NULL,"\t\n"))!=NULL)
	{
		if(temp[0]=='M' && temp[1]=='D')
		{
			strcpy(md,temp);
		}
		else if(temp[0]=='I' && temp[1]=='H')
		{
			strcpy(nhstring,(temp+5));
			nh=atoi(nhstring) ;
		}
        
	}
    
	
	if(nh==1 && md[5]!='^')
	{
		contigNo=getContigNo(rname);
		if(isGaussian==1)
		{
			updateInsertCounts(insertSize);
		}
		else
		{
			if(contigLengths[contigNo]>inputMean)
				updateInsertCounts(insertSize);
		}
		processErrorTypes(cigar,md,readString,strandNo);
		uniqueMappedReads++;
	}
	
    
}

long int getEffectiveLength(int insertSize)
{
	if(insertSize<0)
		return effectiveLengths[0];
    
	if(insertSize>=maxInsertSize)
	{
		long int effectiveLength=0;
		for(int i=0;i<contigLengths.size();i++)
		{
			if(contigLengths[i]>=insertSize)
				effectiveLength+=(contigLengths[i]-insertSize+1);
		}
		return effectiveLength;
        
	}
	if(effectiveLengths[insertSize]==-1)
	{
		long int effectiveLength=0;
		for(int i=0;i<contigLengths.size();i++)
		{
			if(contigLengths[i]>=insertSize)
				effectiveLength+=(contigLengths[i]-insertSize+1);
		}
		effectiveLengths[insertSize]=effectiveLength;
	}
	return effectiveLengths[insertSize];
}

long double computeErrorProb(char *cigar, char *md, char *read, int strandNo)
{
    
    unsigned long readLength=strlen(read);
	
    
	long double errorProb=noErrorProbs[readLength-1];
    
    
    
	if(md[5]=='^')
		return errorProb;
	
	
	char tempMD[1000], tempCigar[1000];
    
	unsigned long mdLength=strlen(md)-5;
	unsigned long tempLength=0;
	
	char *temp;
	int index=0,totalLength=0;
	
    
	int curIndex=0;
	int *inserts=new int[readLength];
	
	for(int i=0;i<readLength;i++)
	{
		inserts[i]=0;
	}
    
    //	int cigarLength=strlen(cigar);
	char cigarChar;
    
	strcpy(tempCigar,cigar);
    
	temp=strtok(tempCigar,"IDM^\t\n ");
    
	while(temp!=NULL)
	{
        
		tempLength=atoi(temp);
		totalLength+=strlen(temp);
		cigarChar=cigar[totalLength];
        
		if(cigarChar=='M')
		{
			index+=tempLength;
			curIndex+=tempLength;
		}
		else if(cigarChar=='I')
		{
			unsigned long i;
			if(strandNo==0)
			{
				//look up insert probs
				i=index;
				
			}
			else
			{
				i=readLength-index-1;
			}
            
			errorProb=errorProb*inPosDist[i]*inLengthDist[tempLength-1]/(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);
            
			inserts[curIndex]=tempLength;
            
			index+=tempLength;
		}
		else if(cigarChar=='D')
		{
			unsigned long i;
			if(strandNo==0)
			{
				i=index;
                //	look up delete probs
			}
			else
			{
				i=readLength-index-1;
			}
            
			errorProb=errorProb*delPosDist[i]*delLengthDist[tempLength-1]/(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);
		}
		totalLength++;
		temp=strtok(NULL,"IDM^\t\n ");
	}
    
    
	strcpy(tempMD,md);
    
	strtok(tempMD,":");
	strtok(NULL,":");
    
	index=0,totalLength=0,tempLength=0;
    
	int f, t;
    
	while((temp=strtok(NULL,"ACGTN^\t\n "))!=NULL)
	{
		tempLength=strlen(temp);
        
		totalLength+=tempLength;
		
		if(totalLength<mdLength)
		{
			char from=md[5+totalLength];
            
			if(from=='^')
			{
				totalLength++;
				index+=atoi(temp);
				for(int i=totalLength;i<mdLength;i++)
				{
					from=md[5+totalLength];
					if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
						totalLength++;
					else
						break;
				}
			}
			else if(from=='A' || from=='C' || from=='G' || from=='T'|| from=='N')
			{
				totalLength++;
				index+=atoi(temp)+1;
				
                
				curIndex=0;
				for(int i=0;i<index;i++)
				{
					curIndex+=inserts[i];
				}
				char to=read[index-1+curIndex];
                
				int i;
				if(strandNo==0)
					i=index-1+curIndex;
				else
					i=readLength-index-curIndex;
                
                
				errorProb=errorProb*errorPosDist[i]/(1-errorPosDist[i]-inPosDist[i]-delPosDist[i]);
                
                
				switch(from)
				{
					case 'A':
						f=0;
						break;
					case 'C':
						f=1;
						break;
					case 'G':
						f=2;
						break;
					case 'T':
						f=3;
						break;
					default:
						f=4;
				}
                
				switch(to)
				{
					case 'A':
						t=0;
						break;
					case 'C':
						t=1;
						break;
					case 'G':
						t=2;
						break;
					case 'T':
						t=3;
						break;
					default:
						t=4;
				}
                
				if(f==t)
				{
					
                    
				}
				else
				{
					//errorTypeProb
					errorProb*=baseErrorRates[f]*errorTypeProbs[f][t];
				}
			}
			else
				break;
            
		}
		
	}
	
	delete []inserts;
	return errorProb;
}


double computeLikelihood(char const *file)
{
    mapFile=fopen(file, "r");
	char line1[MAX_REC_LEN];
	char line2[MAX_REC_LEN];
    
	char *qname1,*qname2,preqname1[500],preqname2[500];
    
	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line1[0]);
    
	long double sum=0.0;
	long double logsum=0.0;
    
	char * temp;
	char *rname1, *rname2;
	int	pos1,pos2,flag,strandNo1, strandNo2, insertSize1, insertSize2;
	char *cigar1, *cigar2, *readString1, *readString2, md1[1000], md2[1000];
    
    
	long double insertSizeProb;
    
	long double errorProb1, errorProb2;
	
	long double gapProb;
    
    int tempInsertSize=0;
    long double tempProb=0;
    
	preqname1[0]=preqname2[0]='*';
	preqname1[1]=preqname2[1]='\0';
    
    
	int it=0;
    
	while(fgets(line1, MAX_FILE_READ, mapFile)!=NULL)
	{
		if(line1[0]=='@')
			continue;
        //????
		if(fgets(line2, MAX_FILE_READ, mapFile)==NULL)
			break;
        
		qname1=strtok(line1,"\t");
		
        temp=strtok(NULL,"\t");
        flag=atoi(temp);
		
        
		strandNo1=(flag&16)>>4;
        
        temp=strtok(NULL,"\t");
		
		temp=strtok(NULL,"\t");
		pos1=atoi(temp);
        
		
		cigar1=strtok(NULL,"\t");
        
		
		temp=strtok(NULL,"\t");
		insertSize1=atoi(temp);
        
        readString1=strtok(NULL,"\t");
        
        
		while((temp=strtok(NULL,"\t\n"))!=NULL)
		{
			if(temp[0]=='M' && temp[1]=='D')
			{
				strcpy(md1,temp);
			}
		}
        
        //        		cout<<insertSize1<<" "<<cigar1<<" "<<md1<<endl;
		
        //second of the pair
        
		qname2=strtok(line2,"\t");
		temp=strtok(NULL,"\t");
		flag=atoi(temp);
		
        
		strandNo2=(flag&16)>>4;
        
        temp=strtok(NULL,"\t");
		
		temp=strtok(NULL,"\t");
		pos2=atoi(temp);
        
		cigar2=strtok(NULL,"\t");
        
		temp=strtok(NULL,"\t");
		
		insertSize2=atoi(temp);
        
		readString2=strtok(NULL,"\t");
        
        
		while((temp=strtok(NULL,"\t\n"))!=NULL)
		{
			if(temp[0]=='M' && temp[1]=='D')
			{
				strcpy(md2,temp);
			}
		}
        
        //cout<<insertSize2<<" "<<cigar2<<" "<<md2<<endl;

		int insertSize=max(insertSize1, insertSize2);
		
        
		insertSizeProb=0;
        
		if(insertSize>=0 && insertSize<maxInsertSize)
		{
			insertSizeProb=insertLengthDist[insertSize];
		}
        
		if(insertSizeProb==0)
		{
			insertSizeProb=1/(double)uniqueMappedReads;
		}
		
        
		errorProb1=computeErrorProb(cigar1,md1,readString1,strandNo1);
		errorProb2=computeErrorProb(cigar2,md2,readString2,strandNo2);

		long int totalEffectiveLength=getEffectiveLength(insertSize);
		//cout<<insertSize<<"\t"<<totalEffectiveLength<<endl;
        
        
		long double prob=(1/(long double)(totalEffectiveLength))*insertSizeProb*errorProb1*errorProb2;

        //        cout<<errorProb1<<" "<<errorProb2<<" "<<insertSizeProb<<" "<<prob<<endl;
        if(strcmp(qname1,preqname1)==0 && strcmp(qname2,preqname2)==0)
		{
            if(tempProb<prob)
            {
                tempProb=prob;
                tempInsertSize=insertSize;
                tempInsertSize=tempInsertSize<0?0:tempInsertSize;
                gapProb=errorProb2;
            }
            
			sum+=prob;
            //cout<<"In if\n";
		}
		else if(strcmp("*",preqname1)!=0 && strcmp("*",preqname2)!=0)
		{
			if(sum<1e-320 || isnan(sum))
			{
				sum=1e-320;
			}
			logsum+=log10(sum);
            
            
            int gapIndex=-log10(gapProb);
            gapIndex++;
            
            if(gapIndex<1000 && gapIndex>=0)
            {
            	gapProbs[gapIndex]++;
            }
            else
            {
            	gapProbs[999]++;
            }
            
            
            if(tempInsertSize>=maxInsertSize)
                insertCountsMapped[maxInsertSize-1]++;
            else
                insertCountsMapped[tempInsertSize]++;
			
			sum=prob;
            
            tempProb=prob;
            tempInsertSize=insertSize;
            tempInsertSize=tempInsertSize<0?0:tempInsertSize;
            
            gapProb=errorProb2;
            //cout<<"In else if\n";
		}
		else
		{
			sum=prob;
			
            tempProb=prob;
            tempInsertSize=insertSize;
            tempInsertSize=tempInsertSize<0?0:tempInsertSize;
            
			gapProb=errorProb2;
			//cout<<"In else\n";
		}
        
		strcpy(preqname1,qname1);
		strcpy(preqname2,qname2);
		it++;
        
		
		if(isinf( logsum ))
		{
			//cout<<it<<endl;
			exit(1);
		}
        
        
	}
	if(sum!=0)
    {
        if(sum<1e-320 || isnan(sum))
        {
            sum=1e-320;
        }
		logsum+=log10(sum);
    }
    
	fclose(mapFile);
	
	return logsum;
}

void printHelp()
{
    
	cout<<"cgal v0.9.5-beta"<<endl;
	cout<<"----------------"<<endl;
	cout<<endl;
	cout<<"cgal - computes likelihood"<<endl;
	cout<<"Usage:"<<endl;
	cout<<"cgal [options] <contigfile>"<<endl;
	cout<<endl;
	cout<<"Required arguments:"<<endl;
	cout<<"<contigfile>\t Assembly file in FASTA format"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
	cout<<"-h [--help]\t\t Prints this message"<<endl;
	cout<<endl;
	cout<<"Output: "<<endl;
	cout<<"(In file out.txt) <numberContigs> <totalLikelihood> <mappedLikelihood> <unmappedLikelihood> <noReads> <noReadsUnmapped>"<<endl;
	cout<<"<numberContigs>\t\t Number of contigs"<<endl;
	cout<<"<totalLikelihood>\t Total log likelihood value"<<endl;
	cout<<"<mappedLikelihood>\t Likelihood value of reads mapped by the mapping tool"<<endl;
	cout<<"<unmappedLikelihood>\t Likelihood value corresponding to reads not mapped by alignment tool"<<endl;
	cout<<"<noReads>\t\t Total number of paired-end reads"<<endl;
	cout<<"<noReadsUnmapped>\t Number of reads not mapped by the alignment tool"<<endl;
	cout<<endl;
	exit(1);
    
}

void reverseStr(char *read)
{
	char ch;
	int start = 0;
	int end = strlen(read)-1;
	//cout<<end<<endl;
	int len = end;

	while (start < end)
    {
        ch = read[start];
        read[start] = read[end];
        read[end] = ch;
        start++;
        end--;
    }
	read[len+1]='\0';
}


void reverse(char *reverse, char *read)
{
    
	char ch='A';
	int readLength=strlen(read);
	//cout<<readLength<<endl;
	for(int i=1; i<=readLength;i++)
	{
		ch=read[i-1];
		if(ch=='A' || ch=='a')
			reverse[readLength-i]='T';
		else if(ch=='C' || ch=='c')
			reverse[readLength-i]='G';
		else if(ch=='G' || ch=='g')
			reverse[readLength-i]='C';
		else if(ch=='T' || ch=='t')
			reverse[readLength-i]='A';
		else
			reverse[readLength-i]='N';
	}
	reverse[readLength]='\0';
    
}


int getDistance(char *s, char * t, int sStart, int sEnd, int tStart, int tEnd, int ** dis)
{
	int m=sEnd-sStart+1;
    int n=tEnd-tStart+1;
	
    
	int val1;
	int val2;
	int val3;
    
	
    int i,j;
    
	for(i=0;i<=m;i++)
        dis[i][0]=i;
    
	for(j=0;j<=n;j++)
        dis[0][j]=j;
    
	for(i=1;i<=m;i++)
	{
		for(j=1;j<=n;j++)
		{
			val2=dis[i-1][j]+1;
			val3=dis[i][j-1]+1;
			val2=val2<val3?val2:val3;
			
            val1=dis[i-1][j-1];
            if(s[i-1+sStart]!=t[j-1]+tStart)
			{
				val1++;
                
			}
			val3=val1<val2?val1:val2;
			
            dis[i][j]=val3;
			
		}
		
	}
    
    
    
    return 	dis[m][n];
    
}

int getMismatch(char *s, char * t, int sStart, int sEnd, int tStart, int tEnd, double error)
{
	int m=sEnd-sStart+1;
    int n=tEnd-tStart+1;
    
	int mismatches=0;
    
	for(int i=0;i<m,i<n;i++)
	{
		if(s[sStart+i]!=t[tStart+i])
		{
			mismatches++;
			if(mismatches>error)
				return mismatches;
		}
	}
	return mismatches;
    
}


int getOverlap(char *s, char *t, double error)
{
    int sLength=strlen(s);
    int tLength=strlen(t);
    int distance=0;
    
    int **dis;
    
    dis=new int*[sLength+1];

    for(int i=0;i<=sLength;i++)
    {
        dis[i]=new int[tLength+1];
    }
    
    
    for(int i=sLength>tLength?sLength-tLength:0;i<sLength;i++)
    {
        distance=getDistance(s,t , i, sLength-1, 0, sLength-i-1,dis);
        
        if(distance<=error*(sLength-i))
        {
            for(int j=0;j<=sLength;j++)
            {
                delete [] dis[j];
            }
            
            delete [] dis;
            
            return (sLength-i);
        }
    }
    
    for(int i=0;i<=sLength;i++)
    {
        delete [] dis[i];
    }
    
    delete [] dis;
    
    return 0;
}

class GapFiller
{
    double **countsGap;
    double **new_counts_gap;
    double **probsGap;
    double **errorProbsGap;
    double **qual_gap;
    double **partial_count_array;

    int *gap_coverage;
    double **partial_prob_arr;
    double** partial_read_quality;


    int gap_no;
    int contigNo;
    long int gapStart;
    int gapLength;
    int valid_count;
    int invalid_count;
    int partial_read_count;
    int p_count;
    double region_perct,region_perct_max_gap_estm;
    long int end_pos_max;


    int originalGap;
    int partial_read_len,unmapped_read_len,max_read_len;

    int fillNotfill,perfectReadGap,perfectReadGaplen;
    int **repeatflag;
    int rep_flag,one_side_repeat_flag,large_gap_flag,fillupflag,part_read_c;
    int fill_original,comp_count;

    long int startPos, endPos;
	char gapFileName[500];
	char partialgapFileName[500];
    char *concensus;
    char *bestString;
	char *secondBestString;
	char *current_str,*previous_str,*current_str2;

	char ** reads_gap;
	int * pos_reads;
	int * isReverse;

	char ** partialreads_gap;
    int * partialpos_reads;
    int * preadflag;
    int * clipped_index;

	int num_reads_in_gap;

    int *mark_accepted_reads;
    double *maxlikelihood_value;
    int **final_readpos;
    int *saved_reads;
    int **partial_read_flag;

    int maxGap;
    int left_max,right_min,discont_or_not;
    char *gap_left,*gap_right;
    char partial_left[100],partial_right[100];
    int partial_saved_read_temp[2],partial_saved_read_final[2];
    int side_limit;
    int **unmapped_read_pos_arr_org;
    int **partial_read_pos_arr_org;
    int maxSize;
    int mid_limitu,mid_limitp;
    float gp_frac1,gp_frac2;
    int umaxleftf,umaxrightf,ucoverf;
    int negoverlap;

public:

   void allocate(int maxGap,int read_count,int read_count2,int len,int stat1,int stat2,int stat3,int neglap,int mp,int up,int lgf,float gpf1,float gpf2)
    {
        this->maxGap=maxGap;
        this-> negoverlap = neglap;
        this->mid_limitu=up;
        this->mid_limitp=mp;
        this->gp_frac1=gpf1;
        this->gp_frac2=gpf2;
        this->large_gap_flag=lgf;
        num_reads_in_gap = read_count;
        partial_read_count = read_count2;
        side_limit=30;
        this->maxSize=maxGap+2*tempmaxDistance;
        concensus=new char[maxGap+1];
        concensus[0]='N';
        end_pos_max=0;
        valid_count=invalid_count=0;
        
        fillupflag=0;
        p_count=0;

        countsGap=new double *[maxSize];
        probsGap=new double *[maxSize];
        errorProbsGap=new double *[maxSize];
        new_counts_gap=new double *[maxSize];

        partial_prob_arr = new double*[maxGap];
        qual_gap = new double*[maxGap];
        gap_coverage = new int[maxGap];
        partial_count_array = new double*[maxGap];

        unmapped_read_len=len;
        partial_read_len = mp/2;;

        fillNotfill = stat1;
        perfectReadGap = stat2;
        perfectReadGaplen = stat3;
        rep_flag=0;
        one_side_repeat_flag=0;
        fill_original=0;

        valid_count=0;
        invalid_count=0;
        region_perct=0;
        part_read_c=0;

        if(read_count>0)
        {
            pos_reads = new int[read_count];
            isReverse = new int[read_count];
            mark_accepted_reads = new int[read_count];
            maxlikelihood_value = new double[read_count];

            saved_reads = new int[read_count];

            unmapped_read_pos_arr_org = new int*[read_count];
            final_readpos = new int*[read_count];


            for(int i=0;i<read_count;i++)
            {

               unmapped_read_pos_arr_org[i] = new int[3];
               final_readpos[i] = new int[3];
            }

            reads_gap = new char*[read_count];

            for(int i=0;i<read_count;i++)
            {
                reads_gap[i] = new char[unmapped_read_len+1];
                saved_reads[i] = 0;

                pos_reads[i] = 0;
                isReverse[i] = 0;
                maxlikelihood_value[i] = 0;
                unmapped_read_pos_arr_org[i][0] = -200;
                unmapped_read_pos_arr_org[i][1] = 0;
                unmapped_read_pos_arr_org[i][2] = 0;

                final_readpos[i][0] = -200;
                final_readpos[i][1] = 0;//length
                final_readpos[i][2] = -1;//read number
            }
        }

        if(partial_read_count > 0)
        {
            maxlikelihood_value = new double[partial_read_count];
            partial_read_pos_arr_org = new int*[partial_read_count];


            for(int i=0;i<partial_read_count;i++)
            {
               partial_read_pos_arr_org[i] = new int[3];
               partial_read_pos_arr_org[i][0] = 0;
               partial_read_pos_arr_org[i][1] = -200;
               partial_read_pos_arr_org[i][2] = 0;

               maxlikelihood_value[i] = 0;
            }
        }

        for(int i=0;i<maxGap;i++)
        {
            partial_prob_arr[i] = new double[5];
            partial_count_array[i] = new double[4];

            qual_gap[i] = new double[5];
            for(int j=0;j<5;j++)
            {
                if(j<4)partial_prob_arr[i][j] = 1;
                else partial_prob_arr[i][j] = 0;
            }
            gap_coverage[i] = 0;
        }

        bestString=new char[maxGap+1];
        secondBestString=new char[maxGap+1];
        current_str=new char[maxGap+1];
        current_str2=new char[maxGap+1];
        previous_str=new char[maxGap+1];
        gap_left = new char[side_limit+1];
        gap_right = new char[side_limit+1];

        gap_left[0]='\0';
        gap_right[0]='\0';

        for(long int i=0;i<maxSize;i++)
        {
            countsGap[i]=new double[5];
            probsGap[i]=new double[5];
            errorProbsGap[i]=new double[5];
            new_counts_gap[i]=new double[5];
        }

        for(int i=0;i<2;i++)
        {
            partial_saved_read_temp[i] = partial_saved_read_final[i] = -1;
        }
    }

    void qualityFilter(char* qual_str,int index,int called)
    {
        int flag=0,Q;
        double expected_num_errors = 0;
        int Q_count=0;
        int Q_threshold = 6,Q_c_thresh=10;
        double temp,E_MAX = 10,val;//0.5,.25
        int read_len=strlen(qual_str);

        for(int i=0;i<read_len;i++)
        {
            Q = qual_str[i]-33;//ASCII-Base-33 for Illumina
            val = pow(10,-Q/10.0);
            expected_num_errors += val;
            if(Q<Q_threshold)Q_count++;
            if(called == 1)partial_read_quality[index][i] = val;//partial
        }
    }

    int findRepeat()
    {
        FILE* partial_gap_file =fopen(partialgapFileName,"r");
        int flag1=0,flag2=0,flag=0;
        char* readstring;
        int clipped_index;
        char *line1= new char[MAX_REC_LEN];
        int MAX_FILE_READ=MAX_REC_LEN/sizeof(line1[0]);
        char* token;
        int match=-1,len,p_count=0;
        vector<size_t> positions;

        while(fgets(line1, MAX_FILE_READ, partial_gap_file)!=NULL)
        {
            p_count++;
            if(p_count>partial_limit)break;
        }

        fclose(partial_gap_file);

        this->part_read_c = p_count;

        if(p_count>0)
        {
            partial_gap_file =fopen(partialgapFileName,"r");

            repeatflag = new int*[p_count];
            for(int i=0;i<p_count;i++)
                repeatflag[i] = new int[3];

            for(int i=0;i<p_count;i++)
                for(int j=0;j<3;j++)
                    this->repeatflag[i][j]=-1;

            p_count=0;
            while(fgets(line1, MAX_FILE_READ, partial_gap_file)!=NULL)
            {
                token = strtok(line1,"\t");
                readstring = token;

                token = strtok(NULL, "\t");
                clipped_index = atoi(token);

                token = strtok(NULL, "\t");
                match = atoi(token);
                len = strlen(readstring);

                int n = 20;//lowest amount of character allowed to match
                std::string s1(this->gap_left),s2(this->gap_right),s3,s4,s5(readstring);

                int lim = strlen(this->gap_left)-n;

                for(int i = 0;i<lim;i++)
                {
                    s3 = s1.substr(i,s1.size());
                    size_t pos = s5.find(s3, 0);

                    while(pos != string::npos)
                    {
                        positions.push_back(pos);
                        pos = s5.find(s3,pos+1);
                    }
                    if(positions.size() > 1)
                    {
                        //cout<<"Found repeat in gap = "<<this->gap_no<<",for left side in read - "<<p_count<<"\tRepeat = "<<s3<<"\tLength = "<<s3.length()<<endl;
                        repeatflag[p_count][0]=1;//left
                        repeatflag[p_count][1] = s3.length();
                        repeatflag[p_count][2] = positions[0];
                        //cout<<repeatflag[p_count][1]<<"\t"<<repeatflag[p_count][2]<<endl;
                        flag1=1;
                        flag2=p_count;
                        one_side_repeat_flag=1;
                        break;
                    }
                    positions.clear();
                }

                positions.clear();
                lim = strlen(this->gap_right)-n;

                for(int i = 0;i<lim;i++)
                {
                    s4 = s2.substr(0,s2.size()-i);
                    size_t pos = s5.find(s4, 0);

                    while(pos != string::npos)
                    {
                        positions.push_back(pos);
                        pos = s5.find(s4,pos+1);
                    }
                    if(positions.size() > 1)
                    {
                        //for(int i=0;i<;i++)
                        //cout<<"Position - "<<positions2[i]<<endl;
                        //cout<<"Found repeat in gap = "<<this->gap_no<<",for right side in read - "<<p_count<<"\tRepeat = "<<s4<<"\tLength = "<<s4.length()<<endl;
                        repeatflag[p_count][0]=2;//left
                        repeatflag[p_count][1] = positions[positions.size()-1];
                        if(flag1 ==1 && p_count == flag2)flag=1;
                        one_side_repeat_flag=1;
                        break;
                    }
                    positions.clear();
                }

                positions.clear();
                p_count++;
                if(p_count>partial_limit)break;
            }
            fclose(partial_gap_file);
        }
        delete[] line1;
        return flag;
    }

    void update_partial_prob(int gaplen)
    {
        FILE* partial_gap_file =fopen(partialgapFileName,"r");

        char* readstring;
        int clipped_index;
        char line1[MAX_REC_LEN];
        int MAX_FILE_READ=MAX_REC_LEN/sizeof(line1[0]);
        char* token;
        int match=-1,len;
        int val=MAX_GAP;
        left_max=-val,right_min=val;

        for(int i=0;i<gaplen;i++)
            for(int j=0;j<=4;j++)
                partial_count_array[i][j] = 1;

        for(int i=0;i<maxGap;i++)
        {
            for(int j=0;j<5;j++)
            {
                if(j<4)partial_prob_arr[i][j] = 1;
                else partial_prob_arr[i][j] = 0;
            }
        }

        int p_count=0;

        while(fgets(line1, MAX_FILE_READ, partial_gap_file)!=NULL)
        {
            token = strtok(line1,"\t");
            readstring = token;

            token = strtok(NULL, "\t");
            clipped_index = atoi(token);

            token = strtok(NULL, "\t");
            match = atoi(token);
            len = strlen(readstring);

            if(len > this->partial_read_len)this->partial_read_len=len;

            if(repeatflag[p_count][0] != -1)
            {
                if(repeatflag[p_count][0] == 1)
                {
                    //cout<<"Found repeat for left side in read - "<<p_count<<endl;
                    //for(int i=0;i<positions1.size();i++)
                    //cout<<"Position - "<<positions1[i]<<endl;
                    clipped_index = repeatflag[p_count][2]+repeatflag[p_count][1]-1;
                }

                if(repeatflag[p_count][0] == 2)
                {
                    //cout<<"Found repeat for right side in read - "<<p_count<<endl;
                    //for(int i=0;i<;i++)
                    //cout<<"Position - "<<positions2[positions2.size()-1]<<endl;
                    clipped_index = repeatflag[p_count][1];
                }
            }

            int stop_index1 = min(len - clipped_index - 1,gaplen);

            int stop_index2 = -1;
            if(clipped_index <= gaplen)stop_index2 = 0;
            else stop_index2= clipped_index - gaplen;

           // cout<<match<<" "<<stop_index1<<" "<<stop_index2<<endl;
            if(match==1 || match ==4)
            {
                int j=0;
                for(int i = clipped_index+1; i<clipped_index+1+stop_index1;i++,j++)
                {
                    if(charCodes[readstring[i]]<4)partial_count_array[j][charCodes[readstring[i]]] += 1;
                    else
                    {
                        for(int h=0;h<4;h++)partial_count_array[j][h] += 1;
                    }
                    //cout<<readstring[i];
                }

                if(j-1 > left_max)left_max = j-1;
                
            }
            else if(match==2 || match ==3)
            {
                int j= gaplen-1;
                for(int i = clipped_index-1; i>=stop_index2;i--,j--)
                {
                    if(charCodes[readstring[i]]<4)partial_count_array[j][charCodes[readstring[i]]] += 1;
                    else
                    {
                        for(int h=0;h<4;h++)partial_count_array[j][h] += 1;
                    }
                    //cout<<readstring[i];
                }
                if(j+1 < right_min)right_min = j+1;
                
            }
           
            p_count++;
            if(p_count>partial_limit)break;
        }

        fclose(partial_gap_file);

        if(partial_flag ==1)this->partial_read_count = p_count;

        //cout<<"max index left- "<<left_max<<"\t"<<"min index right - "<<right_min<<endl;

        if(0)
        {
            for(int i=0;i<gaplen;i++)
            {
                cout<<i<<" - ";
                for(int j=0;j<4;j++)
                {
                    cout<<partial_count_array[i][j]<<" ";
                }
                cout<<endl;
            }
        }

        int total_count_p =0,left_c=0,right_c=0;
        //cout<<left_max<<" "<<right_min<<" "<<partial_read_count<<endl;

        for(int i=0;i<gaplen;i++)
        {
            total_count_p =0;
            int max_index=0,max_val=-1;
            for(int k=0;k<4;k++)
            {
                total_count_p += partial_count_array[i][k];
                if(partial_count_array[i][k] > max_val)
                {
                    max_val = partial_count_array[i][k];
                    max_index = k;
                }
            }
            
            if(i<=left_max-5 || i>=right_min+5)
            {
                char ch;
                if(max_index==0)ch='A';
                else if(max_index==1)ch = 'C';
                else if(max_index==2)ch = 'G';
                else ch = 'T';

                if(i<=left_max-5)
                {
                    partial_left[left_c++] = ch;
                }
                else
                {
                    partial_right[right_c++] = ch;
                }
                //cout<<i<<"\t"<<ch<<endl;
            }

           for(int k=0;k<4;k++)
            {
                probsGap[i+left_maxDistance][k] = partial_count_array[i][k]/total_count_p;
                if(total_count_p>4)
                {
                    partial_prob_arr[i][k] = partial_count_array[i][k];//it's a count
                    partial_prob_arr[i][4] = 1;
                }
            }
        }

        partial_left[left_c]='\0';
        partial_right[right_c]='\0';
        //cout<<left_c<<"\t"<<right_c<<endl;
        //cout<<partial_left<<"\t"<<partial_right<<endl;

    }

	void computeProbsGap(int flag)
	{
	    double total=0;
	    double Ncount=0;
	    for(long int i=0;i<endPos-startPos;i++)
	    {
	        total=0;
	        for(int j=0;j<5;j++)
	        {
	            total=total+countsGap[i][j];
	        }
	        Ncount=countsGap[i][4];

	        for(int j=0;j<4;j++)
	        {
	            if(total)probsGap[i][j]=((countsGap[i][j]+(Ncount/4))/total);
	            else probsGap[i][j] = .25;
	        }
	        probsGap[i][4]=0;
	    }
	    //zx
	    if(flag == 1)
	    {
	        //draw_read(draw,this->gapLength,"S",0,0,0,0,' ');
	        update_partial_prob(gapLength);
        }
	}
	
	void computeErrorProbsGap()
	{
	    double sum=0;
	    
	    for(long int i=0;i<endPos-startPos;i++)
	    {
	        for(int j=0;j<5;j++)
	        {
	            sum=0;
	            for(int k=0;k<4;k++)
	            {
	                if(j==k)
	                    continue;
	                
	                sum+=probsGap[i][k]*(errorTypeProbs[k][j]);
	            }
	            errorProbsGap[i][j]=sum;
	        }
	    }
	}
	
	void initGapFiller(int g,int contigNo,long int gapStart, int originalGap, char * gapFileName, char* partialgapFileName)
	{
	    this->gap_no = g;
	    this->gapStart=gapStart;
	    this->contigNo=contigNo;
	    
	    this->originalGap=originalGap;
	    
	    strcpy(this->gapFileName,gapFileName);
	    strcpy(this->partialgapFileName,partialgapFileName);
	}

    void findGapLeftRight()
    {
        this->gapLength=this->originalGap;
        initialize_start_end();



        int c_j=0;
        for(long int i=0;i<left_maxDistance;i++)
        {
            if(left_maxDistance-i <= this->side_limit)gap_left[c_j++]=contigs[contigNo][i+startPos];
        }

        gap_left[c_j]='\0';
        c_j=0;

        for(long int i=left_maxDistance+gapLength;i<left_maxDistance+gapLength+right_maxDistance;i++)
        {
            if(c_j < this->side_limit)gap_right[c_j++]=contigs[contigNo][i+startPos-gapLength+originalGap];
        }
        gap_right[c_j]='\0';
        //cout<<"Left = "<<gap_left<<"\tRight = "<<gap_right<<endl;
        //cout<<left_maxDistance<<"\t"<<right_maxDistance<<endl;
    }

    int find_contig_match(char* left,char* right)
    {
        //return 0;
        char line1[MAX_REC_LEN];
        if(originalGap > this->negoverlap)return 0;
        int n = 3;//lowest amount of character allowed to match
        std::string s1(left),s2(right),s3,s4;

        for(int i = 0;i<side_limit-n;i++)
        {
            s3 = s1.substr(i,s1.size());
            s4 = s2.substr(0,s2.size()-i);

            int pos = s3.find(s4);
            if(pos != -1)
            {
                std::string rem_str = s2.substr(s4.size(),s2.size());

                FILE* partial_gap_file =fopen(partialgapFileName,"r");
                char* readstring;

                int MAX_FILE_READ=MAX_REC_LEN/sizeof(line1[0]);
                char* token;
                int len,max_match=-1,max_pos=-1,match_count=0,mismatch=0,error_thresh=2;
                int part_count=0;

                while(fgets(line1, MAX_FILE_READ, partial_gap_file)!=NULL)
                {
                    max_match=-1,max_pos=-1;
                    token = strtok(line1,"\t");
                    readstring = token;
                    len = strlen(readstring);

                    for(int j=0;j<len-s1.size();j++)
                    {
                        match_count=0,mismatch=0;
                        for(int k=0,l=0;k<s1.size();k++,l++)
                        {
                            if(readstring[j+l] == s1[k])match_count++;
                            else mismatch++;
                            if(mismatch> error_thresh)break;
                        }
                        //cout<<readstring<<endl;
                        //cout<<"pos = "<<j<<" "<<"\tmatch count = "<<match_count<<"\tcurrent max = "<<max_match<<endl;
                        if(match_count> max_match)
                        {
                            max_match=match_count;
                            max_pos = j;

                        }
                    }

                    if(s1.size() - max_match <= error_thresh)
                    {
                        int newpos = max_pos + s1.size(),match=0;
                        int matchval=s4.size();

                        //for(int g=newpos;g+rem_str.size()<len-1;g++)
                        {
                            //cout<<newpos<<endl;
                            //cout<<s3<<" "<<s4<<endl;
                            match=0;

                             //cout<<readstring<<" "<<max_match<<" "<<max_pos<<" "<<s1.size()<<endl;

                            for(int j=0;j<rem_str.size();j++)
                            {
                                //cout<<rem_str[j]<<"\t"<<readstring[g+j]<<endl;
                                if(rem_str[j] == readstring[newpos+j])match++;
                            }
                            //cout<<rem_str<<endl;
                            //cout<<match<<" "<<endl;

                            if(rem_str.size() - match <= error_thresh)
                            {
                                //cout<<"Match found= "<<matchval<<"\tGaplength = "<<originalGap<<"\tpos = "<<pos<<endl;
                                fclose(partial_gap_file);

                                return matchval;
                            }
                            //matchval--;
                        }
                    }
                    part_count++;
                    if(part_count > partial_limit)break;
                }//while
                fclose(partial_gap_file);
            }
        }

        return 0;
    }

	void initialize_start_end()
	{

	    if(this->gapStart-left_maxDistance < 0)
        {
            startPos=0;
            //cout<<"Left_contig limit passed"<<endl;
            left_maxDistance= this->gapStart;
            if(side_limit>left_maxDistance)side_limit = left_maxDistance;
        }
        else
        {
            startPos = this->gapStart-left_maxDistance;
        }

        if(this->gapStart+this->gapLength+right_maxDistance > contigLengths[this->contigNo])
        {
            endPos=contigLengths[this->contigNo];
            //cout<<"Right_contig limit passed"<<endl;
            right_maxDistance = contigLengths[this->contigNo]-(this->gapStart+this->gapLength);
            if(side_limit>right_maxDistance)side_limit=right_maxDistance;
        }
        else
        {
            endPos=this->gapStart+this->gapLength+right_maxDistance;
        }
        if(endPos > end_pos_max)end_pos_max = endPos;
	}
	
	int initialize(int gapLength,int flag,int negGapCheck)//called using gapestimate
	{
	    this->gapLength=gapLength;

	    initialize_start_end();

	   // cout<<"endpos = "<<endPos<<" ,start = "<<startPos<<endl;
	    for(long int i=0;i<endPos-startPos;i++)
	    {
	        for(int j=0;j<5;j++)
	        {
	            countsGap[i][j]=0;
	            probsGap[i][j]=0;
	            errorProbsGap[i][j]=0;
	            new_counts_gap[i][j]= 0;
	        }
	        //cout<<i<<endl;
	    }

        region_perct=0;

	    for(int i=0;i<num_reads_in_gap;i++)
        {
            //final_readpos[i] =0;
            maxlikelihood_value[i] = 0;
        }
        for(int i=0;i<this->partial_read_count;i++)
        {
            maxlikelihood_value[i] = 0;
        }

        for(int i=0;i<maxGap;i++)
        {
            for(int j=0;j<5;j++)
            {
                partial_prob_arr[i][j] = 1;
                qual_gap[i][j]=0;
            }
            gap_coverage[i] = 0;
        }

	    
	    char ch;

	    for(long int i=0;i<left_maxDistance;i++)
	    {
	        ch=contigs[contigNo][i+startPos];

	        //cout<<"After contig arr "<<i<<endl;
	        if(charCodes[ch] < 0 || charCodes[ch] > 4)
            {
                //cout<<"Contig No. "<<contigNo<<", charcode value = "<<charCodes[ch]<<", character = "<<ch<<endl;
                //countsGap[i][0]++;
            }
            else countsGap[i][charCodes[ch]]++;
            //cout<<"After countsgap arr "<<i<<endl;
	    }

	    //cout<<"B???"<<endl;
	    for(long int i=left_maxDistance;i<left_maxDistance+gapLength;i++)
	    {
	        countsGap[i][4]++;
	    }
        //cout<<"C???"<<endl;
	    for(long int i=left_maxDistance+gapLength;i<left_maxDistance+gapLength+right_maxDistance;i++)
	    {
	        ch=contigs[contigNo][i+startPos-gapLength+originalGap];
	        if(charCodes[ch] < 0 || charCodes[ch] > 4)
	        {
	            //cout<<"Contig No. "<<contigNo<<", charcode value = "<<charCodes[ch]<<", character = "<<ch<<endl;
	            //countsGap[i][0]++;
	        }
	        else countsGap[i][charCodes[ch]]++;
	        //cout<<"After countsgap arr "<<i<<endl;
	    }

	    int gfp=0;
	    if(side_limit>0 && negGapCheck == 0 && flag == 1)gfp = find_contig_match(gap_left,gap_right);
	    //cout<<"HERE???"<<endl;
	    computeProbsGap(flag);
	    //cout<<"HERE???"<<endl;
	    computeErrorProbsGap();
	    //cout<<"HERE???"<<endl;
	    return gfp;
	    
	}

	void draw_read(FILE* fp,int length,const char *s,int readno, double temp_val,int isz,int g,char type)
	{
	    int readlen = unmapped_read_len;//this->unmapped_read_len> this->partial_read_len?this->unmapped_read_len:this->partial_read_len;
	    if(s[0] == 'S')
	    {
	        for(int i=0;i<readlen;i++)
                fprintf(fp,"%s"," ");
	        fprintf(fp,"====================+Gap = %d starting,length = %d===============================\n",g,length);
	        for(int i=0;i<readlen;i++)
                fprintf(fp,"%s"," ");
	        char new_s[length+1];
	        for(int i=0;i<length;i++)
	        {
	            new_s[i]='N';
	        }
	        new_s[length]='\0';
            fprintf(fp,"%s\n",new_s);
        }
        else if(s[0] == 'Z')
        {
            for(int num=0;num<2;num++)
            {
                for(int i=0;i<readlen;i++)
                    fprintf(fp,"%s"," ");
                for(int i=0;i<length;i++)
                {
                    if(num==0)fprintf(fp,"%d\t",i);
                    else fprintf(fp,"%d\t",gap_coverage[i]);
                }
                fprintf(fp,"\n");
             }
        }
        else
        {
            for(int i=0;i<readlen+length;i++)
                fprintf(fp,"%s"," ");
            fprintf(fp,"%s[",s);

            fprintf(fp,"%d %d isz = %d %c]",readno,length,isz,type);
            //else fprintf(fp,"%d %d isz = %d %c]",readno,(int)temp_val,isz,type);
            fprintf(fp,"\n");
        }
	}

	double getDiff(const char * target, const char* ref, int length)
    {
        double diff=0;
        double frac=0;
        for(int i=0;i<length;i++)
        {
            if(toupper(target[i])!=ref[i])
            {
                diff++;
            }
        }
        frac = diff/length;
        //cout<<frac<<endl;
        return frac;
    }

    int find_partial_match(char* ref,char* search,int pos,int c,int len_t)//zxcf
    {
        int len_r = strlen(ref),len_s = strlen(search);
        std::string s1(ref),s2(search),s3,s4;
        int len_thresh;
        if(c==0)len_thresh=len_t;
        else len_thresh=4;

        //cout<<"len_r = "<<len_r<<"\tlen_s = "<<len_s<<endl;
        //cout<<"ref = "<<ref<<"\tread = "<<search<<endl;
        if(len_r >= len_s && len_s >=len_thresh)
        {
            if(pos == 0)//left
            {
                s3 = s1.substr(s1.size()-len_s,len_s);
                s4 = s2.substr(s2.size()-len_s,len_s);//read
                //cout<<"gap left= "<<s3<<"\tread left = "<<s4<<endl;
                double frac=getDiff(s3.c_str(),s4.c_str(),len_s);
                if(frac < 0.2 && c ==1)return 1;
                if(frac <= 0.08 && c ==0)
                {
                    //cout<<"Frac = "<<frac<<"\tleft side len = "<<len_s<<endl;
                    return 1;
                }
            }
            else
            {
                s3 = s1.substr(0,len_s);//ref
                s4 = s2.substr(0,len_s);//read
                //cout<<"gap right = "<<s3<<"\tread right = "<<s4<<endl;
                double frac=getDiff(s3.c_str(),s4.c_str(),len_s);

                if(frac < 0.2 && c ==1)
                {
                    //cout<<"Overlap Frac = "<<frac<<"\tright side len = "<<len_s<<endl;
                    return 1;
                }
                if(frac <= 0.08 && c ==0)
                {
                    //cout<<"Perfect Frac = "<<frac<<"\tright side len = "<<len_s<<endl;
                    return 1;
                }
            }
        }
        return 0;
    }

    char* get_read_frag(char* read,int placed_pos)
    {
        std::string s1(read),s2;
        static char c[30];
        int neg = -placed_pos;

        if(placed_pos < 0)//left side
        {
            if(neg < side_limit)s2 = s1.substr(0,neg);
            else s2 = s1.substr(neg-side_limit,side_limit);
        }
        else
        {
            s2 = s1.substr(this->gapLength-placed_pos,side_limit);
        }

        strcpy(c,s2.c_str());

        return c;
    }

    void detect_overlap_gapestimate(int **pflag,int gaplen,int *ret_v,int len_thresh)
    {
        FILE* partial_gap_file =fopen(partialgapFileName,"r");
        char* readstring;
        char perfstr[MAX_READLENGTH];
        char line1[MAX_REC_LEN];
        int MAX_FILE_READ=MAX_REC_LEN/sizeof(line1[0]);
        char* token;
        int match=-1;
        int l_max=-MAX_GAP,r_min=MAX_GAP;
        int left_cross[1000],right_cross[1000];
        int left_c=0,right_c=0;
        int p_count=0,len,placed_pos=0,overlap_count=0,read_covered=0;
        double mismatch_threshold = .1;

        char partial_reads[this->partial_read_count][MAX_READLENGTH];
        int match_reads[this->partial_read_count];
        int sm_flag[this->partial_read_count];
        
        for(int f=0;f<this->partial_read_count;f++)sm_flag[f]=0;

        while(fgets(line1, MAX_FILE_READ, partial_gap_file)!=NULL)
        {
            token = strtok(line1,"\t");
            readstring = token;
            strtok(NULL, "\t");
            token = strtok(NULL, "\t");
            match = atoi(token);
            match_reads[p_count] = match;
            len = strlen(readstring);

            for(int i=0;i<=len;i++)partial_reads[p_count][i] = readstring[i];
            //cout<<p_count<<" = "<<placed_pos<<", len = "<<len<<" "<<readstring<<endl;
            if(pflag[p_count][0] == 0)
            {
                p_count++;
                if(p_count > partial_limit)break;
                continue;
            }

            int pos = pflag[p_count][1];
            int j,flag=0,start=-1;

            for(j=0;j<len;j++)
            {
                if(pos+j>=0 && pos+j<gaplen)
                {
                    if(flag==0){flag=1;start=pos+j;}
                }
                if(pos+j == gaplen)break;
            }
            //A single read that covers the whole gap, can't affect both of these values
            if(match==1 || match ==4 || pos < 0)
            {
                if(pos+j-1 > l_max)l_max = pos+j-1;
            }
            else
            {
                if(start < r_min)
                {
                    r_min = start;
                    //cout<<p_count<<"\t"<<r_min<<endl;
                }
            }

            p_count++;
            if(p_count > partial_limit)break;
        }

        fclose(partial_gap_file);


        if(l_max == -MAX_GAP)l_max = -1;
        if(r_min == MAX_GAP)r_min = gaplen;

        //cout<<"max index from left side = "<<l_max<<"\t"<<"min index from right side= "<<r_min<<endl;
        //for(int i=0;i<this->partial_read_count;i++)
          //         cout<<pflag[i][0]<<"\t"<<pflag[i][1]<<endl;

        overlap_count = l_max - r_min + 1;

        int ovflag=0;

        for(int k=0;k<this->partial_read_count;k++)
        {
            if(pflag[k][0] == 1)//read used to fill gap
            {
                strcpy(readstring,partial_reads[k]);
                len = strlen(readstring);
                placed_pos = pflag[k][1];

                if(placed_pos <0 && placed_pos + len > this->gapLength)//covers the whole gap
                {
                    read_covered=1;
                    char l[side_limit],r[side_limit];
                    std::string s1(readstring),s2;
                    //cout<<s1<<"\t"<<placed_pos<<endl;
                    if(-placed_pos < side_limit)s2 = s1.substr(0,-placed_pos);
                    else s2 = s1.substr(-placed_pos-side_limit,side_limit);
                    strcpy(l,s2.c_str());

                    s2 = s1.substr(-placed_pos+this->gapLength,side_limit);
                    strcpy(r,s2.c_str());

                    //cout<<"L = "<<l<<"\tR = "<<r<<endl;
                    if(find_partial_match(gap_left,l,0,0,len_thresh) && find_partial_match(gap_right,r,1,0,len_thresh))
                    {
                        //cout<<"Perfect read - "<<k<<endl;
                        ovflag=1;
                        strcpy(perfstr,readstring);
                    }
                }
                if(placed_pos < 0 && placed_pos + len -1 >= r_min && placed_pos + len <= this->gapLength)left_cross[left_c++] = k;
                if(placed_pos > 0 && placed_pos <= l_max)right_cross[right_c++] = k;
                if(placed_pos < 0 && placed_pos+ len > this->gapLength && (match_reads[k] ==2 || match_reads[k] ==3))
                {
                    right_cross[right_c++] = k;
                    sm_flag[k]=1;
                    //cout<<"sm_flag for read = "<<k<<endl;;
                }
            }
        }

        //cout<<ovflag<<"\t"<<read_covered<<endl;
        if(ovflag || (this->perfectReadGap ==1 && originalGap <=20 && gaplen == this->perfectReadGaplen))
        {
            ret_v[0] = 300;
            ret_v[1] = 0;

            //cout<<"Perfect Read found, gaplength = "<<gaplen<<endl;
            //cout<<perfstr<<endl;
            return;
        }

        if(r_min <= l_max)
        {
            int max_overlap=0;
            int max_overlap_index=-1;
            int false_overlap_flag=0;
            int p1,p2;
            char readstring1[MAX_READLENGTH],readstring2[MAX_READLENGTH];
            //cout<<l_max<<"\t"<<r_min<<"\t"<<left_c<<"\t"<<right_c<<endl;

            //for(int i=0;i<left_c;i++)cout<<left_cross[i]<<endl;
            //for(int i=0;i<right_c;i++)cout<<right_cross[i]<<endl;

            for(int i=0;i<left_c;i++)
            {
                for(int j=0;j<right_c;j++)
                {
                    if(left_cross[i] != right_cross[j])
                    {
                        strcpy(readstring1,partial_reads[left_cross[i]]);
                        p1 = pflag[left_cross[i]][1];
                        len = strlen(readstring1);

                        strcpy(readstring2,partial_reads[right_cross[j]]);
                        p2 = pflag[right_cross[j]][1];

                        int diff_gap = p1 + len - gaplen;
                        if(diff_gap > 0)
                        {
                            overlap_count = (p1 + len -1) - p2 + 1 - diff_gap;
                        }
                        else
                        {
                            overlap_count = (p1 + len -1) - p2 + 1;
                            diff_gap=0;
                        }

                        //cout<<overlap_count<<"\t"<<left_cross[i]<<"\t"<<right_cross[j]<<"\t"<<p1<<"\t"<<p2<<endl;
                        if(overlap_count >= overlap_threshold)
                        {
                            //cout<<overlap_count<<"\t"<<left_cross[i]<<"\t"<<right_cross[j]<<"\t"<<p1<<"\t"<<p2<<endl;
                            std::string common_left,common_right;
                            if(sm_flag[right_cross[j]] != 1)
                            {
                                if(find_partial_match(gap_left,get_read_frag(readstring1,p1),0,1,-1))
                                {
                                    std::string s1(readstring1);
                                    common_left = s1.substr(s1.size()-overlap_count-diff_gap,overlap_count);
                                    //cout<<left_cross[i]<<"- cl = "<<common_left<<endl;
                                }
                                if(find_partial_match(gap_right,get_read_frag(readstring2,p2),1,1,-1))
                                {
                                    std::string s2(readstring2);
                                    common_right = s2.substr(0,overlap_count);
                                    //cout<<right_cross[j]<<"- cr = "<<common_right<<endl;
                                }
                            }
                            else
                            {
                                //cout<<"In else\n";

                                int x = pflag[right_cross[j]][1];
                                if(find_partial_match(gap_left,get_read_frag(readstring1,p1),0,1,-1))
                                {
                                    std::string s1(readstring1);
                                    common_left = s1.substr(s1.size()-overlap_count-x,overlap_count-x);

                                }
                                std::string s2(readstring2),s3;
                                s3 = s2.substr(-x+this->gapLength,side_limit);
                                char temp_r[side_limit];
                                strcpy(temp_r,s3.c_str());
                                if(find_partial_match(gap_right,temp_r,1,1,-1))
                                {
                                    common_right = s2.substr(-x,overlap_count+x);
                                }
                            }

                            int len1 = common_left.size(),len2 = common_right.size();

                            //cout<<len1<<"\t"<<len2<<endl;
                            if( len1 > 0 && len2 > 0 && len1 == len2)
                            {
                                double mismatch_frac = getDiff(common_left.c_str(),common_right.c_str(),len1);

                                //cout<<"mismatch frac = "<<mismatch_frac<<endl;
                                if(mismatch_frac <= mismatch_threshold)
                                {
                                    //cout<<"Overlap Length = "<<overlap_count<<"\tMismatch% = "<<mismatch_frac<<"\tOverlapped string = "<<common_left<<"\t"<<common_right<<endl;
                                    if(len1 > max_overlap)
                                    {
                                        max_overlap = len1;
                                        max_overlap_index = p2;
                                        partial_saved_read_temp[0] = left_cross[i];
                                        partial_saved_read_temp[1] = right_cross[j];
                                    }
                                }
                                else
                                {
                                     false_overlap_flag = -1;
                                     //cout<<"In placeread, False Overlap detected..with mismatch% = "<<mismatch_frac<<", The false strings = "<<common_left<<"\t"<<common_right<<endl;
                                }
                            }
                            else
                            {
                                //if(c==1)cout<<"Scn3 found for gap\n";
                            }
                        }
                    }
                }
            }

            //cout<<"max overlap found = "<<max_overlap<<"\tFalse detection = "<<false_overlap_flag<<endl;

            if((false_overlap_flag ==0 && max_overlap >= overlap_threshold) ||
                false_overlap_flag ==-1 && max_overlap >= 2*overlap_threshold)
            {
                ret_v[0] = max_overlap;
                ret_v[1] = 0;
            }
            else if(false_overlap_flag ==-1 || max_overlap < overlap_threshold)
            {
                ret_v[0] = 0;
                ret_v[1] = -1;
                partial_saved_read_temp[0] = partial_saved_read_temp[1] = -1;
            }
            return;
        }

        ret_v[0] = 0;
        ret_v[1] = 0;

        return;
    }

    static bool sortcol( const vector<int>& v1,const vector<int>& v2 )
    {
        return v1[0] < v2[0];
    }
  
  int detectFalseOverlapUnm(int **pflag,int gaplen)
  {
    int l_max=-MAX_GAP,r_min=MAX_GAP;
    int left_cross[1000],right_cross[1000];
    int left_c=0,right_c=0;
    int p_count=0,len,placed_pos=0,overlap_count=0,read_covered=0;
    double mismatch_threshold = .1;
    
    for(int i=0;i<this->num_reads_in_gap;i++)
    {
      if(pflag[i][2] != -1)
      {
        //cout<<pflag[i][0]<<"\t"<<pflag[i][1]<<"\t"<<pflag[i][2]<<endl;
        
        len = pflag[i][1];
        int pos = pflag[i][0];
        int j,flag=0,start=-1;
        
        for(j=0;j<len;j++)
        {
          if(pos+j>=0 && pos+j<gaplen)
          {
            if(flag==0)
            {
              flag=1;
              start=pos+j;
            }
          }
          if(pos+j == gaplen)break;
        }
        //A single read that covers the whole gap, can't affect both of these values
        if(pos < 0 || pos+len < gaplen)
        {
          if(pos+j-1 > l_max)l_max = pos+j-1;
        }
        else
        {
          if(start < r_min)
          {
            r_min = start;
            //cout<<p_count<<"\t"<<r_min<<endl;
          }
        }
      }
    }
    
    if(l_max == -MAX_GAP)l_max = -1;
    if(r_min == MAX_GAP)r_min = gaplen;
    
    //cout<<"max index from left side = "<<l_max<<"\t"<<"min index from right side= "<<r_min<<endl;
    
    overlap_count = l_max - r_min + 1;
    
    //cout<<overlap_count<<endl;
    
    for(int k=0;k<this->num_reads_in_gap;k++)
    {
      if(pflag[k][2] != -1)//read used to fill gap
      {
        len = pflag[k][1];
        placed_pos = pflag[k][0];
        
        if((placed_pos < 0 || placed_pos + len <= gaplen) && placed_pos + len -1 >= r_min)left_cross[left_c++] = k;
        if(placed_pos > 0 && placed_pos + len -1 >= gaplen && placed_pos <= l_max)right_cross[right_c++] = k;
      }
    }
    
    int max_overlap=0;
    int max_overlap_index=-1;
    int false_overlap_flag=1;
    int p1,p2;
    
    if(r_min <= l_max)
    {
      char readstring1[MAX_READLENGTH],readstring2[MAX_READLENGTH];
      //cout<<l_max<<"\t"<<r_min<<"\t"<<left_c<<"\t"<<right_c<<endl;
      
      //for(int i=0;i<left_c;i++)cout<<left_cross[i]<<endl;
      //for(int i=0;i<right_c;i++)cout<<right_cross[i]<<endl;
      
      for(int i=0;i<left_c;i++)
      {
        for(int j=0;j<right_c;j++)
        {
          if(left_cross[i] != right_cross[j])
          {
            strcpy(readstring1,reads_gap[left_cross[i]]);
            p1 = pflag[left_cross[i]][0];
            len = pflag[left_cross[i]][1];
            
            strcpy(readstring2,reads_gap[right_cross[j]]);
            p2 = pflag[right_cross[j]][0];
            
            int diff_gap = p1 + len - gaplen;
            
            if(diff_gap > 0)
            {
              overlap_count = (p1 + len -1) - p2 + 1 - diff_gap;
            }
            else
            {
              overlap_count = (p1 + len -1) - p2 + 1;
              diff_gap=0;
            }
            
            //cout<<overlap_count<<"\t"<<left_cross[i]<<"\t"<<right_cross[j]<<"\t"<<p1<<"\t"<<p2<<endl;
            
            if(overlap_count >= 5)
            {
              //cout<<overlap_count<<"\t"<<left_cross[i]<<"\t"<<right_cross[j]<<"\t"<<p1<<"\t"<<p2<<endl;
              std::string common_left,common_right;
              
              //cout<<readstring1<<endl;
              //cout<<readstring2<<endl;
              
              if((p1 < 0 && find_partial_match(gap_left,get_read_frag(readstring1,p1),0,1,-1)) || (p1 > 0 && p1 + len <= gaplen))
              {
                std::string s1(readstring1);
                common_left = s1.substr(s1.size()-overlap_count-diff_gap,overlap_count);
                //cout<<left_cross[i]<<"- cl = "<<common_left<<endl;
              }
              
              if(find_partial_match(gap_right,get_read_frag(readstring2,p2),1,1,-1))
              {
                std::string s2(readstring2);
                common_right = s2.substr(0,overlap_count);
                //cout<<right_cross[j]<<"- cr = "<<common_right<<endl;
              }
              
              int len1 = common_left.size(),len2 = common_right.size();
              
              //cout<<len1<<"\t"<<len2<<endl;
              
              if( len1 > 0 && len2 > 0 && len1 == len2)
              {
                double mismatch_frac = getDiff(common_left.c_str(),common_right.c_str(),len1);
                
                //cout<<"mismatch frac = "<<mismatch_frac<<endl;
                if(mismatch_frac <= mismatch_threshold)
                {
                  //cout<<"Overlap Length = "<<overlap_count<<"\tMismatch% = "<<mismatch_frac<<"\tOverlapped string = "<<common_left<<"\t"<<common_right<<endl;
                  if(len1 > max_overlap)
                  {
                    max_overlap = len1;
                    max_overlap_index = p2;
                    false_overlap_flag = 0;
                  }
                }
              }
            }
          }
        }
      }
    }
    
    //cout<<"max overlap found = "<<max_overlap<<"\tFalse detection = "<<false_overlap_flag<<endl;
    return false_overlap_flag;
  }

    double findOverlapUnmapped(int **pos,int len)
    {
        int correct_overlap=0,incorrect_penalty=0,gap_penalty=0,diff=0,fullclosed=0;
        vector<vector<int> > vec(this->num_reads_in_gap);
        for (int i = 0; i < this->num_reads_in_gap; i++)
        {
            vec[i] = vector<int>(3);
            for (int j = 0; j < 3; j++)
                vec[i][j] = pos[i][j];
        }

        sort(vec.begin(), vec.end(),sortcol);
        /*for(int i=0;i<this->num_reads_in_gap;i++)
        {
            for(int j=0;j<3;j++)
            {
                if(vec[i][0]!=-200)cout<<vec[i][j]<<" ";
            }
            if(vec[i][0]!=-200)cout<<endl;
        }*/

        for(int i=0;i<this->num_reads_in_gap-1;i++)
        {
            if(vec[i][0] != -200)
            {
                diff = (vec[i][0]+ vec[i][1] - vec[i+1][0]);//the num of matched char

                if(diff >= match_count_discont)//Correct overlaps
                {
                    //correct_overlap += diff;
                }
                else if(diff>=0)//0<=diff<match_count_discont---->Incorrect/No overlap
                {
                    //cout<<"Discontinuity detected - "<<vec[i][0]<<" "<<vec[i][1]<<" "<<vec[i+1][0]<<endl;
                    incorrect_penalty +=-250;
                    this->discont_or_not=1;
                }
                else//Gap found between reads
                {
                    //cout<<"gap between read = "<<vec[i][2]<<" and Read = "<<vec[i+1][2]<<"\tDiscarding both"<<endl;
                    gap_penalty += -4*50;
                    mark_accepted_reads[vec[i][2]] = 0;
                    mark_accepted_reads[vec[i+1][2]] = 0;
                    if(this->gapLength == originalGap)
                    {
                        unmapped_read_pos_arr_org[vec[i][2]][0] = -200;
                        unmapped_read_pos_arr_org[vec[i+1][2]][0] = -200;

                        unmapped_read_pos_arr_org[vec[i][2]][1] = 0;
                        unmapped_read_pos_arr_org[vec[i+1][2]][1] = 0;
                    }
                    valid_count-=2;
                }
            }
        }

        int left_right_match_count=0,left_right_match_adv=0;

        for(int i=0;i<this->num_reads_in_gap;i++)
        {
            if(vec[i][0] != -200)
            {
                if(vec[i][0] < 0 && -vec[i][0] >=3 && vec[i][0] + vec[i][1] > 0)left_right_match_count++;
                if(vec[i][0] < this->gapLength && vec[i][0] + vec[i][1] - this->gapLength >=3)
                    left_right_match_count++;
            }
        }

        left_right_match_adv = left_right_match_count*50;

        double return_value = correct_overlap + incorrect_penalty + gap_penalty + left_right_match_adv;
        //cout<<"Total correct overlap = "<<correct_overlap<<"\tIncorrect Penalty = "
          //  <<incorrect_penalty<<"\tGap penalty = "<<gap_penalty<<"\tLeft_right adv = "<<left_right_match_adv<<endl;
        return return_value;
    }


	double placeReads(int ge,int g,int finalize_flag,FILE* draw,int gapoffset,int updateflag)//ge=iteration no.//g=gapNo,
	{
	    FILE * scaffoldMap2=NULL;
	    if(partial_flag==1)scaffoldMap2=fopen(partialgapFileName, "r");
	    
	    char line1[MAX_REC_LEN];
		char line2[MAX_REC_LEN];
	
		int MAX_FILE_READ=MAX_REC_LEN/sizeof(line1[0]);

		char * temp;

		int	pos1,pos2,flag1, flag2, strandNo1, strandNo2, insertSize1, insertSize2;
		char *readString1;
		char *readString2 = new char[this->unmapped_read_len +1];

	    int readLength1, readLength2;
	    int charCode,index;
	    int fromCharCode, toCharCode, k;

	    double maxProb=0,tempProb,maxLikelihood=0,antilog_value=0,z=0,tlb=0;

	    int insertSize,mleInsertSize,tempInsertSize;

	    int maxMatch, tempMatch, maxPos,read_start,read_end;

	    int isRev=0, readIndex;

	    for(long int i=left_maxDistance;i<endPos-startPos-right_maxDistance;i++)
	    {
	        for(int j=0;j<5;j++)
	        {
	            countsGap[i][j]=0;
	        }
	    }

	    int count_of_read = this->num_reads_in_gap;
        //cout<<count_of_read<<endl;

        for(int i=0;i<count_of_read;i++)
        {
            maxlikelihood_value[i]=0;
            mark_accepted_reads[i] = 0;
            final_readpos[i][0] = -200;
            final_readpos[i][1] = 0;//length
            final_readpos[i][2] = -1;//read number

            if(this->gapLength == originalGap)
            {
                unmapped_read_pos_arr_org[i][0]=-200;
                unmapped_read_pos_arr_org[i][1]=0;
                unmapped_read_pos_arr_org[i][2]=0;
            }
        }
        for(int i=0;i<this->partial_read_count;i++)
        {
            maxlikelihood_value[i] = 0;
        }

        //place partial read
        if(partial_flag)
        {
            //cout<<"Using Partial reads..\n";
            count_of_read = 0;
            char* token;
            int ref_pos;

            while(fgets(line1, MAX_FILE_READ, scaffoldMap2)!=NULL)
            {
                token = strtok(line1,"\t");
                readString1 = token;
                readLength1 = strlen(readString1);
                token = strtok(NULL, "\t");
                int clipped_index = atoi(token);

                token = strtok(NULL, "\t");

                flag1 = atoi(token);

                if(flag1 == 1 || flag1 == 2)strandNo1=0;
                else strandNo1 = 1;

                pos1 = atoi(strtok(NULL, "\t"));
                strtok(NULL, "\t");
                ref_pos = atoi(strtok(NULL, "\t"));

                //zxcpartial
                isRev=0;
                maxProb=-DBL_MAX;
                maxMatch=0;
                maxPos=0;

                maxProb=-DBL_MAX;
                mleInsertSize=0;

                if(flag1 == 1 || flag1 == 4){read_start = clip_thresh;read_end = 0;}
                else {read_start = 0;read_end = clip_thresh;}

                //cout<<"Read = "<<count_of_read<<"\tPos1 = "<<pos1<<"\tGapstart = "<<gapStart<<endl;
                count_of_read++;
                if(count_of_read > partial_limit)break;

                if(pos1<gapStart)
                {
                    insertSize=gapStart- ref_pos +readLength1;

                    for(long int i=gapStart-readLength1+1;i<gapStart;i++)//zxcfp1
                    {
                        tempInsertSize=insertSize+i-gapStart;
                        
                        if((tempInsertSize<insertThresholdMin || tempInsertSize > insertThresholdMax) && ref_pos !=-1)
                        {
                            continue;
                        }

                        //if(ref_pos !=-1)tempProb = insertLengthDistSmoothed[tempInsertSize];
                        //else tempProb = 1;
                        
                        tempProb = 1;
                        
                        for(int j=read_start;j<readLength1-read_end;j++)
                        {
                            charCode=charCodes[readString1[j]];
                            index=i+j-startPos;

                            
                            if(isRev==1)
                            {
                                readIndex=readLength1-1-j;
                            }
                            else
                            {
                                readIndex=j;
                            }
                            
                            if(index<0 || readIndex < 0)continue;
                            if(charCode<4)
                            {
                                tempProb *= (probsGap[index][charCode]*(1-errorPosDist[readIndex])+errorPosDist[readIndex]*errorProbsGap[index][charCode]);
                            }
                            else
                            {
                                tempProb *= (errorPosDist[readIndex]*errorProbsGap[index][4]);
                            }

                        }

                        tempProb = log(tempProb);

                        if(tempProb>maxProb)
                        {
                            maxMatch=tempMatch;//Not used here
                            maxProb=tempProb;
                            mleInsertSize=insertSize+i-gapStart;
                            maxPos=i-startPos;
                        }

                        antilog_value= pow(10,tempProb);

                        for(int j=0;j<readLength1;j++)
                        {
                            if(i-startPos+j>=left_maxDistance && i-startPos+j<left_maxDistance+gapLength)
                            {
                                countsGap[i-startPos+j][charCodes[readString1[j]]] += antilog_value;
                            }
                        }
                    }//for end
                }
                else
                {
                    if(ref_pos != -1)ref_pos += gapoffset;
                    insertSize=ref_pos-(gapStart+gapLength)+readLength1;

                    for(long int i=gapStart+gapLength-readLength1+1;i<gapStart+gapLength;i++)//zxcfp2
                    {
                       
                        tempInsertSize=insertSize+gapStart+gapLength-i;
                        
                        if((tempInsertSize<insertThresholdMin || tempInsertSize > insertThresholdMax) && ref_pos !=-1)
                        {
                            continue;
                        }
                     
                        //if(ref_pos !=-1)tempProb = insertLengthDistSmoothed[tempInsertSize];
                        //else tempProb = 1;
                        
                        tempProb = 1;

                        for(int j=read_start;j<readLength1-read_end;j++)
                        {
                            charCode=charCodes[readString1[j]];
                            index=i+j-startPos;

                            if(isRev==1)
                            {
                                readIndex=readLength1-1-j;
                            }
                            else
                            {
                                readIndex=j;
                            }

                            if(index<0 || readIndex < 0){continue;}

                            if(charCode<4)
                            {
                                tempProb *= (probsGap[index][charCode]*(1-errorPosDist[readIndex])+errorPosDist[readIndex]*errorProbsGap[index][charCode]);
                            }
                            else
                            {
                                tempProb *= (errorPosDist[readIndex]*errorProbsGap[index][4]);
                            }
                        }

                        tempProb = log(tempProb);

                        if(tempProb>maxProb)
                        {
                            maxMatch=tempMatch;
                            maxProb=tempProb;
                            mleInsertSize=insertSize+gapStart+gapLength-i;
                            maxPos=i-startPos;
                        }

                        antilog_value = pow(10,tempProb);

                        for(int j=0;j<readLength1;j++)
                        {
                            if(i-startPos+j>=left_maxDistance && i-startPos+j<left_maxDistance+gapLength)
                            {
                                countsGap[i-startPos+j][charCodes[readString1[j]]]+= antilog_value;
                            }
                        }                 
                    }//for end
                }

                if(maxProb != -DBL_MAX)
                {
                    maxLikelihood += (maxProb);
                }

            }
            fclose(scaffoldMap2);
            
            //==========================Calculate maxlikelihood based on finalize?=====================
            if(1)
            {
                int **pflag;
                pflag = new int*[this->partial_read_count];

                for(int i=0;i<this->partial_read_count;i++)
                {
                    pflag[i] = new int[2];
                    pflag[i][0] = 1;
                    pflag[i][1] = 0;

                    if(this->gapLength == originalGap)
                    {
                        partial_read_pos_arr_org[i][0] = 0;
                        partial_read_pos_arr_org[i][1] = -200;
                        partial_read_pos_arr_org[i][2] = 0;
                    }
                }

                char* new_bestString = new char[gapLength+1];
                new_bestString[0]='\0';
                computeSequence(0,0);
                strcpy(new_bestString,concensus);

                //cout<<new_bestString<<endl;

                int read_start,read_end;
                //cout<<left_maxDistance<<" "<<right_maxDistance<<endl;

                for(long int i=left_maxDistance;i<endPos-startPos-right_maxDistance;i++)
                {
                    for(int j=0;j<5;j++)
                    {
                        new_counts_gap[i][j]=0;
                    }
                }

                int maxSize=gapLength+(left_maxDistance+right_maxDistance);

                char *gapString=new char[maxSize+1];

                for(long int i=0;i<left_maxDistance;i++)
                {
                    gapString[i]=contigs[contigNo][i+startPos];

                }
                for(long int i=left_maxDistance;i<left_maxDistance+gapLength;i++)
                {
                    gapString[i]=new_bestString[i-left_maxDistance];

                }
                for(long int i=left_maxDistance+gapLength;i<left_maxDistance+gapLength+right_maxDistance;i++)
                {
                    gapString[i]=contigs[contigNo][i+startPos-gapLength+originalGap];

                }
                gapString[maxSize]='\0';

                //cout<<gapString<<endl;

                scaffoldMap2=fopen(partialgapFileName, "r");

                count_of_read = -1;

                while(fgets(line1, MAX_FILE_READ, scaffoldMap2)!=NULL)
                {
                    isRev=0;
                    maxProb=-DBL_MAX;
                    maxMatch=0;
                    maxPos=0;
                    mleInsertSize=0;

                    token = strtok(line1,"\t");
                    readString1 = token;
                    readLength1 = strlen(readString1);
                    int clipped_index = atoi(strtok(NULL, "\t"));

                    token = strtok(NULL, "\t");
                    flag1 = atoi(token);

                    if(flag1 == 1 || flag1 == 2)strandNo1=0;
                    else strandNo1 = 1;

                    pos1 = atoi(strtok(NULL, "\t"));

                    strtok(NULL, "\t");
                    ref_pos = atoi(strtok(NULL, "\t"));
                    //cout<<readString1<<" "<<flag1<<" "<<pos1<<" "<<strandNo1<<" "<<readLength1<<endl;
                    //cout<<strandNo1<<endl;

                    strandNo1=0;

                    if(flag1 == 1 || flag1 == 4){read_start = clip_thresh;read_end = 0;}
                    else {read_start = 0;read_end = clip_thresh;}

                    count_of_read++;
                    if(count_of_read>(partial_limit-1))break;
                    //cout<<"============>Read = "<<count_of_read<<endl;

                    if(pos1<gapStart)
                    {
                        insertSize = gapStart- ref_pos +readLength1;
                        
                        for(long int i=gapStart-readLength1+1;i<gapStart;i++)//zxcfp3
                        {
                            tempInsertSize=insertSize+i-gapStart;
                          
                            if((tempInsertSize<insertThresholdMin || tempInsertSize > insertThresholdMax) && ref_pos !=-1)
                            {
                                continue;
                            }
                     
                            tempProb = 1;
                            
                            tempMatch=0;

                            for(int j=read_start;j<readLength1-read_end;j++)
                            {
                                toCharCode=charCodes[readString1[j]];
                                index=i+j-startPos;
                                fromCharCode=charCodes[gapString[index]];

                                if(strandNo1==0)
                                {
                                    k=j;
                                }
                                else
                                {
                                    k=readLength1-j-1;
                                }

                                if(fromCharCode==toCharCode)
                                {
                                    tempProb*=(1-errorPosDist[k]-inPosDist[k]-delPosDist[k]);
                                }
                                else
                                {
                                    if(fromCharCode >=0 && fromCharCode<=4)
                                        tempProb*=errorPosDist[k]*errorTypeProbs[fromCharCode][toCharCode];
                                }
                            }

                            if(tempProb>maxProb)
                            {
                                maxMatch=tempMatch;
                                maxProb=tempProb;
                                mleInsertSize=insertSize+i-gapStart;
                                maxPos=i-startPos;
                            }
                       }
                    }
                    else
                    {
                        if(ref_pos != -1)ref_pos += gapoffset;
                        insertSize=ref_pos-(gapStart+gapLength)+readLength1;

                        for(long int i=gapStart+gapLength-readLength1+1;i<gapStart+gapLength;i++)//zxcfp4
                        {
                            tempInsertSize=insertSize+gapStart+gapLength-i;

                            if((tempInsertSize<insertThresholdMin || tempInsertSize > insertThresholdMax) && ref_pos !=-1)
                            {
                                continue;
                            }
                         
                            tempProb = 1;
                           
                            tempMatch=0;

                            for(int j=read_start;j<readLength1-read_end;j++)
                            {
                                toCharCode=charCodes[readString1[j]];
                                index=i+j-startPos;
                                fromCharCode=charCodes[gapString[index]];

                                if(strandNo1==0)
                                {
                                    k=j;
                                }
                                else
                                {
                                    k=readLength1-j-1;
                                }

                                if(fromCharCode==toCharCode)
                                {
                                    tempProb*=(1-errorPosDist[k]-inPosDist[k]-delPosDist[k]);
                                    
                                }
                                else
                                {
                                    if(fromCharCode >=0 && fromCharCode<=4)
                                        tempProb*=errorPosDist[k]*errorTypeProbs[fromCharCode][toCharCode]; 
                                }
                            }
                            if(tempProb>maxProb)
                            {
                                maxMatch=tempMatch;
                                maxProb=tempProb;
                                mleInsertSize=insertSize+gapStart+gapLength-i;
                                maxPos=i-startPos;
                            }
                        }
                    }

                    double temp_log_val = -log10(maxProb);
                    
                    if(temp_log_val < gapProbCutOff)
                    {
                        valid_count++;
                        pflag[count_of_read][1] = maxPos - left_maxDistance;
                        if(this->gapLength == originalGap)
                        {
                            partial_read_pos_arr_org[count_of_read][0] = 1;
                            partial_read_pos_arr_org[count_of_read][1] = maxPos - left_maxDistance;
                            partial_read_pos_arr_org[count_of_read][2] = readLength1;
                        }
                    }
                    else
                    {
                        invalid_count++;
                        pflag[count_of_read][0] = 0;
                    }
                }

                fclose(scaffoldMap2);

                delete [] gapString;
                delete [] new_bestString;

                int ret_val[2];
                ret_val[0] = 0;
                ret_val[1] = 0;

                detect_overlap_gapestimate(pflag,gapLength,ret_val,8);
                //cout<<ret_val[0]<<"\t"<<ret_val[1]<<endl;
                
                if(ret_val[0]==300)
                {
                    //cout<<"Perfect read returned\n";
                    maxLikelihood += ret_val[0];
                }
                else if(ret_val[0] >=1 && ret_val[0] < MAX_READLENGTH)
                {
                    //cout<<"Overlap returned = "<<ret_val[0]<<endl;
                    maxLikelihood +=30*ret_val[0];

                }
                else if(ret_val[1] == -1)
                {
                    //cout<<"Wrong overlap returned"<<endl;
                    maxLikelihood +=-100;//wrong overlap
                }

                for(int i=0;i<this->partial_read_count;i++)delete[] pflag[i];
                delete[] pflag;
            }
        }
        //==================//partial done============================

        //Placement of unmapped begins==================================
        count_of_read = this->num_reads_in_gap;

        if(unmapped==1)//zxcfu
        {
            for(int read_itr= 0;read_itr<count_of_read;read_itr++)
            {
                maxProb=-DBL_MAX;
                tlb=0;
                tempProb=0;
                maxMatch=0;
                maxPos=0;
                mleInsertSize=0;
                z=0;

                pos1 = pos_reads[read_itr];
                readLength2 = strlen(reads_gap[read_itr]);
                strcpy(readString2,reads_gap[read_itr]);

                if(pos1<gapStart)
                {
                    insertSize=gapStart-pos1+readLength2;

                    for(long int i=gapStart-readLength2+1;i<gapStart+gapLength;i++)
                    {
                        tempInsertSize=insertSize+i-gapStart;

                        if(tempInsertSize<(insertThresholdMin) || tempInsertSize > (insertThresholdMax))
                        {
                            continue;
                        }

                        tempProb = insertLengthDistSmoothed[tempInsertSize];

                        tempMatch=0;

                        for(int j=0;j<readLength2;j++)
                        {
                            charCode=charCodes[readString2[j]];
                            index=i+j-startPos;


                            if(isReverse[read_itr]==1)
                            {
                                readIndex=readLength2-1-j;
                            }
                            else
                            {
                                readIndex=j;
                            }

                            if(index<0 || readIndex < 0)continue;

                            if(charCode<4)
                            {
                                tempProb *=(probsGap[index][charCode]*(1-errorPosDist[readIndex])+errorPosDist[readIndex]*errorProbsGap[index][charCode]);
                            }
                            else
                            {
                                tempProb *= (errorPosDist[readIndex]*errorProbsGap[index][4]);
                            }
                            
                        }

                        tempProb = log10(tempProb);

                        if(tempProb > maxProb)
                        {
                            maxMatch=tempMatch;
                            maxProb=tempProb;
                            mleInsertSize=insertSize+i-gapStart;
                            maxPos=i-startPos;
                        }
                        
                        antilog_value = exp(0.5*tempProb);

                        for(int j=0;j<readLength2;j++)
                        {
                            if(i-startPos+j>=left_maxDistance && i-startPos+j<left_maxDistance+gapLength)
                            {
                                
                                int index= charCodes[readString2[j]];
                                countsGap[i-startPos+j][index] += antilog_value;
                            }
                        }

                    }//for end
                }
                else
                {
                    pos1 += gapoffset;
                    insertSize=pos1-(gapStart+gapLength)+readLength2;

                    for(long int i=gapStart-readLength2+1;i<gapStart+gapLength;i++)
                    {
                        tempInsertSize=insertSize+gapStart+gapLength-i;

                        if(tempInsertSize<insertThresholdMin || tempInsertSize > insertThresholdMax)
                        {
                            continue;
                        }

                        tempProb = insertLengthDistSmoothed[tempInsertSize];
                        tempMatch=0;


                        for(int j=0;j<readLength2;j++)
                        {
                            charCode=charCodes[readString2[j]];
                            index=i+j-startPos;

                            if(isReverse[read_itr]==1)
                            {
                                readIndex=readLength2-1-j;
                            }
                            else
                            {
                                readIndex=j;
                            }

                            if(index<0 || readIndex < 0)continue;
                            if(charCode<4)
                            {
                                tempProb *=(probsGap[index][charCode]*(1-errorPosDist[readIndex])+errorPosDist[readIndex]*errorProbsGap[index][charCode]);
                            }
                            else
                            {
                                tempProb *= (errorPosDist[readIndex]*errorProbsGap[index][4]);
                            }
                        }

                        tempProb = log10(tempProb);

                        if(tempProb >maxProb)
                        {
                            maxMatch=tempMatch;
                            maxProb=tempProb;
                            mleInsertSize=insertSize+gapStart+gapLength-i;
                            maxPos=i-startPos;
                        }

                        antilog_value = exp(0.5*tempProb);
                        for(int j=0;j<readLength2;j++)
                        {
                            if(i-startPos+j>=left_maxDistance && i-startPos+j<left_maxDistance+gapLength)
                            {
                                int index= charCodes[readString2[j]];
                                countsGap[i-startPos+j][index] += antilog_value;//antilog
                            }
                        }
                    }//for end
                }

                if(maxProb != -DBL_MAX)
                {
                    maxlikelihood_value[read_itr] = maxProb;
                }
                else
                {
                    maxlikelihood_value[read_itr] = 0;
                    invalid_count++;
                }
            }

            //============================Unmapped Placement done==========================================

            
            char* new_bestString = new char[gapLength+1];//cfu
            new_bestString[0]='\0';

            computeSequence(0,0);
            strcpy(new_bestString,concensus);

            //cout<<ge<<" - "<<new_bestString<<endl;

            for(long int i=left_maxDistance;i<endPos-startPos-right_maxDistance;i++)
            {
                for(int j=0;j<5;j++)
                {
                    new_counts_gap[i][j]=0;
                }
            }

            int maxSize=gapLength+(left_maxDistance+right_maxDistance);

            char *gapString=new char[maxSize+1];

            for(long int i=0;i<left_maxDistance;i++)
            {
                gapString[i]=contigs[contigNo][i+startPos];

            }
            for(long int i=left_maxDistance;i<left_maxDistance+gapLength;i++)
            {
                gapString[i]=new_bestString[i-left_maxDistance];

            }
            for(long int i=left_maxDistance+gapLength;i<left_maxDistance+gapLength+right_maxDistance;i++)
            {
                gapString[i]=contigs[contigNo][i+startPos-gapLength+originalGap];

            }
            gapString[maxSize]='\0';
            //if(ge == num_itr -1)cout<<gapString<<endl;

            for(int read_itr= 0;read_itr<count_of_read;read_itr++)
            {
                //if(-maxlikelihood_value[read_itr] < 15 * gapProbCutOff)
                {
                    maxProb=0;
                    maxMatch=0;
                    maxPos=0;
                    maxProb=-DBL_MAX;
                    mleInsertSize = 0;

                    pos1 = pos_reads[read_itr];
                    strcpy(readString2,reads_gap[read_itr]);
                    readLength2 = strlen(readString2);

                    if(pos1<gapStart)
                    {
                        //if(gapno == 0)cout<<"Pos1 = "<<pos1<<", Gapstart = "<<gapStart<<endl;
                        insertSize = gapStart-pos1+readLength2;
                        for(long int i=gapStart-readLength2+1;i<gapStart+gapLength;i++)
                        {
                            tempInsertSize=insertSize+i-gapStart;
                            if(tempInsertSize<insertThresholdMin || tempInsertSize > insertThresholdMax)
                             {
                                continue;
                             }

                            tempProb=1;
                            tempMatch=0;

                            for(int j=0;j<readLength2;j++)
                            {
                                toCharCode=charCodes[readString2[j]];
                                index=i+j-startPos;
                                fromCharCode=charCodes[gapString[index]];


                                if(isReverse[read_itr]==0)
                                {
                                    k=j;
                                }
                                else
                                {
                                    k=readLength2-j-1;
                                }

                                if(fromCharCode==toCharCode)
                                {
                                    tempProb*=(1-errorPosDist[k]-inPosDist[k]-delPosDist[k]);
                                }
                                else
                                {
                                    if(fromCharCode >=0 && fromCharCode<=4)tempProb*=errorPosDist[k]*errorTypeProbs[fromCharCode][toCharCode];
                                }
                            }

                            if(tempProb>maxProb)
                            {
                                maxMatch=tempMatch;
                                maxProb=tempProb;
                                mleInsertSize=insertSize+i-gapStart;
                                maxPos=i-startPos;
                            }
                        }
                    }
                    else
                    {
                        pos1 += gapoffset;
                        insertSize=pos1-(gapStart+gapLength)+readLength2;

                        for(long int i=gapStart-readLength2+1;i<gapStart+gapLength;i++)
                        {
                            tempInsertSize=insertSize+gapStart+gapLength-i;
                            if(tempInsertSize<insertThresholdMin || tempInsertSize > insertThresholdMax)
                            {
                                continue;
                            }

                            tempProb=1;

                            tempMatch=0;

                            for(int j=0;j<readLength2;j++)
                            {
                                toCharCode=charCodes[readString2[j]];
                                index=i+j-startPos;
                                fromCharCode=charCodes[gapString[index]];

                                if(isReverse[read_itr]==0)
                                {
                                    k=j;
                                }
                                else
                                {
                                    k=readLength2-j-1;
                                }

                                if(fromCharCode==toCharCode)
                                {
                                    tempProb*=(1-errorPosDist[k]-inPosDist[k]-delPosDist[k]);
                                }
                                else
                                {
                                    if(fromCharCode >=0 && fromCharCode<=4)tempProb*=errorPosDist[k]*errorTypeProbs[fromCharCode][toCharCode];
                                }
                            }

                            if(tempProb>maxProb)
                            {
                                maxMatch=tempMatch;
                                maxProb=tempProb;
                                mleInsertSize=insertSize+gapStart+gapLength-i;
                                maxPos=i-startPos;
                            }
                        }
                    }

                    double temp_log_val = -log10(maxProb);

                    //if(read_itr==79)cout<<maxPos<<"\t"<<temp_log_val<<endl;

                   if(temp_log_val < gapProbCutOff)//zxcfu
                   {
                       //maxLikelihood += -temp_log_val;
                        maxlikelihood_value[read_itr] = -temp_log_val;
                        maxLikelihood += maxlikelihood_value[read_itr];
                        valid_count++;
                   }
                   else
                   {
                        maxLikelihood += -50;//-50;
                   }


                   if(temp_log_val < gapProbCutOff)
                   {
                       int flag=0;
                       for(int j=0;j<readLength2;j++)
                       {
                           if(maxPos+j>=left_maxDistance && maxPos+j<left_maxDistance+gapLength)
                           {
                              new_counts_gap[maxPos+j][charCodes[readString2[j]]]+=1;
                              if(flag==0)
                              {
                                  flag=1;

                              }
                           }
                       }
                       mark_accepted_reads[read_itr]=1;
                       final_readpos[read_itr][0] = maxPos-left_maxDistance;
                       final_readpos[read_itr][1] = readLength2;
                       final_readpos[read_itr][2] = read_itr;//length

                       if(this->gapLength == originalGap)
                       {
                            unmapped_read_pos_arr_org[read_itr][0] = maxPos-left_maxDistance;
                            unmapped_read_pos_arr_org[read_itr][1] = readLength2;
                            unmapped_read_pos_arr_org[read_itr][2] = read_itr;
                            //cout<<read_itr<<" - "<<(maxPos-left_maxDistance)<<endl;
                            this->fill_original=1;
                       }
                       
                       if(originalGap <=30)
                       {
                           int val = final_readpos[read_itr][0] + final_readpos[read_itr][1] -this->gapLength;
                           if(final_readpos[read_itr][0] < 0 && val > 0)
                           {
                                if(-final_readpos[read_itr][0] > 3 && val > 3)ucoverf=1;
                           }
                           
                           if(final_readpos[read_itr][0] < 0 && final_readpos[read_itr][0] + final_readpos[read_itr][1] > 0)
                           {
                               if(-final_readpos[read_itr][0] > 3)umaxleftf=1;
                           }
                           
                           if(final_readpos[read_itr][0] > 0 && final_readpos[read_itr][0] < this->gapLength && val > 0)
                           {
                               if(val > 3)umaxrightf=1;
                           }
                       }
                   }
                }
            }

            computeSequence(1,1);
            strcpy(current_str,concensus);
            
            if(strcmp(current_str,previous_str)==0)
            {
                comp_count++;
            }
            else
            {
                strcpy(previous_str,current_str);
                comp_count=0;
            }
            //cout<<concensus<<endl;

            if(finalize_flag)//zxcfu
            {
                //Modify likelihood by ignoring reads based on coverage
                // ==================================================

                int region[50000];
                int region_start=0;
                int region_count=0;
                int region_left_right[2]={MAX_GAP,-MAX_GAP};

                for(int i=0;i<gapLength;i++)
                {
                    if(gap_coverage[i]< covergae_threshold2 && region_start ==0)
                    {
                        region[region_count] = i;
                        region_start=1;
                        //cout<<"starting here for "<<i<<endl;
                    }
                    else if(gap_coverage[i] >= covergae_threshold2 && region_start ==1)
                    {
                        if(i - 1 - region[region_count] >= 10)
                        {
                            region_start = 0;
                            region[region_count+1] = i-1;
                            region_count += 2;
                        }
                        //cout<<"ending here for "<<i<<endl;
                    }

                    if(i == gapLength-1 &&  region_start == 1)
                    {
                        if(i - region[region_count] >= 10)
                        {
                            region_start = 0;
                            region[region_count+1] = i;
                            region_count += 2;
                        }
                    }
                    //if(ge == num_itr-1)cout<<i<<" = "<<gap_coverage[i]<<endl;
                }

                if(region_count)
                {
                    region_left_right[0] = region[0];
                    region_left_right[1] = region[region_count-1];
                    region_perct = (region_left_right[1] - region_left_right[0]*1.0)/gapLength;
                }
                else region_perct = 0;

                if(0)
                {
                    //for(int i=0;i<gapLength;i++)
                      //  cout<<gap_coverage[i]<<" ";
                    //cout<<endl;
                    cout<<region_count<<endl;
                    for(int j=0;j<region_count;j+=2)
                    {
                        cout<<"gaplen = "<<gapLength<<"\tregion start = "<<region[j]<<"\tregion end = "<<region[j+1]<<endl;
                    }
                    cout<<"Region left min = "<<region_left_right[0]<<"\t"<<"region right max = "<<region_left_right[1]<<endl;
                    cout<<"Diff = "<<(region_left_right[1] - region_left_right[0])<<"\tperct. = "<<(region_left_right[1] - region_left_right[0]*1.0)/gapLength<<endl;
                }

                for(int i=0;i<count_of_read;i++)
                {
                    if(mark_accepted_reads[i] == 1)
                    {
                        if(final_readpos[i][0] >= region_left_right[0] &&
                            final_readpos[i][0] + final_readpos[i][1] -1 < region_left_right[1])
                        {
                            maxLikelihood += -50;
                            //cout<<"Discarding read - "<<i<<"\t- val - "<<maxlikelihood_value[i]<<endl;
                            mark_accepted_reads[i] = 0;
                            if(this->gapLength == originalGap)
                            {
                                unmapped_read_pos_arr_org[i][0] = -200;
                                unmapped_read_pos_arr_org[i][1] = 0;
                            }
                            valid_count--;
                            final_readpos[i][0] = -200;
                            final_readpos[i][1] = 0;
                        }
                    }
                }

                

                //=====================Modify Likelihood based on overlap in unmapped too=========================

                maxLikelihood += findOverlapUnmapped(final_readpos,gapLength);

                //===============Update counsgap based on live partial reads for subsequent EM iterations==============
                
                //cout<<concensus<<endl;
                //cout<<comp_count<<endl;

                ///zxcfc
                int condition;

                condition = unmapped && comp_count >= 1 && region_perct !=0 && ge != num_itr-1;
                

                if(condition && updateflag)
                {
                    //cout<<"In update\t"<<comp_count<<endl;
                    char text_left[MAX_READLENGTH],text_right[MAX_READLENGTH],pattern[MAX_READLENGTH];
                    int index_pair[1000];
                    for(int i=0;i<1000;i++)index_pair[i]=-1;;
                    int start_n=0,pair_count=0,start_index1,end_index1;
                    char **str_arr;
                    int N_count=0,num_match_count0=0,num_match_count1=1;
                    char temp_read_str[MAX_READLENGTH];
                    char temp_text_right[MAX_READLENGTH];

                    int match_threshold=this->unmapped_read_len * 0.25;
                    int max_segment_length=this->unmapped_read_len * 0.67;
                    //will generate a max_segment_length-1 sized string from both sides of gap

                    int min_gap_len = this->unmapped_read_len/2 +1;

                    int checkinsertsize=1;

                    double count_pos[gapLength][4];
                    for(int j=0;j<gapLength;j++)
                    {
                        for(int k=0;k<4;k++)
                            count_pos[j][k] = 0;
                    }

                    for(int i=0;i<strlen(concensus);i++)
                    {
                        if(concensus[i] == 'N' && start_n==0)
                        {
                            start_n=1;
                            if(i>0)index_pair[pair_count++] = i-1;
                            else index_pair[pair_count++] = i;
                            N_count++;

                        }
                        else if(concensus[i] != 'N' && start_n==1)
                        {
                            start_n =0;
                            index_pair[pair_count++] = i;
                            if(N_count<min_gap_len)
                            {
                                pair_count -= 2;
                            }
                            N_count=0;
                        }
                        else if(concensus[i] == 'N' && start_n==1)N_count++;
                        if(i == strlen(concensus)-1 && start_n ==1)
                        {
                            index_pair[pair_count++] = i;
                            if(N_count<min_gap_len)
                            {
                                pair_count -= 2;
                            }
                        }
                    }

                    
                    if(pair_count > 2)
                    {
                        index_pair[1] = index_pair[pair_count-1];
                    }
                    
                    //cout<<pair_count<<"\t"<<index_pair[0]<<"\t"<<index_pair[1]<<endl;

                    int flag1=1,flag2=1;
                    if(pair_count<2){flag1=0;flag2=0;}

                    if(flag1==1 || flag2==1)
                    {
                        int i=0;
                        str_arr = new char*[2];
                        //cout<<"pair count = "<<pair_count<<endl;
                        end_index1 = index_pair[i++];
                        start_index1 = index_pair[i];

                        int index_s = -1,stop_index = -1,m=0,n=0;
                        if(end_index1 >= max_segment_length)index_s = end_index1 - max_segment_length+1;
                        else index_s=0;
                        
                        //cout<<"start point = "<<index_s<<endl;

                        for(int j=index_s;concensus[j] != 'N';j++)
                        {
                            text_left[j-index_s] = concensus[j];
                            m++;
                        }
                        text_left[m] = '\0';

                        if(start_index1+max_segment_length <= strlen(concensus))stop_index = start_index1+max_segment_length-1;
                        else stop_index = strlen(concensus)-1;

                        for(int j=start_index1;j<=stop_index;j++)
                        {
                            text_right[j-start_index1] = concensus[j];
                            n++;
                        }
                        text_right[n] = '\0';
                        
                        //cout<<"end point = "<<stop_index<<endl;

                        str_arr[i-1] = new char[strlen(text_left)+1];
                        strcpy(str_arr[i-1],text_left);
                        
                        //cout<<"left = "<<str_arr[i-1]<<"\tlen = "<<strlen(str_arr[i-1])<<endl;
                                                
                        str_arr[i] = new char[strlen(text_right)+1];
                        strcpy(str_arr[i],text_right);
                        
                        //cout<<"right = "<<str_arr[i]<<"\tlen = "<<strlen(str_arr[i])<<endl;

                        int temp_flag1 = flag1,temp_flag2 = flag2;

                        for(int read_itr= 0;read_itr<count_of_read;read_itr++)
                        {
                            //if(read_itr ==371)cout<<maxlikelihood_value[read_itr]<<"\t"<<mark_accepted_reads[read_itr]<<"\t"<<strand_flag[read_itr]<<endl;
                            if(mark_accepted_reads[read_itr]==0 && maxlikelihood_value[read_itr]!=0)
                            {
                                int R2 = strlen(reads_gap[read_itr]);
                                int placed_pos;
                                int mapped_pos= pos_reads[read_itr];
                                int insertsize1=-1,insertsize2=-1,insertsize=-1;

                                flag1 = temp_flag1;
                                flag2 = temp_flag2;
                                char ifelse=-1;

                                if(checkinsertsize)
                                {
                                    for(int g=0;g<2;g++)
                                    {
                                        placed_pos = gapStart+ index_pair[g];
                                        if(mapped_pos < gapStart)
                                        {
                                            if(g==0)insertsize1 = placed_pos + R2 - mapped_pos;
                                            else insertsize2 = placed_pos + R2 - mapped_pos;
                                            insertsize = placed_pos + R2 - mapped_pos;
                                            ifelse = 'I';
                                        }
                                        else
                                        {
                                            
                                            if(g==0)insertsize1 = mapped_pos - placed_pos + R2;
                                            else insertsize2 = mapped_pos - placed_pos + R2;
                                            insertsize = mapped_pos - placed_pos + R2;
                                            ifelse = 'E';
                                        }

                                        if(insertsize <(insertThresholdMin)+100 || insertsize > (insertThresholdMax-100))
                                        {
                                            if(mapped_pos < gapStart)
                                            {
                                                //cout<<"If\t"<<"Mapped = "<<mapped_pos<<"\tPlaced pos= "<<placed_pos;
                                                //cout<<"\tRejecting read - "<<read_itr<<", due to wrong insertsize = "<<insertsize<<" in update\n";
                                            }
                                            else
                                            {
                                                //cout<<"Else\t"<<"Mapped = "<<mapped_pos<<"\tPlaced pos= "<<placed_pos;
                                                //cout<<"\tRejecting read - "<<read_itr<<", due to wrong insertsize = "<<insertsize<<" in update\n";
                                            }

                                            if(g==0)flag1 = 0;
                                            else flag2 = 0;
                                        }
                                        //else cout<<"ACC = "<<insertsize<<"\tRead = "<<read_itr<<"\tSide = "<<g<<endl;
                                    }
                                }

                                if(flag1 ==1)
                                {
                                    i=0;
                                    //cout<<"Updating for left side"<<endl;
                                    for(int j=0;j<strlen(str_arr[i]);j++)//position of left side
                                    {
                                        int match = 0;
                                        int match_flag=1;
                                        for(int k=0;k<strlen(str_arr[i])-j;k++)
                                        {
                                            if(str_arr[i][j+k] != reads_gap[read_itr][k])
                                            {
                                                match_flag=0;
                                                break;
                                            }
                                            else match++;
                                        }

                                        if(match_flag==1 && match > match_threshold)
                                        {
                                            //cout<<insertsize1<<" "<<insertsize2<<" "<<insertsize<<endl;
                                            //cout<<ifelse<<" "<<i<<" "<<read_itr<<"\t"<<"\tmatch count = "<<match<<"\tlikelihood = "<<maxlikelihood_value[read_itr]<<", Insertsize = "<<insertsize1<<endl;
                                            num_match_count0++;
                                            int len = strlen(reads_gap[read_itr]);
                                            for(int r=0;r<len;r++)
                                            {
                                                int ind_read = charCodes[reads_gap[read_itr][r]];
                                                int ind_ref = index_pair[i]-match+1+r;

                                                if(ind_ref == gapLength)break;
                                                if(ind_ref > index_pair[i])
                                                {
                                                    if(ind_read < 4)count_pos[ind_ref][ind_read] +=match;
                                                    else
                                                    {
                                                        for(int q=0;q<4;q++)count_pos[ind_ref][q] +=match;
                                                    }
                                                }
                                            }
                                            break;
                                        }//if end
                                    }//for end
                                }

                                if(flag2==1)
                                {
                                    i=1;
                                    
                                    //cout<<"Updating for right side"<<endl;
                                    
                                    strcpy(temp_read_str,reads_gap[read_itr]);
                                    reverseStr(temp_read_str);//so reverse one is in temp

                                    strcpy(temp_text_right,str_arr[i]);//original is in str
                                    reverseStr(temp_text_right);//so reverse one is in temp

                                    for(int j=0;j<strlen(temp_text_right);j++)
                                    {
                                        int match = 0;
                                        int match_flag=1;
                                        for(int k=0;k<strlen(temp_text_right)-j;k++)
                                        {
                                            if(temp_text_right[j+k] != temp_read_str[k])
                                            {
                                                match_flag=0;
                                                break;
                                            }
                                            else match++;
                                        }

                                        if(match_flag==1 && match > match_threshold)
                                        {
                                            if(0)
                                            {
                                                cout<<"Original = "<<reads_gap[read_itr]<<endl;
                                                cout<<"Reversed = "<<temp_read_str<<endl;
                                                cout<<"right segment org = "<<str_arr[i]<<endl;
                                                cout<<"right segment reversed = "<<temp_text_right<<endl;
                                            }
                                            //cout<<ifelse<<" "<<i<<" "<<read_itr<<"\t"<<"\tmatch count = "<<match<<"\tlikelihood = "<<maxlikelihood_value[read_itr]<<", Insertsize = "<<insertsize2<<endl;
                                            //cout<<"Updating for right side"<<endl;
                                            num_match_count1++;
                                            int len = strlen(reads_gap[read_itr]);

                                            for(int r=0;r<len;r++)
                                            {
                                                int ind_read = charCodes[reads_gap[read_itr][len-r-1]];
                                                int ind_ref = index_pair[i]+match-1-r;

                                               if(ind_ref <0)break;
                                               if(ind_ref < index_pair[i])
                                               {
                                                    if(ind_read < 4)count_pos[ind_ref][ind_read] +=match;
                                                    else
                                                    {
                                                        for(int q=0;q<4;q++)count_pos[ind_ref][q] +=match;
                                                    }
                                                }
                                            }
                                            break;
                                        }
                                    }
                                }//else end
                            }//if end
                        }//for end of reads

                        //cout<<"Iteration = "<<ge<<"\tLeft part match = "<<num_match_count0<<
                        //"\tRight part match = "<<num_match_count1<<"\tpair count = "<<pair_count<<endl;

                        if(flag1 ==1 && num_match_count0 ==0 && index_pair[0] < strlen(partial_left))
                        {
                            //cout<<"Re-Updating using partial left"<<endl;
                            //cout<<partial_left<<endl;
                            //cout<<index_pair[0]<<" "<<strlen(partial_left)<<endl;
                            for(int f=index_pair[0]+1;f<strlen(partial_left);f++)
                            {
                                count_pos[f][charCodes[partial_left[f]]] += 1;
                                //cout<<charCodes[partial_left[f]]<<endl;
                            }
                        }
                        if(flag2 ==1 && num_match_count1 ==0 && index_pair[1] > gapLength - strlen(partial_right))
                        {
                            //cout<<"Re-Updating using partial right"<<endl;
                            //cout<<partial_right<<endl;
                            //cout<<index_pair[2]<<" "<<strlen(partial_right)<<endl;
                            for(int f=gapLength - strlen(partial_right),c=0;f<index_pair[1];f++)
                            {
                                count_pos[f][charCodes[partial_right[c++]]] += 1;
                                //cout<<charCodes[partial_left[f]]<<endl;
                            }
                        }

                        //update countsgap
                        int total_count_p = 0;
                        //cout<<"Update of countsgap done\n";
                        for(int j=0;j<gapLength;j++)
                        {
                            total_count_p =0;

                            for(int k=0;k<4;k++)
                            {
                                total_count_p += count_pos[j][k];
                            }
                            if(total_count_p>0)
                            {

                                for(int k=0;k<4;k++)
                                {
                                    countsGap[left_maxDistance+j][k] = (count_pos[j][k]/total_count_p);
                                }
                            }
                        }

                        if(0)
                        {
                            for(int j=0;j<gapLength;j++)
                            {
                                cout<<"Gap index = "<<j<<" ==> ";
                                for(int k=0;k<4;k++)
                                {
                                    cout<<count_pos[j][k]<<" ";
                                }
                                cout<<endl;
                            }
                        }


                        for(int v=0;v<2;v++)
                        {
                            delete[] str_arr[v];
                        }
                        delete [] str_arr;

                    }

                }
            }//finalize end
            delete [] gapString;
                delete [] new_bestString;
        }
        //========Unmapped done======================================


        delete [] readString2;

	    return maxLikelihood;
	}

	int checkDuplicate(double *arr,int maxval)
	{
	    int flag=0;

	    int *count = new int[maxval+1];
        for(int i=0;i<maxval+1;i++)
	    {
	        count[i] = 0;
	    }
        for(int i = 0; i < 4; i++)
         {
             if((int)arr[i] == maxval)
             {
                 if(count[(int)arr[i]] == 1)
                 {
                    flag=1;
                    break;
                 }
                 else count[(int)arr[i]]++;
             }
         }


        delete[] count;
        //cout<<"flag = "<<flag<<endl;
        return flag;
	}

	void computeSequence(int check,int choice_array)//based on the 2Darray countsGap
	{
	    //printArr(left_maxDistance-1,right_maxDistance+this-> gapLength);
	    double max=0;
	    int maxIndex=0;

	    int index=0;
        index = endPos-startPos-(left_maxDistance+right_maxDistance);

        //cout<<left_maxDistance<<"\t"<<endPos-startPos-right_maxDistance<<"\t"<<index<<endl;

        for(long int i=left_maxDistance;i<endPos-startPos-right_maxDistance;i++)
	    {
	        max=0;
	        maxIndex=-1;

	        int flag=0;//all zero
	        for(int j=0;j<=4;j++)
	        {
	            if(choice_array ==0)
	            {
                    if(countsGap[i][j]>max)
                    {
                        max=countsGap[i][j];
                        maxIndex=j;
                        flag=1;//all elements in the array is not 0
                    }
                }
                else
                {
                    if(new_counts_gap[i][j]>max)
                    {
                        max=new_counts_gap[i][j];
                        maxIndex=j;
                        flag=1;//all elements in the array is not 0
                    }
                }
             }
             //cout<<"Here2"<<endl;
	        ///Check =0 when, this fn is called from  func to compute best string, we don't
	                            //want to put N for that call
	        ///Check =1 when, this fn is called from finalize() func to compute final string, we
            	                            //want to put N for that call

	        int flag2=0;
	        int coverage_flag=1;

	        if(check == 1)
	        {
	            gap_coverage[i-left_maxDistance] = (int)max;
	            if(gap_coverage[i-left_maxDistance] <= covergae_threshold1)coverage_flag=0;
	        }
	        //cout<<"Here3, flag2= "<<flag2<<endl;

	        if(((flag == 1 && flag2 == 0) ||(!check)) && coverage_flag)
            {
                if(maxIndex==0)
                {
                    concensus[i-left_maxDistance]='A';
                }
                else if(maxIndex==1)
                {
                    concensus[i-left_maxDistance]='C';

                }
                else if(maxIndex==2)
                {
                    concensus[i-left_maxDistance]='G';

                }
                else if(maxIndex==3)
                {
                    concensus[i-left_maxDistance]='T';

                }
                else
                {
                    concensus[i-left_maxDistance]='N';
                }
            }
            else
            {
                //cout<<"Countsgap row index = "<<i<<" , flag = "<<flag<<" , flag2 = "<<flag2<<endl;
                concensus[i-left_maxDistance]='N';
            }
	    }

	    //cout<<"Before NULL"<<endl;
	    concensus[index]='\0';

	   //cout<<"Putting Null in "<<(index)<<endl;
	}

	int getConcensus(char *s,int gaplen)
	{
	    int len;
	    len = this-> gapLength;
	    int char_count=0;

	    for(int i=0;i<len;i++)
        {
            s[i]=concensus[i];
        }
        s[len]='\0';
        return len;
    }

    int checkAllZero(double* arr,int flag)
    {
        int non_zero = 0;
        for(int k=0;k<4;k++)
        {
            if(arr[k]>flag)non_zero++;
        }
        if(non_zero == 0)return 0;
        else return 1;
    }

    int check_update(double* arr,int g,int j)
    {
        double max_val=-DBL_MAX,second_max = -DBL_MAX;
        int maxp=-MAX_GAP,sec_p=-MAX_GAP;

        for(int k=0;k<4;k++)
        {
            if(arr[k] > max_val)
            {
                second_max = max_val;
                sec_p = maxp;
                max_val = arr[k];
                maxp = k;
            }
            else if (arr[k] >= second_max)
            {
                second_max = arr[k];
                sec_p=k;
            }
        }

        int diff = max_val - second_max;
        //if(g>=0)cout<<j<<" -\t"<<arr[0]<<"\t"<<arr[1]<<"\t"<<arr[2]<<"\t"<<arr[3]<<"\t"<<maxp<<"\t"<<sec_p<<endl;
        if(diff >= partial_threshold)
        {
            if(max_val > 3 && second_max > 3)
            {
                if(maxp==-MAX_GAP || sec_p==-MAX_GAP)cout<<"Error- "<<j<<"\t"<<arr[0]<<"\t"<<arr[1]<<"\t"<<arr[2]<<"\t"<<arr[3]<<"\t"<<maxp<<"\t"<<sec_p<<endl;

                if(qual_gap[j][maxp] <= qual_gap[j][sec_p])return maxp;
                else return sec_p;
            }
            return 50;
        }
        else
        {
            if(max_val >= 1 && second_max >= 1)
            {
                //if(g==8)cout<<j<<" -\t"<<arr[0]<<"\t"<<arr[1]<<"\t"<<arr[2]<<"\t"<<arr[3]<<"\t"<<qual_gap[j][maxp]<<"\t"<<qual_gap[j][sec_p]<<endl;
              if(maxp==-MAX_GAP || sec_p==-MAX_GAP)cout<<"Error2- "<<j<<"\t"<<arr[0]<<"\t"<<arr[1]<<"\t"<<arr[2]<<"\t"<<arr[3]<<"\t"<<maxp<<"\t"<<sec_p<<endl;

                if(qual_gap[j][maxp] <= qual_gap[j][sec_p])return maxp;
                else return sec_p;
            }
            return -1;
        }
    }

    void clear_countsGap(int a)
    {
        for(int j=0;j<a;j++)
        {
            for(int k=0;k<4;k++)
            {
                countsGap[j+left_maxDistance][k] = 0;
            }
        }
    }

    int findRegion(int region[5000])
    {
        int Nstart=0,region_count=0;

        int len = this->gapLength;
        Nstart=0,region_count=0;

        for(int i=0;i<len;i++)
        {
            if(concensus[i] =='N' && Nstart ==0)
            {
                region[2*region_count] = i;
                Nstart = 1;
            }
            else if(concensus[i] != 'N' && Nstart ==1)
            {
                Nstart=0;
                region[2*region_count+1] = i-1;
                region_count++;
            }
            if(i==len-1 && Nstart==1)
            {
                region[2*region_count+1] = i;
                region_count++;
            }
        }
        return region_count;
    }

    int findDiscontinous(int ** pos,int rv[])
    {
        vector<vector<int> > vec(this->num_reads_in_gap);
        for (int i = 0; i < this->num_reads_in_gap; i++)
        {
            vec[i] = vector<int>(3);
            for (int j = 0; j < 3; j++)
            {
                vec[i][j] = pos[i][j];
            }
        }

        sort(vec.begin(), vec.end(),sortcol);
/*
        cout<<"The vector\n";
        for(int i=0;i<this->num_reads_in_gap;i++)
        {
            for(int j=0;j<3;j++)
            {
                if(vec[i][0]!=-200)cout<<vec[i][j]<<" ";
            }
            if(vec[i][0]!=-200)cout<<endl;
        }
*/
        int diff,dis_pos=-1,discount=0;

        for(int i=0;i<this->num_reads_in_gap-1;i++)
        {
            if(vec[i][0] != -200)
            {
                diff = (vec[i][0]+ vec[i][1] - vec[i+1][0]);//the num of matched char
                if(diff >=0 && diff <= match_count_discont/2)
                {
                    dis_pos =vec[i][0]+ vec[i][1];
                    //if(diff == match_count-1)dis_pos--;
                    rv[discount++] = dis_pos;
                    //rv[1] = diff;
                    //return 1;
                }
                //else cout<<pos[i]<<endl;
            }
        }
        return discount;
    }

    int recheck_sequence(int ** pos)
    {
        int region[5000];
        int region_count = findRegion(region);
        int flag = 0;
        int len = this->gapLength;

        int nec_ret[2000];
        flag = findDiscontinous(pos,nec_ret);

        if(flag > 0)
        {
            //cout<<"Filled discontionously"<<"\tTotal point of discontinuity = "<<flag<<endl;;

            for(int k=0;k<flag;k++)
            {
                //cout<<"Filled discontionously,putting N, at index = "<<nec_ret[k]<<endl;
                concensus[nec_ret[k]]='N';
            }
            region_count = findRegion(region);
        }

        //cout<<region_count<<" "<<region[0]<<" "<<region[1]<<endl;
        double reduction_factor;
        
        if(originalGap<400)reduction_factor=1;
        else if(originalGap<1200)reduction_factor = 1.5;
        else reduction_factor = 2;
        
        int readchar = 30;

        if(region_count<=1)
        {
            if(region_count==1)
            {
                //cout<<region_count<<" "<<region[0]<<" "<<region[1]<<" "<<region_perct_max_gap_estm<<endl;
                //if(fillupflag !=1)
                if(region_perct_max_gap_estm < .75 || flag >0)//so more than 25% is covered
                {
                    //cout<<"Reducting\n";
                    int i,j;
                    for(i=region[0]-1;i>=region[0]-reduction_factor*readchar && i>=0;i--)concensus[i]='N';
                    for(j=region[1]+1;j<=region[1]+reduction_factor*readchar && j<len;j++)concensus[j]='N';
                    if(i<0 && j==len)
                    {
                        return 1;
                    }
                }
            }
        }
        else
        {
            int start = region[0],end = region[2*region_count-1];
            //cout<<region_count<<" "<<start<<" "<<end<<endl;
            //cout<<concensus<<endl;
            for(int j=start;j<end;j++)
            {
                concensus[j]='N';
            }

            //if(region_perct_max_gap_estm < .75 || flag > 0)
            {
                int i,j;
                for(i=start-1;i>start-1-reduction_factor*readchar && i>=0;i--)concensus[i]='N';
                for(j=end+1;j<end+1+reduction_factor*readchar && j<len;j++)concensus[j]='N';
                if(i<0 && j==len)
                {
                    this->gapLength=originalGap;
                    initialize_start_end();
                    
                    return 1;
                }
            }
        }
        return 0;
    }

    static bool sortcol2( const vector<int>& v1,const vector<int>& v2 )
    {
        return v1[1] < v2[1];
    }

    int recheck_partial(int ** pos,int gc)
    {
        int region[5000];
        int discontinous_flag = 0,len = this->gapLength,diff,dis_pos=-1,start=-1,end=-1,rcountleft=0,rcountright=0,readcover=0;
        int region_count = findRegion(region);

        for(int i=0;i<this->partial_read_count;i++)
        {
            //cout<<pos[i][0]<<"\t"<<pos[i][1]<<endl;
            if(pos[i][0] == 1 && pos[i][1] < 0)rcountleft++;
            else if(pos[i][0] == 1 && pos[i][1] > 0) rcountright++;
            if(pos[i][0] == 1 && pos[i][1] < 0 && pos[i][1] + pos[i][2] >= this->gapLength)readcover++;
        }
        //cout<<readcover<<"\t"<<rcountleft<<" "<<rcountright<<"\tregioncount = "<<region_count<<endl;

        if(region_count == 0 && rcountleft > 0 && rcountright > 0 && readcover ==0)
        {
            vector<vector<int> > vec(this->partial_read_count);
            for (int i = 0; i < this->partial_read_count; i++)
            {
                vec[i] = vector<int>(3);
                for (int j = 0; j < 3; j++)
                {
                    vec[i][j] = pos[i][j];
                }
            }

            sort(vec.begin(), vec.end(),sortcol2);
            /*cout<<"The vector\n";
            for(int i=0;i<this->partial_read_count;i++)
            {
                cout<<i<<" - ";
                for(int j=0;j<3;j++)
                {
                    cout<<vec[i][j]<<" ";
                }
                cout<<endl;
            }*/
            for(int i=0;i<this->partial_read_count-1;i++)
            {
                if(vec[i][1] != -200 && vec[i][1] <0 && vec[i+1][1] >0)
                {
                    diff = (vec[i][1]+ vec[i][2] - vec[i+1][1]);//the num of matched char

                    if(diff==0)
                    {
                        discontinous_flag = 1;
                        dis_pos =  vec[i+1][1];
                        //cout<<"Found Disocontinous Seq"<<endl;
                        //cout<<vec[i][1]<<"\t"<<vec[i+1][1]<<"\ndispos = "<<dis_pos<<endl;
                        concensus[dis_pos]='N';
                        region_count = findRegion(region);
                        break;
                    }
                    else if(diff < 0)
                    {
                        //cout<<"\nImpossssssssssssssssssssssssssible = "<<diff<<endl;
                    }
                }
            }
        }

        if(discontinous_flag ==0 && region_count == 0 && (readcover > 0 || (rcountleft==0 || rcountright ==0)))
        {
            //cout<<"Extreme Problem\t No read on one side or diff<0 or reads covering whole gap wrongly\n";
            return -1;
        }

        if(region_count > 1)
        {
            start = region[0],end = region[2*region_count-1];
            //cout<<"\n\n multi region beacuse of N in read\n\n";
            for(int j=start;j<= end;j++)
            {
                for(int k=0;k<4;k++)
                {
                    countsGap[j+left_maxDistance][k] = 0;
                }
            }
        }
        else if(region_count ==1)
        {
            start = region[0],end = region[1];
        }
        else return 1;
        //cout<<concensus<<endl;
        //cout<<start<<" "<<end<<" "<<rcountleft<<" "<<rcountright<<"\tFlag = "<<discontinous_flag<<endl;

        int remove_char=10;//so total 10*2+1 = 21, which will cross the threshold of 20 of perfectGap...
        int min_N = 21;
        if(discontinous_flag ==0 && end - start >= min_N && rcountright >0 && rcountleft >0)return 1;
        else
        {
            //if you are here, it means either discontinous_flag=1...discontinous
            //or region<8 or leftcount==0 or rightcount==0
            //cout<<start<<" "<<end<<" "<<rcountleft<<" "<<rcountright<<"\tFlag = "<<discontinous_flag<<endl;

            if(discontinous_flag==1 || end - start < min_N)
            {
                if(discontinous_flag ==0)
                {
                    int rem_char = min_N - end + start;
                    remove_char = rem_char/2 + rem_char%2;
                }
                for(int j=start;j>start-1-remove_char && j>=0;j--)
                {
                    for(int k=0;k<4;k++)
                    {
                        countsGap[j+left_maxDistance][k] = 0;
                    }
                }

                for(int j=end+1;j<end+1+remove_char && j<len;j++)
                {
                    for(int k=0;k<4;k++)
                    {
                        countsGap[j+left_maxDistance][k] = 0;
                    }
                }
                //cout<<"Removing "<<remove_char<<" characters from each side|| Gap case = "<<gc<<endl;
            }
            return 0;
        }
    }

    void recompute1(int ** position)
    {
        FILE* fp = fopen(partialgapFileName,"r");
        char line1[MAX_REC_LEN];
        char * token,* readString1;
        int MAX_FILE_READ=MAX_REC_LEN/sizeof(line1[0]);
        int readLength1;

        int read_itr = 0;
        while(fgets(line1, MAX_FILE_READ, fp)!=NULL)
        {
            token = strtok(line1,"\t");
            readString1 = token;
            readLength1 = strlen(readString1);

            if(position[read_itr][0] == 1)
            {
                int pos = position[read_itr][1];
                for(int j=0;j<readLength1;j++)
                {
                    if(pos+j>=0 && pos+j<this->gapLength)
                    {
                        countsGap[pos+j+left_maxDistance][charCodes[readString1[j]]] +=1;
                    }
                }
            }
            read_itr++;
            if(read_itr>partial_limit)break;
        }

        fclose(fp);
    }

    void recompute2(int ** position)
    {
        //printArr(left_maxDistance-1,right_maxDistance+this-> gapLength);
        for(int read_itr= 0;read_itr<this->num_reads_in_gap;read_itr++)
        {
            if(position[read_itr][1] > 0)
            {
                int pos = position[read_itr][0];
                //cout<<read_itr<<" "<<pos<<endl;

                for(int j=0;j<strlen(reads_gap[read_itr]);j++)
                {
                    if(pos+j>=0 && pos+j<this->gapLength)
                    {
                        countsGap[pos+j+left_maxDistance][charCodes[reads_gap[read_itr][j]]]+=1;
                    }
                }
            }
        }
    }

	void finalize(int gapLength,int gapno)
    {
        int gapoffset = gapLength - originalGap;
        FILE * scaffoldMap2=NULL;
        if(partial_flag==1)scaffoldMap2=fopen(partialgapFileName, "r");

        partial_read_flag = new int*[this->partial_read_count];

        for(int i=0;i<this->partial_read_count;i++)
        {
            partial_read_flag[i] = new int[3];
            partial_read_flag[i][0] = 0;
            partial_read_flag[i][1] = -200;
            partial_read_flag[i][2] = this->partial_read_len;
        }

        char line1[MAX_REC_LEN];
        char line2[MAX_REC_LEN];
        int MAX_FILE_READ=MAX_REC_LEN/sizeof(line1[0]);

        char * temp,*readString2,*readString1;

        int	pos1,pos2,flag1, flag2, strandNo1, strandNo2, insertSize1, insertSize2;
        int readLength1, readLength2;
        int fromCharCode, toCharCode, index, k;

        double maxProb=0,tempProb,maxLikelihood=0;

        int **unmapped_read_pos_arr;

        unmapped_read_pos_arr = new int*[this->num_reads_in_gap];

        for(int i=0;i<this->num_reads_in_gap;i++)
        {
            unmapped_read_pos_arr[i] = new int[3];
            unmapped_read_pos_arr[i][0]= -200;
            unmapped_read_pos_arr[i][1]=0;
            unmapped_read_pos_arr[i][2]=-1;
        }

        int insertSize,mleInsertSize,tempInsertSize;
        int left_right_check[2]={0,0};
        int leftcount=0,rightcount=0,left_start_zero=0,right_fin_glen=0;
        int maxMatch, tempMatch, maxPos;
        int totalCount=0, discardedCount=0;
        int maxSize=gapLength+(left_maxDistance+right_maxDistance);
        int unmapped_max_left=0,unmapped_max_right=0;
        int newleft=0,newright=0,newmaxleft=0,newmaxright=0;

        char *gapString=new char[maxSize+1];

        int read_itr;
        //cout<<left_maxDistance<<" "<<right_maxDistance<<endl;
        for(long int i=0;i<left_maxDistance;i++)
        {
            gapString[i]=contigs[contigNo][i+startPos];

        }
        for(long int i=left_maxDistance;i<left_maxDistance+gapLength;i++)
        {
            gapString[i]=bestString[i-left_maxDistance];

        }
        for(long int i=left_maxDistance+gapLength;i<left_maxDistance+gapLength+right_maxDistance;i++)
        {
            gapString[i]=contigs[contigNo][i+startPos-gapLength+originalGap];
        }
        gapString[maxSize]='\0';

        int end_lim;
        end_lim=left_maxDistance+alloc_arg;

        for(long int i=left_maxDistance;i<end_lim && i < end_pos_max;i++)
        {
            for(int j=0;j<5;j++)
            {
                countsGap[i][j]=0;
            }
        }

        //Don;tmove this up, otherwise the counsgap won't be cleared to the maximum gapestimate

        this->gapLength=gapLength;
        initialize_start_end();

        if(unmapped)//zxcf fin unm
        {
            draw_read(draw,this->gapLength,"S",0,0,0,gapno,' ');

            for(read_itr= 0;read_itr<this->num_reads_in_gap;read_itr++)
            {
                maxProb=0;
                maxMatch=0;
                maxPos=0;
                totalCount++;

                pos1 = pos_reads[read_itr];
                readLength2 = strlen(reads_gap[read_itr]);

                if(pos1<gapStart)
                {
                    insertSize = gapStart-pos1+readLength2;
                    for(long int i=gapStart-readLength2+1;i<gapStart+gapLength;i++)
                    {
                        tempInsertSize=insertSize+i-gapStart;
                        if(tempInsertSize<insertThresholdMin || tempInsertSize > insertThresholdMax)
                         {
                            continue;
                         }

                        tempProb=1;//insertLengthDistSmoothed[tempInsertSize];
                        tempMatch=0;

                        for(int j=0;j<readLength2;j++)
                        {
                            toCharCode=charCodes[reads_gap[read_itr][j]];
                            index=i+j-startPos;
                            fromCharCode=charCodes[gapString[index]];

                            if(isReverse[read_itr]==0)
                            {
                                k=j;
                            }
                            else
                            {
                                k=readLength2-j-1;
                            }

                            if(fromCharCode==toCharCode)
                            {
                                tempProb*=(1-errorPosDist[k]-inPosDist[k]-delPosDist[k]);
                            }
                            else
                            {
                                if(fromCharCode >=0 && fromCharCode<=4)tempProb*=errorPosDist[k]*errorTypeProbs[fromCharCode][toCharCode];
                            }
                        }
                        if(tempProb>maxProb)
                        {
                            maxMatch=tempMatch;
                            maxProb=tempProb;
                            mleInsertSize=insertSize+i-gapStart;
                            maxPos=i-startPos;
                        }
                   }
                }
                else
                {
                    pos1 += gapoffset;
                    insertSize=pos1-(gapStart+gapLength)+readLength2;

                    for(long int i=gapStart-readLength2+1;i<gapStart+gapLength;i++)
                    {
                        tempInsertSize=insertSize+gapStart+gapLength-i;
                        if(tempInsertSize<insertThresholdMin || tempInsertSize > insertThresholdMax)
                        {
                            continue;
                        }

                        tempProb=1;//insertLengthDistSmoothed[tempInsertSize];
                        tempMatch=0;

                        for(int j=0;j<readLength2;j++)
                        {

                            toCharCode=charCodes[reads_gap[read_itr][j]];
                            index=i+j-startPos;
                            fromCharCode=charCodes[gapString[index]];

                            if(isReverse[read_itr]==0)
                            {
                                k=j;
                            }
                            else
                            {
                                k=readLength2-j-1;
                            }

                            if(fromCharCode==toCharCode)
                            {
                                tempProb*=(1-errorPosDist[k]-inPosDist[k]-delPosDist[k]);

                            }
                            else
                            {
                                if(fromCharCode >=0 && fromCharCode<=4)tempProb*=errorPosDist[k]*errorTypeProbs[fromCharCode][toCharCode];

                            }
                        }
                        if(tempProb>maxProb)
                        {
                            maxMatch=tempMatch;
                            maxProb=tempProb;
                            mleInsertSize=insertSize+gapStart+gapLength-i;
                            maxPos=i-startPos;
                        }
                    }
                }

                if(-log10(maxProb) < gapProbCutOff && saved_reads[read_itr] ==1)
                {
                    //if(pos1<gapStart)cout<<"If ACC\t"<<mleInsertSize<<"\t";
                    //else cout<<"Else ACC\t"<<mleInsertSize<<"\t";
                    //cout<<"Final Read = "<<read_itr<<", pos = "<<(maxPos-left_maxDistance)<<endl;
                    //cout<<"maxprob = "<<maxProb<<", pos = "<<(maxPos-left_maxDistance)<<", insertsize = "<<mleInsertSize<<"log = "<<-log10(maxProb)<<endl;
                    //cout<<"Read = "<<read_itr<<", prob. = "<<-log10(maxProb)<<endl;
                    //cout<<read_itr<<" "<<reads_gap[read_itr]<<endl;
                    //cout<<read_itr<<endl;
                    char ifelse;
                    if(pos1<gapStart)ifelse='I';
                    else ifelse='E';
                    draw_read(draw,maxPos-left_maxDistance,reads_gap[read_itr],read_itr,-log10(maxProb),mleInsertSize,gapno,ifelse);

                    //====================================================
                    unmapped_read_pos_arr[read_itr][0] = maxPos-left_maxDistance;
                    unmapped_read_pos_arr[read_itr][1] = readLength2;
                    unmapped_read_pos_arr[read_itr][2] = read_itr;

                    if(maxPos-left_maxDistance ==0)left_start_zero=1;
                    if(maxPos-left_maxDistance + readLength2==this->gapLength)right_fin_glen=1;

                    if(maxPos-left_maxDistance < 0 && maxPos-left_maxDistance+readLength2 > 0)
                    {
                        left_right_check[0]=1;
                        if(-(maxPos-left_maxDistance) > unmapped_max_left)unmapped_max_left = -(maxPos-left_maxDistance);
                    }

                    int val = maxPos-left_maxDistance + readLength2 - this-> gapLength ;
                    if(maxPos-left_maxDistance < this-> gapLength && val > 0)
                    {
                        left_right_check[1]=1;
                        if(val > unmapped_max_right)unmapped_max_right = val;
                    }
                    
                    
                    {
                         if(maxPos-left_maxDistance < 0 && maxPos-left_maxDistance+readLength2 > 0)
                         {
                            newleft=1;
                            if(-(maxPos-left_maxDistance) > newmaxleft)newmaxleft = -(maxPos-left_maxDistance);
                         }
                         if(maxPos-left_maxDistance > 0 && maxPos-left_maxDistance < this-> gapLength && val > 0)
                         {
                            newright=1;
                            if(val > newmaxright)newmaxright = val;
                         }
                    }

                    for(int j=0;j<readLength2;j++)
                    {
                        if(maxPos+j>=left_maxDistance && maxPos+j<left_maxDistance+gapLength)
                        {
                            countsGap[maxPos+j][charCodes[reads_gap[read_itr][j]]]+=1;
                        }
                    }
                }
                else
                {
                   //cout<<"Final Read = "<<read_itr<<", Discarded unmapped "<<-log10(maxProb)<<endl;
                    //cout<<read_itr<<endl;
                    //cout<<read_itr<<"\t"<<-log10(maxProb)<<"\t"<<reads_gap[read_itr]<<endl;
                    discardedCount++;
                }
            }//for end
        }
        //======================finalize of unmapped done====================

        if(partial_flag)//zxcf fin part
        {
            char* token;
            int ref_pos,newInsertSize=-1,read_start,read_end;
            read_itr = -1;
            draw_read(draw,this->gapLength,"S",0,0,0,gapno,' ');
            //================================================
            while(fgets(line1, MAX_FILE_READ, scaffoldMap2)!=NULL)
            {
                newInsertSize=-1;
                token = strtok(line1,"\t");
                readString1 = token;
                readLength1 = strlen(readString1);

                int clipped_index = atoi(strtok(NULL, "\t"));

                token = strtok(NULL, "\t");
                flag1 = atoi(token);

                if(flag1 == 1 || flag1 == 2)strandNo1=0;
                else strandNo1 = 1;

                pos1 = atoi(strtok(NULL, "\t"));
                //cout<<readString1<<" "<<flag1<<" "<<pos1<<" "<<strandNo1<<" "<<readLength1<<endl;
                //cout<<strandNo1<<endl;
                strtok(NULL, "\t");
                ref_pos = atoi(strtok(NULL, "\t"));

               strandNo1=0;


               ///This function also helps to place the read in correct position
               //if there are repeats in read and bowtie2 fails to give correct alignment

                maxProb=0;
                maxMatch=0;
                maxPos=0;

                totalCount++;

                if(flag1 == 1 || flag1 == 4){read_start = clip_thresh;read_end = 0;}
                else {read_start = 0;read_end = clip_thresh;}

                read_itr++;
                if(read_itr > (partial_limit-1))break;
                //cout<<"Read = "<<read_itr<<"\tPos1 = "<<pos1<<"\tGapstart = "<<gapStart<<endl;

                if(pos1<gapStart)
                {
                    insertSize = gapStart-ref_pos+readLength1;
                    for(long int i=gapStart-readLength1+1;i<gapStart;i++)
                    {
                        tempInsertSize=insertSize+i-gapStart;

                        if((tempInsertSize<insertThresholdMin || tempInsertSize > insertThresholdMax) && ref_pos !=-1)
                        {
                            continue;
                        }
                     
                        tempProb = 1;//(insertLengthDistSmoothed[tempInsertSize]);
                        
                        tempMatch=0;

                        for(int j=read_start;j<readLength1-read_end;j++)
                        {
                            toCharCode=charCodes[readString1[j]];
                            index=i+j-startPos;
                            fromCharCode=charCodes[gapString[index]];

                            if(strandNo1==0)
                            {
                                k=j;
                            }
                            else
                            {
                                k=readLength1-j-1;
                            }

                            if(fromCharCode==toCharCode)
                            {
                                tempProb*=(1-errorPosDist[k]-inPosDist[k]-delPosDist[k]);
                            }
                            else
                            {
                                if(fromCharCode >=0 && fromCharCode<=4)tempProb*=errorPosDist[k]*errorTypeProbs[fromCharCode][toCharCode];
                            }
                        }

                        if(tempProb>maxProb)
                        {
                            maxMatch=tempMatch;
                            maxProb=tempProb;
                            mleInsertSize=insertSize+i-gapStart;
                            maxPos=i-startPos;
                        }
                   }
                }
                else
                {
                    ref_pos += gapoffset;
                    insertSize=ref_pos-(gapStart+gapLength)+readLength1;

                    for(long int i=gapStart+gapLength-readLength1+1;i<gapStart+gapLength;i++)
                    {
                        tempInsertSize=insertSize+gapStart+gapLength-i;

                        if((tempInsertSize<insertThresholdMin || tempInsertSize > insertThresholdMax) && ref_pos !=-1)
                        {
                            continue;
                        }
                     
                        tempProb = 1;//(insertLengthDistSmoothed[tempInsertSize]);
                        
                        tempMatch=0;

                        for(int j=read_start;j<readLength1-read_end;j++)
                        {
                            toCharCode=charCodes[readString1[j]];
                            index=i+j-startPos;
                            fromCharCode=charCodes[gapString[index]];

                            if(strandNo1==0)
                            {
                                k=j;
                            }
                            else
                            {
                                k=readLength1-j-1;
                            }

                            if(fromCharCode==toCharCode)
                            {
                                tempProb*=(1-errorPosDist[k]-inPosDist[k]-delPosDist[k]);
                            }
                            else
                            {
                                if(fromCharCode >=0 && fromCharCode<=4)tempProb*=errorPosDist[k]*errorTypeProbs[fromCharCode][toCharCode];
                            }
                        }

                        if(tempProb>maxProb)
                        {
                            maxMatch=tempMatch;
                            maxProb=tempProb;
                            mleInsertSize=insertSize+gapStart+gapLength-i;
                            maxPos=i-startPos;
                        }
                    }
                }
                if(-log10(maxProb) < gapProbCutOff || partial_saved_read_final[0] == read_itr || partial_saved_read_final[1] == read_itr)
                {

                    //if(pos1<gapStart)cout<<"If ACC\t";
                    //else cout<<"Else ACC\t";
                    //cout<<"Read = "<<read_itr<<", maxProb_final = "<<-log10(maxProb)<<endl;
                    //cout<<"Read = "<<read_itr<<"\tmaxprob = "<<maxProb<<", pos = "<<(maxPos-left_maxDistance)<<", insertsize = "<<mleInsertSize<<endl;
                    //cout<<"Read = "<<read_itr<<", Else Accepted unmapped"<<endl;
                    //cout<<readString1<<endl;
                    if(maxPos-left_maxDistance <0)leftcount++;
                    else rightcount++;
                    if(maxPos-left_maxDistance <0 && maxPos-left_maxDistance + readLength1 >= this->gapLength)
                    {
                        if(-(maxPos-left_maxDistance) >= 3 && (maxPos-left_maxDistance + readLength1 - this->gapLength) >=3)
                            leftcount = rightcount = 10;
                    }
                    partial_read_flag[read_itr][0] = 1;
                    partial_read_flag[read_itr][1] = maxPos-left_maxDistance;
                    partial_read_flag[read_itr][2] = readLength1;

                    int placed_pos = maxPos-left_maxDistance + gapStart;
                    if(maxPos-left_maxDistance <0)
                    {
                        if(ref_pos != -1)
                        {
                            newInsertSize = placed_pos - ref_pos + readLength1;
                        }
                    }
                    else
                    {
                        if(ref_pos != -1)
                        {
                            newInsertSize = ref_pos + readLength1 - placed_pos ;
                        }
                    }

                    draw_read(draw,maxPos-left_maxDistance,readString1,read_itr,flag1,newInsertSize,gapno,'P');

                    int flag=0;
                    for(int j=0;j<readLength1;j++)
                    {
                        if(maxPos+j>=left_maxDistance && maxPos+j<left_maxDistance+gapLength)
                        {
                            countsGap[maxPos+j][charCodes[readString1[j]]]+=1;
                            qual_gap[maxPos+j-left_maxDistance][charCodes[readString1[j]]] += partial_read_quality[read_itr][j];
                        }
                    }
               }
               else
               {

                    discardedCount++;
               }
            }
        }
        //=========================================================================
        total_read += totalCount;
        discarded_read += discardedCount;

        //======================New Modifications==================================
        int used_read = totalCount - discardedCount,u_flag=1,recompute_flag=0;

        int Nflag[2],lflag[2];
        for(int i=0;i<2;i++)
        {
            Nflag[i] = lflag[i] = -1;
        }
        //cout<<used_read<<endl;
        //cout<<leftcount<<"\t"<<rightcount<<"\t"<<this->gapLength<<"\t"<<partial_flag<<endl;

        int unmapped_sidethresh;

        if(unmapped == 1)
        {
            //cout<<unmapped_max_left<<"\t"<<unmapped_max_right<<endl;
            //cout<<left_right_check[0]<<"\t"<<left_right_check[1]<<endl;
            //cout<<left_start_zero<<"\t"<<right_fin_glen<<endl;

            int new_flag=0;
            unmapped_sidethresh=4;
            

            if((unmapped_max_left <2*unmapped_sidethresh && unmapped_max_left >0 )||
                (unmapped_max_right < 2*unmapped_sidethresh && unmapped_max_right >0))
            {
                if(region_perct_max_gap_estm > .75)
                {
                    
                        used_read=0;
                }
            }

            if((unmapped_max_left  < unmapped_sidethresh && unmapped_max_left > 0) ||
                (unmapped_max_right < unmapped_sidethresh && unmapped_max_right > 0))
            {
                int region[5000];
                computeSequence(1,0);
                int rc = findRegion(region);

                if(rc>=1)
                {
                    if(unmapped_max_left  < unmapped_sidethresh && unmapped_max_left > 0)
                    {
                        lflag[0]=1;
                        //cout<<"Left side cleared\n";
                    }

                    if(unmapped_max_right < unmapped_sidethresh && unmapped_max_right > 0)
                    {
                        lflag[1]=1;
                        //cout<<"Right right side cleared\n";
                    }
                }
                else if(rc==0)
                {
                    used_read = 0;
                    unmapped_max_left = unmapped_max_right = -1;
                    //cout<<"Gap gone2\n";
                }
            }

            if(left_right_check[0]== 0 && left_right_check[1] == 0 && used_read !=0)
            {
                 used_read = 0;
                 unmapped_max_left = unmapped_max_right = -1;
                 //cout<<"Gap gone3\n";
            }
            if((left_right_check[0] == 0 && left_start_zero != 0) || (left_right_check[1] == 0 && right_fin_glen != 0))
            {
                int region[5000];
                computeSequence(1,0);
                int rc = findRegion(region);

                if(rc>=1)
                {
                    if(left_right_check[0] == 0)Nflag[0]=1;
                    if(left_right_check[1] == 0)Nflag[1] = 1;
                }
           }

            int jump_condition = 1;

            if(used_read == 0 || new_flag==1 || (!(left_right_check[0] ==1 && left_right_check[1] == 1) && jump_condition))
            {
                this->gapLength=originalGap;
                initialize_start_end();

                int offset = this->gapLength > gapLength?0:(gapLength-this->gapLength);
                clear_countsGap(this->gapLength+offset);

                if(!left_right_check[0] && left_right_check[1] && unmapped_max_right>= unmapped_sidethresh)
                {
                    //cout<<"Recomputing for right\n";
                    recompute2(unmapped_read_pos_arr_org);
                    recompute_flag=1;
                }
                else if(left_right_check[0] && !left_right_check[1] && unmapped_max_left >= unmapped_sidethresh)
                {
                    //cout<<"Recomputing for left\n";
                    recompute2(unmapped_read_pos_arr_org);
                    recompute_flag=1;
                }
            }
            
            
        }
//zxcf
        if(partial_flag==1)
        {
            int ret_val[2]={0,0};
            int fcount=0,f_flag=0,u_flag=1;

            detect_overlap_gapestimate(partial_read_flag,this->gapLength,ret_val,8);
            //cout<<"Final overlap count = "<<ret_val[0]<<"\tCorrect/False = "<<ret_val[1]<<endl;
            //cout<<"Max overlap found between read "<<partial_saved_read_final[0]<<" and read "<<partial_saved_read_final[1]<<endl;

            int gap_case=0;
            
            if((originalGap - gapLength) > 0 && ret_val[0] > 0)
            {
                gap_case=1;
                //fprintf(cases,"Gap = %d\tCase = 1\tretv[0] = %d\tretv[1] = %d\n",gapno,ret_val[0],ret_val[1]);
            }
            else if((originalGap - gapLength) > 0 && ret_val[0] == 0)
            {
                gap_case=2;
                //fprintf(cases,"Gap = %d\tCase = 2\tretv[0] = %d\tretv[1] = %d\n",gapno,ret_val[0],ret_val[1]);
            }
            else if((originalGap - gapLength) < 0 && ret_val[0] > 0)
            {
                gap_case=3;
                //fprintf(cases,"Gap = %d\tCase = 3\tretv[0] = %d\tretv[1] = %d\n",gapno,ret_val[0],ret_val[1]);
            }
            else if((originalGap - gapLength) < 0 && ret_val[0] == 0)
            {
                gap_case=4;
                //fprintf(cases,"Gap = %d\tCase = 4\tretv[0] = %d\tretv[1] = %d\n",gapno,ret_val[0],ret_val[1]);
            }
            else
            {
                gap_case=5;
                //fprintf(cases,"Gap = %d\tCase = 5\tretv[0] = %d\tretv[1] = %d\n",gapno,ret_val[0],ret_val[1]);
            }

            //cout<<"Gap case = "<<gap_case<<"\tIn script itr = "<<script_itr<<"\tUsed reads = "<<used_read<<endl;

            if(used_read < partial_threshold || gap_case == 2 || gap_case == 4)
            {
                this->gapLength=originalGap;
                initialize_start_end();
                int offset = this->gapLength > gapLength?0:(gapLength-this->gapLength);
                clear_countsGap(this->gapLength+offset);

                if(used_read < partial_threshold || gap_case == 4)
                {
                    u_flag=0;
                    //cout<<"One/zero Read in gap\n";
                    //cout<<"Clearing gap case -4\n";
                }
                else //gap_case==2
                {
                    recompute1(partial_read_pos_arr_org);
                    for(int i=0;i<this->partial_read_count;i++)
                    {
                        partial_read_flag[i][0] = partial_read_pos_arr_org[i][0];
                        partial_read_flag[i][1] = partial_read_pos_arr_org[i][1];
                        partial_read_flag[i][2] = partial_read_pos_arr_org[i][2];
                    }

                    detect_overlap_gapestimate(partial_read_flag,this->gapLength,ret_val,8);

                    if(ret_val[1] ==-1)
                    {
                        clear_countsGap(this->gapLength+offset);
                        u_flag=0;
                        //cout<<"Clearing gap case -2, due to error in original too\n";
                    }
                }
            }

            if(u_flag ==1 && ret_val[0] == 0 && ret_val[1]==0)
            {
                for(int j=0;j<this->gapLength;j++)
                {
                    if(checkAllZero(countsGap[j+left_maxDistance],0) == 1)//there are some non zero indices
                    {
                        int update_flag = check_update(countsGap[j+left_maxDistance],gapno,j);

                        if(update_flag != -1)
                        {
                            if(update_flag != 50)
                            {
                                countsGap[j+left_maxDistance][update_flag] += 10;
                            }
                        }
                        else
                        {
                            for(int k=0;k<4;k++)
                            {
                                countsGap[j+left_maxDistance][k] = 0;
                            }
                        }
                    }
                }
            }
        }
        //====================================================================
        //fprintf(minmax2,"%d\n",this->gapLength);
        finalUsedReads[2*this->gap_no] = totalCount;
        finalUsedReads[2*this->gap_no+1] = used_read;
        //printArr(left_maxDistance-1,left_maxDistance+this-> gapLength);

        computeSequence(1,0);

        if((unmapped && (left_right_check[0] || left_right_check[1] || used_read !=0)))//nc
        {
            //cout<<recompute_flag<<endl;
            //cout<<Nflag<<endl;
            if(Nflag[0] ==1)concensus[0]='N';
            if(Nflag[1] ==1)concensus[this->gapLength-1]='N';

            if(lflag[0] ==1)concensus[0]='N';
            if(lflag[1] ==1)concensus[this->gapLength-1]='N';

            int clear_val;
            if(recompute_flag ==0)clear_val=recheck_sequence(unmapped_read_pos_arr);
            else clear_val=recheck_sequence(unmapped_read_pos_arr_org);
            
            if(clear_val ==1)
            {
                this->gapLength=originalGap;
                initialize_start_end();
                
                int offset = this->gapLength > gapLength?0:(gapLength-this->gapLength);
                clear_countsGap(this->gapLength+offset);
                computeSequence(1,0);
            }
        }

        delete [] gapString;

        for(int i=0;i<this->num_reads_in_gap;i++)
        {
            delete [] unmapped_read_pos_arr[i];
        }
        delete[] unmapped_read_pos_arr;

        for(int i=0;i<this->partial_read_count;i++)
        {
            delete [] partial_read_flag[i];
        }
        delete[] partial_read_flag;

        if(partial_flag)fclose(scaffoldMap2);
    }

    void parseUnmapped(int g)
    {
        FILE * scaffoldMap1=NULL;
        scaffoldMap1=fopen(gapFileName, "r");

        char line1[MAX_REC_LEN];
        char line2[MAX_REC_LEN];

        int MAX_FILE_READ=MAX_REC_LEN/sizeof(line1[0]);

        char * temp,*readString1,*readString2,*quality_str;

        int readLength1, readLength2, strandNo1,strandNo2;

        char tempReadString[MAX_READLENGTH];

        int readIndex,flag;

        int r_c=0,rev_flag;
        this->max_read_len=0;
        int limit=50000;//change allocate if you change here
        int total_read=0;
        
        while(fgets(line1, MAX_FILE_READ, scaffoldMap1)!=NULL)
        {
            if(line1[0]=='@')
                continue;
            
            if(strlen(line1) < 60)
            {
                //fgets(line2, MAX_FILE_READ, scaffoldMap1);
                cout<<"unknown problem in sam, skipping read no - "<<r_c<<" for gap - "<<g<<endl;
                
                continue;
            }
            
            if(fgets(line2, MAX_FILE_READ, scaffoldMap1)==NULL)
                break;

            //1st of the pair
            strtok(line1,"\t");
            temp=strtok(NULL,"\t");
            flag=atoi(temp);

            strandNo1=(flag&16)>>4;
            strtok(NULL,"\t");
            temp=strtok(NULL,"\t");
            pos_reads[r_c]=atoi(temp);

            //second of the pair
            strtok(line2,"\t");
            temp=strtok(NULL,"\t");
            flag=atoi(temp);

            for(int j=0;j<4;j++)strtok(NULL,"\t");
            readString2=strtok(NULL,"\t");
            readLength2=strlen(readString2);

            if(readLength2 > this->max_read_len)this->max_read_len = readLength2;

            quality_str = strtok(NULL,"\t");
            quality_str[this->unmapped_read_len]='\0';
            //qualityFilter(quality_str,r_c,0);
            
            

            if(strandNo1==0)
            {
                reverse(tempReadString, readString2);
                strcpy(readString2,tempReadString);
                rev_flag=1;
            }
            else
            {
                rev_flag=0;
            }
            
            
            if(readLength2 > limit)
            {
                string s1(readString2);
                string s2 = s1.substr(0,limit);
                string s3 = s1.substr(limit,readLength2-limit);
                
                strcpy(reads_gap[r_c],s2.c_str());
                strcpy(reads_gap[r_c+1],s3.c_str());

                pos_reads[r_c+1]=pos_reads[r_c];
                
                isReverse[r_c]=rev_flag;
                isReverse[r_c+1]=rev_flag;
                r_c += 2;
            }
            else
            {
                strcpy(reads_gap[r_c],readString2);
                isReverse[r_c]=rev_flag;
                //cout<<reads_gap[r_c]<<endl;
                r_c++;
            }
            
            total_read++;
            if(total_read == unmapped_limit)break;
        }
        fclose(scaffoldMap1);

    }

    void parsePartial()
    {
        FILE* partial_gap_file =fopen(partialgapFileName,"r");

        char line1[MAX_REC_LEN];
        char * token,*readstring;
        int MAX_FILE_READ=MAX_REC_LEN/sizeof(line1[0]);

        int countp=0;
        if(partial_flag)countp = this->partial_read_count;
        else
        {
            while(fgets(line1, MAX_FILE_READ, partial_gap_file)!=NULL)
            {
                countp++;
                if(countp>partial_limit)break;
            }

            fclose(partial_gap_file);
            partial_gap_file =fopen(partialgapFileName,"r");
        }

        this->p_count = countp;


        if(this->p_count > 0)
        {
            if(partial_flag)partial_read_quality = new double*[this->partial_read_count];

            partialreads_gap = new char*[this->p_count];

            for(int i=0;i<this->p_count;i++)
            {
                partialpos_reads = new int[this->p_count];
                preadflag = new int[this->p_count];
                clipped_index = new int[this->p_count];
            }

            countp=0;
            int len=0;

            while(fgets(line1, MAX_FILE_READ, partial_gap_file)!=NULL)
            {
                token = strtok(line1,"\t");
                readstring = token;
                len = strlen(readstring);

                partialreads_gap[countp] = new char[partial_read_len];

                strcpy(partialreads_gap[countp],readstring);

                token = strtok(NULL, "\t");
                clipped_index[countp] = atoi(token);

                token = strtok(NULL, "\t");
                preadflag[countp] = atoi(token);

                token=strtok(NULL,"\t");
                partialpos_reads[countp] = atoi(token);

                if(partial_flag)
                {
                    for(int i=0;i<3;i++)token = strtok(NULL,"\t");

                    token[partial_read_len]='\0';
                    partial_read_quality[countp] = new double[partial_read_len];
                    qualityFilter(token,countp,1);
                }
                countp++;
                if(countp>partial_limit)break;
            }
        }
        fclose(partial_gap_file);
    }

	void printArr(int a,int b)
	{

	    printf("Printing probsgap array\n\n");
        for(int i=a;i<=b;i++)
        {
            cout<<"i = "<<i-left_maxDistance<<"==>";
            for(int j=0;j<4;j++)
                cout<<probsGap[i][j]<<" ";
            cout<<endl;
        }

/*
	    printf("Printing probsgap array\n\n");
	    for(int i=a;i<=b;i++)
	    {
	        cout<<"i = "<<i<<"==>";
	        for(int j=0;j<4;j++)
	            cout<<probsGap[i][j]<<" ";
	        cout<<endl;
	    }


	    printf("Printing errorprobsgap array\n\n");
	    for(int i=a;i<b;i++)
        {
            cout<<"i = "<<i<<"==>";
            for(int j=0;j<4;j++)
                cout<<errorProbsGap[i][j]<<"    ";
            cout<<endl;
        }

        printf("Printing errorposdist array\n\n");
        for(int i=0;i<maxReadLength;i++)
        {
            cout<< errorPosDist[i]<<endl;
        }
        cout<<endl;

*/
	}

	int check_change(int arr[],int n)
	{
	    //for(int i=0;i<n;i++)cout<<arr[i]<<endl;
	    if(n==1)return 0;
	    for(int i=1;i<n;i++)
	    {
	        if(arr[0] != arr[i])return 1;
	    }
	    return 0;
	}

	int findmaxgap(double arr[][2],int n)
	{
	    double min = -1000000;
	    int min_index = 100000;
	    for(int i=0;i<n;i++)
	    {
	        if(arr[i][1] > min || (arr[i][1] == min && arr[i][0] < min_index))
	        {
	            min_index = arr[i][0];
	            min = arr[i][1];
	            //cout<<min_index<<endl;
	        }
	    }
	    return min_index;
	}

	double run(int gaplen,FILE* fp,int g,int finalize_flag,int c)///c =called from
	{
	    double likelihood=0;
	    left_maxDistance = tempmaxDistance;
        right_maxDistance = tempmaxDistance;
        initialize(gaplen,1,1);

        if(side_limit <10){return 10000;}

        int p=0;
        
        step2match=0;

        int preset_unfilled_len = 2*this->unmapped_read_len;
        comp_count=0;
        
        for(p=0;p<num_itr;p++)
        {
            valid_count=0;
            invalid_count= 0;

            likelihood = placeReads(p,g,finalize_flag,fp,gaplen-originalGap,0);
            computeProbsGap(0);
            computeErrorProbsGap();

            if(unmapped)
            {
                if(comp_count >=5)
                {
                    //cout<<"EM done, Final itr = "<<p<<endl;
                    break;
                }
                
            }
            else
            {
                if(p == 2)break;
            }
            //computeSequence(0,0);
            //cout<<concensus<<endl;

        }
        //cout<<"ran with "<<gaplen<<", valid count = "<<valid_count<<endl;
        //for(int k=0;k<this->num_reads_in_gap;k++)saved_reads[k] = mark_accepted_reads[k];
        region_perct_max_gap_estm = region_perct;
        
        double retvall;
        
        if(c==0)retvall=valid_count;
        if(c==1)retvall= likelihood;
        
        return retvall;
	}

	void findRange(int len,double org_read,FILE* fp,int g,int finalize_flag,double lf,int sample_count,int *val_ret)
	{
	    int verbose =0;
	    if(verbose ==1)cout<<"Original Gap = "<<len<<"\tLikelihood = "<<org_read<<endl;
	    int range = 40;
	    double maxlikelihood=-1000000,likelihood=0;

	    int rc_gapl[100],rc_gapr[100];
	    double rc_left[100],rc_right[100];
	    int min_gap = lf*originalGap;
	    int skip = (originalGap - min_gap)/sample_count,skip_count=0;
	    int new_gap = originalGap;
	    int c1=0,c2=0,c=0,cl=0,cr=0,r;
        if(verbose ==1)cout<<"Sampling gapestimates\n";

	    while(skip_count<sample_count)
	    {
	        new_gap -= skip;
	        likelihood = run(new_gap,fp,g,finalize_flag,0);
	        if(likelihood == 10000){val_ret[0] = -3;return;}
	        rc_left[c1++] = likelihood;
	        rc_gapl[c1-1] = new_gap;
	        if(verbose ==1)cout<<"Gap = "<<new_gap<<"\tReads = "<<valid_count<<"\tLikelihood = "<<rc_left[c1-1]<<endl;
	        skip_count++;
	    }

        new_gap = originalGap;
        skip_count=0;
       if(verbose ==1) cout<<"=================="<<endl;
        double temp_valid=0;
        int right_same=0,rsameflag=0,fixedgap=3*len;

	    while(skip_count<sample_count)
        {
            new_gap += skip;
            likelihood = run(new_gap,fp,g,finalize_flag,0);
            if(likelihood == 10000){val_ret[0] = -3;return;}
            rc_right[c2++] = likelihood;
            if(verbose ==1)cout<<"Gap = "<<new_gap<<"\tReads = "<<valid_count<<"\tLikelihood = "<<rc_right[c2-1]<<endl;
            rc_gapr[c2-1] = new_gap;
            skip_count++;
            if(likelihood==temp_valid)
            {
                if(right_same==0)fixedgap = new_gap;
                right_same++;
            }
            else
            {
                temp_valid = likelihood;
                right_same=0;
            }

            if(right_same > .4*sample_count)
            {
                for(int d=0;d<sample_count-skip_count;d++)
                {
                    rc_right[c2++]=likelihood;
                    rc_gapr[c2-1] = fixedgap;
                }
                rsameflag=1;
                break;
            }
        }
        c = min(c1,c2);

        //for(int i=0;i<c;i++)
          //  cout<<rc_left[i]<<" "<<rc_gapl[i]<<" "<<rc_right[i]<<" "<<rc_gapr[i]<<endl;

        int max_index=-1,temp_max=0;
        double curr_max_rc=org_read;

        for(int i=0;i<c;i++)
        {
            if(rc_left[i] > org_read)cl++;
            if(rc_left[i] >= org_read)
            {
                if(rc_left[i] > curr_max_rc){max_index = i;curr_max_rc = rc_left[max_index];}
                //cout<<cl<<" "<<rc_left[i]<<" "<<curr_max_rc<<" "<<max_index<<" "<<rc_left[max_index]<<endl;
            }
        }

        if(verbose ==1)cout<<"Left side best gap estimate = "<<rc_gapl[max_index]<<"\tUsed reads = "<<rc_left[max_index]<<endl;
        if(max_index != -1)temp_max = rc_gapl[max_index];

        /*if(rsameflag==1)
        {
            //cout<<"Returning as right side became constant\n";
            val_ret[0] = -5;
            val_ret[1] = (temp_max-50)>=0?(temp_max-50):0;
            val_ret[2] = fixedgap;
            return ;
        }*/

        max_index=-1,curr_max_rc=org_read;

        for(int i=0;i<c;i++)
        {
            if(rc_right[i] > org_read)cr++;
            if(rc_right[i] >= org_read)
            {
                if(rc_right[i] > curr_max_rc){max_index = i;curr_max_rc = rc_right[max_index];}
                //cout<<cr<<" "<<rc_right[i]<<" "<<curr_max_rc<<" "<<max_index<<" "<<rc_gapr[max_index]<<endl;
            }
        }
        if(verbose ==1)cout<<"Right side best gap estimate = "<<rc_gapr[max_index]<<"\tUsed reads = "<<rc_right[max_index]<<endl;
        if(verbose ==1)cout<<"Left greater = "<<cl<<endl;
        if(verbose ==1)cout<<"Right greater = "<<cr<<endl;

        int diff = cl-cr < 0?cr-cl:cl-cr;

        if(diff <=4)
        {
            val_ret[0]=0;
            //cout<<"\nLarge gap Unsure,so returning without...\n"<<endl;
            return;
        }
        else if(cl > cr)r=-1;
        else if(cl < cr) r = 1;

        int chosen_gap=0,best_gap=-1,valid_max=-1;
        if(r ==-1)chosen_gap = temp_max;
        else if(r==1)chosen_gap = rc_gapr[max_index];
        else chosen_gap = len;

        //cout<<r<<" "<<"\tChosen gap center = "<<chosen_gap<<endl;
         maxlikelihood=-1000000,likelihood=0;
         int same_countg=0;
         /*
        for(int i=chosen_gap -range;i<chosen_gap+range;i+=1)
        {
            likelihood = run(i,fp,g,finalize_flag,1,1);
            
            if(likelihood > maxlikelihood)
            {
                maxlikelihood = likelihood;
                best_gap = i;
                same_countg=0;
            }
            else
            {
                same_countg++;
            }
            
            //if(same_countg>(2*range/3))break;
            //cout<<"GapLength = "<<i<<"\tUsed read = "<<valid_count<<"\tLikelihood = "<<likelihood<<"\tCurrent best = "<<best_gap<<endl;
        }
        */
        //cout<<"Returning r = "<<r<<", cl = "<<cl<<", cr = "<<cr<<", best = "<<best_gap<<endl;
	    val_ret[0]=1;
	    val_ret[1] = chosen_gap;//best_gap;

	    return;
	}

	int checkGapReads(int org,FILE* fp,int g)
	{
        int step;
        int readcountThresh=3;

        if(org<30)
        {
            if(org<15)step = 10;
            else step = 20;
            for(int i=0;i<80;i+=step)
            {
                double val = run(i,fp,g,1,1);
                if(val == 10000)return -2;
                //cout<<"Gapestimate = "<<i<<"\tUsed read = "<<valid_count<<endl;
                if(valid_count> readcountThresh)return -1;
            }
        }
        else
        {
            for(int k=0;k<4;k++)
            {
                int gap;
                if(k==0)gap = org/2;
                else gap = org*k;

                double val = run(gap,fp,g,1,1);
                if(val == 10000)return -2;
                //cout<<"Gapestimate = "<<gap<<"\tUsed read = "<<valid_count<<endl;
                if(valid_count >= readcountThresh)return -1;
            }
        }
        return 1;
	}

	
    
    void setParameters()
    {
        covergae_threshold1 = 0;//for final read placement
        match_count_discont=4;//means >= match_count_discont char must overlap betwn 2 reads
        covergae_threshold2 = 1;//for likelihood modification-3
        
        partial_threshold=2;
        clip_thresh= 2;
    }


    void analyzeGap(int * g,int fillflag)
    {
        findGapLeftRight();

        this->rep_flag = findRepeat();//it will be 1 if both sides are found in the read

        if(this->rep_flag == 1 && partial_flag)
        {
            //cout<<"Gap = "<<g<<" not filled due to repeat\n";
            g[0]=-1;g[1]=-1;
        }
        else if(this->one_side_repeat_flag ==1 && partial_flag && originalGap > 3*mid_limitp)
        {
            //cout<<"In case of frag, > 600 Gaplen gap = "<<g<<" not filled due to repeat problem\n";
            g[0]=-1;g[1]=-1;
        }
        else if(fillflag == -1){g[0]=-1;g[1]=-1;}
        else
        {
            if(this->partial_read_count==0 && this->num_reads_in_gap==0)
            {//cout<<"No partial read in gap\n";
                g[0]=-1;g[1]=-1;
            }
            else
            {
                int q=-1;
                if(script_itr ==1){g[0]=0;g[1]=765466;}///zxcfg
                //if(script_itr ==1){g[0]=q;g[1]=q;}
                else{g[0]=0;g[1]=965757;}
            }
        }
    }

	void fillGap(int g,FILE * fp,int fillflag)
	{
        float gp_frac1,gp_frac2;

        int finalize_flag=1,gnum1,gnum2,sampleflag=0;
        int gapnum[2];
        
        if(unmapped)
        {
            if(large_gap_flag ==0)finalize_flag=0;
        }

	    if(unmapped)parseUnmapped(g);
	    parsePartial();

	    setParameters();

        analyzeGap(gapnum,fillflag);

        gnum1=gapnum[0];
        gnum2=gapnum[1];

	    if(g>=gnum1 && g <=gnum2)
        {
            num_itr=200;
        }
        else num_itr=0;

        if(!(g>=gnum1 && g <=gnum2))
        {
            this->gp_frac1=1;
            this->gp_frac2=1;
        }

        //cout<<g<<"\t"<<gnum1<<"\t"<<gnum2<<"\t"<<this->gp_frac1<<"\t"<<this->gp_frac2<<endl;

        int gapMin= originalGap*this->gp_frac1;
        int gapMax=originalGap*this->gp_frac2;
        
		int gapEstimate=gapMin;

        int maxGapEstimate = gapMin;
        int secondMaxGapEstimate = gapMin;

        double maxLikelihood=-DBL_MAX;
        double secondMaxLikelihood=-DBL_MAX;

        double likelihood=0,prevlikelihood=0;

        bestString[0]='\0';
		secondBestString[0]='\0';
		current_str[0]='\0';
		current_str2[0]='\0';
        previous_str[0]='\0';
        char original_str[originalGap+1];

        int fill_or_not = 0;
        int same_count=0,partial_same=4,jump_same = 50,same_thresh,gaplimit=900;
        int stuckCount=0;

        if(unmapped)same_thresh=jump_same;
        else same_thresh=partial_same;

        int range=gapMax-gapMin+1;
        int used_read_arr[range];
        double diff1=0,diff2=0;
        int side_flag=0;
        int j=0,k=range-1;

        double likelihood_arr[range][2];
        int left_right_maxDistance[range][2];


        for(int i=0;i<range;i++)
        {
            for(int j=0;j<2;j++)
            {

                left_right_maxDistance[i][j]=0;
                likelihood_arr[i][j]=0;
            }
        }


        int less_read_flag=0;
        //if(g>=gnum1 && g <=gnum2)cout<<"Gap = "<<g<<"\tOriginal Gapestimate = "<<originalGap<<endl;

       // if(g>=gnum1 && g <=gnum2)cout<<"Gap "<<g<<"- Gap left = "<<gap_left<<"\tGap right = "<<gap_right<<endl;

        if(sampleflag == 0 && unmapped && originalGap <=mid_limitu && (g>=gnum1 && g <=gnum2))
            less_read_flag= checkGapReads(originalGap,fp,g);

        if(less_read_flag == 1)range=0;
        if(less_read_flag == -2){side_flag =1;range=0;}
        
        int prev_best=-1,curr_best=0,prev_u=-1,curr_u=0,sec_same=0,sec_same2=0;
        
        for(;j<range;j++)
        {
            umaxleftf=umaxrightf=ucoverf=0;
            fill_or_not = initialize(gapEstimate,1,j);
            //fill_or_not=0;
            if(side_limit <10){side_flag =1;break;}

            if(this->one_side_repeat_flag == 1)fill_or_not=0;
            if(fill_or_not !=0 && (g>=gnum1 && g <=gnum2))break;

            int i=0;
            this->discont_or_not=0;

            time_t now,end;
            time(&now);
            
            int preset_unfilled_len = 2*this->unmapped_read_len;
            comp_count=0;
            
            overlap_threshold=5;

            step2match=0;

            int breakflag=0;

            for(i=0;i<num_itr;i++)//zxcf
            {
                valid_count=0;
                invalid_count= 0;
                likelihood = placeReads(i,g,finalize_flag,fp,gapEstimate-originalGap,large_gap_flag);
                computeProbsGap(0);
                computeErrorProbsGap();
                
                if(unmapped)
                {
                    if(comp_count >=5)
                    {
                        //cout<<"EM done, Final itr = "<<i<<endl;
                        break;
                    }
                    if(large_gap_flag == 1 && region_perct*gapEstimate < preset_unfilled_len){breakflag=1;break;}   
                }
                else
                {
                    if(i == 2)break;
                }
            }
            
            //cout<<i<<endl;
            
            if(unmapped && !finalize_flag && (g>=gnum1 && g <=gnum2))
            {
                valid_count=0;
                likelihood = placeReads(i,g,1,fp,gapEstimate-originalGap,0);
            }

            time(&end);

            computeSequence(0,0);
            
            if(likelihood>maxLikelihood)
            {
                secondMaxLikelihood=maxLikelihood;
                secondMaxGapEstimate=maxGapEstimate;
                strcpy(secondBestString,bestString);
                maxLikelihood=likelihood;
                maxGapEstimate=gapEstimate;
                //computeSequence(0,0);

                strcpy(bestString,concensus);
                //cout<<"CS called from if- gapl = "<<gapEstimate<<" "<<bestString<<endl;
                for(int k=0;k<this->num_reads_in_gap;k++)saved_reads[k] = mark_accepted_reads[k];
                region_perct_max_gap_estm = region_perct;

                for(int k=0;k<2;k++)partial_saved_read_final[k] = partial_saved_read_temp[k];
                
                curr_best = j;
                prev_u = valid_count;
            }
            else if(likelihood>secondMaxLikelihood)
            {
                secondMaxLikelihood=likelihood;
                secondMaxGapEstimate=gapEstimate;
                //cout<<"CS called from ielse"<<endl;
                //computeSequence(0,0);
                strcpy(secondBestString,concensus);
                //cout<<"CS called from else if- gapl = "<<gapEstimate<<" "<<secondBestString<<endl;
            }
            if(g>=gnum1 && g <=gnum2)
            {    
                //cout<<"GapEstimate = "<<gapEstimate<<"\t"<<"Likelihood = "<<likelihood<<"\t"<<"Used Reads = "<<valid_count<<"\tTime taken = "<<difftime(end,now)<<"\tCurrent best gapest. = "<<maxGapEstimate<<endl;
            }
            likelihood_arr[j][0] = gapEstimate;
            likelihood_arr[j][1] = likelihood;

            left_right_maxDistance[j][0] = left_maxDistance;
            left_right_maxDistance[j][1] = right_maxDistance;


            if(gapEstimate == originalGap)strcpy(original_str,concensus);

            if(partial_flag || unmapped)
            {
                used_read_arr[j] = valid_count;

                diff1 = abs(prevlikelihood - likelihood);
                if(diff1 <=0.9)same_count++;
                else same_count=0;
                prevlikelihood = likelihood;
                
                if(same_count == same_thresh && this->gapLength >= originalGap)break;
                else if(same_count == same_thresh && this->gapLength < originalGap)
                {
                    run(originalGap,fp,g,1,0);
                    computeSequence(0,0);
                    strcpy(original_str,concensus);

                    break;
                }
                
                if(unmapped)
                {
                    curr_u = valid_count;
                    
                    //if(g>=gnum1 && g <=gnum2)cout<<sec_same<<"\t"<<curr_best<<"\t"<<prev_best<<"\t"<<curr_u<<"\t"<<prev_u<<endl;
                    
                    if(curr_best == prev_best && abs(curr_u - prev_u) <=2)sec_same++;
                    else
                    {
                        prev_best = curr_best;
                        sec_same =0;
                        //prev_u = curr_u;
                    }
                    
                    if(sec_same >= 2*same_thresh)
                    {
                        if(this->gapLength >= originalGap)break;
                        else if(this->gapLength < originalGap)
                        {
                            run(originalGap,fp,g,1,0);
                            computeSequence(0,0);
                            strcpy(original_str,concensus);
                            break;
                        }
                    }
                    
                    if(originalGap <=30)
                    {
                        if(!(umaxleftf ==1 || umaxrightf ==1 || ucoverf ==1))sec_same2++;
                        else sec_same2=0;
                        
                        //cout<<sec_same2<<endl;
                        
                        if(sec_same2 >= 1.5*same_thresh)
                        {
                            if(this->gapLength >= originalGap)break;
                            else if(this->gapLength < originalGap)
                            {
                                run(originalGap,fp,g,1,0);
                                computeSequence(0,0);
                                strcpy(original_str,concensus);
                                break;
                            }
                        }
                    }
                
                    if(this->discont_or_not==1 && valid_count < 5)stuckCount++;
                    else stuckCount=0;
                    if(stuckCount > 3* same_thresh)
                    {
                        if(this->gapLength >= originalGap)break;
                        else if(this->gapLength < originalGap)
                        {
                            run(originalGap,fp,g,1,0);
                            computeSequence(0,0);
                            strcpy(original_str,concensus);
                            break;
                        }
                    }
                }
            }

            gapEstimate++;
            
        }

        //cout<<"Gap "<<g<<"- Partial Gap left = "<<partial_left<<"\tPartial Gap right = "<<partial_right<<endl;

        if(fill_or_not !=0 && (g>=gnum1 && g <=gnum2))
        {
            //cout<<"Gap Not filled= "<<g<<", ----#MaxGapEstimate = 0"<<endl;
            //fprintf(minmax1,"%d\n",0);
            this->gapLength=0;
            gaptofill[g]=fill_or_not;
        }
        else
        {
            if(g>=gnum1 && g <=gnum2)
            {
                //cout<<"Gap = "<<g<<", ----#MaxGapEstimate = "<<maxGapEstimate<<"\t"<<"MaxLikelihood = "<<maxLikelihood<<"\tperct. = "<<region_perct_max_gap_estm<<endl;

            }

            //fprintf(minmax1,"%d\n",originalGap);

            if(g>=gnum1 && g <=gnum2)
            {
                if(unmapped)
                {
                    if(less_read_flag==1)
                    {
                        //cout<<"Gap = "<<g<<"\tOriginalGap = "<<originalGap<<"\tLess reads, filling with original length\n";
                        run(originalGap,fp,g,1,0);
                        computeSequence(0,0);
                        strcpy(original_str,concensus);
                        finalize(originalGap,g);
                    }
                    else
                    {
                        if(side_flag)
                        {
                            //cout<<"Side break- Unmapped\n";
                            run(originalGap,fp,g,1,0);
                            computeSequence(0,0);
                            strcpy(bestString,concensus);
                            finalize(originalGap,g);
                        }
                        else if(check_change(used_read_arr,j) == 1)
                        {
                            //cout<<"here2 - "<<j<<endl;
                            finalize(maxGapEstimate,g);
                        }
                        else
                        {
                            //cout<<"here3"<<endl;
                            strcpy(bestString,original_str);
                            finalize(originalGap,g);
                        }
                    }
                }
                else
                {
                    if(maxGapEstimate == 0)
                    {
                        //strcpy(bestString,original_str);
                        if(used_read_arr[maxGapEstimate] != 0)finalize(maxGapEstimate,g);
                        else
                        {
                            if(maxGapEstimate < originalGap)
                            {
                                run(originalGap,fp,g,1,0);
                                computeSequence(0,0);
                                strcpy(original_str,concensus);
                            }
                            finalize(originalGap,g);
                        }
                    }
                    else
                    {
                        if(side_flag)
                        {
                            //cout<<"Side break- Partial\n";
                            left_maxDistance = left_right_maxDistance[maxGapEstimate-gapMin][0];
                            right_maxDistance = left_right_maxDistance[maxGapEstimate-gapMin][1];
                        }
                        finalize(maxGapEstimate,g);
                    }

                }
            }
            //if(g>=gnum1 && g <=gnum2)cout<<concensus<<endl;
	    }
	}
    
    void freeGapFiller()
    {
        for(long int i=0;i<this->maxSize;i++)
        {
            delete [] countsGap[i];
            delete [] probsGap[i];
            delete [] errorProbsGap[i];
            delete [] new_counts_gap[i];
        }

        delete [] countsGap;

        delete [] probsGap;

        delete [] errorProbsGap;

        delete [] new_counts_gap;

        delete [] concensus;

        delete [] bestString;

        delete [] secondBestString;

        delete [] current_str;

        delete[] current_str2;

        delete [] previous_str;

        delete[] gap_left;

        delete[] gap_right;

        if(part_read_c>0)
        {
            for(int i=0;i<this->part_read_c;i++)
            {
                delete[] repeatflag[i];
            }
            delete[] repeatflag;
        }

        if(partial_read_count>0)
        {
            for(int i=0;i<partial_read_count;i++)
            {

                delete[] partial_read_pos_arr_org[i];

            }

            delete[] partial_read_pos_arr_org;
            delete[] maxlikelihood_value;

        }

        if(this->num_reads_in_gap>0)
        {
            for(int i=0;i<this->num_reads_in_gap;i++)
            {
                delete [] reads_gap[i];

                delete[] unmapped_read_pos_arr_org[i];
                delete[] final_readpos[i];
            }


            delete[] unmapped_read_pos_arr_org;
            delete[] reads_gap;

            delete[] pos_reads;

            delete[] isReverse;

            delete[] mark_accepted_reads;

            delete[] maxlikelihood_value;

            delete[] final_readpos;

            delete[] saved_reads;
        }

        delete[] gap_coverage;

        for(int i=0;i<maxGap;i++)
        {
            delete [] partial_prob_arr[i];
            delete [] qual_gap[i];
            delete [] partial_count_array[i];
        }

        delete[] partial_count_array;
        delete[] partial_prob_arr;
        delete[] qual_gap;

        if(this->p_count > 0)
        {
            for(int i=0;i<this->p_count;i++)
            {
                delete[] partialreads_gap[i];
                if(partial_flag)delete[] partial_read_quality[i];
            }
            delete[] partialpos_reads;
            delete[] preadflag;
            delete[] clipped_index;
            delete[] partialreads_gap;
            if(partial_flag)delete[] partial_read_quality;
        }
        
    }
};

int findcount_file(char *fname,int flag)
{
    FILE * scaffoldMap1=NULL;
    scaffoldMap1=fopen(fname, "r");

    char line1[MAX_REC_LEN];

    int MAX_FILE_READ=MAX_REC_LEN/sizeof(line1[0]);

    int r_c=0;

    while(fgets(line1, MAX_FILE_READ, scaffoldMap1)!=NULL)
    {
        if(flag == 0)
        {
            if(fgets(line1, MAX_FILE_READ, scaffoldMap1)==NULL)
                break;
        }
        r_c++;
        if(flag==1 && r_c > partial_limit)break;
    }
    

    fclose(scaffoldMap1);
    return r_c;
}

void checkComplete(char * gapstr,int r_c[4])
{
    int Nstart=0,region_count=0,c_left=0,c_right=0,lc=0,rc=0;

    int len = strlen(gapstr);
    Nstart=0,region_count=0;

    for(int i=0;i<len;i++)
    {
        if(gapstr[i] =='N' && Nstart ==0)
        {
            Nstart = 1;
            c_left=1;
        }
        else if(gapstr[i] != 'N' && Nstart ==1)
        {
            Nstart=0;
            region_count++;
            c_right=1;
        }
        if(i==len-1 && Nstart==1)
        {
            region_count++;
        }
        if(c_left==0 && gapstr[i] !='N')lc++;
        if(c_right==1 && gapstr[i] !='N')rc++;
    }
    r_c[0] = region_count;
    r_c[1] = lc;
    r_c[2] = rc;

    if(region_count != 0)r_c[3] = lc+rc;
    else r_c[3] = len;
}

void storeModel()
{
    FILE* out = fopen("model.txt","w");

    for(int i=0;i<maxReadLength;i++)
        fprintf(out,"%lf\n",inPosDist[i]);
    fprintf(out,">\n");

    for(int i=0;i<maxReadLength;i++)
        fprintf(out,"%lf\n",delPosDist[i]);
    fprintf(out,">\n");

    for(int i=0;i<maxReadLength;i++)
        fprintf(out,"%lf\n",errorPosDist[i]);
    fprintf(out,">\n");

    for(int i=0;i<5;i++)
    {
        for(int j=0;j<5;j++)
        {
            fprintf(out,"%lf\t",errorTypeProbs[i][j]);
        }
        fprintf(out,"\n");
    }
    fprintf(out,">\n");

    fprintf(out,"%lf\n",insertSizeMean);
    fprintf(out,"%lf\n",leftSD);
    fprintf(out,"%lf\n",rightSD);
    fprintf(out,"%d\n",gapProbCutOff);

    fclose(out);
}

void printModel()
{
    for(int i=0;i<maxReadLength;i++)
            printf("%lf\n",inPosDist[i]);
        printf("\n\n");

    for(int i=0;i<maxReadLength;i++)
            printf("%lf\n",delPosDist[i]);
        printf("\n\n");

    for(int i=0;i<maxReadLength;i++)
            printf("%lf\n",errorPosDist[i]);
        printf("\n\n");

    for(int i=0;i<5;i++)
        {
            for(int j=0;j<5;j++)
            {
                printf("%lf\t",errorTypeProbs[i][j]);
            }
            printf("\n");
        }
        printf("\n\n");

    printf("%lf\n",insertSizeMean);
    printf("%lf\n",leftSD);
    printf("%lf\n",rightSD);
    printf("%d\n",gapProbCutOff);
}

void loadModel()
{
    //zzzz
    inPosDist=new double[maxReadLength];
    delPosDist=new double[maxReadLength];
    errorPosDist=new double[maxReadLength];

    FILE* in = fopen("model.txt","r");
    char *line= new char[MAX_REC_LEN];
    int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);

    int linecount=0;

    while(fgets(line, MAX_FILE_READ, in)!=NULL)
    {
        if(line[0]=='>')break;
        line[strlen(line)] = '\0';
        inPosDist[linecount++] = atof(line);
    }

    linecount = 0;

    while(fgets(line, MAX_FILE_READ, in)!=NULL)
    {
        if(line[0]=='>')break;
        line[strlen(line)] = '\0';
        delPosDist[linecount++] = atof(line);
    }

    linecount = 0;

    while(fgets(line, MAX_FILE_READ, in)!=NULL)
    {
        if(line[0]=='>')break;
        line[strlen(line)] = '\0';
        errorPosDist[linecount++] = atof(line);
    }

    linecount = 0;
    char *temp = new char[100];

    while(fgets(line, MAX_FILE_READ, in)!=NULL)
    {
        if(line[0]=='>')break;
        line[strlen(line)] = '\0';

        temp=strtok(line,"\t");
        int colcount=0;

        while(temp != NULL)
        {
            errorTypeProbs[linecount][colcount++]= atof(temp);
            temp=strtok(NULL,"\t");
        }
        linecount++;
    }

    fscanf(in,"%lf\n",&insertSizeMean);
    fscanf(in,"%lf\n",&leftSD);
    fscanf(in,"%lf\n",&rightSD);
    fscanf(in,"%d\n",&gapProbCutOff);

    fclose(in);

    //printModel();
}

int findFrac(int originalGap,float *info,int mid_limitp,int mid_limitu)
{
    int ret_val;
    if(partial_flag)
    {
        if(originalGap<=mid_limitp/2){info[0]=.00001;info[1]=(float)factor/originalGap;ret_val=-1;}
        else if(originalGap<=mid_limitp){info[0]=.00001;info[1]=5.0;ret_val=5;}
        else
        {
            info[0]=1;info[1]=1;
            ret_val=3;
        }
    }
    else
    {
        if(originalGap<=mid_limitu/3){info[0]=.3;info[1]=(float)factor/originalGap;ret_val=-1;}
        else if(originalGap <= mid_limitu){info[0]=.5;info[1]=2.5;ret_val=3;}
        else
        {
            info[0]=1;info[1]=1;
            info[2]=1;
            ret_val=1;
        }
        //cout<<"If = "<<mid_limitu*0.15<<endl;
        //cout<<"Else If = "<<mid_limitu<<endl;
    }
    return ret_val;
}


int main(int argc, char *argv[])
{
	//==========================================================================
    //Declaration of all variables
    
    time_t nowb, endb;
    
    time(&nowb);
    
    char *line= new char[MAX_REC_LEN];
   	char *line1= new char[MAX_REC_LEN];
   	char *line2= new char[MAX_REC_LEN];
    char *line3= new char[MAX_REC_LEN];

   	int read;
   	
   	char dir[1000];
   	
   	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);
   	long int bufferLength=1024;

   	char *contig=new char[bufferLength];
   	contig[0]='\0';
   	char *newcontig;
   	char *contigName;
   	contigLength=0;
   	char* asm_name;

    long int contigNo;
    double readCount;
    char *temp;

   	long int tempContigLength=0;

   	long int totalInserts=0;
    long int insertSum=0;

    int nonzeroCount=0,read_length=0;
    int gapstartno,gapendno,thread_id,neg_lap,partial_len,totalgapstofill,unm_limit,setinputmeanflag;

    char *gapped_read_path;
    char *temp_path;

	//======================Declarations complete====================================================

    //======================Start reading the ref genome with gaps file==============================

    //Open gapped contig file for reading
	contigFileName=argv[1];
	tempmaxDistance=atoi(argv[2]);
    read_length = atoi(argv[3]);
	script_itr= atoi(argv[4]);
	partial_flag=atoi(argv[5]);
	unmapped=atoi(argv[6]);
	thread_id = atoi(argv[7]);
	totalgapstofill=atoi(argv[8]);
	mapFileName = argv[9];
	temp_path=argv[10];
	gapped_read_path=argv[11];
	neg_lap=atoi(argv[12]);
	partial_len=atoi(argv[13]);
    unm_limit = atoi(argv[14]);
    setinputmeanflag=atoi(argv[15]);
    
    if(setinputmeanflag == 1)inputMean = atoi(argv[16]);
	//cout<<gapstartno<<"\t"<<gapendno<<"\t"<<thread_id<<endl;


    //cases = fopen("Final_Case.txt","w");

	contigFile=fopen(contigFileName, "r");

	if (contigFile == NULL)
	{
		printf("Can't open contig file\n");
		exit(1);
	}
	while(fgets(line, MAX_FILE_READ, contigFile)!=NULL)
	{
		if(line[0]==';')
		{
			continue;
		}
		else if(line[0]=='>')
		{
			contigName=new char[strlen(line)];
			strcpy(contigName,line+1);
			contigName[strlen(contigName)-1]='\0';
			contigNames.push_back(strtok(contigName," \t\n"));

			if(contigLength>0)
			{
				noContigs++;
				contigs.push_back(contig);
				contigLengths.push_back(contigLength);
				totalContigLength+=contigLength;
				contigLength=0;
				bufferLength=1024;
				contig=new char[bufferLength];
				contig[0]='\0';
			}
		}
		else
		{
			read=strlen(line);
			tempContigLength=contigLength;
			if(read<MAX_FILE_READ-1)
			{
				contigLength+=(read-1);
			}
			else
			{
				contigLength+=MAX_FILE_READ-1;
				read++;

			}
			if(contigLength>bufferLength)
			{
				bufferLength=max(bufferLength*2,contigLength+1);
				newcontig=new char[bufferLength];
				strcpy(newcontig,contig);
				line[read-1]='\0';
				strcat(newcontig, line);
				delete []contig;
				contig=newcontig;
			}
			else
			{
				line[read-1]='\0';
				strcpy(contig+tempContigLength, line);
			}
		}
	}

	noContigs++;
	contigs.push_back(contig);
	contigLengths.push_back(contigLength);
	totalContigLength+=contigLength;

    fclose(contigFile);

    //=====================Finished with reading the contig file===========================
    
	for(long int i=0;i<noContigs;i++)
	{
		for(long int j=0;j<contigLengths[i];j++)
		{
			contigs[i][j]=toupper(contigs[i][j]);
		}
	}

	for(int i=0;i<256;i++)
    {
        if(i=='A')
        {
            charCodes[i]=0;
        }
        else if(i=='C')
        {
            charCodes[i]=1;
        }
        else if(i=='G')
        {
            charCodes[i]=2;
        }
        else if(i=='T')
        {
            charCodes[i]=3;
        }
        else
        {
            charCodes[i]=4;
        }
    }

    strcpy(dir,temp_path);
    strcat(dir,"stat.txt");

    summaryFile=fopen(dir,"r");
    fscanf(summaryFile,"%ld %ld %d %d",&totalCount, &unCount, &maxReadLength, &MAX_INSERT_SIZE);
    fclose(summaryFile);

    MAX_INSERT_SIZE=MAX_INSERT_SIZE>20000?MAX_INSERT_SIZE:20000;

    insertCutoffMin=MAX_INSERT_SIZE;//20k

    //There are two variables,MAX_INSERT_SIZE and maxInsertSize,
    //at this fn, initInsertCounts, 2nd one is initialized by 1st one .i.e. 20k
	initInsertCounts(MAX_INSERT_SIZE);
    initErrorTypes(maxReadLength);
    
	sprintf(noErrorCigar,"%d",maxReadLength);
    strcpy(noErrorMD,"MD:Z:");
    strcat(noErrorMD,noErrorCigar);
    strcat(noErrorCigar,"M");

    //So, noErrorMD = MD:Z:100                   noErrorCigar = 100M

    //============================Start building/loading the model=======================

	
    mapFile=fopen(mapFileName, "r");
    //cout<<"out.sam filename = "<<mapFileName<<endl;
    
    if (mapFile == NULL)
    {
	    printf("Can't open map file\n");
	    exit(1);
    }
    while(fgets(line1, MAX_FILE_READ, mapFile)!=NULL)
    {
	
	    if(line1[0]=='@')
		    continue;
        
	    processMapping(line1);
	    
    }
    //cout<<"Myout.sam file has "<<count<<" lines"<<endl;
    fclose(mapFile);

    //===========================Done with myout.sam========================================

    computeProbabilites();
    double val1 = computeLikelihood(mapFileName);
   // cout<<"likelihood //val 1  = "<<val1<<endl;
    
    //FILE * mappedInsertFile=fopen("mappedInserts.txt","w");
    
    for(long int i=0;i<maxInsertSize;i++)
    {
        if(insertCountsMapped[i]>0)
        {
            //    cout<<i<<","<<insertCountsAll[i]<<endl;
            //fprintf(mappedInsertFile, "%ld,%ld\n",i,insertCountsMapped[i]);
            totalInserts+=insertCountsMapped[i];
            insertSum+=insertCountsMapped[i]*i;
            nonzeroCount++;
        }
    }

    double insertMean=insertSum/(double)totalInserts;
    //cout<<"Total Inserts = "<<totalInserts<<endl;
    
    //fclose(mappedInsertFile);
    
    long int gapProbSum=0;

    for(int i=0;i<1000;i++)
    {
        gapProbSum+=gapProbs[i];
        //cout<<gapProbs[i]<<endl;
    }

    long int gapProbCount=0;
    double value;

    value = .8;

    //cout<<"gap_prob_cutoff fraction = "<<value<<endl;
    //zxv
    for(int i=0;i<1000;i++)
    {
        gapProbCount+=gapProbs[i];
        if(gapProbCount>=value*gapProbSum)
        {
            gapProbCutOff=i;
            break;
        }
    }

    //cout<<"gap_prob_cutoff in main = "<<gapProbCutOff<<endl;

    time(&endb);
    
    //cout<<"Time to build the model = "<<difftime(endb,nowb)<<" seconds"<<endl;
    
    

    left_coeff = 3;
    right_coeff =3;
    

    //zx
    insertThresholdMin=max((int)(insertSizeMean-left_coeff*leftSD),1);
    insertThresholdMax=min((int)(insertSizeMean+right_coeff*rightSD),maxInsertSize);
    
    if(partial_flag)
    {
        insertThresholdMin -= partial_len;
        insertThresholdMax += partial_len;
    }

    //cout<<"isz mean = "<<insertSizeMean<<"\tinsertThresholdMin = "<<insertThresholdMin<<"\tinsertThresholdMax = "<<insertThresholdMax<<endl;

    //==============Done with insertthteshold========================================

    /*gap_prob_dist = fopen("gap_prob_dist.txt","w");

    for(int i=0;i<1000;i++)
    {
        fprintf(gap_prob_dist,"%d\t%ld\n",i,gapProbs[i]);
    }

    fclose(gap_prob_dist);
    */

    time_t now, end;
    time(&now);

    //====================================Start filling the gaps from here===================
    int gapContigNo;
    long int gapStart;
    int gapLength;
    int gapStringLength;
    int diff;

    char gapFileName[500];
    char partialgapFileName[500];
    char *gapString=new char[MAX_GAP];

    temp=new char[100];
    char *temp2=new char[100];
    
    int gapNo=0;

    strcpy(dir,temp_path);
    strcat(dir,"gapInfo.txt");

    gapInfoFile=fopen(dir,"r");

    char utilfilename[100];
    sprintf(temp2,"%d",thread_id);

    strcpy(utilfilename,temp_path);
    strcat(utilfilename,"draw");
    strcat(utilfilename,temp2);
    strcat(utilfilename,".txt");
    draw = fopen(utilfilename,"w");

    strcpy(utilfilename,temp_path);
    strcat(utilfilename,"gapout");
    strcat(utilfilename,temp2);
    strcat(utilfilename,".txt");
    gapOutFile=fopen(utilfilename,"w");


    strcpy(dir,temp_path);
    strcat(dir,"stat2.txt");
    stat2 = fopen(dir,"r");
    
    int tot_gl=0,tot_filled_gl=0,r_count1=0,r_count2=0,mainstat1,mainstat2,mainstat3,tot_gaps=0;
    while(fgets(line1, MAX_FILE_READ, gapInfoFile)!=NULL)tot_gaps++;
    
    gaptofill = new int[tot_gaps];
    finalUsedReads = new int[2*tot_gaps];
    timereq = new double[tot_gaps];

    for(int i=0;i<tot_gaps;i++)gaptofill[i] = 0;
    for(int i=0;i<2*tot_gaps;i++)finalUsedReads[i] = 0;
    for(int i=0;i<tot_gaps;i++)timereq[i] = 0;

    fclose(gapInfoFile);
    
    strcpy(dir,temp_path);
    strcat(dir,"gapInfo.txt");
    gapInfoFile=fopen(dir,"r");

    int gappedthread[totalgapstofill];
    
    strcpy(dir,temp_path);
    strcat(dir,"gaploads.txt");
        
    
    FILE* loadgap = fopen(dir,"r");
    
    int linecountload=0,gcount=0;
    
    while(1)
    {
        if(linecountload == thread_id)
        {
            while(gcount< totalgapstofill)
            {
                fscanf(loadgap,"%d\t",&gappedthread[gcount]);
                gcount++;
            }
            //fscanf(loadgap,"\n");
        }
        else
        {
            while(1)
            {
                fgets(line, MAX_FILE_READ, loadgap);
                string str(line);
                char c = '\n';

                // Find first occurrence of 'g'
                size_t found = str.find(c);
                if (found != string::npos)
                    break;
            }
        }
        linecountload++;
        
        if(linecountload > thread_id)break;
    }
    
    fclose(loadgap);
    
    //cout<<totalgapstofill<<"\t"<<gcount<<endl;
    //for(int i=0;i<totalgapstofill;i++)cout<<gappedthread[i]<<"\t";
    //cout<<thread_id<<endl;
    
    
    int gapsdone=0,fillperct=0;
    
    float info[3];
    int trackgap=0,perctstep=20;

    while(fgets(line1, MAX_FILE_READ, gapInfoFile)!=NULL && fgets(line3, MAX_FILE_READ, stat2)!=NULL)
    {
        for(int i=0;i<3;i++)info[i]=0;
        
        if(gapNo == gappedthread[trackgap])
        {
            time_t nowt,endt;
            time(&nowt);
            
            trackgap++;
             
            int fillflag=1;

            left_maxDistance = tempmaxDistance;
            right_maxDistance = tempmaxDistance;

        	strcpy(gapFileName,gapped_read_path);
            strcpy(partialgapFileName,gapped_read_path);

            strcat(gapFileName,"gaps_");
            strcat(partialgapFileName,"partial_gaps_");
            //itoa(gapNo,temp,10);
            sprintf(temp2,"%d",gapNo);
            //cout<<"Check"<<endl;
            strcat(gapFileName,temp2);
            strcat(partialgapFileName,temp2);

            strcat(gapFileName,".sam");
            strcat(partialgapFileName,".sam");

            temp=strtok(line1,"\t");
            gapContigNo=atoi(temp);

            temp=strtok(NULL,"\t");
            gapStart=atol(temp);

            temp=strtok(NULL,"\t\n");
            gapLength=atoi(temp);

            //=============new from main============
            temp=strtok(line3,"\t");
            mainstat1 = atoi(temp);
            temp=strtok(NULL,"\t");
            mainstat2 = atoi(temp);
            temp=strtok(NULL,"\t");
            mainstat3 = atoi(temp);
            //===================================
            //cout<<partial_flag<<" "<<unmapped<<endl;
            if(unmapped==1)
            {
                r_count1 = findcount_file(gapFileName,0);
                if(r_count1 > unmapped_limit)
                {
                    cout<<"Gap = "<<gapNo<<"\tReads = "<<r_count1<<endl;
                    r_count1 = unmapped_limit;
                    fillflag = -1;
                }
            }
            if(partial_flag==1)r_count2 = findcount_file(partialgapFileName,1);
            //cout<<"unampped =  "<<r_count1<<"\tpartial = "<<r_count2<<endl;
            //cout<<"Allocating length of size = "<<gapLength<<endl;
            
            GapFiller gf;
            
            factor=3*partial_len;
            
            int allocation_factor = findFrac(gapLength,info,2*partial_len,unm_limit);

            if(allocation_factor==-1)alloc_arg = factor*3;
            else alloc_arg = gapLength*allocation_factor;

            gf.allocate(alloc_arg,r_count1,r_count2,read_length,mainstat1,mainstat2,mainstat3,neg_lap,2*partial_len,unm_limit,info[2],info[0],info[1]);
            
            //cout<<"After allocate"<<endl;
		    gf.initGapFiller(gapNo,gapContigNo,gapStart,gapLength, gapFileName,partialgapFileName);
		    //cout<<"After init- "<<gapNo<<endl;

            
            gf.fillGap(gapNo,draw,fillflag);
    //zxcf
            //cout<<"After fillgap"<<endl;

		    gapStringLength = gf.getConcensus(gapString,gapLength);

		    fprintf(gapOutFile,"%d\t%d\t%ld\t%d\t%d\t%s\n",gapNo,gapContigNo,gapStart,gapLength,gapStringLength,gapString);
            
            tot_gl += gapLength;
            tot_filled_gl += gapStringLength;
            //cout<<"Before free\n";
		
		    gf.freeGapFiller();
		    //cout<<"After free\n";

		    time(&endt);
		    timereq[gapNo] = difftime(endt,nowt);
		    //cout<<"Gap = "<<gapNo<<"\ttime = "<<timereq[gapNo]<<endl;
		    gapsdone++;
            
            fillperct = 100*gapsdone/totalgapstofill;
		    
		    if(fillperct >= perctstep)
		    {   
		        cout<<"Thread - "<<thread_id<<"\tGaps Completed = "<<perctstep<<"%"<<endl;
		        perctstep += 20;
		    }
		}
		else
		{
		    gaptofill[gapNo] = -1;
		}
		gapNo++;
    }

    //fprintf(reads_per_gap,"Total Gaps = %d ,Total Reads = %d, Total Discarded = %d\n\n",gapNo,total_read,discarded_read);
    //fprintf(reads_per_gap,"Total Gap Length = %d, Total_filled Gap Length = %d",tot_gl,tot_filled_gl);

    fclose(gapInfoFile);
    fclose(gapOutFile);
    
    //fclose(minmax1);
    //fclose(minmax2);

    fclose(draw);
    //fclose(cases);
    

    time(&end);
    
    sprintf(temp2,"%d",thread_id);
    strcpy(utilfilename,temp_path);
    strcat(utilfilename,"gaptofill");
    strcat(utilfilename,temp2);
    strcat(utilfilename,".txt");

    FILE* gaptf = fopen(utilfilename,"w");
    for(int i=0;i<tot_gaps;i++)
        fprintf(gaptf,"%d\n",gaptofill[i]);
    fclose(gaptf);

    double seconds = difftime(end,now);



    delete[] gaptofill;
    delete[] finalUsedReads;
    delete[] timereq;

    delete[] line;
    delete[] line1;
    delete[] line2;
    delete[] line3;

    delete[] gapString;

    delete[] temp2;

    delete[] insertCounts;

    delete[] errorPos;
    delete[] inPos;
    delete[] inLengths;
    delete[] delPos;
    delete[] delLengths;
    delete[] readLengths;

    delete[] errorPosDist;
    delete[] inPosDist;
    delete[] delPosDist;
    delete[] inLengthDist;
    delete[] delLengthDist;
    delete[] insertLengthDist;
    delete[] insertLengthDistSmoothed;

    delete[] noErrorProbs;
    delete[] effectiveLengths;
    delete[] insertCountsMapped;
    delete[] insertCountsUnmapped;

	return 0;
}

