#include <iostream>
#include <bits/stdc++.h>
#include <string>
#include <limits>
#include <vector>
#include <cstring>

using namespace std;

#include <math.h>

#define MAX_REC_LEN 1024
#define MAX_GAPLEN 100000

vector<char*> contigs;

vector<unsigned long> contigLengths;
vector<char*> contigNames;

///=================================
int main(int argc, char *argv[])
{
    FILE * contigFile=NULL,*trimmedContigFile=NULL;
    FILE* flankfile = NULL;
    char* contigFileName = argv[1];
    
    int trimsize=atoi(argv[2]);
    int readlen = atoi(argv[3]);
    char* trimmedfile=argv[4];
    //char* flankfilename=argv[5];
    //cout<<contigFileName<<endl;
    //cout<<readlen<<endl;
    
    contigFile=fopen(contigFileName, "r");

    if (contigFile == NULL)
    {
        printf("Can't open gapped genome file\n");
        exit(1);
    }

    trimmedContigFile = fopen(trimmedfile,"w");
    
    //flankfile = fopen(flankfilename,"w");

    char *line= new char[MAX_REC_LEN];

    int MAX_FILE_READ = 0;
	int noContigs=0,gappedcontigs=0;

	long int contigLength=0;
    unsigned long read;
    long int bufferLength=1024;
    long int tempContigLength=0;

    char *newcontig;
    char *contigName;

    char *contig;
    contig=new char[bufferLength];
    contig[0]='\0';

	MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);
    int Nflag=0;

    char *temp=new char[100];
    
    
//==================================================================
    while(fgets(line, MAX_FILE_READ, contigFile)!=NULL)
    {
        if(line[0]==';')//There is no ';' in contig file, why this?
		{
			continue;
		}
		else if(line[0]=='>')//fasta format name of the contig
		{
			contigName = new char[strlen(line)];
			strcpy(contigName,line+1);
			contigName[strlen(contigName)-1]='\0';
			contigNames.push_back(strtok(contigName," \t\n"));
			
			if(contigLength>0)//This block is accessed after completing the full reading of each contig, not at the start
            {
                noContigs++;
                contig=(char *)realloc(contig, (contigLength+1)*sizeof(char));
                contigs.push_back(contig);
				contigLengths.push_back(contigLength);

                contigLength=0;
                bufferLength=1024;
                contig=new char[bufferLength];
                contig[0]='\0';
            }

            contigName = new char[strlen(line)];
            strcpy(contigName,line+1);
            contigName[strlen(contigName)-1]='\0';

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
                strcpy(newcontig+tempContigLength, line);
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
	contig=(char *)realloc(contig, (contigLength+1)*sizeof(char));
	contigs.push_back(contig);
	contigLengths.push_back(contigLength);

    fclose(contigFile);

    int gapcount=0,gapStart,gapLength,gapContigNo;
    int nStart=0;//How many gaps are there?
    unsigned long nStartPos=0,nCount;//What are those gaps start positions?

    int offset,index;
    int dummy=-1;
    
    char leftflank[trimsize+1],rightflank[trimsize+1];
    
    int trimmed=0,condition=0;
    
    for(int i=0;i<noContigs;i++)//For each contig
    {
        fprintf(trimmedContigFile,">%s\n",contigNames[i]);
    	
        for(long int j=0;j<contigLengths[i];j++)//For the entire contig length
        {
            condition = (contigs[i][j]=='N' || contigs[i][j]=='n');
            
            if(condition)//If there is a gap at that position
            {
                if(nStart==0)
                {
                    nStart=1;
                    nCount=1;
                    nStartPos=j;
                    gapcount++;
                }
                else
                {
                    nCount++;
                }
            }
            
            if((!condition && nStart==1) || (condition && j == contigLengths[i] - 1))//if there was ever any gap
            {
                if(trimsize > 0 && nCount > 1 && (int)nCount < readlen && (int)nStartPos - trimsize > 2*trimsize 
                    && (contigLengths[i] - nStartPos - nCount) > 2*trimsize)
                {
                    for(int t=0;t<trimsize;t++)
                    {
                        leftflank[trimsize-t-1] = contigs[i][nStartPos-t-1];
                        //contigs[i][nStartPos-t-1]='N';
                    }
                    for(int t=0;t<trimsize;t++)
                    {
                        rightflank[t] = contigs[i][nStartPos+nCount+t];
                        //contigs[i][nStartPos+nCount+t]='N';
                    }
                    
                    leftflank[trimsize]='\0';
                    rightflank[trimsize]='\0';
                    
                    if(!strpbrk(leftflank,"N") && !strpbrk(rightflank,"N"))//No N in flanks
                    {  
                        for(int t=0;t<trimsize;t++)
                        {
                          contigs[i][nStartPos-t-1]='N';
                          contigs[i][nStartPos+nCount+t]='N';
                        }
                        
                        j += trimsize;
                    
                        //fprintf(flankfile,"%d\n",gapcount);
                        //fprintf(flankfile,"%s\n",leftflank);
                        //fprintf(flankfile,"%s\n",rightflank);
                        
                        trimmed++;
                    }
                    else 
                    {
                        //fprintf(flankfile,"%d\n",dummy);
                    }
                }
                else 
                {
                    //fprintf(flankfile,"%d\n",dummy);
                }
                nStart=0;
            }
        }
        
        fprintf(trimmedContigFile,"%s\n",contigs[i]);
    }
    
    fclose(trimmedContigFile);
    //fclose(flankfile);
    
    contigLengths.clear();
    contigs.clear();
    contigNames.clear();

}
