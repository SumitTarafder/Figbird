#include <iostream>
#include <bits/stdc++.h>
#include <string>
#include <limits>
#include <vector>
#include <cstring>

using namespace std;
#include <math.h>

#define MAX_REC_LEN 1024

///=================================
int main(int argc, char *argv[])
{
    FILE * contigFile=NULL,*mergedContigFile=NULL;
    char* contigFileName = argv[1];
    char * temp_path;
    temp_path=argv[2];
    
    char dir[1000];

    //cout<<contigFileName<<endl;

    contigFile=fopen(contigFileName, "r");

    if (contigFile == NULL)
    {
        cerr<<"Can't open gapped genome file during reduction\n";
        exit(1);
    }
    
    strcpy(dir,temp_path);
    strcat(dir,"newgenome.fa");

    mergedContigFile = fopen(dir,"w");

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
			if(contigLength>0)//This block is accessed after completing the full reading of each contig, not at the start
            {
                noContigs++;
                contig=(char *)realloc(contig, (contigLength+1)*sizeof(char));

                if(Nflag ==1)
                {
                    Nflag=0;
                    fprintf(mergedContigFile,">%s\n",contigName);
                    fprintf(mergedContigFile,"%s\n",contig);
                    gappedcontigs++;
                }
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

			if(Nflag ==0)
			{
                for(int z=0;z<read;z++)
                {
                    if(line[z] == 'N' || line[z] == 'n')
                    {
                        Nflag = 1;
                        break;
                    }
                }
            }

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

    if(Nflag ==1)
    {
        Nflag=0;
        fprintf(mergedContigFile,">%s\n",contigName);
        fprintf(mergedContigFile,"%s\n",contig);
        gappedcontigs++;
    }
    //cout<<"total = "<<noContigs<<"\tGapped contigs = "<<gappedcontigs<<endl;


    fclose(contigFile);
    fclose(mergedContigFile);

}
