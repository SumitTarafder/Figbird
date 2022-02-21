#include <stdio.h>
#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <cmath>
#include <cfloat>
#include<cstdlib>
#include <assert.h>
#include <bits/stdc++.h>
#include<unistd.h>
using namespace std;

#define _USE_MATH_DEFINES
#include <thread>
#include <math.h>
#include <time.h>

#define MAX_REC_LEN 1024
#define MAX_READLENGTH 200

int gapthresh=400;

int noContigs=0;
int noReads=0;
long int contigLength=0;
int maxReadLength=0;

long int totalContigLength=0;
vector<char*> contigs;
vector<char*> contigNames;
vector<long int> contigLengths;
double *contigReadCounts;
const int MAX_GAP_LEN=100000;
int *gaptofill;
char command[10000];

long int gapStart;
int gapNo,gapContigNo;

char *contigFileName;
char* samoutfilename;
char* tp,*grp;

int tempmaxDistance,gtf,read_length,script_itr,partial_flag,unmapped,neg_overlap,partial_readlen,trim;

int **gapallocate;
int g2fcount=0;


void thread_fillgap(int tid,int read_length,int load)
{
    char *temp2=new char[100];
    char exefile[1000];

    sprintf(temp2,"%d",tid);

    //strcpy(exefile,tp);

    strcpy(exefile,"a");
    strcat(exefile,temp2);
    strcat(exefile,".out");

    strcpy(command,"g++ Figbird.cpp -o ");
    strcat(command,exefile);
	strcat(command," && ./");

    strcat(command,exefile);
    strcat(command," ");

	strcat(command,contigFileName);
	strcat(command," ");

	sprintf(temp2,"%d",tempmaxDistance);
	strcat(command,temp2);
	strcat(command," ");

	sprintf(temp2,"%d",read_length);
	strcat(command,temp2);
	strcat(command," ");

	sprintf(temp2,"%d",script_itr);
	strcat(command,temp2);
	strcat(command," ");

	sprintf(temp2,"%d",partial_flag);
	strcat(command,temp2);
	strcat(command," ");

	sprintf(temp2,"%d",unmapped);
	strcat(command,temp2);
    strcat(command," ");

	sprintf(temp2,"%d",tid);
    strcat(command,temp2);
	strcat(command," ");
	
	sprintf(temp2,"%d",load);
    strcat(command,temp2);
	strcat(command," ");

    strcat(command,samoutfilename);
    strcat(command," ");

    strcat(command,tp);
    strcat(command," ");

    strcat(command,grp);
    strcat(command," ");

    sprintf(temp2,"%d",neg_overlap);
    strcat(command,temp2);
    strcat(command," ");
    
    sprintf(temp2,"%d",partial_readlen);
    strcat(command,temp2);
    strcat(command," ");
    
    sprintf(temp2,"%d",gapthresh);
    strcat(command,temp2);
    
    //cout<<"In thread - "<<tid<<endl;
    //cout<<command<<endl;

    system(command);
    remove(exefile);
    
    delete[] temp2;

}

void readgapouts(char * fname, int num_threads,int totalgaps)
{
    char filename[100];
    char *temp=new char[100];
    
    typedef struct Gapstruct
    {
        int gapno;
        int gapcontigno;
        long int gapstart;
        int gaplength;
        int gapstringlength;
        char *gapstr;
    }Gap;
    
    Gap g[totalgaps];
    
    int gapNo,gapContigNo,gapLength,gapStringLength;
    long int gapStart;
    
    
    for(int i=0;i<num_threads;i++)
    {
        strcpy(filename,fname);
        sprintf(temp,"%d",i);

        strcat(filename,temp);
        strcat(filename,".txt");
        
        FILE* fp = fopen(filename,"r");
        
        while(fscanf(fp,"%d\t%d\t%ld\t%d\t%d\t",&gapNo,&gapContigNo,&gapStart,&gapLength,&gapStringLength) != -1)
        {
       		//cout<<gapNo<<" "<<gapContigNo<<" "<<gapStringLength<<" "<<gapLength<<endl;
       		
       		g[gapNo].gapno = gapNo;
       		g[gapNo].gapcontigno = gapContigNo;
       	    g[gapNo].gapstart = gapStart;
            g[gapNo].gaplength = gapLength;
            g[gapNo].gapstringlength = gapStringLength;		

       		if(gapStringLength > 0)
       		{
       		    g[gapNo].gapstr = new char[gapStringLength+1];
       		    fscanf(fp,"%s\n",g[gapNo].gapstr);
       		}
       		else
       		{
       		    fscanf(fp,"\n");
       		}
        }
        fclose(fp);
        
        remove(filename);
    }
    
    
    strcpy(filename,fname);
    strcat(filename,".txt");
    
    FILE* fp = fopen(filename,"w");
    
    for(int i=0;i<totalgaps;i++)
    {
        fprintf(fp,"%d\t%d\t%ld\t%d\t%d\t",g[i].gapno,g[i].gapcontigno,g[i].gapstart,g[i].gaplength,g[i].gapstringlength);
        
        if(g[i].gapstringlength > 0)fprintf(fp,"%s\n",g[i].gapstr);
        else fprintf(fp,"\n");
    }        
    
    fclose(fp);
    
    for(int i=0;i<totalgaps;i++)
    {
        if(g[i].gapstringlength > 0)
            delete[] g[i].gapstr;
    }
    
    delete[] temp;
}


void mergeFiles(char * fname, int num_threads)
{
    char filename[100],inifilename[100],finalfilename[100],Gname[100];

    strcpy(finalfilename,fname);
    strcat(finalfilename,".txt");

    strcpy(Gname,tp);
    strcat(Gname,"G.txt");

    char *temp2=new char[100];

    strcpy(inifilename,fname);
    strcat(inifilename,"0.txt");

    for(int i=1;i<num_threads;i++)
    {
        strcpy(filename,fname);
        sprintf(temp2,"%d",i);

        strcat(filename,temp2);
        strcat(filename,".txt");

        std::ifstream if_a(inifilename, std::ios_base::binary);
        std::ifstream if_b(filename, std::ios_base::binary);
        std::ofstream of_c(Gname, std::ios_base::binary);

        of_c << if_a.rdbuf() << if_b.rdbuf();
        remove(filename);

        rename(Gname,inifilename);
    }

    rename(inifilename,finalfilename);
     delete[] temp2;
     
}


void readGtf(int tid)
{
    char filename[100];
    char *temp2 = new char[100];

    sprintf(temp2,"%d",tid);
    strcpy(filename,tp);
    strcat(filename,"gaptofill");
    strcat(filename,temp2);
    strcat(filename,".txt");

    char *line= new char[MAX_REC_LEN];

   	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);

    //cout<<filename<<endl;

    FILE* fp = fopen(filename,"r");

    if(fp==NULL)
    {
        cerr<<"Couldn't find threaded gap filled files to merge in Fillgaps.cpp\n";
        exit(1);
    }

    int gapno=0,val;
    
    while(fscanf(fp,"%d\n",&val) != -1)
    {
        if(val != -1)
        {
            gaptofill[gapno]= val;
            //cout<<val<<endl;
            g2fcount++;
        }
        gapno++;
    }

    fclose(fp);

    remove(filename);
    
    delete[] temp2;
    delete[] line;
    
}

static bool sortcol( const vector<int>& v1,const vector<int>& v2 )
{
    return v1[1] > v2[1];
}

void writeGapLoad(int t, int* arr)
{
    char loadname[1000];

    strcpy(loadname,tp);
    strcat(loadname,"gaploads.txt");
    
    FILE* fp = fopen(loadname,"w");
    
    for(int i = 0; i < t; i++)
        sort(gapallocate[i],gapallocate[i]+ arr[i]);
        
    for(int i = 0; i < t; i++)
    {
        for(int j = 0;j<arr[i];j++)
        {
            fprintf(fp,"%d\t",gapallocate[i][j]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}

void printAllocate(int row,int col)
{
    
    for(int i=0;i< row;i++)
    {
        cout<<"Thread - "<<i<<"\t";
        for(int j=0;j< col;j++)
            cout<<gapallocate[i][j]<<"\t";
        cout<<endl;
    }
    cout<<"=========================="<<endl;
}

void printloadrem(vector<vector <int>> arr ,int t)
{
    for(int i=0;i< t;i++)
    {
        for(int j=0;j< 2;j++)
            cout<<arr[i][j]<<"\t";
        cout<<endl;
    }
}

int findNcount(char * str)
{
    int Ncount=0;
    
    for(int i=0;i<strlen(str);i++)
    {
        if(str[i] == 'N')Ncount++;
    }
    return Ncount;
}


int main(int argc, char *argv[])
{
    time_t nowt,endt;
    time(&nowt);

    ///==========================================================================
    ///Declaration of all variables

    char *line= new char[MAX_REC_LEN];
    char *line1= new char[MAX_REC_LEN];
   	char *line2= new char[MAX_REC_LEN];
    char *line3= new char[MAX_REC_LEN];
    
    char dir[1000];

   	int read;
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

   	long int tempContigLength=0;

   	long int totalInserts=0;
    long int insertSum=0;

    int nonzeroCount=0;


   	int gapLength;
    int gapStringLength;

    char gapString[MAX_GAP_LEN];
    char *temp = new char[100];

    int tot_gaps=0;
    int num_threads;
    
    //===========================================================================
    //Read command line arguments
	contigFileName=argv[1];
	tempmaxDistance=atoi(argv[2]);
	read_length = atoi(argv[3]);
	script_itr= atoi(argv[4]);
	partial_flag=atoi(argv[5]);
	unmapped=atoi(argv[6]);
	num_threads = atoi(argv[7]);
	samoutfilename=argv[8];
	tp=argv[9];
	grp=argv[10];
    neg_overlap=atoi(argv[11]);
    partial_readlen=atoi(argv[12]);
    trim=atoi(argv[13]);

    //=================More processing...===================
    
    strcpy(dir,tp);
    strcat(dir,"gapInfo.txt");

    FILE* gapInfoFile=fopen(dir,"r");

    if(gapInfoFile==NULL)
    {
        cerr<<"Couldn't open gapinfo in Fillgaps.cpp\n";
        exit(1);
    }
    
    while(fgets(line1, MAX_FILE_READ, gapInfoFile)!=NULL)tot_gaps++;
    fclose(gapInfoFile);

    gaptofill = new int[tot_gaps];
    for(int i=0;i<tot_gaps;i++)gaptofill[i] = 0;

    ///=========================Start Multi-threaded Gapfilling================================
  
    int gap_per_thread=1;
    int sumalloc=0;
    cout<<"Total # of gaps = "<<tot_gaps<<endl;
    
    if(tot_gaps <= num_threads)
    {
        num_threads = tot_gaps;
        gap_per_thread = 1;
    }
    else
    {
        float t_f = (tot_gaps * 1.0 / num_threads);
        gap_per_thread = t_f;

        if(t_f - float(gap_per_thread) > 0)
            gap_per_thread++ ;
        //gap_per_thread = tot_gaps / num_threads + 1 ;
    }
    
    //cout<<"Thread count = "<<num_threads<<endl;
    //cout<<"Per thread gap count = "<<gap_per_thread<<endl;
    
    int threadstat[num_threads][2];
    
    gapallocate = new int*[num_threads];
        
    for(int i=0;i< num_threads;i++)
    {
        gapallocate[i] = new int[gap_per_thread];
        for(int j=0;j< gap_per_thread;j++)gapallocate[i][j] = -1;
        for(int j=0;j < 2;j++)threadstat[i][j] = 0;
    }
    
    thread t[num_threads];
    int loadperthread[num_threads];
    for(int i=0;i<num_threads;i++)loadperthread[i] = 0;
    
    if(tot_gaps <= num_threads)
    {    
        for(int i=0;i<num_threads;i++)
        {
            gapallocate[i][0] = i;
            loadperthread[i] = 1;
            sumalloc++;
        }
    }
    else
    {
        int cnt=0;
        int smallgap=0,largegap=0;
        long int smallarr[tot_gaps],largearr[tot_gaps];
        
        strcpy(dir,tp);
        strcat(dir,"gapInfo.txt");
        
        gapInfoFile=fopen(dir,"r");
        
        int gaplen;
        
        while(fgets(line1, MAX_FILE_READ, gapInfoFile)!=NULL)
        {
            temp=strtok(line1,"\t");
            temp=strtok(NULL,"\t");
            
            temp=strtok(NULL,"\t\n");
            gaplen=atoi(temp);
            
            if(gaplen > gapthresh)
            {
               largearr[largegap++] = cnt;//Y
            }
            else//<=300
            {
                smallarr[smallgap++] = cnt;//X
            }
            cnt++;
        }
        
        fclose(gapInfoFile);
        
        //cout<<"Critical gap count\t"<<smallgap<<endl;
        
        //for(int i=0;i<smallgap;i++)cout<<smallarr[i]<<"\t";
        //cout<<endl;
        
        //cout<<"Easy gap count = \t"<<largegap<<endl;
        
        //for(int i=0;i<largegap;i++)cout<<largearr[i]<<"\t";
        //cout<<endl;
        
        int totalalloc=0,tempsmall=smallgap;
        
        if(smallgap <= num_threads)
        {
            for(int i=0;i<num_threads;i++)
            {
                if(totalalloc == smallgap)break;
                gapallocate[i][loadperthread[i]++] = smallarr[i];
                totalalloc++;
                sumalloc++;
                threadstat[i][0]++;
            }
        }
        else
        {
            while(1)
            {
                for(int i=0;i<num_threads;i++)
                {
                    gapallocate[i][loadperthread[i]++] = smallarr[totalalloc++];
                    sumalloc++;
                    threadstat[i][0]++;
                    
                    if(totalalloc == smallgap)
                    {
                        break;
                    }
                }
                
                if(totalalloc == smallgap)
                {
                    break;
                }
            }
        }
        
        //for (int i = 0; i < num_threads; i++)
            //cout<<"Thread = "<<i<<"\tLoad = "<<loadperthread[i]<<endl;
        
        vector<vector<int>> loadrem(num_threads);
        
        for (int i = 0; i < num_threads; i++)
        {
            loadrem[i] = vector<int>(2);
            for (int j = 0; j < 2; j++)
                loadrem[i][j] = 0;
        }

        for(int i=0;i<num_threads;i++)
        {
            loadrem[i][0] = i;
            loadrem[i][1] = gap_per_thread - loadperthread[i];
        }
        
        
        //printloadrem(loadrem,num_threads);
       
        sort(loadrem.begin(), loadrem.end(),sortcol);
        
        //cout<<"Load rem after sort\n"<<endl;
        
        //printloadrem(loadrem,num_threads);
       
        totalalloc = 0;
        
        if(largegap <= num_threads)
        {
            while(1)
            {
                for(int i=0;i<num_threads;i++)
                {
                    if(totalalloc == largegap)break;
                    int threadno = loadrem[i][0];
                    
                    gapallocate[threadno][loadperthread[threadno]++] = largearr[totalalloc++];

                    sumalloc++;
                    threadstat[threadno][1]++;
                        
                }
                if(totalalloc == largegap)break;
            }
        }
        else
        {
            int flag=0;
            for(int i=0;i<num_threads;i++)
            {
                for(int j=0;j<loadrem[i][1];j++)
                {   
                    int threadno = loadrem[i][0]; 
                    gapallocate[threadno][loadperthread[threadno]++] = largearr[totalalloc++];
                    sumalloc++;
                    threadstat[threadno][1]++;
                    if(totalalloc == largegap)
                    {
                        flag =1;
                        break;
                    }
                }
                if(flag==1)break;
            }
        }
    }
    
    /*cout<<"=========Stats=============="<<endl;
    for(int i=0;i< num_threads;i++)
    {
        int sum=0;
        cout<<"Thread - "<<i<<"\t";
        for(int j=0;j< 2;j++)
        {
            //cout<<threadstat[i][j]<<"\t";
            sum += threadstat[i][j];
        }
        cout<<"Total = "<<sum<<endl;
    }*/
    
    writeGapLoad(num_threads,loadperthread);
    
    //cout<<"Total allocated gaps = "<<sumalloc<<endl;
    
    if(1)
    {
        for(int i=0;i<num_threads;i++)
        {
            t[i] = thread(thread_fillgap,i,read_length,loadperthread[i]);
            
            unsigned int microsecond = 1000000;
            usleep(1 * microsecond);
        }
    
        for(int i=0;i<num_threads;i++)t[i].join();
    }
    
    
    ///=========================================================================
    
    //cout<<"All joined\n";


    ///=================Join Threaded output files=================
    char top[1000];

    strcpy(top,tp);
    strcat(top,"draw");
    mergeFiles(top,num_threads);

    strcpy(top,tp);
    strcat(top,"gapout");
    readgapouts(top,num_threads,tot_gaps);

    ///===============================Construct gap filled scaffold file================================

    ///collect gaptofill

    for(int i=0;i<num_threads;i++)
    {
        readGtf(i);
    }
    //cout<<"G2f count = "<<g2fcount<<endl;
    
    FILE* contigFile=fopen(contigFileName, "r");

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


    int nStart=0;
    long int nStartPos=0;
    int nCount=0,afterCount=0;

    long int offset=0, index=0;

    strcpy(dir,tp);
    strcat(dir,"gapout.txt");

    FILE* gapOutFile=fopen(dir,"r");

    strcpy(dir,tp);
    strcat(dir,"filledContigs.fa");

    FILE* filledContigFile=fopen(dir,"w");

    //cout<<"# of Contig = "<<noContigs<<endl;
    int gap_count=-1;
    
    strcpy(dir,tp);
    strcat(dir,"Ncount.txt");

    FILE* Nfile = fopen(dir,"w");
    
    int prevNcount=0,newNcount=0;

    //FILE* gapstat = fopen("gapstat.txt","w");
    //fprintf(gapstat,"Gap\t\tResult\t\tTotal_filled_char\tLeft side\tRight Side\t\tTotal Reads\t\tUsed Reads\t\tTime taken(s)\n\n");

    for(int i=0;i<noContigs;i++)
    {
        //cout<<"Contig len = "<<contigLengths[i]<<endl;
    	contig=new char[contigLengths[i]+1];
    	//cout<<"Allocation Done"<<endl;
    	offset=0;
    	index=0;

    	fprintf(filledContigFile,">%s\n",contigNames[i]);

        for(long int j=0;j<contigLengths[i];j++)
        {
            if(contigs[i][j]=='N' || contigs[i][j]=='n')
            {
                //cout<<"found N"<<endl;
                if(nStart==0)
                {
                    nStart=1;
                    nCount=1;
                    nStartPos=j;
                    gap_count++;
                }
                else
                {
                    nCount++;
                }
            }
            else
            {
            	if(nStart==1)
                {
                	//cout<<"Gap found in contig - "<<i<<endl;
                	if(nCount>=1)
                	{
                   		prevNcount += nCount;
                   		fscanf(gapOutFile,"%d\t%d\t%ld\t%d\t%d\t",&gapNo,&gapContigNo,&gapStart,&gapLength,&gapStringLength);
                   		//cout<<gapNo<<" "<<gapContigNo<<" "<<gapStringLength<<" "<<gapLength<<endl;

                   		if(gapStringLength >0)
                   		{
                   		    fscanf(gapOutFile,"%s\n",gapString);
                   		}
                   		else
                   		{
                   		    fscanf(gapOutFile,"\n");
                   		}
                   		
                   		newNcount += findNcount(gapString);
                   		//==================================

                    	contig[index]='\0';
						fprintf(filledContigFile,"%s",contig);
                    	if(gapStringLength >0)fprintf(filledContigFile,"%s",gapString);

                    	offset=offset+(gapStringLength-gapLength);
                    	//cout<<"Offset = "<<offset<<endl;
                    	contig=(char *)realloc(contig, (contigLengths[i]+1+offset-index)*sizeof(char));
                    	//cout<<"Realloc size = "<<(contigLengths[i]+1+offset-index)*sizeof(char)<<endl;
                    	index=0;
                    	//gapString="";
                	}
                	else
                	{

                		contig[index]='\0';
                		fprintf(filledContigFile,"%s\n",contig);
                		index=0;
                		for(int k=0;k<nCount;k++)
                		{
                			contig[k]='N';
                		}
                		contig[nCount]='\0';
                		fprintf(filledContigFile,"%s\n",contig);
                	}
                	nStart=0;
                }
                if(gap_count >= 0 && gaptofill[gap_count] > 0)
                {
                    //cout<<"Skipping character "<<contigs[i][j]<<", for gap "<<(gap_count-1)<<endl;
                    gaptofill[gap_count]--;
                }
                else
                {
                    contig[index]=contigs[i][j];
                    index++;
                }
            }
        }
        contig[index]='\0';
    	fprintf(filledContigFile,"%s\n",contig);
        //cout<<"Done with Contig"<<i<< "successfully"<<endl;
    }
    
    int diff = prevNcount - newNcount;
    int fillmore=-1;
    
    if(diff == 0 || newNcount == 0)fillmore = 0;
    if(newNcount == 0)fillmore = 0;
    else fillmore =1;
    
    //cout<<prevNcount<<"\t"<<newNcount<<"\t"<<fillmore<<endl;
        
    fprintf(Nfile,"%d",fillmore);

	fclose(gapOutFile);
	fclose(filledContigFile);
	fclose(Nfile);
	//fclose(gapstat);
	
	delete[] gaptofill;
	
	time(&endt);
    
    double seconds = difftime(endt,nowt);
    cout<<"Time taken = "<<seconds<<" seconds"<<endl;

    cout<<"======================================"<<endl;
    cout<<"Iteration "<<script_itr<<" ends successfully"<<endl;
    cout<<"======================================"<<endl;
    //cout<<"Effective length for mean insertsize = "<<getEffectiveLength(insertSizeMean)<<endl;


	return 0;
}
