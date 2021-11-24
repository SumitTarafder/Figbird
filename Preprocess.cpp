//
//  main.cpp
//  gap
//
//  Created by Atif Rahman on 2/4/16.
//  Copyright (c) 2016 Atif Rahman. All rights reserved.
//
#include <bits/stdc++.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <iomanip>
#include <unistd.h>
#include <limits>
#include <vector>
#include <cstring>

using namespace std;
#include <math.h>
#include <map>

#define MAX_REC_LEN 1024
#define MAX_READLENGTH 200
#define MAX_READ_PER_GAP 5000
#define MAX_NAMELENGTH 200
#define HASH_TABLE_SIZE 10001

//===================Defined Structs========================

struct SAM
{
	char qname[MAX_NAMELENGTH];
	int flag;
	char rname[MAX_NAMELENGTH];
	int pos;
	int mapq;
	char cigar[MAX_READLENGTH];
	char rnext[MAX_NAMELENGTH];
	int pnext;
	int tlen;
	char seq[MAX_NAMELENGTH];
	char qual[MAX_NAMELENGTH];
	char md[MAX_NAMELENGTH];
	unsigned long ih;
	int nm,as;
	long int contigNo;
};

struct Gap
{
    int contigNo;
    long int gapStart;
    int gapLength;
    int read_count;
    int partial_read_count;
    char **unmapped_jump_reads;
    char **partial_frag_reads;
    int total_allocated_rows;
};


struct Contig
{
	char contigName[1000];
	long int contigNo;
};

//==================================================

//==================File variables==============
FILE *outFile;
FILE *linkingFile;
FILE *singletonFile;
FILE *mapFile;
FILE *longmapFile;
FILE *contigFile;
FILE *unFile;
FILE *statFile;
FILE *statFile2;
FILE *gapInfoFile;
FILE* output1,*output2;
//=================Int variables================
long int noGaps=0;
int **readCounts;
long int unCount=0;
long int totalCount=0;
unsigned long maxReadLength=0;
long int MAX_FRAGMENT_SIZE=5000;
long int noContigs;
int maxDistance;
int r_c_p=0,rcp=0;
int case_count[4]={0};
int check_pos_count=0;
int partial_check_pos_count=0;
int small_gap_count=0,large_gap_count=0,huge_gap=0;
int gapFileNo=-1;
int partial_gapFileNo=-1;
int pos_strand=0,neg_strand=0;
int *read_in_gap_f;
int *read_in_gap_r;
int samflag=-1;
int read_mean;

int cigar_val[3]={0,0,0};
int script_itr=-1;
int partial_case_count[2]={ 0 };
int * gaptofill;
int * perfectread_gap;
int * perfectread_gaplen;
char *tempwrite,*cig1,*cig2;
int* gappedcontigs;
int maxInsertSize=0;
long int *insertCounts;
long int discardedReads=0;//zxcf
//=================Double vars===================

double *contigReadCounts;
double insertSizeMean=0;
double insertSizeVar=0;
double squaredError=0;

//==================All vectors==================
vector <SAM *> reads1;
vector <SAM *> reads2;

vector <SAM *> mixedReads1;
vector <SAM *> mixedReads2;

vector<char*> contigs;
vector<char*> contigNames;
vector<unsigned long> contigLengths;

vector<Contig *> hashTable[HASH_TABLE_SIZE];
vector<Gap *> gaps;

//================Char arrays====================

char line1[MAX_REC_LEN];
char line2[MAX_REC_LEN];

char **gapFiles;
char **partial_gapFiles;


//=================Declaration Done================================
void reverse(char *reverse, char *read)
{
	char ch='A';
	unsigned long readLength=strlen(read);
	//cout<<readLength<<endl;
	for(unsigned long i=1; i<=readLength;i++)
	{
		ch=read[i-1];
		if(ch=='A')
			reverse[readLength-i]='T';
		else if(ch=='C')
			reverse[readLength-i]='G';
		else if(ch=='G')
			reverse[readLength-i]='C';
		else if(ch=='T')
			reverse[readLength-i]='A';
		else
			reverse[readLength-i]='N';
	}
	reverse[readLength]='\0';
}

int parseDel(char* cigar)
{
    //For reads that match with right side, cigar will be like this 4S97M
    //It will return del= 4, amount of soft clipped part
    //otherwise it will return 0 for all other case of reads
    //cout<<"parseDel"<<endl;
    int del=0;
    char*p1,*p2;

    p1 = strpbrk(cigar,"S");
    p2 = strpbrk(cigar,"M");

    std::string s="";

    int index1=0,index2=0;

    //cout<<p1<<" "<<p2<<endl;

    if(p1 && p2)
    {
        index1=(int)(p1-cigar);
        index2=(int)(p2-cigar);

        if(index1 < index2)
        {
            for(int k=0;k<index1;k++)
            {
                s+= cigar[k];
            }
            const char * c = s.c_str();
            del = atoi(c);
        }
    }
    return del;
}

void parse_Cigar(char* cigar,int readlen)
{
    //For cases like 101M or 70M31S, it will do nothing, 3 cigar values will be 0
    //For case 1 and 4- if S comes before M- 10S80M11S, then this code works
    //For case 2 and 3- S comes naturally before M- 10S91M, then this code works although we don't call it from there!
    //cigar[0] = S cigar[1] = M cigar[2] = S

    //zx
    char*p1,*p2;

    p1 = strpbrk(cigar,"S");
    p2 = strpbrk(cigar,"M");

    std::string s="";

    int index1=0,index2=0;

    if(p1 && p2)
    {
        index1=(int)(p1-cigar);
        index2=(int)(p2-cigar);

        if(index1 < index2)
        {
            for(int k=0;k<index1;k++)
            {
                s+= cigar[k];
            }
            const char * c = s.c_str();
            cigar_val[0] = atoi(c);

            s="";
            for(int k=index1+1;k<index2;k++)
            {
                s+= cigar[k];
            }
            c = s.c_str();
            cigar_val[1] = atoi(c);
            s="";
            if((cigar_val[0] + cigar_val[1]) != readlen)
            {
                for(int k=index2+1;k<strlen(cigar);k++)
                {
                    s+= cigar[k];
                }

                int n = s.length();
                char new_s[n + 1];

                strcpy(new_s, s.c_str());

                int s_index[5];
                std::string s1="S";
                int s_count=0;

                for (int i = 0; i < s.length(); i++) {
                    if (s.substr(i, s1.length()) == s1) {
                        s_index[s_count++] = i;
                    }
                }

               //if(s_index[s_count-1] == (s.length()-1))
               if(s_count>0)
               {
                    if(strpbrk(new_s,"D") || strpbrk(new_s,"I")|| strpbrk(new_s,"M")||strpbrk(new_s,"X")||strpbrk(new_s,"="))
                    {
                        s="";
                        for(int k = s_index[s_count-1]-2;k<s_index[s_count-1];k++)
                        {
                            s+= new_s[k];
                        }

                        strcpy(new_s, s.c_str());

                        if(strpbrk(new_s,"D") || strpbrk(new_s,"I")|| strpbrk(new_s,"M")||strpbrk(new_s,"X")||strpbrk(new_s,"="))
                        {
                            s="";
                            s += new_s[s_index[s_count-1]-1];
                        }
                        c = s.c_str();
                        cigar_val[2] = atoi(c);

                    }
                    else
                    {
                        cigar_val[2] = readlen - cigar_val[0] - cigar_val[1];
                    }
                }
            }
        }
    }
}


unsigned long getHash(char *str)
{
    unsigned long hash = 5381;
    int c;
    
    while (c = *str++)
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
    
    return hash;
}

void initHashTable()
{
	unsigned long index;
	for(long int i=0;i<noContigs;i++)
	{
		index=getHash(contigNames[i]) % HASH_TABLE_SIZE;
		Contig *c=new Contig();
		strcpy(c->contigName,contigNames[i]);
		c->contigNo=i;
		hashTable[index].push_back(c);
	}
}

int getContigNo(char *contigName)
{
	unsigned long index=getHash(contigName) % HASH_TABLE_SIZE;
    
	for(int i=0;i<hashTable[index].size();i++)
	{
		if(strcmp(hashTable[index][i]->contigName,contigName)==0)
			return hashTable[index][i]->contigNo;
	}
	return -1;
    
}

int qualityFilter(char* qual_str)
{
    int flag=0,Q;

    if(1)return flag;

    double expected_num_errors = 0;
    int Q_count=0;

    int Q_threshold = 6,Q_c_thresh=10;
    double temp,E_MAX = 7,val;//0.5,.25
    int read_len=strlen(qual_str);

    for(int i=0;i<read_len;i++)
    {
        Q = qual_str[i]-33;//ASCII-Base-33 for Illumina
        val = pow(10,-Q/10.0);
        expected_num_errors += val;
        if(Q<Q_threshold)Q_count++;
    }
    //cout<<expected_num_errors<<endl;

    if(expected_num_errors > E_MAX)flag=1;
    //if(Q_count > Q_c_thresh)flag=1;

    return flag;//YOU MUST remove the return above
}

int check_duplicate(char* read,int gapno)
{
    int flag = 0;

    std::string s2(read);
    std::string s3;

    if(samflag == 2)
    {
        for(int i=0;i<gaps[gapno]->read_count;i++)
        {
            if(strcmp(read,gaps[gapno]->unmapped_jump_reads[i])==0)
            {
                return 1;
            }

            std::string s1(gaps[gapno]->unmapped_jump_reads[i]);
            s3 = s2.substr (2,s2.size()-4);//clip 2 from begining and 2 from ending
            //cout<<s2<<endl;

            if (s1.find(s3) != std::string::npos)
            {
                return 1;
            }
        }
        return 0;
    }
    else
    {
        for(int i=0;i<gaps[gapno]->partial_read_count;i++)
        {
            if(strcmp(read,gaps[gapno]->partial_frag_reads[i])==0)
            {
                //cout<<"Found duplicate.. read = "<<i<<" "<<gapno<<" "<<gaps[gapno]->partial_read_count<<endl;
                return 1;
            }
        }
        return 0;
    }
     return 0;
}

void writeSam(SAM* read, int gapFileNo)
{
  FILE* out = fopen(gapFiles[gapFileNo],"a");
	fprintf(out,"%s\t%d\t%ld\t%d\t%s\t%d\t%s\t%s\t%s\tIH:i:%ld\n",read->qname,read->flag,read->contigNo,read->pos,
            read->cigar,read->tlen,read->seq,read->qual,read->md,read->ih);
  fclose(out);  
}

void writeSam2(SAM* read, FILE* out)
{
  fprintf(out,"%s\t%d\t%ld\t%d\t%s\t%d\t%s\t%s\t%s\tIH:i:%ld\n",read->qname,read->flag,read->contigNo,read->pos,
          read->cigar,read->tlen,read->seq,read->qual,read->md,read->ih);
}

void printCigar(int gapNo,char* cigar,int match)
{
    //cout<<"Gap = "<<gapNo<<"\tCigar = "<<cigar<<"\tmatch = "<<match<<endl;

}


int writePartialSam(SAM* read, int gapNo,int strandNo,int del,int pos2)
{
    FILE *out = fopen(partial_gapFiles[gapNo],"a");
    int clipped_index = 0;
    int match=-1;
    int gap_s = gaps[gapNo]->gapStart;
    int gap_e = gaps[gapNo]->gapStart + gaps[gapNo]->gapLength;
    int contiglen = contigLengths[gaps[gapNo]->contigNo];
    int readlength = strlen(read->seq);

    int w_flag=0;
    char tempReadString[500];

    if(read->pos < gap_s)
    {
        if(strandNo == 0)match = 1;
        else match = 4;

        //cout<<"Going to parse_cigar for "<<read->cigar<<" match = "<<match<<" pos = "<<read->pos<<" gaps = "<<gap_s<<endl;
        //cout<<read->seq<<endl;

        for(int i=0;i<3;i++)cigar_val[i] = 0;
        parse_Cigar(read->cigar,readlength);

        if(cigar_val[0])//In case S comes before M,then extra check
        {
            if(cigar_val[2])//S--M--S
            {
                clipped_index = readlength - cigar_val[2] - 1;
                fprintf(out,"%s\t%d\t%d\t%d\t%s\t%d\t%s\n",read->seq,clipped_index,match,read->pos,read->cigar,pos2,read->qual);
                w_flag=1;
                //printCigar(gapNo,read->cigar,match);
            }
            else//Only S--M, so discard them
            {

            }
        }
        else//M comes before S or no s at all like 101M,so write it//MS or M or MIM
        {
            clipped_index = gap_s -read->pos;
            fprintf(out,"%s\t%d\t%d\t%d\t%s\t%d\t%s\n",read->seq,clipped_index,match,read->pos,read->cigar,pos2,read->qual);
             w_flag=1;
             //printCigar(gapNo,read->cigar,match);
        }
    }
    else if(read->pos > gap_s)
    {
        parse_Cigar(read->cigar,readlength);
        if(strandNo == 0)match = 2;
        else match = 3;

        clipped_index = gap_e - 1 -read->pos + del + 2;//+2 because of 1-based coordinate system in sam reporting
        fprintf(out,"%s\t%d\t%d\t%d\t%s\t%d\t%s\n",read->seq,clipped_index,match,read->pos,read->cigar,pos2,read->qual);
        //cout<<read->cigar<<endl;
        //if(cigar_val[2])printCigar(gapNo,read->cigar,match);
        w_flag=1;
    }
    
    fclose(out);

    if(w_flag)
    {
        partial_check_pos_count++;
        //cout<<"Read1 = / "<<"Gap = "<<partial_gapFileNo<<" "<<read1->cigar<<" "<<read1->pos<<" "<<strandNo1<<" "<<del<<endl;
        if(strandNo==0)
        {
            read_in_gap_f[partial_gapFileNo]++;
            pos_strand++;
        }
        else
        {
            read_in_gap_r[partial_gapFileNo]++;
            neg_strand++;
        }
    }
    return match;
}

void printSam(SAM* read)
{
	//printf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tIH:i:%ld\tAS:i:%d\n",read->qname,read->flag,read->rname,read->pos,
      //     read->mapq,read->cigar,read->rnext,read->pnext,read->tlen,read->seq,read->qual,read->md,read->ih,read->as);
    
    
    /*	printf("%s\t%d\t%d\t%s\t%d\t%s\t\%s\tIH:i:%d\n",read->qname,read->flag,read->pos,
     read->cigar,read->tlen,read->seq,read->md,read->ih);
     */
}

int checkInsert(int a,int b,double mean,int gap)
{
    double diff1,diff2;
    diff1 = abs(mean-a);
    diff2 = abs(mean-b);

    if(diff1 < diff2)return a;
    else return b;
}

int checkRange(int a, int b, int mean)
{
    int insertsizemin = mean - 1000,insertsizemax=mean + 1000;

    if(a > insertsizemin && a < insertsizemax)return 1;
    if(b > insertsizemin && b < insertsizemax)return 1;
    if((a < insertsizemin && b > insertsizemax) || (b < insertsizemin && a > insertsizemax))return 1;

    return 0;
}

int checkPos(int contigNo, long int pos, int strandNo,int readlength)
{
    int num_gaps= gaps.size();

    int tempinsertsize[num_gaps],tempinsertsize_abs[num_gaps];
    int val[2];

    int map = -1;

    int flag = 0,gap_index=0;
    int min_val=1000000,min_index=-1;

    for(int i=0;i<num_gaps;i++)
    {
        tempinsertsize[i] = 0;

        if(contigNo==gaps[i]->contigNo && ((strandNo==0 && pos>(gaps[i]->gapStart-maxDistance) &&
        pos<(gaps[i]->gapStart)) || (strandNo==1 &&
        pos>(gaps[i]->gapStart+gaps[i]->gapLength) && pos<(gaps[i]->gapStart+gaps[i]->gapLength+maxDistance))))
        {
            if(maxDistance<= 250)return i;//Should be commented out if you want to execute following code below
            if(pos < gaps[i]->gapStart)
            {
                //Calculated assuming the other end is placed at the end of the gap
                //pos----------------->
                //======================================NNNNNNNNNNNNNNNNNNNNN=======================
                //                                                          ---------------->
                val[0] = (gaps[i]->gapStart+gaps[i]->gapLength - pos + readlength);//based on gap end
                val[1] = (gaps[i]->gapStart - pos + 1);//based on start
                map=0;

                if(checkRange(val[0],val[1],read_mean))
                    tempinsertsize[i] = checkInsert(val[0],val[1],read_mean,i);//returns closer to mean
            }
            else
            {
                //Calculated assuming the other end is placed at gapstart - readlen + 1
                //                                                                      pos----------------->
                //======================================NNNNNNNNNNNNNNNNNNNNN==================================
                //                     ----------------->
                val[0] = (pos - gaps[i]->gapStart + 2*readlength - 1);//based on gap start
                val[1] = (pos - gaps[i]->gapStart - gaps[i]->gapLength + readlength + 1);//based on gap end
                map=1;

                if(checkRange(val[0],val[1],read_mean))
                    tempinsertsize[i] = checkInsert(val[0],val[1],read_mean,i);//returns closer to mean
            }

            if(tempinsertsize[i] !=0)
            {
                flag++;
                gap_index=i;
            }

            tempinsertsize_abs[i] = abs((read_mean) - tempinsertsize[i]);

            if(tempinsertsize_abs[i] < min_val)
            {
                min_val = tempinsertsize_abs[i];
                min_index = i;
            }
        }
        if(gaps[i]->contigNo > contigNo)break;
    }

    if(maxDistance <= 250)return -1;//Should be commented out if you want to execute following code below

    if(flag == 0)return -1;

    int min_thresh = read_mean - read_mean* 0.6;

    int c_index;

    if(flag == 1)c_index=gap_index;
    else c_index=min_index;

    if(tempinsertsize[c_index] < min_thresh)return -1;
    else return c_index;
}

int checkPos2(int contigNo, long int pos, int strandNo,int readlength, int del)
{
    for(int i=0;i<gaps.size();i++)
    {
        int gapEnd = gaps[i]->gapStart+gaps[i]->gapLength;//143

        if(contigNo==gaps[i]->contigNo &&
        ((strandNo==0 &&
        (
            (pos>(gaps[i]->gapStart-readlength+1) && pos <= gaps[i]->gapStart) ||
            (pos > gapEnd && del && (pos - del) <= gapEnd))
        )||
        (strandNo==1 &&
        (
            (pos>(gaps[i]->gapStart-readlength+1) && pos <= gaps[i]->gapStart) ||
            (pos > gapEnd && del && (pos - del) <= gapEnd))
        )))
        {
            return i;
        }
        //if(gaps[i]->contigNo > contigNo)break;
    }
    return -1;
}

void printVectors(FILE * out)
{
	SAM *read1, *read2;
	unsigned long readLength1,readLength2;
	int pos1, pos2;
    //	int insertSize;
    
	char temp[MAX_READLENGTH];
    
	unsigned long ih=reads1.size();
    
	int ih1_0=0;
	int ih2_0=0;
    
    
	for(int i=0;i<ih;i++)
	{
		if(strcmp(reads1[i]->rname,"*")!=0 && reads1[i]->nm==0)
		{
			ih1_0++;

		}

	}
    
	for(int i=0;i<ih;i++)
	{
		if(strcmp(reads2[i]->rname,"*")!=0  && reads2[i]->nm==0)
		{
			ih2_0++;

		}

	}
	if(ih1_0>0)
	{
        for(int i=0;i<ih;i++)
        {
            if(strcmp(reads1[i]->rname,"*")!=0  && reads1[i]->nm==0)
            {
                //contigReadCounts[reads1[i]->contigNo]+=1/(double)ih1_0;
            }
        }
	}
    
	if(ih2_0>0)
	{
        for(int i=0;i<ih;i++)
        {
            if(strcmp(reads2[i]->rname,"*")!=0  && reads2[i]->nm==0)
            {
                //contigReadCounts[reads2[i]->contigNo]+=1/(double)ih2_0;
            }
        }
	}

	for(unsigned long i=0;i<ih;i++)
	{
		read1=reads1[i];
		read2=reads2[i];

		if(strcmp(read1->rname,"*")==0 || strcmp(read2->rname,"*")==0)
		{
		    //Never here
		    //		int i=0;
            //		int j=0;
			int ncount1=0;
			int ncount2=0;
            
			readLength1=strlen(read1->seq);
			if(readLength1>maxReadLength)
				maxReadLength=readLength1;
            
			readLength2=strlen(read2->seq);
			if(readLength2>maxReadLength)
				maxReadLength=readLength2;
            
            
			for(int i=0;i<readLength1;i++)
			{
				if(read1->seq[i]=='N'||read1->seq[i]=='n')
				{
					ncount1++;
				}
			}
            
			for(int i=0;i<readLength2;i++)
			{
                
				if(read2->seq[i]=='N'||read2->seq[i]=='n')
				{
					ncount2++;
				}
			}
            
			if(ncount1/(double)readLength1 < 0.8 && ncount2/(double)readLength2 < 0.8)
			{
			    //cout<<"PrintVectors*************"<<ncount1<<" $$$$"<<ncount2<<endl;
				unCount++;
				totalCount++;
                
                /*
				fputs("@",u);
				fputs(read1->qname,u);
				fputs("\n",u);
                
                
				int strandNo=(read1->flag&16)>>4;
				if(strandNo==1)
				{
					reverse(temp,read1->seq);
					fputs(temp,u);
				}
				else
				{
					fputs(read1->seq,u);
				}
				fputs("\n",u);
                
				fputs("+",u);
				fputs(read1->qname,u);
				fputs("\n",u);
                
				fputs(read1->qual,u);
				fputs("\n",u);
                
				fputs("@",u);
				fputs(read2->qname,u);
				fputs("\n",u);
                
				strandNo=(read2->flag&16)>>4;
				if(strandNo==1)
				{
					reverse(temp,read2->seq);
					fputs(temp,u);
				}
				else
				{
					fputs(read2->seq,u);
                    
				}
				fputs("\n",u);
                
				fputs("+",u);
				fputs(read2->qname,u);
				fputs("\n",u);
                
				fputs(read2->qual,u);
				fputs("\n",u);
				*/
				
				for(int i=0;i<reads1.size();i++)
					delete reads1[i];
                
				for(int i=0;i<reads2.size();i++)
					delete reads2[i];
                
                
				reads1.clear();
				reads2.clear();
				
				return;
			}
            
		}
		else if(strcmp(read1->rname,read2->rname)!=0)
		{
            //Never here
			//readCounts[read1->contigNo][read2->contigNo]++;
			//printSam(read1);
			//printSam(read2);
			//cout<<"Warning: contig names different"<<endl;
            
            //		getchar();
            
		}
		else
		{
			//readCounts[read1->contigNo][read2->contigNo]++;
            
			pos1=read1->pos;
			readLength1=strlen(read1->seq);
			if(readLength1>maxReadLength)
				maxReadLength=readLength1;
			read1->ih=ih;
            
            
			pos2=read2->pos;
			readLength2=strlen(read2->seq);
			if(readLength2>maxReadLength)
				maxReadLength=readLength2;
			read2->ih=ih;
            
			writeSam2(read1,out);
			writeSam2(read2,out);
			r_c_p++;

			//if(strcmp(read1->cigar,cig1)==0 && strcmp(read2->cigar,cig2)==0)rcp++;

            //cout<<"In else--->>>>>>>>>>>>>>>>>>>. \n";
		}
	}
	totalCount++;

	for(int i=0;i<reads1.size();i++)
		delete reads1[i];
    
	for(int i=0;i<reads2.size();i++)
		delete reads2[i];
	
    
	reads1.clear();
	reads2.clear();
    
}

int check_Ncount_partial(char* read)
{
    int count=0;
    for(int i=0;i<strlen(read);i++)
    {
        if(read[i] == 'N')count++;
    }
    if(count<=3)return 1;
    else return 0;
}

int checkChar(char* read)
{
    int flag=0;

    for(int i=0;i<strlen(read);i++)
    {
        if(!(read[i] == 'A' || read[i] == 'C' || read[i] == 'G'|| read[i] == 'T' || read[i] == 'N'
            || read[i] == 'a' || read[i] == 'c' || read[i] == 'g'|| read[i] == 't' || read[i] == 'n'))
        {
           flag=1;
           break;
        }
    }

    return flag;
}

void checkMIM(char* cigar,int g)
{

    int index1=0,index2=0,index3=0,m_count=0,i_count=0,i_len=0;
    std::string s="";

    for(int i=0;i<strlen(cigar);i++)
    {
        if(cigar[i] == 'S' || cigar[i] == 'D' || cigar[i] == '=' || cigar[i] == 'X')return;
        if(cigar[i] == 'M')
        {
            if(m_count == 0)index1 = i;
            else if(m_count==1)index3=i;
            else return;
            m_count++;
        }
        else if(cigar[i] == 'I')
        {
            if(i_count==1)return;
            index2 = i;
            i_count++;
        }
    }

    if(index1 && index2 && index3 && index1 < index2 && index2 < index3)
    {
        //cout<<index1 <<"\t"<<index2<<"\t"<<index3<<endl;
        //cout<<"MIM found, cigar = "<<cigar<<"\tGap = "<<g<<endl;

        for(int k=index1+1;k<index2;k++)
        {
            s+= cigar[k];

        }
        i_len = atoi(s.c_str());

        perfectread_gap[g] = 1;
        perfectread_gaplen[g] = i_len+1;//gaps[g]->gapLength;
        //cout<<"MIM found, cigar = "<<cigar<<"\tGap = "<<g<<"\tgaplen = "<<perfectread_gaplen[g]<<endl;
    }
}

static void free_numbers(char **array, size_t size)
{
    for (size_t i = 0; i < size; i++)
        free(array[i]);
    free(array);
}

void allocateJ(int gapFileNo)
{
    int offset=10;
    //cout<<"In allocate J\n";
    if (gaps[gapFileNo]->read_count % offset == 0)
    {
        char **newptr = (char **)realloc(gaps[gapFileNo]->unmapped_jump_reads, (gaps[gapFileNo]->read_count + offset) * sizeof(*gaps[gapFileNo]->unmapped_jump_reads));
        if (newptr == NULL)
        {
            free_numbers(gaps[gapFileNo]->unmapped_jump_reads, gaps[gapFileNo]->total_allocated_rows);
            cerr<<"Allocation problem in main, exiting...\n";
            exit(1);
        }

        gaps[gapFileNo]->unmapped_jump_reads = newptr;
        gaps[gapFileNo]->total_allocated_rows= gaps[gapFileNo]->read_count;

        for(int i=gaps[gapFileNo]->read_count;i<gaps[gapFileNo]->read_count + offset;i++)
        {
            gaps[gapFileNo]->unmapped_jump_reads[i] = (char *)malloc(MAX_READLENGTH * sizeof(char));

            if (gaps[gapFileNo]->unmapped_jump_reads[i] == 0)
            {
                free_numbers(gaps[gapFileNo]->unmapped_jump_reads, gaps[gapFileNo]->total_allocated_rows);
                cerr<<"Allocation problem in main, exiting...\n";
                exit(1);
            }
            gaps[gapFileNo]->total_allocated_rows++;
        }
    }
}

void allocateP(int gapFileNo)
{
    int offset = 10;
    //cout<<"In allocate P\n";

    if(gaps[gapFileNo]->partial_read_count % offset == 0)
    {
        char **newptr = (char **)realloc(gaps[gapFileNo]->partial_frag_reads, (gaps[gapFileNo]->partial_read_count + offset)* sizeof(*gaps[gapFileNo]->partial_frag_reads));
        if (newptr == NULL)
        {
            free_numbers(gaps[gapFileNo]->partial_frag_reads, gaps[gapFileNo]->partial_read_count);
            cerr<<"Allocation problem in main, exiting...\n";
            exit(1);
        }

        gaps[gapFileNo]->partial_frag_reads = newptr;
        gaps[gapFileNo]->total_allocated_rows=gaps[gapFileNo]->partial_read_count;

        for(int i=gaps[gapFileNo]->partial_read_count;i<gaps[gapFileNo]->partial_read_count + offset;i++)
        {
            gaps[gapFileNo]->partial_frag_reads[i] = (char *)malloc(MAX_READLENGTH * sizeof(char));

            if(gaps[gapFileNo]->partial_frag_reads[i] == 0)
            {
                free_numbers(gaps[gapFileNo]->partial_frag_reads, gaps[gapFileNo]->total_allocated_rows);
                cerr<<"Allocation problem in main, exiting...\n";
                exit(1);
            }
            gaps[gapFileNo]->total_allocated_rows++;
        }
    }
}

void printMixedVectors()
{
	SAM *read1, *read2;
	unsigned long readLength1,readLength2;
    long int pos1, pos2;
    int contigNo1, contigNo2;
    int strandNo1, strandNo2;

    char temp[MAX_READLENGTH];
    
	int ih1=mixedReads1.size();
    int ih2=mixedReads2.size();
    
	int ih1_0=0;
	int ih2_0=0;
    
    
	for(int i=0;i<ih1;i++)
	{
		if(strcmp(mixedReads1[i]->rname,"*")!=0 && mixedReads1[i]->nm==0)
		{
			ih1_0++;
		}
	}
    
	for(int i=0;i<ih2;i++)
	{
		if(strcmp(mixedReads2[i]->rname,"*")!=0  && mixedReads2[i]->nm==0)
		{
			ih2_0++;
		}
	}
	
	if(ih1_0>0)
	{
        for(int i=0;i<ih1;i++)
        {
            if(strcmp(mixedReads1[i]->rname,"*")!=0 && mixedReads1[i]->nm==0)
            {
                //contigReadCounts[mixedReads1[i]->contigNo]+=1/(double)ih1_0;
            }
        }
	}
	
	if(ih2_0>0)
	{
        for(int i=0;i<ih2;i++)
        {
            if(strcmp(mixedReads2[i]->rname,"*")!=0 && mixedReads2[i]->nm==0)
            {
                //contigReadCounts[mixedReads2[i]->contigNo]+=1/(double)ih2_0;
            }
        }
	}
	
	if(ih1>1 || ih2>1)
    {
        for(int i=0;i<mixedReads1.size();i++)
        {
            read1=mixedReads1[i];
        }
        for(int j=0;j<mixedReads2.size();j++)
        {
            read1=mixedReads2[j];
        }
    }

	for(int i=0;i<mixedReads1.size();i++)
	{
		read1=mixedReads1[i];
        
		for(int j=0;j<mixedReads2.size();j++)
		{
			read2=mixedReads2[j];
            if(i==0 && j==0)//
            {
                readLength1=strlen(read1->seq);
                if(readLength1>maxReadLength)
                    maxReadLength=readLength1;
                
                readLength2=strlen(read2->seq);
                if(readLength2>maxReadLength)
                    maxReadLength=readLength2;
                
                
                int ncount1=0;
                int ncount2=0;
                
                
                for(int i=0;i<readLength1;i++)
                {
                    if(read1->seq[i]=='N'||read1->seq[i]=='n')
                    {
                        ncount1++;
                    }
                }
                
                for(int i=0;i<readLength2;i++)
                {
                    
                    if(read2->seq[i]=='N'||read2->seq[i]=='n')
                    {
                        ncount2++;
                    }
                }

                //IF both the mates in the read pair contains more than 80% Ns, we put them in unmapped

                if(ncount1/(double)readLength1 < 0.8 && ncount2/(double)readLength2 < 0.8)
                {
                    unCount++;
                    //out<<"Yes, this pair is unmapped\n";
                    totalCount++;

                    /*fputs("@",u);
                    fputs(read1->qname,u);
                    fputs("\n",u);

                    int strandNo=(read1->flag&16)>>4;
                    if(strandNo==1)// If the seq is reverse complemented, then reverse it again
                    {
                        reverse(temp,read1->seq);
                        fputs(temp,u);
                    }
                    else
                    {
                        fputs(read1->seq,u);
                    }
                    fputs("\n",u);
                    
                    fputs("+",u);
                    fputs(read1->qname,u);
                    fputs("\n",u);
                    
                    fputs(read1->qual,u);
                    fputs("\n",u);
                    
                    fputs("@",u);
                    fputs(read2->qname,u);
                    fputs("\n",u);
                    
                    strandNo=(read2->flag&16)>>4;
                    if(strandNo==1)
                    {
                        reverse(temp,read2->seq);
                        fputs(temp,u);
                    }
                    else
                    {
                        fputs(read2->seq,u);
                    }
                    fputs("\n",u);
                    
                    fputs("+",u);
                    fputs(read2->qname,u);
                    fputs("\n",u);
                    
                    fputs(read2->qual,u);
                    fputs("\n",u);

                    */
                }
                else
                {
                    //Otherwise they are cleared from the mixed reads

                    for(int i=0;i<mixedReads1.size();i++)
                        delete mixedReads1[i];
                    
                    for(int i=0;i<mixedReads2.size();i++)
                        delete mixedReads2[i];
                    
                    
                    mixedReads1.clear();
                    mixedReads2.clear();
                    
                    return;
                    
                }
            }

            
            //=====================For rest of the i and j combinations including 0 and 0==========================================
            //4 cases
            //zxc
            //Case1: If both are unmapped,just clear,Do nothing
            if(((read1->flag) & 4) != 0 && ((read2->flag) & 4) != 0)
			{
			    case_count[0]++;
                for(int i=0;i<mixedReads1.size();i++)
                    delete mixedReads1[i];
                for(int i=0;i<mixedReads2.size();i++)
                    delete mixedReads2[i];
                
                //cout<<"Case -1 #################### Unmapped\n";
				mixedReads1.clear();
				mixedReads2.clear();
				
				return;
			}
            //Case2: If one of them i.e. read 1 is mapped and other one(read 2) is unmapped, write in singleton file
            //Because read2 is unmapped, we can use it to fill gap, find which gap read2 falls in
            //zxc
            else if(((read1->flag & 4) == 0 && (read2->flag & 4) != 0) ||
                       ((read1->flag & 4) == 0 && (read2->flag & 4) == 0 && maxDistance > 250))
            {
                case_count[1]++;
                //cout<<"Case -2 ####################\n";
                for(int i=0;i<mixedReads1.size();i++)
                {
                    //writeSam2(mixedReads1[i], single);
                }
                for(int i=0;i<mixedReads2.size();i++)
                {
                    //writeSam2(mixedReads2[i], single);
                }
                //cout<<"Yes, this pair is singleton\n";
                
                for(int i=0;i<mixedReads1.size();i++)
                {
                    contigNo1=mixedReads1[i]->contigNo;
                    pos1=mixedReads1[i]->pos;
                    strandNo1=(mixedReads1[i]->flag&16)>>4;

                    if(samflag==2 && !checkChar(mixedReads2[0]->seq) && !qualityFilter(mixedReads2[0]->qual))
                    {
                        if((read2->flag & 4) != 0)//unmapped
                        {
                            gapFileNo=checkPos(contigNo1,pos1, strandNo1,strlen(mixedReads2[0]->seq));

                            if(gapFileNo>=0 && gaps[gapFileNo]->read_count <=3000)
                            {
                                char temp[MAX_READLENGTH];
                                reverse(temp,mixedReads2[0]->seq);

                                if((strandNo1 ==1 && !check_duplicate(mixedReads2[0]->seq,gapFileNo)) ||
                                    (strandNo1 == 0 && !check_duplicate(temp,gapFileNo)))
                                {
                                    check_pos_count++;
                                    writeSam(mixedReads1[i], gapFileNo);
                                    writeSam(mixedReads2[0], gapFileNo);

                                    if(strandNo1 ==0)strcpy(mixedReads2[0]->seq,temp);

                                    allocateJ(gapFileNo);

                                    strcpy(gaps[gapFileNo]->unmapped_jump_reads[gaps[gapFileNo]->read_count],mixedReads2[0]->seq);

                                    gaps[gapFileNo]->read_count++;
                                }
                            }
                        }
                        else
                        {
                            char original1[MAX_READLENGTH],original2[MAX_READLENGTH];
                            strcpy(original2,mixedReads2[0]->seq);
                            strcpy(original1,mixedReads1[i]->seq);

                            contigNo2=mixedReads2[0]->contigNo;
                            pos2=mixedReads2[0]->pos;
                            strandNo2=(mixedReads2[0]->flag&16)>>4;

                            contigNo1=mixedReads1[i]->contigNo;
                            pos1=mixedReads1[i]->pos;
                            strandNo1=(mixedReads1[i]->flag&16)>>4;

                            gapFileNo=checkPos(contigNo1,pos1, strandNo1,strlen(mixedReads2[0]->seq));

                            if(gapFileNo>=0 && gaps[gapFileNo]->read_count <=3000)
                            {
                                char temp[MAX_READLENGTH];
                                char rev[MAX_READLENGTH];
                                reverse(temp,mixedReads2[0]->seq);

                                if(strandNo2==1)
                                {
                                    strcpy(mixedReads2[0]->seq,temp);
                                }

                                reverse(rev,mixedReads2[0]->seq);

                                if((strandNo1 == 1 && !check_duplicate(mixedReads2[0]->seq,gapFileNo)) ||
                                   (strandNo1 ==0 && !check_duplicate(rev,gapFileNo)))
                                {
                                    check_pos_count++;

                                    writeSam(mixedReads1[i], gapFileNo);
                                    writeSam(mixedReads2[0], gapFileNo);

                                    if(strandNo1 ==0)
                                    {
                                        strcpy(mixedReads2[0]->seq,rev);
                                    }

                                    allocateJ(gapFileNo);

                                    strcpy(gaps[gapFileNo]->unmapped_jump_reads[gaps[gapFileNo]->read_count],mixedReads2[0]->seq);
                                    gaps[gapFileNo]->read_count++;

                                }
                            }

                            strcpy(mixedReads2[0]->seq,original2);
                            strcpy(mixedReads1[i]->seq,original1);

                            gapFileNo=checkPos(contigNo2,pos2, strandNo2,strlen(mixedReads1[i]->seq));


                            if(gapFileNo>=0 && gaps[gapFileNo]->read_count <=3000)
                            {
                                char temp[MAX_READLENGTH];
                                char rev[MAX_READLENGTH];
                                reverse(temp,mixedReads1[i]->seq);

                                if(strandNo1==1)
                                {
                                    //for linking this is extra, if it is form reverse strand, it has to be reversed
                                    //Now we got the read, we will send it to figbird,it will reverse or not based on 1st of the pair
                                    strcpy(mixedReads1[i]->seq,temp);
                                }

                                //From this step, it matches with above unmapped conditions steps
                                reverse(rev,mixedReads1[i]->seq);

                                if((strandNo2 == 1 && !check_duplicate(mixedReads1[i]->seq,gapFileNo)) ||
                                    (strandNo2 == 0 && !check_duplicate(rev,gapFileNo)))
                                {
                                    check_pos_count++;

                                    writeSam(mixedReads2[0], gapFileNo);
                                    writeSam(mixedReads1[i], gapFileNo);

                                    if(strandNo2 ==0)
                                    {
                                        strcpy(mixedReads1[i]->seq,rev);
                                    }

                                    allocateJ(gapFileNo);

                                    strcpy(gaps[gapFileNo]->unmapped_jump_reads[gaps[gapFileNo]->read_count],mixedReads1[i]->seq);
                                    gaps[gapFileNo]->read_count++;

                                }
                            }
                        }
                    }

                    if(samflag==1)//case2
                    {
                        int del = parseDel(mixedReads1[i]->cigar);
                        partial_gapFileNo = checkPos2(contigNo1,pos1, strandNo1,strlen(mixedReads1[i]->seq),del);
                        if(partial_gapFileNo>=0 && gaps[partial_gapFileNo]->partial_read_count <= 3000)
                        {
                            if(check_Ncount_partial(mixedReads1[i]->seq) && !check_duplicate(mixedReads1[i]->seq,partial_gapFileNo))
                            {
                                int match = writePartialSam(mixedReads1[i],partial_gapFileNo,strandNo1,del,-1);

                                allocateP(partial_gapFileNo);
                                strcpy(gaps[partial_gapFileNo]->partial_frag_reads[gaps[partial_gapFileNo]->partial_read_count],mixedReads1[i]->seq);
                                gaps[partial_gapFileNo]->partial_read_count++;
                                partial_case_count[1]++;
                                checkMIM(mixedReads1[i]->cigar,partial_gapFileNo);
                                //cout<<"Chosen read1, Gap = "<<partial_gapFileNo<<"\tMatch = "<<match<<"\tgapstart = "<<gaps[partial_gapFileNo]->gapStart<<", starnd = "<<strandNo1<<"\tpos = "<<pos1<<endl;

                            }
                        }
                    }
                }
                
                for(int i=0;i<mixedReads1.size();i++)
                    delete mixedReads1[i];
                for(int i=0;i<mixedReads2.size();i++)
                    delete mixedReads2[i];
                
                
                mixedReads1.clear();
                mixedReads2.clear();
                
                return;
                
            }
            ///zxc
            /*
            //Case3: If read 1 is unmapped and other one(read 2) is mapped, write in singleton file
            //Because read1 is unmapped, find which gap read1 falls in
            else if(((read1->flag & 4) != 0 && (read2->flag & 4) == 0) ||
                    ((read1->flag & 4) == 0 && (read2->flag & 4) == 0 && maxDistance > 250))
            {
                case_count[2]++;
                //cout<<"Case -3 ####################\n";
                for(int i=0;i<mixedReads1.size();i++)
                {
                    //writeSam2(mixedReads1[i], single);
                }
                for(int i=0;i<mixedReads2.size();i++)
                {
                    //writeSam2(mixedReads2[i], single);
                }
                //cout<<"Yes, this pair is singleton\n";

                for(int i=0;i<mixedReads2.size();i++)
                {
                    contigNo2=mixedReads2[i]->contigNo;
                    pos2=mixedReads2[i]->pos;
                    strandNo2=(mixedReads2[i]->flag&16)>>4;

                    if(samflag==2 && !checkChar(mixedReads1[0]->seq) && !qualityFilter(mixedReads1[0]->qual))
                    {
                        gapFileNo=checkPos(contigNo2,pos2, strandNo2,strlen(mixedReads1[0]->seq));

                        char temp[MAX_READLENGTH];
                        reverse(temp,mixedReads1[0]->seq);
                        if(gapFileNo>=0)
                        {
                            if((strandNo2 == 1  && !check_duplicate(mixedReads1[0]->seq,gapFileNo)) ||
                                (strandNo2 == 0 && !check_duplicate(temp,gapFileNo)))
                            {
                                check_pos_count++;

                                writeSam(mixedReads2[i], gapFileNo);
                                writeSam(mixedReads1[0], gapFileNo);
                            //cout<<"Case-3: "<<", Gap = "<<gapFileNo<<", Pos = "<<mixedReads2[i]->pos<<" , "<<mixedReads2[i]->cigar<<endl;
                                if(strandNo2 == 0)strcpy(mixedReads1[0]->seq,temp);
                                strcpy(gaps[gapFileNo]->unmapped_jump_reads[gaps[gapFileNo]->read_count],mixedReads1[0]->seq);
                                gaps[gapFileNo]->read_count++;

                            }
                        }
                    }

                    if(samflag==1)//case3
                    {
                        int del = -1;

                        del = parseDel(mixedReads2[i]->cigar);
                        partial_gapFileNo = checkPos2(contigNo2,pos2, strandNo2,strlen(mixedReads2[i]->seq),del);
                        if(partial_gapFileNo>=0 && gaps[partial_gapFileNo]->partial_read_count <= 3000)
                        {
                            if(check_Ncount_partial(mixedReads2[i]->seq) && !check_duplicate(mixedReads2[i]->seq,partial_gapFileNo))
                            {
                                int match = writePartialSam(mixedReads2[i],partial_gapFileNo,strandNo2,del,-1);

                                allocateP(partial_gapFileNo);

                                strcpy(gaps[partial_gapFileNo]->partial_frag_reads[gaps[partial_gapFileNo]->partial_read_count],mixedReads2[i]->seq);
                                gaps[partial_gapFileNo]->partial_read_count++;
                                partial_case_count[1]++;
                                checkMIM(mixedReads2[i]->cigar,partial_gapFileNo);
                                //cout<<"Chosen read2, Gap = "<<partial_gapFileNo<<"\t Match = "<<match<<"\tgapstart = "<<gaps[partial_gapFileNo]->gapStart<<", starnd = "<<strandNo2<<"\tpos = "<<pos2<<endl;

                            }
                        }
                    }
                }

                for(int i=0;i<mixedReads1.size();i++)
                    delete mixedReads1[i];
                for(int i=0;i<mixedReads2.size();i++)
                    delete mixedReads2[i];
                
                mixedReads1.clear();
                mixedReads2.clear();
                return;
            }
            */
            ///zxc
            //case4: Those that are both mapped but from different contigs
            else if(strcmp(read1->rname,read2->rname)!=0)
			{
			    //cout<<"Case -4 ####################Linker Read\n";
			    case_count[3]++;
                read1->ih=ih1;
                read2->ih=ih2;
                //readCounts[read1->contigNo][read2->contigNo]++;
                //writeSam2(read1,linking);
                //writeSam2(read2, linking);

			}
		}
	}
    
	for(int i=0;i<mixedReads1.size();i++)
		delete mixedReads1[i];
	for(int i=0;i<mixedReads2.size();i++)
		delete mixedReads2[i];
	
    
	mixedReads1.clear();
	mixedReads2.clear();
    
}

SAM *getSAM(char *line)
{
	SAM *sam=new SAM;
	char *temp;
    
	temp=strtok(line,"\t\n ");
	strcpy(sam->qname,temp);
    
	temp=strtok(NULL,"\t\n ");
	sam->flag=atoi(temp);
    
	temp=strtok(NULL,"\t\n ");
	strcpy(sam->rname,temp);
    
	temp=strtok(NULL,"\t\n ");
	sam->pos=atoi(temp);


	temp=strtok(NULL,"\t\n ");
	sam->mapq=atoi(temp);
    
	temp=strtok(NULL,"\t\n ");
	strcpy(sam->cigar,temp);
	
	temp=strtok(NULL,"\t\n ");
	strcpy(sam->rnext,temp);
    
	temp=strtok(NULL,"\t\n ");
	sam->pnext=atoi(temp);
    
	temp=strtok(NULL,"\t\n ");
	sam->tlen=atoi(temp);
    
	temp=strtok(NULL,"\t\n ");
	strcpy(sam->seq,temp);
	
	temp=strtok(NULL,"\t\n ");
	strcpy(sam->qual,temp);
    
	sam->nm=-1;
    
	while((temp=strtok(NULL,"\t\n "))!=NULL)
	{
		if(temp[0]=='M' && temp[1]=='D')
		{
			strcpy(sam->md,temp);
		}
		if(temp[0]=='N' && temp[1]=='M')
		{
			sam->nm=atoi(&temp[5]);
		}
		if(temp[0]=='A' && temp[1]=='S')
        {
            sam->as=atoi(&temp[5]);
        }
	}
    
	sam->contigNo=getContigNo(sam->rname);
    
	return sam;
}

void printHelp()
{
    
	cerr<<"Invalid parameters"<<endl;
	exit(1);
    
}

void printStat(int numgap_contig[])
{
    for (int i = 0; i<noContigs; ++i)
    {
        int flag=0,index_count=0;;
        for(int j=0;j<gaps.size();j++)
        {
            if(gaps[j]->contigNo == i)
            {
                if(flag==0)
                {
                    //cout<<"Contig "<<i<<"/"<<contigNames[i]<<",\tlen =  "<<contigLengths[i]<<" :\t";flag=1;
                }
                //if(samflag==2)cout<<" Gap = "<<j<<"\tstart = "<<gaps[j]->gapStart<<"\tgaplen = "<<gaps[j]->gapLength<<"\t#Reads = "<<gaps[j]->read_count<<"||";
                //else cout<<" Gap = "<<j<<"\tstart = "<<gaps[j]->gapStart<<"\tgaplen = "<<gaps[j]->gapLength<<"\t#Reads = "<<gaps[j]->partial_read_count<<"||";

                numgap_contig[i]++;
            }
        }
        //if(flag)cout<<endl;
    }

    int gap_index=0;
    for(int i=0;i<noContigs;i++)
    {
        if(numgap_contig[i])
        {
            //cout<<"Contig "<<i<<"===>Num gaps = "<<numgap_contig[i];

       }
        for(int j=0;j<numgap_contig[i];j++)
        {
            double frac= gaps[gap_index]->gapLength*100.0/contigLengths[i];
            if(numgap_contig[i]==1)
            {
                //cout<<"\tGap = \t"<<gap_index<<"\tStart = "<<gaps[gap_index]->gapStart<<setw(16)<<"\tGaplength = "<<gaps[gap_index]->gapLength<<setw(20)<<"\tContig Length = "<<contigLengths[i]<<" - "<<frac<<"%";
            }

            else
            {
                //cout<<"\n\t\t\t\t\t\t\t\tGap = \t"<<gap_index<<"\tStart = "<<gaps[gap_index]->gapStart<<setw(16)<<"\tGaplength = "<<gaps[gap_index]->gapLength<<setw(20)<<"\tContig Length = "<<contigLengths[i]<<" - "<<frac<<"%";
            }
            //if(frac <50)cout<<"\tAvailable";
            gap_index++;
        }

        //if(numgap_contig[i]>0)cout<<endl;
    }
}

void findOverlappedGap(int numgap_contig[],int mean)
{
    int gap_index=0,pg=0;
    int maxDistance_thresh=mean;

    for(int i=0;i<noContigs;i++)
    {
        int ng = numgap_contig[i];
        if(ng>1)
        {
            int left_diff=0,right_diff=0;
            for(int j=0;j<ng;j++)//j = gap index
            {
                if(j==0)
                {
                    left_diff = gaps[gap_index]->gapStart;
                    right_diff = gaps[gap_index+1]->gapStart - (gaps[gap_index]->gapStart + gaps[gap_index]->gapLength);
                    if(right_diff < maxDistance_thresh)
                    {
                        //cout<<"Gap Jam  Found, gap = "<<gap_index<<"\tLength = "<<gaps[gap_index]->gapLength<<", in contig = "<<i<<", left side = "<<left_diff<<", right side = "<<right_diff<<endl;
                    }
                }
                else if(j == ng-1)
                {
                    left_diff = gaps[gap_index]->gapStart - (gaps[gap_index-1]->gapStart + gaps[gap_index-1]->gapLength);
                    right_diff = contigLengths[i]- (gaps[gap_index]->gapStart + gaps[gap_index]->gapLength);
                    if(left_diff < maxDistance_thresh)
                    {
                        //cout<<"Gap Jam  Found, gap = "<<gap_index<<"\tLength = "<<gaps[gap_index]->gapLength<<", in contig = "<<i<<", left side = "<<left_diff<<", right side = "<<right_diff<<endl;
                    }
                }
                else
                {
                    left_diff = gaps[gap_index]->gapStart - (gaps[gap_index-1]->gapStart + gaps[gap_index-1]->gapLength);
                    right_diff = gaps[gap_index+1]->gapStart - (gaps[gap_index]->gapStart + gaps[gap_index]->gapLength);

                    if(left_diff < maxDistance_thresh || right_diff < maxDistance_thresh)
                    {
                        //cout<<"Gap Jam  Found, gap = "<<gap_index<<"\tLength = "<<gaps[gap_index]->gapLength<<", in contig = "<<i<<", left side = "<<left_diff<<", right side = "<<right_diff<<endl;
                    }
                }
                if(left_diff < 2000 && right_diff < 2000)
                {
                    //cout<<"Count = "<<pg<<"\tGap = "<<gap_index<<"\tLength = "<<gaps[gap_index]->gapLength<<endl;
                    pg++;
                    //gaptofill[gap_index] = 0;
                }
                gap_index++;
            }
        }
        else if(ng==1)gap_index++;
    }
    //cout<<"Total = "<<pg<<endl;
}


void collectPartialSAM(SAM* read,int pos2)
{
    int contigNo=read->contigNo;
    int pos=read->pos;
    int strandNo=(read->flag&16)>>4;
    
    int del = -1;
    del = parseDel(read->cigar);
    
    partial_gapFileNo = checkPos2(contigNo,pos, strandNo,strlen(read->seq),del);

    if(partial_gapFileNo>=0 && gaps[partial_gapFileNo]->partial_read_count <= 3000)
    {
        //cout<<"Before W"<<endl;
        if(check_Ncount_partial(read->seq) && !check_duplicate(read->seq,partial_gapFileNo))
        {
            int match = writePartialSam(read,partial_gapFileNo,strandNo,del,pos2);

            allocateP(partial_gapFileNo);

            strcpy(gaps[partial_gapFileNo]->partial_frag_reads[gaps[partial_gapFileNo]->partial_read_count],read->seq);
            gaps[partial_gapFileNo]->partial_read_count++;
            partial_case_count[0]++;
            
            checkMIM(read->cigar,partial_gapFileNo);
        }
    }
}

void reWriteReadset(FILE* output1,FILE* output2, SAM* read1,SAM* read2)
{
    char temp1[MAX_READLENGTH];
    char swap;
    int strandNo=(read1->flag&16)>>4;
    //cout<<"read- "<<(linecount/2+1)<<"\t"<<strandNo;
    if(strandNo==1)
    {
        reverse(temp1,read1->seq);
        int len = strlen(read1->qual);
        for(int i=0;i<len/2;i++)
        {
            swap = read1->qual[i];
            read1->qual[i] = read1->qual[len - i - 1];
            read1->qual[len - i - 1] = swap;
        }
        fprintf(output1,"@%s\n%s\n+\n%s\n",read1->qname,temp1,read1->qual);
    }
    else fprintf(output1,"@%s\n%s\n+\n%s\n",read1->qname,read1->seq,read1->qual);

    strandNo=(read2->flag&16)>>4;
    //cout<<"read- "<<(linecount/2+1)<<"\t"<<strandNo;
    if(strandNo==1)
    {
        reverse(temp1,read2->seq);
        int len = strlen(read2->qual);
        for(int i=0;i<len/2;i++)
        {
            swap = read2->qual[i];
            read2->qual[i] = read2->qual[len - i - 1];
            read2->qual[len - i - 1] = swap;
        }
        fprintf(output2,"@%s\n%s\n+\n%s\n",read2->qname,temp1,read2->qual);
    }
    else fprintf(output2,"@%s\n%s\n+\n%s\n",read2->qname,read2->seq,read2->qual);
}

//zxcf
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
	    if(index>MAX_FRAGMENT_SIZE)
		{
			discardedReads++;
			return;
		}
	}
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
		contigNo=atol(rname);
		if(contigLengths[contigNo]>0)
			updateInsertCounts(insertSize);
	}
}



int main(int argc, char *argv[])
{
	time_t nowt,endt;
	time(&nowt);
	
	int MAX_FILE_READ = 0;
	noContigs=0;
	long int contigLength=0;
	unsigned long read;
	long int bufferLength=1024;
    long int tempContigLength=0;
    int nStart=0;//How many gaps are there?
    long int nStartPos=0;//What are those gaps start positions?
    int nCount=0,afterCount=0;//ncount is number of NNNN in a particular gap

    int r_c=0, fullmap=0,partmap=0,partial_uncount=0,myoutflag=1;
    int pm_count=0,lenr;
    int totlinecounter=0,myoutcounter=0,partialcounter=0,writeflag=0,reWritecall=0;
	
    char *line= new char[MAX_REC_LEN];
    tempwrite=new char[10];
    cig1 = new char[10];
    cig2= new char[10];
    strcpy(cig1,"101M");
    strcpy(cig2,"101M");
    cig1[4]='\0';
    cig2[4]='\0';
    //char *templine= new char[MAX_REC_LEN];

    char *contigFileName=argv[1];
    maxDistance=atoi(argv[2]);
    samflag = atoi(argv[3]);
    char * mapFileName = argv[4];
    char* outfilename = argv[5];

	char * filled_contigfilename = argv[6];

	int genome_reduction=atoi(argv[12]);
    int read_reduction=atoi(argv[13]);
    

	char *newcontig;
	char *contigName;

    char *contig;
	contig=new char[bufferLength];
	contig[0]='\0';

	MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);
    //======================================Start reading the contig file=======================================
    
    map<int, int> contignums;
	int gapcount=0;

	//enn

    if(genome_reduction==1)
    {

    contigFile=fopen(filled_contigfilename, "r");

    if(contigFile == NULL)
    {
        cerr<<"Can't open contig file\n";
        exit(1);
    }

    //It contains x scaffolds in total with gaps

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
                contigs.push_back(contig);
                contigLengths.push_back(contigLength);
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

    for(int i=0;i<noContigs;i++)//For each contig
    {
        for(long int j=0;j<contigLengths[i];j++)//For the entire contig length
        {
            if(contigs[i][j]=='N' || contigs[i][j]=='n')//If there is a gap at that position
            {
                if(nStart==0)
                {
                    nStart=1;
                    nCount=1;
                    nStartPos=j;
                }
                else
                {
                    nCount++;
                }
            }
            else if(nStart==1)//if there was ever any gap
            {
                if(nCount>=1)
                {
                    contignums.insert(std::pair<int, int>(gapcount, i));
                    gapcount++;
                }
                nStart=0;
            }
        }
    }

    contigs.clear();
    contigLengths.clear();

    nStart=0;
    nStartPos=0;
    bufferLength=1024;
    tempContigLength=0;
    noContigs=0;
    contigLength=0;
    contig=new char[bufferLength];
    contig[0]='\0';

    }
	///===================================================================================================
	contigFile=fopen(contigFileName, "r");

    if (contigFile == NULL)
	{
		cerr<<"Can't open contig file\n";
		exit(1);
	}

    //It contains x scaffolds in total with gaps

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
			contigNames.push_back(strtok(contigName," \t\n"));//In case of multiple contigs, 
			                                                    //we are pushing each name into a vector of char array
			
			if(contigLength>0)//This block is accessed after completing the full reading of each contig, not at the start
			{
				//cout<<"Found contig - "<<noContigs<<"\t"<<contigNames[noContigs]<<endl;
				noContigs++;
				contig=(char *)realloc(contig, (contigLength+1)*sizeof(char));
				contigs.push_back(contig);
				contigLengths.push_back(contigLength);
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

    //=================================End of reading the contig file==========================================//

    //=================================Gather Info about gaps in the contig====================================//

    char gif[1000];
    strcpy(gif,argv[10]);
    strcat(gif,"gapInfo.txt");

    gapInfoFile=fopen(gif,"w");
    gapcount=0;

    for(int i=0;i<noContigs;i++)//For each contig
    {
        for(long int j=0;j<contigLengths[i];j++)//For the entire contig length
        {
            if(contigs[i][j]=='N' || contigs[i][j]=='n')//If there is a gap at that position
            {
                if(nStart==0)
                {
                    nStart=1;
                    nCount=1;
                    nStartPos=j;
                }
                else
                {
                    nCount++;
                }
            }
            else if(nStart==1)//if there was ever any gap
            {
                if(nCount>=1)
                {
                    //cout<<"Found new gap\n";
                    int contigNotowrite;
                    Gap * g=new Gap();

                    g->contigNo=i;
                    contigNotowrite = i;

                    g->gapStart=nStartPos;
                    g->gapLength=nCount;
                    g->read_count=0;
                    g->partial_read_count=0;
                    g->total_allocated_rows = 0;

                    gaps.push_back(g);

                    if(genome_reduction==1)
                    {
                        std::map<int, int>::iterator it = contignums.find(gapcount);
                        if (it != contignums.end())
                        {

                            contigNotowrite = it->second;
                        }
                    }
                    
                    fprintf(gapInfoFile,"%d\t%ld\t%d\n",contigNotowrite,g->gapStart,g->gapLength);
                    if(g->gapLength <=100)small_gap_count++;
                    else large_gap_count++;
                    if(g->gapLength >900)huge_gap++;
                    gapcount++;

                }
                nStart=0;
            }
        }
    }
    //enn


    noGaps = small_gap_count + large_gap_count;

    gaptofill = new int[noGaps];
    perfectread_gap = new int[noGaps];
    perfectread_gaplen = new int[noGaps];

    read_in_gap_f = new int[noGaps];
    read_in_gap_r = new int[noGaps];

    for(int i=0;i<noGaps;i++)
    {
        gaptofill[i] = 1;
        perfectread_gap[i] = 0;
        perfectread_gaplen[i] = 0;
        read_in_gap_r[i] = read_in_gap_f[i] =0;
    }

    //===============================================Done with gaps=================================//
    
    initHashTable();

	strcpy(gif,argv[10]);
    strcat(gif,"stat.txt");

	statFile=fopen(gif,"w");

	strcpy(gif,argv[10]);
    strcat(gif,"stat2.txt");

	statFile2=fopen(gif,"w");
    //cout<<"No. contigs = "<<noContigs<<"\tNo Gaps = "<<noGaps<<endl;
    
	//=======Start working with the SAM file output by Bowtie2====================
	
	//=================================Preprocessing works========================

	

	//cerr<<outfilename<<endl;

	outFile=fopen(outfilename,"w");

	if (outFile == NULL)
    	{
    		cerr<<"Can't create myout file\n";
    		exit(1);
    	}

    gapFiles=new char*[gaps.size()];
    partial_gapFiles = new char*[gaps.size()];
    
    for(int i=0;i<gaps.size();i++)
    {
      gapFiles[i] = new char[400];
      partial_gapFiles[i] = new char[400];
    }

    char *gapFileName=new char[1000];
    char *partial_gapFileName=new char[1000];

    char *temp=new char[2000];

    if(samflag==2){
    for(int i=0;i<gaps.size();i++)
	{
        strcpy(gapFileName,argv[9]);
        strcat(gapFileName,"gaps_");
        sprintf(temp,"%d",i);
        strcat(gapFileName,temp);
        strcat(gapFileName,".sam");
        FILE* fp = fopen(gapFileName,"w");
        fclose(fp);
        gapFileName[strlen(gapFileName)]='\0';
        strcpy(gapFiles[i],gapFileName);
    }}
    else{
    for(int i=0;i<gaps.size();i++)
    {
        strcpy(partial_gapFileName,argv[9]);
        strcat(partial_gapFileName,"partial_gaps_");
        sprintf(temp,"%d",i);
        strcat(partial_gapFileName,temp);
        strcat(partial_gapFileName,".sam");
        FILE* fp = fopen(partial_gapFileName,"w");
        fclose(fp);
        partial_gapFileName[strlen(partial_gapFileName)]='\0';
        strcpy(partial_gapFiles[i],partial_gapFileName);
    }}

	int it=0; 
    int end=0;
    
	char preqname1[100];
	char preqname2[100];
	
	strcpy(preqname1,"*");
	strcpy(preqname2,"*");
    
    char qname[200];
    int flag_segment;

    //===============================Starting the read of SAM file result.sam line by line================

    int def = atoi(argv[11]);

    if(def == 1)
    {
        if(read_reduction == 1)
            writeflag=1;
    }
    else
    {
        if(read_reduction == 1 && samflag ==1)writeflag=1;
    }

    if(writeflag)
    {

        string s1(argv[7]);
        string s2(argv[8]);

        size_t found1 = s1.find_last_of(".");
        size_t found2 = s2.find_last_of(".");

        string s3 = s1.substr(0,found1);
        string s4 = s2.substr(0,found2);

        string ext = s1.substr(found1,s3.size());
        string r = "_reduced";

        s3 += r + ext;
        s4 += r + ext;

        output1 = fopen(s3.c_str(),"w");
        output2 = fopen(s4.c_str(),"w");

        if(output1 == NULL || output2 == NULL)
        {
            cerr<<"Can't create reduced read pair during preproscessing...exiting.\n";
            exit(1);
        }

        cout<<s3<<endl;
        cout<<s4<<endl;
    }
//enn

        mapFile=fopen(mapFileName, "r");

        if (mapFile == NULL)
	    {
		    cerr<<"Can't open alignment file\n";
		    exit(1);
	    }
	
	    if(samflag==2 && maxDistance > 250)
	    {
	        while(fgets(line, MAX_FILE_READ, mapFile)!=NULL)
	        {
		        if(line[0]=='@')//Skip headers
			        continue;
               
		        SAM *read1=getSAM(line);
		       
                while((read1->flag & 2) == 0)
                {
                    strcpy(qname,read1->qname);
                    flag_segment=read1->flag & 192;//Only interested in MSBth and (MSB-1)th position of the flag according to the AND op.
                    
                    while(strcmp(qname,read1->qname)==0 && (read1->flag & 192) == flag_segment)
                    {
                        delete read1;
                        if(fgets(line, MAX_FILE_READ, mapFile)!=NULL)
                        {
                            read1=getSAM(line);          		
                        }
                    }
                    strcpy(qname,read1->qname);
                    flag_segment=read1->flag & 192;

                    while(strcmp(qname,read1->qname)==0 && (read1->flag & 192) == flag_segment)
                    {
                        delete read1;

                        if(fgets(line, MAX_FILE_READ, mapFile)!=NULL)
                        {
                            read1=getSAM(line);
                        }
                        else
                        {                  
                            end=1;
                            break;
                        }
                    }
                   
                    if(end==1)
                    {
                        break;
                    }
                }

                if(end==1)
                {
                    break;
                }

		        fgets(line, MAX_FILE_READ, mapFile);

		        SAM *read2=getSAM(line);
				       
		        if(strcmp(read1->qname,preqname1)!=0 || strcmp(read2->qname,preqname2)!=0)
		        {
			        strcpy(preqname1,read1->qname);
			        strcpy(preqname2,read2->qname);

			        printVectors(outFile);

                    reads1.push_back(read1);
			        reads2.push_back(read2);
		        }
		        else
		        {
			        reads1.push_back(read1);
			        reads2.push_back(read2);
		        }
	        }

	        fclose(mapFile);
            printVectors(outFile);
            
            fclose(outFile);
	        
            
            //================================
            
            initInsertCounts(MAX_FRAGMENT_SIZE);
            
            outFile=fopen(outfilename, "r");
            
	        if (outFile == NULL)
            {
	            printf("Can't open new out file\n");
	            exit(1);
            }
            while(fgets(line1, MAX_FILE_READ, outFile)!=NULL)
            {
	
	            if(line1[0]=='@')
		            continue;
                
	            processMapping(line1);
	            
            }
            //cout<<"Myout.sam file has "<<count<<" lines"<<endl;
            fclose(outFile);
            
            //====================
            //zxcf
            
            long int insCount=discardedReads;
        
	        double sum=0;
            
	        for(int i=0;i<maxInsertSize;i++)
	        {
		        insCount+=(insertCounts[i]-1);
		        sum+=i*(insertCounts[i]-1);
                
	        }
	        
	        read_mean=sum/insCount;
	        //cerr<<"RM = "<<read_mean<<endl;
	        //================================
	        
	        mapFile=fopen(mapFileName, "r");

            if (mapFile == NULL)
	        {
		        printf("Can't open map file\n");
		        exit(1);
	        }
	        
	        it=0; 
            end=0;
    
	        strcpy(preqname1,"*");
	        strcpy(preqname2,"*");    
	    }

        while(fgets(line, MAX_FILE_READ, mapFile)!=NULL)
    	{
    		if(line[0]=='@')//Skip headers
    			continue;

    		it++;

    		SAM *read1=getSAM(line);
    		totlinecounter++;


            while((read1->flag & 2) == 0)//will enter here if this read is "NOT" part of a pair that aligned in a pair end  fashion
            {
                strcpy(qname,read1->qname);
                flag_segment=read1->flag & 192;//Only interested in MSBth and (MSB-1)th position of the flag according to the AND op.
                //192 = 1100 0000
                //MSB = Indicates if this read is mate 2 in the pair
                //MSB-1 = Indicates if this read is mate 1 in the pair

                while(strcmp(qname,read1->qname)==0 && (read1->flag & 192) == flag_segment)
                {
                    mixedReads1.push_back(read1);


                    if(fgets(line, MAX_FILE_READ, mapFile)!=NULL)
                    {

                        read1=getSAM(line);
                        		totlinecounter++;
                        //if(totlinecounter < 1000)cout<<"In 1st while - "<<totlinecounter<<"\tN = "<<read1->qname<<endl;
                    }
                }
                strcpy(qname,read1->qname);
                flag_segment=read1->flag & 192;

                while(strcmp(qname,read1->qname)==0 && (read1->flag & 192) == flag_segment)
                {
                    mixedReads2.push_back(read1);


                    if(fgets(line, MAX_FILE_READ, mapFile)!=NULL)
                    {

                        read1=getSAM(line);
                        		totlinecounter++;
                        //if(totlinecounter < 1000)cout<<"In 2nd while - "<<totlinecounter<<"\tN = "<<read1->qname<<endl;
                    }
                    else
                    {
                        //cout<<"Never in this block\n";
                        end=1;
                        break;
                    }
                }
                //They are called mixed because this pair did not map in pair end fashion

                if(writeflag)
                {
                    reWriteReadset(output1,output2,mixedReads1[0],mixedReads2[0]);
                    reWritecall++;

                }

                printMixedVectors();

                pm_count++;
                if(end==1)
                {
                   //cout<<"Never in this block\n";
                    break;
                }
            }

            if(end==1)
            {
                //Never in this block
                break;
            }


    		fgets(line, MAX_FILE_READ, mapFile);

    		SAM *read2=getSAM(line);
    				totlinecounter++;


    		//if(totlinecounter < 1000)cout<<"2nd read - "<<totlinecounter<<"\tN = "<<read2->qname<<endl;

    		if(strcmp(read1->qname,preqname1)!=0 || strcmp(read2->qname,preqname2)!=0)
    		{
    			strcpy(preqname1,read1->qname);
    			strcpy(preqname2,read2->qname);

    			
    			if(!(samflag==2 && maxDistance > 250))printVectors(outFile);
                //if(totlinecounter < 1000)cout<<"PV ret-> "<<"\tN1 = "<<read1->qname<<"\tN2 = "<<read2->qname<<endl;

                if(samflag==1)//partial
                {
                    int f1 = strcmp(read1->cigar,cig1);
                    int f2 = strcmp(read2->cigar,cig2);

                    if(f1 == 0 && f2 == 0)
                    {
                        fullmap++;
                    }
                    else
                    {
                        collectPartialSAM(read1,read2->pos);
                        collectPartialSAM(read2,read1->pos);
                        partialcounter++;

                        if(writeflag)
                        {
                            reWriteReadset(output1,output2,read1,read2);
                            reWritecall++;
                        }
                    }

                    if(!(f1==1 && f2==1))partmap++;
                    else partial_uncount++;
                }

    			if(!(samflag==2 && maxDistance > 250))
    			{
    			    reads1.push_back(read1);
    			    reads2.push_back(read2);
                    r_c++;
    		    }
    		}
    		else
    		{
    			if(!(samflag==2 && maxDistance > 250))
    			{
    			    reads1.push_back(read1);
    			    reads2.push_back(read2);
    			}
    			    //cout<<"Never Here, else of printvector\n";
    		}
    		
    		if(samflag==2 && maxDistance > 250)
		    {
                delete read1;
                delete read2;
		    }
    		//if(totlinecounter % 5000000 ==0)cout<<"Processed upto - "<<totlinecounter<<endl;
    	}


        if(!(samflag==2 && maxDistance > 250))printVectors(outFile);

    /// =========================Print Stats====================================
    //cout<<"\nTotal lines in sam file = "<<totlinecounter<<endl;
    //cout<<"Rewrite called = "<<reWritecall<<endl;
    //cout<<"\nTotal pair end mapped count = "<<r_c<<"\tPerct. = "<<(r_c * 100.0/(r_c+unCount))<<"%\n";
    //cout<<"Total unmapped read count = "<<unCount<<"\tPerct. = "<<(unCount * 100.0/(r_c+unCount))<<"%"<<endl;

    //if(maxDistance <= 250)
    {
        //cerr<<"\tPerfect full map 101M = "<<fullmap<<"\tPartial map = "<<partmap<<"\tUnmapped? = "<<partial_uncount<<endl;
        //cerr<<"\tPerct. of full map = "<<(fullmap * 100.0/(fullmap+partmap))<<"%"<<endl;
        //cerr<<"\tPerct. of part map = "<<(partmap * 100.0/(fullmap+partmap))<<"%"<<endl;
    }

    //cerr<<"Total partial reads in gap = "<<partial_check_pos_count<<endl;

    //cerr<<"Total reads in gap = "<<check_pos_count<<endl;

    //cout<<"printmixedvector is called = "<<pm_count<<" times"<<endl;

    //cout<<"No. of small gaps(<=100) = "<<small_gap_count<<", no. of large gaps = "<<large_gap_count<<"\tGreater than 900 length = "<<huge_gap<<endl;
    //cout<<"Partial reads in forward strand = "<<pos_strand<<", reverse strand = "<<neg_strand<<endl;

    //cout<<"Total Partial  from both end mapped = "<<partial_case_count[0]<<"\tFrom unmapped = "<<partial_case_count[1]<<endl;
    //cout<<"Total count of reads = "<<totalCount<<"\tmyout.sam = "<<r_c_p<<"\t101M in myout = "<<rcp<<endl;


    fprintf(statFile,"%ld %ld %ld %ld",totalCount, unCount, maxReadLength, MAX_FRAGMENT_SIZE);

    for(int i=0;i<noGaps;i++)
    {
        fprintf(statFile2,"%d\t%d\t%d\n",gaptofill[i],perfectread_gap[i],perfectread_gaplen[i]);
    }

    for( int i=0;i<gaps.size();i++)
    {
         //cout<<"Gap = "<<i<<",# partial read(f) = "<<read_in_gap_f[i]<<", read(r) = "<<read_in_gap_r[i]<<endl;
         if(samflag==2)
         {
            for(int j=0;j<gaps[i]->total_allocated_rows;j++)
                free(gaps[i]->unmapped_jump_reads[j]);

            if(gaps[i]->total_allocated_rows !=0)free(gaps[i]->unmapped_jump_reads);
         }
         else
         {
           //cout<<i<<"\t"<<gaps[i]->total_allocated_rows<<endl;
            for(int j=0;j<gaps[i]->total_allocated_rows;j++)
                free(gaps[i]->partial_frag_reads[j]);

            if(gaps[i]->total_allocated_rows !=0)free(gaps[i]->partial_frag_reads);

         }
    }

	delete perfectread_gap;
    delete perfectread_gaplen;
    delete gaptofill;
    delete read_in_gap_f;
    delete read_in_gap_r;

	fclose(gapInfoFile);
    fclose(mapFile);
	fclose(statFile);
    fclose(statFile2);
    if(!(samflag==2 && maxDistance > 250))fclose(outFile);

    if(writeflag)
    {
        fclose(output1);
        fclose(output2);
    }

//enn
    time(&endt);

    //cerr<<"Time to preprocess the reads = "<<difftime(endt,nowt)<<" seconds"<<endl;

    return 0;
}
