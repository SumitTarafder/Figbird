#include <iostream>
#include <bits/stdc++.h>
#include <string>
#include <limits>
#include <vector>
#include <cstring>

using namespace std;

#include <math.h>

#define MAX_REC_LEN 10024

int newlength=30;

int pos[3];
int num_itr;
int totalgaps = 0;

char* gappedpath;

typedef struct Gaps{

    char *gapstring,*gapleft,*gapright;
    int left_start_N,right_end_N;
    int fully_closed,originalGap,finalGapLen,r_size;
}Gap;

Gap * gaps;
Gap **temp_gaps;

int checkComplete(char * gapstr)
{
    int Nstart=0,region_count=0,Ncount=0;

    int len = strlen(gapstr);
    Nstart=0,region_count=0;

    for(int i=0;i<len;i++)
    {
        if(gapstr[i] =='N' && Nstart ==0)
        {
            Nstart = 1;
            pos[0] = i;
            Ncount++;
        }
        else if(gapstr[i] != 'N' && Nstart ==1)
        {
            Nstart=0;
            region_count++;
            pos[1] = i-1;
            pos[2] = Ncount;
        }
        if(i==len-1 && Nstart==1)
        {
            region_count++;
            pos[1] = i;
            pos[2] = Ncount;
        }
    }

    return region_count;
}

Gap combine(Gap g,int org,int gaplen,char *s,int rc,int itr)
{
    if(itr==1)g.originalGap=org;
    g.fully_closed = 1 - rc;///1- regioncount, so when count is 0 ,fullyclosed =1

    if(gaplen == 0)
    {
        g.gapleft = g.gapright = NULL;
        g.gapstring = new char[1];
        g.gapstring[0]='\0';
        g.finalGapLen=0;
        //cout<<org<<" "<<gaplen<<" "<<g.left_start_N<<" "<<g.right_end_N<<" "<<g.r_size<<endl;
        return g;
    }

    if(itr==1)
    {
        //cout<<gaplen<<endl;
        //cout<<endl;
        //cout<<s<<endl;
        g.gapstring= new char[gaplen+1];
        strcpy(g.gapstring,s);
        g.gapleft=NULL;
        g.gapright=NULL;
        g.finalGapLen = gaplen;


        if(g.fully_closed != 1)
        {
            checkComplete(g.gapstring);
            g.left_start_N = pos[0];
            g.right_end_N = pos[1];
            g.r_size = gaplen - g.right_end_N;
        }
        //cout<<org<<" "<<gaplen<<" "<<g.left_start_N<<" "<<g.right_end_N<<" "<<g.r_size<<endl;
    }
    else
    {
        int newlen = g.left_start_N + gaplen + g.r_size;
        //cout<<org<<" "<<gaplen<<" "<<newlen<<" "<<g.left_start_N<<" "<<g.right_end_N<<" "<<g.r_size<<endl;
        char* newgap_str = new char[newlen];
        newgap_str[newlen-1]='\0';

        int index_count=0;
        for(int i=0;i<g.left_start_N;i++)newgap_str[index_count++] = g.gapstring[i];
        for(int i=0;i<gaplen;i++)newgap_str[index_count++] = s[i];
        for(int i=0;i<g.r_size-1;i++)newgap_str[index_count++] = g.gapstring[i+1+g.right_end_N];

        delete[] g.gapstring;

        g.gapstring = new char[newlen];
        strcpy(g.gapstring,newgap_str);
        checkComplete(g.gapstring);
        g.left_start_N = pos[0];
        g.right_end_N = pos[1];
        g.finalGapLen = newlen-1;
        g.r_size = g.finalGapLen - g.right_end_N;
    }
    return g;
}

void outFile(Gap *g,int n)
{
    char dir[1000];
    strcpy(dir,gappedpath);
    strcat(dir,"combined_gapstring.txt");
    FILE* fp = fopen(dir,"w");
    for(int i=0;i<n;i++)
    {
        //fprintf(fp,"Gap = %d\n",i);
        //fprintf(fp,"%s\n",g[i].gapleft);
        //fprintf(fp,"%s\n",g[i].gapright);
        fprintf(fp,"%s\n",g[i].gapstring);
    }
    fclose(fp);

}

void initialize(Gap *g,int n)
{
    for(int i=0;i<n;i++)
    {
        //fprintf(fp,"Gap = %d\n",i);
        //fprintf(fp,"%s\n",g[i].gapleft);
        //fprintf(fp,"%s\n",g[i].gapright);
        g[i].gapstring = g[i].gapleft = g[i].gapright = NULL;
        g[i].fully_closed = 0;
        g[i].left_start_N = g[i].right_end_N = -1;
    }
}

void saveGapInfo(int itr)
{
    for(int i=0;i<totalgaps;i++)
    {
        temp_gaps[itr][i].originalGap = gaps[i].originalGap;
        temp_gaps[itr][i].finalGapLen = gaps[i].finalGapLen;
        temp_gaps[itr][i].gapstring = new char[temp_gaps[itr][i].finalGapLen+1];
        strcpy(temp_gaps[itr][i].gapstring,gaps[i].gapstring);


    }
}

int main(int argc, char *argv[])
{
    FILE * contigFile=NULL,*filledContigFile=NULL,*new_gap_file,*gapstrfile,*gapoutfile,*gap_cases;

	int MAX_FILE_READ = 0;
	char *line= new char[MAX_REC_LEN];
    char *line2= new char[MAX_REC_LEN];

    MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);

    char gapoutFilename[100];
    char *temp=new char[100];
    
    num_itr=atoi(argv[1]);
    gappedpath=argv[2];
    char dummy[10];
    strcpy(dummy,"");
   
    for(int itr=1;itr<=num_itr;itr++)
    {
        strcpy(gapoutFilename,gappedpath);
        strcat(gapoutFilename,"gapout_");
        sprintf(temp,"%d",itr);
        strcat(gapoutFilename,temp);
        strcat(gapoutFilename,".txt");

        FILE * gapoutfile = fopen(gapoutFilename,"r");

        if(gapoutfile == NULL)
        {
            printf("Can't open gapout txt file\n");
            exit(1);
        }

        if(itr ==1)
        {
            while(fgets(line, MAX_FILE_READ, gapoutfile)!=NULL)
            {
                totalgaps++;
            }
            fclose(gapoutfile);

            //cout<<"tot = "<<totalgaps<<endl;

            gaps = new Gap[totalgaps];
            temp_gaps = new Gap*[num_itr];

            for(int i=0;i<num_itr;i++)
            {
                temp_gaps[i] = new Gap[totalgaps];

            }

            initialize(gaps,totalgaps);
            for(int i=0;i<num_itr;i++)initialize(temp_gaps[i],totalgaps);

            gapoutfile = fopen(gapoutFilename,"r");
        }

        //cout<<gapoutFilename<<endl;
        int gap_count=0,gc2=0;

        while(gap_count < totalgaps)
        {
            char *s;
            int gapLength,gapStringLength;

            if(gaps[gap_count].fully_closed == 0)
            {
                //cout<<"Starting gap "<<gap_count<<"\twith "<<gc2<<endl;
                fscanf(gapoutfile,"%d\t%d\t%d\t%d\t%d\t",&gapLength,&gapLength,&gapLength,&gapLength,&gapStringLength);
                //cout<<gapLength<<"\t"<<gapStringLength<<endl;
                if(gapStringLength>0)
                {
                    s = new char[gapStringLength+1];
                    fscanf(gapoutfile,"%s\n",s);
                    int rc= checkComplete(s);

                    if(rc>1)
                    {
                        //cout<<"Error!! Multiple fragment in a single gapstr in process_scf.cpp\n\n";
                        exit(0);
                    }

                    gaps[gap_count] = combine(gaps[gap_count],gapLength,gapStringLength,s,rc,itr);
                    delete[] s;
                }
                else
                {
                    gaps[gap_count] = combine(gaps[gap_count],gapLength,gapStringLength,dummy,0,1);
                    fscanf(gapoutfile,"\n");
                }
                gc2++;
            }
            //else cout<<"Skipping gap "<<gap_count<<endl;
            gap_count++;

        }

        fclose(gapoutfile);

        outFile(gaps,gap_count);///Don't comment
        saveGapInfo(itr-1);
    }

    delete[] temp;
    
    char dir[1000];
    strcpy(dir,gappedpath);
    strcat(dir,"Individual_gaps.txt");

    gapstrfile = fopen(dir,"w");
    
    strcpy(dir,gappedpath);
    strcat(dir,"combined_gapstring.txt");
    
    FILE* input = fopen(dir,"r");

    int gcount=0,f=-1;

    line[0]= '\0';

    fprintf(gapstrfile,"GapNo\tOriginal_Length\tFilled_Length\n\n");
    while(gcount < totalgaps)
    {
        if(gaps[gcount].finalGapLen>0)
        {
            fscanf(input,"%s\n",line);
        }
        else strcpy(line,"");

        fprintf(gapstrfile,"%d\t%d\t%d\t%s\n",gcount,gaps[gcount].originalGap,gaps[gcount].finalGapLen,line);
        //cout<<"Gap = "<<gcount<<endl;
        //cout<<line<<endl;
        gcount++;

    }
    fclose(gapstrfile);
    fclose(input);

    delete[] gaps;

    for(int i=0;i<num_itr;i++)delete[] temp_gaps[i];
    delete[] temp_gaps;
}
