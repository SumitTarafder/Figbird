#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <cstring>

using namespace std;
#include <math.h>

#define MAX_REC_LEN 1024

void reverse(char *reverse, char *read)
{

	char ch='A';
	unsigned long readLength=strlen(read);
	readLength--;
    //cout<<readLength;
    //cout<<read[readLength-1];
    //cout<<read<<endl;
	for(unsigned long i=1; i<=readLength;i++)
	{
		ch=read[i-1];
		if(ch=='A'|| ch=='a')
			reverse[readLength-i]='T';
		else if(ch=='C'|| ch=='c')
			reverse[readLength-i]='G';
		else if(ch=='G'|| ch=='g')
			reverse[readLength-i]='C';
		else if(ch=='T'|| ch=='t')
			reverse[readLength-i]='A';
		else
			reverse[readLength-i]='N';
	}
	reverse[readLength]='\0';

}


int main(int argc, char *argv[])
{

    FILE* fp1 = fopen(argv[1],"r");
    FILE* fp2 = fopen(argv[2],"r");

    if(fp1 == NULL || fp2 == NULL)
    {
        cerr<<"Can't open read pair files during reversing...exiting.\n";
        exit(1);
    }

    string s1(argv[1]);
    string s2(argv[2]);

    size_t found1 = s1.find_last_of(".");
    size_t found2 = s2.find_last_of(".");

    string s3 = s1.substr(0,found1);
    string s4 = s2.substr(0,found2);

    string ext = s1.substr(found1,s3.size());
    string r = "_reversed";

    s3 += r + ext;
    s4 += r + ext;

    FILE* output1= NULL,*output2=NULL;
    
    output1 = fopen(s3.c_str(),"w");
    output2 = fopen(s4.c_str(),"w");

    if(output1 == NULL || output2 == NULL)
    {
        cerr<<"Can't create new read pair files during reversing...exiting.\n";
        exit(1);
    }


    int MAX_FILE_READ = 0;

    char *line1= new char[MAX_REC_LEN];
    char *line2= new char[MAX_REC_LEN];

    MAX_FILE_READ=MAX_REC_LEN/sizeof(line1[0]);

    char *rev = new char[MAX_REC_LEN];

    int count=0;

    while(fgets(line1, MAX_FILE_READ, fp1)!=NULL && fgets(line2, MAX_FILE_READ, fp2)!=NULL)
    {
        if(count% 4 !=1)//Skip headers
        {
            fprintf(output1,"%s",line1);
            fprintf(output2,"%s",line2);
        }
        else
        {
            reverse(rev,line1);
            fprintf(output1,"%s\n",rev);
            reverse(rev,line2);
            fprintf(output2,"%s\n",rev);

        }
        count++;
        //if(count==2)break;
    }

    fclose(output1);
    fclose(output2);
    

    fclose(fp1);
    fclose(fp2);

    cout<<s3<<endl;
    cout<<s4<<endl;
}
