#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define  PICSIZE 256
#define  MAXMASK 100

int    pic[PICSIZE][PICSIZE];
int    histogram[256];
double outpicy[PICSIZE][PICSIZE];
double outpicx[PICSIZE][PICSIZE];
double gradientpic[PICSIZE][PICSIZE];
double candidatepic[PICSIZE][PICSIZE];
double thresholdedpic[PICSIZE][PICSIZE];
int    edgeflag[PICSIZE][PICSIZE];
double xmask[MAXMASK][MAXMASK];
double ymask[MAXMASK][MAXMASK];
double conv[PICSIZE][PICSIZE];
bool   filter;

// input file, sigma value, percent
int main(argc,argv)
int argc;
char **argv;
{
    int i,ii,j,jj,p,q,s,t,mr,centx,centy,LO,HI,target,sum;
    double maskval,xsum,ysum,sig,maxival,minival,maxval,slope,percent,ZEROTOL;
    FILE *fo1, *fo2, *fo3, *fp1, *fopen();
    char *foobar;

    argc--; argv++;
    foobar = *argv;
    fp1=fopen(foobar,"rb");

    fo1=fopen("gradient.pgm","wb");
    fo2=fopen("candidates.pgm","wb");
    fo3=fopen("final.pgm","wb");
    
    fprintf(fo1, "P5\n%d %d\n255\n", PICSIZE, PICSIZE);
    fprintf(fo2, "P5\n%d %d\n255\n", PICSIZE, PICSIZE);
    fprintf(fo3, "P5\n%d %d\n255\n", PICSIZE, PICSIZE);

    argc--; argv++;
	foobar = *argv;
	sig = atof(foobar);

    argc--; argv++;
	foobar = *argv;
	percent = atof(foobar);

    mr = (int)(sig * 3);
    centx = (MAXMASK / 2);
    centy = (MAXMASK / 2);

    for (i=0;i<PICSIZE;i++)
    { 
        for (j=0;j<PICSIZE;j++)
        {
            pic[i][j]  =  getc (fp1);
            pic[i][j]  &= 0377;
        }
    }

    for (p=-mr;p<=mr;p++)
    {  
        for (q=-mr;q<=mr;q++)
        {
            maskval = p*exp(-1*(((p*p)+(q*q))/(2*(sig*sig))));
            (ymask[p+centy][q+centx]) = maskval;

            maskval = q*exp(-1*(((p*p)+(q*q))/(2*(sig*sig))));
            (xmask[p+centy][q+centx]) = maskval;
        }
    }

    for (i=mr;i<=255-mr;i++)
    { 
        for (j=mr;j<=255-mr;j++)
        {
            ysum = 0;
            xsum = 0;
            for (p=-mr;p<=mr;p++)
            {
                for (q=-mr;q<=mr;q++)
                {
                    ysum += pic[i+p][j+q] * ymask[p+centy][q+centx];
                    xsum += pic[i+p][j+q] * xmask[p+centy][q+centx];
                }
            }
            outpicy[i][j] = ysum;
            outpicx[i][j] = xsum;
        }
    }

    maxival = 0;
    for (i=mr;i<PICSIZE-mr;i++)
    { 
        for (j=mr;j<PICSIZE-mr;j++)
        {
            gradientpic[i][j]=sqrt((double)((outpicx[i][j]*outpicx[i][j]) +
                                    (outpicy[i][j]*outpicy[i][j])));
            if (gradientpic[i][j] > maxival)
                maxival = gradientpic[i][j];
        }
    }

    for (i=0;i<PICSIZE;i++)
    { 
        for (j=0;j<PICSIZE;j++)
        {
            gradientpic[i][j] = (gradientpic[i][j] / maxival) * 255;       
            fprintf(fo1,"%c",(char)((int)(gradientpic[i][j])));
            histogram[((int)gradientpic[i][j])%255]++;
        }
    }

    for(i=mr;i<256-mr;i++)
    {
        for(j=mr;j<256-mr;j++)
        {
            if((outpicx[i][j]) == 0.0) 
            {
                outpicx[i][j] = .00001;
            }

            slope = outpicy[i][j]/outpicx[i][j];

            if( (slope <= .4142)&&(slope > -.4142))
            {
                if((gradientpic[i][j] > gradientpic[i][j-1])&&(gradientpic[i][j] > gradientpic[i][j+1]))
                {
                    candidatepic[i][j] = 255;
                }
            }
            else if( (slope <= 2.4142)&&(slope > .4142))
            {
                if((gradientpic[i][j] > gradientpic[i-1][j-1])&&(gradientpic[i][j] > gradientpic[i+1][j+1]))
                {
                    candidatepic[i][j] = 255;
                }
            }
            else if( (slope <= -.4142)&&(slope > -2.4142))
            {
                if((gradientpic[i][j] > gradientpic[i+1][j-1])&&(gradientpic[i][j] > gradientpic[i-1][j+1]))
                {
                    candidatepic[i][j] = 255;
                }
            }
            else
            {
                if((gradientpic[i][j] > gradientpic[i-1][j])&&(gradientpic[i][j] > gradientpic[i+1][j]))
                {
                    candidatepic[i][j] = 255;
                }
            }
        }
    }

    for (i=0;i<PICSIZE;i++)
    { 
        for (j=0;j<PICSIZE;j++)
        {
            fprintf(fo2,"%c",(char)((int)(candidatepic[i][j])));
        }
    }
    
    LO = HI = sum = 0;
    target = (int)(PICSIZE*PICSIZE*percent);
    for (i=255;i>0;--i)
    { 
        sum += histogram[i];
        if(sum >= target)
        {
            HI = i;
            LO = (int)(0.35*i);
            break;
        }
    }

    for (i=0;i<PICSIZE;i++)
    { 
        for (j=0;j<PICSIZE;j++)
        {
            if(candidatepic[i][j] == 255)
            {
                if(gradientpic[i][j] > HI)
                {
                    thresholdedpic[i][j] = 255;
                    candidatepic[i][i] = 0;
                }
                else if(gradientpic[i][i] < LO)
                {
                    thresholdedpic[i][i] = 0;
                    candidatepic[i][i] = 0;
                }
            }
        }
    }

    filter = true;
    while(filter)
    {
        filter = false;
        for (i=0;i<PICSIZE;i++)
        { 
            for (j=0;j<PICSIZE;j++)
            {
                if(candidatepic[i][j] == 255)
                {
                    for (p=-1;p<=1;p++)
                    { 
                        for (q=-1;q<=1;q++)
                        {
                            ii = i+p;
                            jj = j+q;
                            if(ii<0 || jj<0 || ii==PICSIZE || jj == PICSIZE)
                                continue;
                            if(thresholdedpic[i+p][j+q] == 255)
                            {
                                candidatepic[i][j] = 0;
                                thresholdedpic[i][j] = 255;
                                filter = true;
                            }
                        }
                    }
                }
            }
        }
    }

    for (i=0;i<PICSIZE;i++)
    { 
        for (j=0;j<PICSIZE;j++)
        {
            fprintf(fo3,"%c",(char)((int)(thresholdedpic[i][j])));
        }
    }

    fclose(fo1);
    fclose(fo2);
    fclose(fo3);
    return 0;
}
