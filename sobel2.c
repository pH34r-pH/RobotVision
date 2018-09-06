#include <stdio.h>                          /* Sobel.c */
#include <stdlib.h>
#include <math.h>

    int pic[256][256];
    int outpicx[256][256];
    int outpicy[256][256];
    int maskx[3][3] = {{-1,0,1},{-2,0,2},{-1,0,1}};
    int masky[3][3] = {{1,2,1},{0,0,0},{-1,-2,-1}};
    double ival[256][256],maxival;
    double low[256][256],high[256][256],lowThresh,highThresh;

int main(argc,argv)
int argc;
char **argv;
{
    int i,j,p,q,mr,sum1,sum2;
    double threshold;
    FILE *fo1, *fo2, *fo3, *fp1, *fopen();
    char *foobar;

    argc--; argv++;
    foobar = *argv;
    fp1=fopen(foobar,"rb");

    argc--; argv++;
    foobar = *argv;
    fo1=fopen(foobar,"wb");

    fo2=fopen("low.pgm","wb");
    fo3=fopen("high.pgm","wb");
    
    argc--; argv++;
	foobar = *argv;
	threshold = atof(foobar);

    for (i=0;i<256;i++)
    { 
        for (j=0;j<256;j++)
        {
            pic[i][j]  =  getc (fp1);
            pic[i][j]  &= 0377;
        }
    }

    mr = 1;
    for (i=mr;i<256-mr;i++)
    { 
        for (j=mr;j<256-mr;j++)
        {
            sum1 = 0;
            sum2 = 0;
            for (p=-mr;p<=mr;p++)
            {
                for (q=-mr;q<=mr;q++)
                {
                    sum1 += pic[i+p][j+q] * maskx[p+mr][q+mr];
                    sum2 += pic[i+p][j+q] * masky[p+mr][q+mr];
                }
            }
            outpicx[i][j] = sum1;
            outpicy[i][j] = sum2;
        }
    }

    maxival = 0;
    lowThresh = 75;
    highThresh = 150;
    for (i=mr;i<256-mr;i++)
    { for (j=mr;j<256-mr;j++)
        {
            ival[i][j]=sqrt((double)((outpicx[i][j]*outpicx[i][j]) +
                                    (outpicy[i][j]*outpicy[i][j])));
            if (ival[i][j] > maxival)
                maxival = ival[i][j];

            low[i][j] = ival[i][j] > lowThresh ? 255 : 0;
            high[i][j] = ival[i][j] > highThresh ? 255 : 0;
        }
    }

    fprintf(fo1, "P5\n");
    fprintf(fo1, "%d %d\n", 256, 256);
    fprintf(fo1, "255\n");
    fprintf(fo2, "P5\n");
    fprintf(fo2, "%d %d\n", 256, 256);
    fprintf(fo2, "255\n");
    fprintf(fo3, "P5\n");
    fprintf(fo3, "%d %d\n", 256, 256);
    fprintf(fo3, "255\n");

    for (i=0;i<256;i++)
    { 
        for (j=0;j<256;j++)
        {
            ival[i][j] = (ival[i][j] / maxival) * 255;            
            fprintf(fo1,"%c",(char)((int)(ival[i][j])));
            fprintf(fo2,"%c",(char)((int)(low[i][j])));
            fprintf(fo3,"%c",(char)((int)(high[i][j])));
        }
    }
    
    fclose(fo1);
    fclose(fo2);
    fclose(fo3);
    return 0;
}
