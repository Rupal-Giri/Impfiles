
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>

#define GREEN   "\x1b[32m"
#define RESET   "\x1b[0m"

typedef struct{
  char source[15];
  double mjd;
  double year;
  double frq;
  double elv;
  double flx;
  double err;
  char telescope[15];
} dataPoint;

int dim_file(FILE *);

int compare (const void * a, const void * b)
{
  return ( *(double*)a - *(double*)b );
}

int main(int argc, char **argv)
{
  int i, j, k, iFrq, nDat, nDatPerFrq, nFrq, nFrqOut, nFile;
  int dumInt, srcFlag=0, frqFlag, itrFlag=0;
  double availableFrq[20], mjdMin=-1., mjdMax=100000., yrMin, yrMax, frqMin=-1., frqMax=100000.;
  dataPoint *datum;
  char line[500], source[15], source1[15], band;
  FILE *fpIn, *fpOut;

  // Checking for options
  for(i=0; i<argc; i++)
    {
      if(!strcmp(argv[i], "-src"))
	{
	  srcFlag=1;
	  strcpy(source, argv[i+1]);
	}

      if(!strcmp(argv[i], "-frq"))
	{
	    frqMin=atof(argv[i+1]);
	    frqMax=atof(argv[i+2]);
	}
      if(!strcmp(argv[i], "-bandC"))
	{
	    frqMin=4000.;
	    frqMax=7000.;
	}
      if(!strcmp(argv[i], "-bandX"))
	{
	    frqMin=7500.;
	    frqMax=9000.;
	}
      if(!strcmp(argv[i], "-bandK"))
	{
	    frqMin=19000.;
	    frqMax=25000.;
	}
      if(!strcmp(argv[i], "-bandV"))
	{
	    frqMin=39000.;
	    frqMax=45000.;
	}
      if(!strcmp(argv[i], "-time"))
	{
	    mjdMin=atof(argv[i+1]);
	    mjdMax=atof(argv[i+2]);
	}
      if(!strcmp(argv[i], "-yr"))
	{
	    yrMin=atof(argv[i+1]);
	    yrMax=atof(argv[i+2]);
	    mjdMin=(yrMin-2005.)*365.25+53371.;
	    mjdMax=(yrMax-2005.)*365.25+53371.;
	}
      if(!strcmp(argv[i], "-iterative"))
	{
	  itrFlag=1;
	}

    }
  if(srcFlag==0)
    {
      printf("Please type the name of the source you are interested in: ");
      scanf("%s", (char *)&source);
    }

  // Aliases!!!
  if(!strcmp(source, "0219+428"))
    strcpy(source1, "3C66A");
  if(!strcmp(source, "3C66A"))
    strcpy(source1, "0219+428");

  if(!strcmp(source, "0316+41"))
    strcpy(source1, "3C84");
  if(!strcmp(source, "3C84"))
    strcpy(source1, "0316+41");

  if(!strcmp(source, "0355+508"))
    strcpy(source1, "NRAO150");
  if(!strcmp(source, "NRAO150"))
    strcpy(source1, "0355+508");

  if(!strcmp(source, "0851+202"))
    strcpy(source1, "OJ287");
  if(!strcmp(source, "OJ287"))
    strcpy(source1, "0851+202");

  if(!strcmp(source, "1226+023"))
    strcpy(source1, "3C273");
  if(!strcmp(source, "3C273"))
    strcpy(source1, "1226+023");

  if(!strcmp(source, "1253-055"))
    strcpy(source1, "3C279");
  if(!strcmp(source, "3C279"))
    strcpy(source1, "1253-055");

  if(!strcmp(source, "1641+399"))
    strcpy(source1, "3C345");
  if(!strcmp(source, "3C345"))
    strcpy(source1, "1641+399");
  
  if(!strcmp(source, "2200+420"))
    strcpy(source1, "BLLAC");
  if(!strcmp(source, "BLLAC"))
    strcpy(source1, "2200+420");

  if(!strcmp(source, "2251+158"))
    strcpy(source1, "3C454.3");
  if(!strcmp(source, "3C454.3"))
    strcpy(source1, "2251+158");

  if(!strcmp(source, "Mrk501"))
    strcpy(source1, "MRK501");
   if(!strcmp(source, "MRK501"))
    strcpy(source1, "Mrk501");
 
  if(!strcmp(source, "3c207"))
    strcpy(source1, "3C207");
   if(!strcmp(source, "3C207"))
    strcpy(source1, "3c207");

   if(!strcmp(source, "BS0846+513"))
    strcpy(source1, "0846+513");
   if(!strcmp(source, "0846+513"))
    strcpy(source1, "BS0846+513");

   if(!strcmp(source, "CTA102"))
    strcpy(source1, "2230+114");
   if(!strcmp(source, "2230+114"))
    strcpy(source1, "CTA102");


   //// The name 1221+21 seems wrong! Is it indeed PKS1222+216? In case, we should correct the Unicum file.
   if(!strcmp(source, "PKS1222+216"))
    strcpy(source1, "1221+21");
   if(!strcmp(source, "1221+21"))
    strcpy(source1, "PKS1222+216");

  fpIn=fopen("Unicum_Bach_pre-CAP_measurements.txt", "r");
  if(fpIn==NULL)
    {
      printf("Error reading Unicum_Bach_pre-CAP_measurements.txt\n");
      exit(1);
    }
  nDat=dim_file(fpIn);
  //printf("%d\n", nDat);
  fclose(fpIn);

  fpIn=fopen("Unicum_last.dat", "r");
  if(fpIn==NULL)
    {
      printf("Error reading Unicum_last.dat\n");
      exit(1);
    }
  nDat=nDat+dim_file(fpIn);
  //printf("%d\n", nDat);
 
  fclose(fpIn);
  datum=(dataPoint *)malloc(nDat*sizeof(dataPoint));

  i=0;
  for(nFile=0; nFile<2; nFile++)
    {
      if(nFile==0)
	{
	  fpIn=fopen("Unicum_Bach_pre-CAP_measurements.txt", "r");
	  printf("#\n# Reading Unicum_Bach_pre-CAP_measurements.txt\n");
	}
      if(nFile==1)
	{
	  fpIn=fopen("Unicum_last.dat", "r");
	  printf("#\n# Reading Unicum_last.dat\n#\n");
	}

      while (!feof(fpIn))
	{
	  fgets(line, 600, fpIn);
	  sscanf(line, "%s", (char *)&datum[i].source);
	  if(strcmp(datum[i].source, "FLAGGED") && strcmp(datum[i].source, "Name") && strcmp(datum[i].source, "MHz") && strcmp(datum[i].source, "-99.000"))
	    {
	      if(!strcmp(datum[i].source, source) || !strcmp(datum[i].source, source1))
		{
		  strcpy(datum[i].telescope, "----");
		  sscanf(line, "%s%d%lf%lf%lf%lf%lf%s", 
			 (char *)&datum[i].source, &dumInt, &datum[i].mjd, &datum[i].frq, &datum[i].elv, &datum[i].flx, &datum[i].err, (char *)&datum[i].telescope);
		  
		  datum[i].year=(datum[i].mjd-53371)/365.25+2005.;
		  
		  i=i+1;
		}
	    }
	}
      fclose(fpIn);
    }

  nDat=i;
  availableFrq[0]=datum[0].frq;
  k=1;
  for(i=1; i<nDat-2; i++)
    {
      frqFlag=0;
      for(j=0; j<k; j++)
	{
	  if(fabs(datum[i].frq-availableFrq[j])<DBL_MIN)
	    frqFlag=1;
	}
      if(frqFlag==0)
	{
	  availableFrq[k]=datum[i].frq;
	  k=k+1;
	}
    }
  nFrq=k;
  qsort (availableFrq, nFrq, sizeof(double), compare);
  printf("# List of available frequencies (in MHz):\n");
  for(k=0; k<nFrq; k++)
    {
      printf("#  %8.1f\n", availableFrq[k]);
    }
  printf("# Time interval (in MJD):\n#  %.5f - %.5f\n", datum[0].mjd, datum[nDat-1].mjd);

  if(itrFlag==0)
    nFrqOut=1;
  else
    nFrqOut=4;

  for(iFrq=0; iFrq<nFrqOut; iFrq++)
    {

      if(iFrq==0 && itrFlag==1)
	{
	  frqMin=4000.;
	  frqMax=7000.;
	  band='C';
	}
      if(iFrq==1 && itrFlag==1)
	{
	  frqMin=7500.;
	  frqMax=9000.;
	  band='X';
	}
      if(iFrq==2 && itrFlag==1)
	{
	  frqMin=19000.;
	  frqMax=25000.;
	  band='K';
	}
      if(iFrq==3 && itrFlag==1)
	{
	  frqMin=39000.;
	  frqMax=45000.;
	  band='V';
	}

      fpOut=fopen("readData.dat", "w+");
      if(fpOut==NULL)
	{
	  printf("Error writing readData.dat\n");
	  exit(1);
	}
      
      printf("#\n# Data format:\n");
      printf("#   Src        MJD          year         Frq      Elv       Flx        Ferr    Tsc\n");
      printf("#              (d)           (y)        (MHz)    (deg)     (Jy)        (Jy)\n");
      fprintf(fpOut, "#   Src        MJD          year         Frq      Elv       Flx        Ferr    Tsc\n");
      fprintf(fpOut, "#              (d)           (y)        (MHz)    (deg)     (Jy)        (Jy)\n");
      
      nDatPerFrq=0;
      for(k=0; k<nFrq; k++)
	{
	  for(i=0; i<nDat; i++)
	    {
	      if(fabs(datum[i].frq-availableFrq[k])<DBL_MIN)
		{
		  if(datum[i].frq>frqMin-0.01 && datum[i].frq<frqMax+0.01 && datum[i].mjd>mjdMin-0.0001 && datum[i].mjd<mjdMax+0.0001)
		    {
		      printf("%8s   %.5f   %.5f   %8.2f   %5.2f   %8.5f   %8.5f   %s\n", 
			     source, datum[i].mjd, datum[i].year, datum[i].frq, datum[i].elv, datum[i].flx, datum[i].err, datum[i].telescope);
		      fprintf(fpOut, "%8s   %.5f   %.5f   %8.2f   %5.2f   %8.5f   %8.5f   %s\n", 
			      source, datum[i].mjd, datum[i].year, datum[i].frq, datum[i].elv, datum[i].flx, datum[i].err, datum[i].telescope);
		      nDatPerFrq=nDatPerFrq+1;
		    }
		}
	    }
	}      
      fclose(fpOut);
      if(itrFlag==1 && nDatPerFrq>0)
	{
	  fpOut=fopen("srcTmp.csh", "w+");
	  if(fpOut==NULL)
	    {
	      printf("Error writing temporary script\n");
	      exit(1);
	    }
	  fprintf(fpOut, "cp readData.dat Sources/%s_%c.dat\n", source, band);
	  fclose(fpOut);
	  system("source srcTmp.csh");
	}
    }

  printf("#\n# If you wish to select some frequencies and/or a time range for your data, please use one of the following commands:\n");
  printf("# >" GREEN "./readData -src %s -frq minFrequency maxFrequency -time minMjd maxMjd\n", datum[i-1].source);
  printf(RESET "# >" GREEN "./readData -src %s -frq minFrequency maxFrequency -yr minYear maxYear" RESET "\n#\n", datum[i-1].source);


  return 1;
}

int dim_file(FILE *fpin)
 {
   int n;
   char line[500];

   n=0;
   while (!feof(fpin)) {
      fgets(line, 300, fpin);
      //if(line[0]!='#' && line[0]!='\n')
	{
	  n=n+1;
	}
  }

  return n-1;
}
