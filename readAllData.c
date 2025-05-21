
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
  char band[6];
  char telescope[15];
} dataPoint;

typedef struct{
  char name[20];
  char nameStd[20];
  char type[20];
  char cat[20];
  char catSrc[20];
  char ra[20];
  char dec[20];
  double z;
} sourceObj;

typedef struct{
  double mjd;
  double mjdMin;
  double mjdMax;
  double bC;
  double errC;
  double frqC;
  int nC;
  double bX;
  double errX;
  double frqX;
  int nX;
  double bK;
  double errK;
  double frqK;
  int nK;
  double bV;
  double errV;
  double frqV;
  int nV;
} epoch;

typedef struct{
  int ndat;
  double mind;
  double modsf3;
  double aversf3;
  double modsf6;
  double aversf6;
} result;

int dim_file(FILE *);
result mind(int, epoch*, char);

int compare (const void * a, const void * b)
{
  return ( *(double*)a - *(double*)b );
}

int main(int argc, char **argv)
{
  int i, j, k, iSrc, itSrc, iFrq, iEpc, nDat, nDatPerFrq, nFrqPerEpc, nFrq, nFrqOut, nFile, nSrc, nEpc, nEpcTot, nSI, nFlxTot;
  int dumInt, srcFlag=0, frqFlag, itrFlag=1, epcFlag;
  int nC, nX, nK, nV;
  double availableFrq[20], mjdMin=-1., mjdMax=100000., yrMin, yrMax, frqMin=-1., frqMax=100000.;
  double xmedia, ymedia, Ar, Br, Cr, Dr, axr[4], ayr[4], a, b, c, r, delta, coeffA, coeffB, betaSF;
  double SI, SIavr, SImin, SImax, avrC, avrX, avrK, avrV, sf3_6Cpar, sf3_6Xpar, sf3_6Kpar, sf3_6Vpar;
  dataPoint *datum;
  char line[500], source[15], source1[15], sourceO[15], band, dumC[3], srcType[10], srcTmp[20];
  sourceObj *src;
  epoch *epc;
  result resC, resX, resK, resV;
  FILE *fpIn, *fpOut, *fpList, *fpEpc, *fpChar, *fpTab1, *fpTab2;
  FILE *fpTmp1, *fpTmp2, *fpTmp3, *fpTmp4;

  strcpy(srcType, "------");
  
  // Checking for options
  for(i=0; i<argc; i++)
    {
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
      if(!strcmp(argv[i], "-type"))
	{
	  strcpy(srcType, argv[i+1]);
	}
      if(!strcmp(argv[i], "-source"))
	{
	  strcpy(srcTmp, argv[i+1]);
	}
    }

  system("rm epochs.dat");    
  fpList = fopen("allSources_list_new.txt", "r");
  if(fpList==NULL)
    {
      printf("Error reading allSources_list_new.txt\n");
      exit(1);
    }
  nSrc=dim_file(fpList);
  src=(sourceObj *)malloc(nSrc*sizeof(sourceObj));
  rewind(fpList);

  iSrc=0;
  for(itSrc=0; itSrc<nSrc; )
    {
      fgets(line, 300, fpList);
      if(line[0]!='#' && line[0]!='\n')
      //if(line[0]=='#')
	{
	  sscanf(line, "%s%s%s%s%s%s%s%s%s%s%s%s%lf", 
		 (char *)&src[iSrc].name, (char *)&dumC, (char *)&src[iSrc].nameStd, (char *)&dumC, (char *)&src[iSrc].type, (char *)&dumC, (char *)&src[iSrc].cat, (char *)&src[iSrc].catSrc,\
		 (char *)&dumC, (char *)&src[iSrc].ra, (char *)&src[iSrc].dec, (char *)&dumC, &src[iSrc].z);
	  //printf("%13s %13s    %10s    %5s %10s    %11s %12s\n", src[iSrc].name, src[iSrc].nameStd, src[iSrc].type, src[iSrc].cat, src[iSrc].catSrc, src[iSrc].ra, src[iSrc].dec);  
	  if((!strcmp(srcType, "bll") && !strcmp(src[iSrc].type, "BLL")) || (!strcmp(srcType, "fsrq") && !strcmp(src[iSrc].type, "FSRQ")) || !strcmp(srcType, "------"))
	    iSrc=iSrc+1;

	  itSrc=itSrc+1;	  
	}
    }
  nSrc=iSrc;
  
  fpChar = fopen("sources_char.dat", "w+");
  if(fpChar==NULL)
    {
      printf("Error writing sources_char.dat\n");
      exit(1);
    }
  fpTab1 = fopen("table1.dat", "w+");
  if(fpTab1==NULL)
    {
      printf("Error writing table1.dat\n");
      exit(1);
    }
  fpTab2 = fopen("table2.dat", "w+");
  if(fpTab2==NULL)
    {
      printf("Error writing table2.dat\n");
      exit(1);
    }
  
  nEpcTot=0;
  nFlxTot=0;
  for(iSrc=0; iSrc<nSrc; iSrc++)
    {
      SIavr=0.;
      SImax=-DBL_MAX;
      SImin=DBL_MAX;
      nSI=0;

      strcpy(source, src[iSrc].name);
      strcpy(source1, "-----");
      //printf("\n Source: %s\n", source);
      
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

      fpEpc=fopen("epochsTmp.dat", "w+");
      if(fpEpc==NULL)
	{
	  printf("Error writing epochsTmp.dat\n");
	  exit(1);
	}
            fpTmp1=fopen("bandC.dat", "w+");
      if(fpTmp1==NULL)
	{
	  printf("Error writing bandC.dat\n");
	  exit(1);
	}
      fpTmp2=fopen("bandX.dat", "w+");
      if(fpTmp2==NULL)
	{
	  printf("Error writing bandX.dat\n");
	  exit(1);
	}
      fpTmp3=fopen("bandK.dat", "w+");
      if(fpTmp3==NULL)
	{
	  printf("Error writing bandK.dat\n");
	  exit(1);
	}
      fpTmp4=fopen("bandV.dat", "w+");
      if(fpTmp4==NULL)
	{
	  printf("Error writing bandV.dat\n");
	  exit(1);
	}
      
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
	      //printf("#\n# Reading Unicum_Bach_pre-CAP_measurements.txt\n");
	    }
	  if(nFile==1)
	    {
	      fpIn=fopen("Unicum_last.dat", "r");
	      //printf("#\n# Reading data.dat\n#\n");
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
		      
		      if(datum[i].frq>4000 && datum[i].frq<7000)
			strcpy(datum[i].band, "bandC");
		      if(datum[i].frq>7500 && datum[i].frq<9000)
			strcpy(datum[i].band, "bandX");
		      if(datum[i].frq>19000 && datum[i].frq<25000)
			strcpy(datum[i].band, "bandK");
		      if(datum[i].frq>39000 && datum[i].frq<45000)
			strcpy(datum[i].band, "bandV");

		      i=i+1;
		      nFlxTot=nFlxTot+1;
		    }
		}
	    }
	  fclose(fpIn);
	}
      nDat=i;

      /*
      for(i=1; i<nDat; i++)
	{
	  if(!strcmp(source, "1219+285"))
	    printf("...%s   %f   %f\n", datum[i].source, datum[i].mjd, datum[i].frq);
	}
      */
      
      epc=(epoch *)malloc(nDat*sizeof(epoch));
      
      availableFrq[0]=datum[0].frq;
      epc[0].mjd=datum[0].mjd;
      epc[0].mjdMin=datum[0].mjd-DBL_MIN;
      epc[0].mjdMax=datum[0].mjd+DBL_MIN;
      iEpc=0;
      nEpc=1;
      k=1;
      for(i=1; i<nDat; i++)
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

	  epcFlag=0;
	  //if(!strcmp(source, "TXS0536+145"))
	  //printf("%f\n", datum[i].mjd);
	  for(iEpc=0; iEpc<nEpc; iEpc++)
	    {	      
	      if(fabs(datum[i].mjd-epc[iEpc].mjd)<7.)
		{
		  if(datum[i].mjd<epc[iEpc].mjdMin)
		    epc[iEpc].mjdMin=datum[i].mjd;
		  if(datum[i].mjd>epc[iEpc].mjdMax)
		    epc[iEpc].mjdMax=datum[i].mjd;
		  epcFlag=1;
		  //if(!strcmp(source, "1219+285"))
		  //printf("%f   %f   %f\n", datum[i].mjd, epc[iEpc].mjdMin, epc[iEpc].mjdMax);
		  //if(!strcmp(datum[i].source, "1219+285") && !strcmp(datum[i].band, "bandC"))
		  //	printf("%s   %f   %f\n", datum[i].source, datum[i].mjd, datum[i].flx);
		      				  
		  break;
		}
	    }
	  if(epcFlag==0)
	    {
	      epc[nEpc].mjd=datum[i].mjd;
	      epc[nEpc].mjdMin=datum[i].mjd-DBL_MIN;
	      epc[nEpc].mjdMax=datum[i].mjd+DBL_MIN;
	      //if(!strcmp(source, "1219+285"))
	      //printf("...%f   %f   %f\n", datum[i].mjd, epc[nEpc].mjdMin, epc[nEpc].mjdMax);

	      nEpc=nEpc+1;
	    }
	}
      nFrq=k;
      qsort (availableFrq, nFrq, sizeof(double), compare);
      //printf("# List of available frequencies (in MHz):\n");
      for(k=0; k<nFrq; k++)
	{
	  //printf("#  %8.1f\n", availableFrq[k]);
	}
      //printf("# Time interval (in MJD):\n#  %.5f - %.5f\n", datum[0].mjd, datum[nDat-1].mjd);
      
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
	  
	  //printf("#\n# Data format:\n");
	  //printf("#   Src        MJD          year         Frq      Elv       Flx        Ferr    Tsc\n");
	  //printf("#              (d)           (y)        (MHz)    (deg)     (Jy)        (Jy)\n");
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
			  //printf("%8s   %.5f   %.5f   %8.2f   %5.2f   %8.5f   %8.5f   %s\n", 
			  //	 source, datum[i].mjd, datum[i].year, datum[i].frq, datum[i].elv, datum[i].flx, datum[i].err, datum[i].telescope);
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
      
      //printf("#\n# If you wish to select some frequencies and/or a time range for your data, please use one of the following commands:\n");
      //printf("# >" GREEN "./readData -src %s -frq minFrequency maxFrequency -time minMjd maxMjd\n", datum[i-1].source);
      //printf(RESET "# >" GREEN "./readData -src %s -frq minFrequency maxFrequency -yr minYear maxYear" RESET "\n#\n", datum[i-1].source);
      
      avrC=0.;
      avrX=0.;
      avrK=0.;
      avrV=0.;
      nC=0;
      nX=0;
      nK=0;
      nV=0;
      for(iEpc=0; iEpc<nEpc; iEpc++)
	{
	  epc[iEpc].bC=0.;
	  epc[iEpc].bX=0.;
	  epc[iEpc].bK=0.;
	  epc[iEpc].bV=0.;
	  epc[iEpc].errC=0.;
	  epc[iEpc].errX=0.;
	  epc[iEpc].errK=0.;
	  epc[iEpc].errV=0.;
	  epc[iEpc].nC=0;
	  epc[iEpc].nX=0;
	  epc[iEpc].nK=0;
	  epc[iEpc].nV=0;
	  for(i=0; i<nDat; i++)
	    {
	      //printf("%f   %f   %f\n", datum[i].mjd, epc[iEpc].mjdMin, epc[iEpc].mjdMax);
	      //if(i>20)
	      //exit(1);
	      if(datum[i].mjd>epc[iEpc].mjdMin-0.00001 && datum[i].mjd<epc[iEpc].mjdMax+0.00001)
		{
		  if(!strcmp(datum[i].band, "bandC"))
		    {
		      epc[iEpc].bC=epc[iEpc].bC+datum[i].flx;
		      epc[iEpc].errC=epc[iEpc].errC+datum[i].err;
		      epc[iEpc].frqC=epc[iEpc].frqC+datum[i].frq;
		      epc[iEpc].nC=epc[iEpc].nC+1;
		    }
		  if(!strcmp(datum[i].band, "bandX"))
		    {
		      epc[iEpc].bX=epc[iEpc].bX+datum[i].flx;
		      epc[iEpc].errX=epc[iEpc].errX+datum[i].err;
		      epc[iEpc].frqX=epc[iEpc].frqX+datum[i].frq;

		      epc[iEpc].nX=epc[iEpc].nX+1;
		    }
		  if(!strcmp(datum[i].band, "bandK"))
		    {
		      epc[iEpc].bK=epc[iEpc].bK+datum[i].flx;
		      epc[iEpc].errK=epc[iEpc].errK+datum[i].err;
		      epc[iEpc].frqK=epc[iEpc].frqK+datum[i].frq;
		      epc[iEpc].nK=epc[iEpc].nK+1;
		    }
		  if(!strcmp(datum[i].band, "bandV"))
		    {
		      epc[iEpc].bV=epc[iEpc].bV+datum[i].flx;
		      epc[iEpc].errV=epc[iEpc].errV+datum[i].err;
		      epc[iEpc].frqV=epc[iEpc].frqV+datum[i].frq;
		      epc[iEpc].nV=epc[iEpc].nV+1;
		    }		  
		}
	    }
	  nFrqPerEpc=0;
	  iFrq=0;
	  if(epc[iEpc].nC>0)
	    {
	      epc[iEpc].bC=epc[iEpc].bC/(double)epc[iEpc].nC;
	      epc[iEpc].errC=epc[iEpc].errC/(double)epc[iEpc].nC;
	      epc[iEpc].frqC=epc[iEpc].frqC/(double)epc[iEpc].nC;
	      avrC=avrC+epc[iEpc].bC;
	      nC=nC+1;
	      
	      // !!! Removing C and V band from the calculation of the spectral index, for consistency reasons !!!
	      //axr[iFrq]=log(epc[iEpc].frqC);
	      //ayr[iFrq]=log(epc[iEpc].bC);
	      //iFrq=iFrq+1;
	      //nFrqPerEpc=nFrqPerEpc+1;
	    }
	  else
	    epc[iEpc].bC=-1.;

	  if(epc[iEpc].nX>0)
	    {
	      epc[iEpc].bX=epc[iEpc].bX/(double)epc[iEpc].nX;
	      epc[iEpc].errX=epc[iEpc].errX/(double)epc[iEpc].nX;
	      epc[iEpc].frqX=epc[iEpc].frqX/(double)epc[iEpc].nX;
	      avrX=avrX+epc[iEpc].bX;
	      nX=nX+1;

	      axr[iFrq]=log(epc[iEpc].frqX);
	      ayr[iFrq]=log(epc[iEpc].bX);
	      iFrq=iFrq+1;
	      nFrqPerEpc=nFrqPerEpc+1;
	    }
	  else
	    epc[iEpc].bX=-1.;
	  if(epc[iEpc].nK>0)
	    {
	      epc[iEpc].bK=epc[iEpc].bK/(double)epc[iEpc].nK;
	      epc[iEpc].errK=epc[iEpc].errK/(double)epc[iEpc].nK;
	      epc[iEpc].frqK=epc[iEpc].frqK/(double)epc[iEpc].nK;
	      avrK=avrK+epc[iEpc].bK;
	      nK=nK+1;

	      axr[iFrq]=log(epc[iEpc].frqK);
	      ayr[iFrq]=log(epc[iEpc].bK);
	      iFrq=iFrq+1;
	      nFrqPerEpc=nFrqPerEpc+1;
	    }
	  else
	    epc[iEpc].bK=-1.;
	  if(epc[iEpc].nV>0)
	    {
	      epc[iEpc].bV=epc[iEpc].bV/(double)epc[iEpc].nV;
	      epc[iEpc].errV=epc[iEpc].errV/(double)epc[iEpc].nV;
	      epc[iEpc].frqV=epc[iEpc].frqV/(double)epc[iEpc].nV;
	      avrV=avrV+epc[iEpc].bV;
	      nV=nV+1;

	      // !!! Removing C and V band from the calculation of the spectral index, for consistency reasons !!!
	      //axr[iFrq]=log(epc[iEpc].frqV);
	      //ayr[iFrq]=log(epc[iEpc].bV);
	      //iFrq=iFrq+1;
	      //nFrqPerEpc=nFrqPerEpc+1;
	    }
	  else
	    epc[iEpc].bV=-1.;

	  /* Regression */
	  if(nFrqPerEpc>1)
	    {
	      //printf("%f %f %f %f %f %f %f %f\n", epc[iEpc].frqC, epc[iEpc].bC, epc[iEpc].frqX, epc[iEpc].bX, epc[iEpc].frqK, epc[iEpc].bK, epc[iEpc].frqV, epc[iEpc].bV);
	      xmedia=0.;
	      ymedia=0.;
	      Ar=0.;
	      Br=0.;
	      Cr=0.;
	      Dr=0;
  
	      for(i=0; i<nFrqPerEpc; i++)
		{
		  xmedia=xmedia+axr[i];
		  ymedia=ymedia+ayr[i];
		}
	      xmedia=xmedia/(double)nFrqPerEpc;
	      ymedia=ymedia/(double)nFrqPerEpc;
	      for(i=0; i<nFrqPerEpc; i++)
		{
		  a=(axr[i]-xmedia)*(ayr[i]-ymedia);
		  b=pow((axr[i]-xmedia), 2.);
		  c=pow((ayr[i]-ymedia), 2.);
		  Ar=Ar+a;
		  Br=Br+b;
		  Cr=Cr+c;
		}
	      r=Ar/(pow((Br*Cr), 0.5));
	      if(r<0)
		r=-r;
	      
	      Ar=0.;
	      Br=0.;
	      Cr=0.;
	      
	      for (i=0; i<nFrqPerEpc; i++)
		{
		  a=pow(axr[i], 2.);
		  c=axr[i]*ayr[i];
		  Ar=Ar+a;
		  Br=Br+ayr[i];
		  Cr=Cr+c;
		  Dr=Dr+axr[i];
		}
	      
	      delta=nFrqPerEpc*Ar-pow(Dr, 2.);
	      ///////////////////////////
	      coeffA=(Ar*Br-Dr*Cr)/delta;
	      coeffB=(nFrqPerEpc*Cr-Br*Dr)/delta;
	      ///////////////////////////
	      
	      /* End Regression */
	    }
	  if(nFrqPerEpc<2)
	    coeffB=-9.99;

	  SI=coeffB;
	  if(SI>-9.99)
	    {
	      SIavr=SIavr+SI;
	      if(SImin>SI)
		SImin=SI;
	      if(SImax<SI)
		SImax=SI;

	      nSI=nSI+1;
	    }

	  fprintf(fpEpc, "%14s   %10.4f   %10.4f   %10.4f   %8.4f  %7.4f  %2d    %8.4f  %7.4f  %2d    %8.4f  %7.4f  %2d    %8.4f  %7.4f  %2d   %2d %+4f\n",
		  source, (epc[iEpc].mjdMin+epc[iEpc].mjdMax)/2., epc[iEpc].mjdMin, epc[iEpc].mjdMax, epc[iEpc].bC, \
		  epc[iEpc].errC, epc[iEpc].nC, epc[iEpc].bX, epc[iEpc].errX, epc[iEpc].nX, epc[iEpc].bK, epc[iEpc].errK, \
		  epc[iEpc].nK, epc[iEpc].bV, epc[iEpc].errV, epc[iEpc].nV, nFrqPerEpc, SI);

	  if(!strcmp(source, srcTmp))
	    {
	      if(epc[iEpc].bC>0)
		{
		  fprintf(fpTmp1, "%10.4f   %8.4f  %7.4f\n",
			  (epc[iEpc].mjdMin+epc[iEpc].mjdMax)/2., epc[iEpc].bC, epc[iEpc].errC);
		}
	      if(epc[iEpc].bX>0)
		{
		  fprintf(fpTmp2, "%10.4f   %8.4f  %7.4f\n",
			  (epc[iEpc].mjdMin+epc[iEpc].mjdMax)/2., epc[iEpc].bX, epc[iEpc].errX);
		}
	      if(epc[iEpc].bK>0)
		{
		  fprintf(fpTmp3, "%10.4f   %8.4f  %7.4f\n",
			  (epc[iEpc].mjdMin+epc[iEpc].mjdMax)/2., epc[iEpc].bK, epc[iEpc].errK);
		}
	      if(epc[iEpc].bV>0)
		{
		  fprintf(fpTmp4, "%10.4f   %8.4f  %7.4f\n",
			  (epc[iEpc].mjdMin+epc[iEpc].mjdMax)/2., epc[iEpc].bV, epc[iEpc].errV);
		}
	    }
	}
      fclose(fpEpc);
      fclose(fpTmp1);
      fclose(fpTmp2);
      fclose(fpTmp3);
      fclose(fpTmp4);
      if(!strcmp(source, srcTmp))
	exit(1);
      
      system("sort -k2 epochsTmp.dat >> epochs.dat");

      if(nC>0)
	avrC=avrC/(double)nC;
      else
	avrC=-9.99;
      if(nX>0)
	avrX=avrX/(double)nX;
       else
	avrX=-9.99;
      if(nK>0)
	avrK=avrK/(double)nK;
       else
	avrK=-9.99;
      if(nV>0)
	avrV=avrV/(double)nV;
       else
	avrV=-9.99;

      ///*
      resC=mind(nEpc, epc, 'C');
      resX=mind(nEpc, epc, 'X');
      resK=mind(nEpc, epc, 'K');
      resV=mind(nEpc, epc, 'V');
      //*/
      if(nSI>0)
	SIavr=SIavr/(double)nSI;
      else
	{
	  SImin=-9.99;
	  SIavr=-9.99;
	  SImax=-9.99;
	}
      betaSF=2*log(resX.modsf6/resX.modsf3)/log(2.);
      if(resC.modsf3>0. && resC.modsf6>0.)
	sf3_6Cpar=resC.modsf3/resC.modsf6;
      else
	sf3_6Cpar=-9.99;
      if(resX.modsf3>0. && resX.modsf6>0.)
	sf3_6Xpar=resX.modsf3/resX.modsf6;
      else
	sf3_6Xpar=-9.99;
      if(resK.modsf3>0. && resK.modsf6>0.)
	sf3_6Kpar=resK.modsf3/resK.modsf6;
      else
	sf3_6Kpar=-9.99;
      if(resV.modsf3>0. && resV.modsf6>0.)
	sf3_6Vpar=resV.modsf3/resV.modsf6;
      else
	sf3_6Vpar=-9.99;
      //if(!strcmp(src[iSrc].type, "BLL"))
	{
	  fprintf(fpChar, "%14s   %3d %7.2f %7.2f %+7.2f %+7.2f    %3d %7.2f %7.2f %+7.2f %+7.2f    %3d %7.2f %7.2f %+7.2f %+7.2f    %3d %7.2f %7.2f %+7.2f %+7.2f    %+5.2f  %+5.2f  %+5.2f   %4.2f\n",
		  src[iSrc].nameStd, resC.ndat, avrC, resC.mind*100., resC.modsf3*100., sf3_6Cpar, resX.ndat, avrX, resX.mind*100., resX.modsf3*100., sf3_6Xpar, resK.ndat, avrK, \
		  resK.mind*100., resK.modsf3*100., sf3_6Kpar, resV.ndat, avrV, resV.mind*100., resV.modsf3*100., sf3_6Vpar, SIavr, SImin, SImax, src[iSrc].z);
	  //fprintf(fpChar, "%3d %7.2f %7.2f %+7.2f %+7.2f    %3d %7.2f %7.2f %+7.2f %+7.2f    %3d %7.2f %7.2f %+7.2f %+7.2f    %3d %7.2f %7.2f %+7.2f %+7.2f    %+5.2f  %+5.2f  %+5.2f\n", \
	  resC.ndat, avrC, resC.mind*100., resC.modsf3*100., sf3_6Cpar, resX.ndat, avrX, resX.mind*100., resX.modsf3*100., sf3_6Xpar, resK.ndat, avrK, \
	  resK.mind*100., resK.modsf3*100., sf3_6Kpar, resV.ndat, avrV, resV.mind*100., resV.modsf3*100., sf3_6Vpar, SIavr, SImin, SImax);

	  printf("%14s    %3d %7.2f %7.2f %+7.2f %+7.2f    %3d %7.2f %7.2f %+7.2f %+7.2f    %3d %7.2f %7.2f %+7.2f %+7.2f    %3d %7.2f %7.2f %+7.2f %+7.2f    %+5.2f  %+5.2f  %+5.2f\n",
		 src[iSrc].nameStd, resC.ndat, avrC, resC.mind*100., resC.modsf3*100., sf3_6Cpar, resX.ndat, avrX, resX.mind*100., resX.modsf3*100., sf3_6Xpar, resK.ndat, avrK, resK.mind*100., resK.modsf3*100., sf3_6Kpar, resV.ndat, avrV, resV.mind*100., resV.modsf3*100., sf3_6Vpar, SIavr, SImin, SImax);

	  fprintf(fpTab1, "%s  &  %s %s  &  %s  &  %s & %s & %4.2f & %3d & %7.2f & %3d & %7.2f & %3d & %7.2f & %3d & %7.2f\\\\\n",
		  src[iSrc].nameStd, src[iSrc].cat, src[iSrc].catSrc, src[iSrc].type, src[iSrc].ra, src[iSrc].dec, src[iSrc].z, resC.ndat, avrC, resX.ndat, avrX, resK.ndat, avrK, resV.ndat, avrV);
	  fprintf(fpTab2,"%s  &  %7.2f & %7.2f & %7.2f &  %7.2f & %7.2f & %7.2f &  %7.2f & %7.2f & %7.2f &  %7.2f & %7.2f & %7.2f & %+5.2f\\\\\n",
		  src[iSrc].nameStd, resC.mind*100., resC.modsf3*100., sf3_6Cpar, resX.mind*100., resX.modsf3*100., sf3_6Xpar, resK.mind*100., resK.modsf3*100., sf3_6Kpar, resV.mind*100., resV.modsf3*100., sf3_6Vpar, SIavr);
	}
	nEpcTot=nEpcTot+resC.ndat+resX.ndat+resK.ndat+resV.ndat;
      /*
      if(!strcmp(source, "J0902+0443"))
	{
	  exit(1);
	}
      */
    }
  fclose(fpChar);
  fclose(fpTab1);
  fclose(fpTab2);

  system("rm epochsTmp.dat");
  system("rm srcTmp.csh");
  
  printf("Total number of datapoints: %d\nTotal number of epochs: %d\n", nFlxTot, nEpcTot);
  
  return 1;
}


int dim_file(FILE *fpin)
 {
   int n;
   char line[500];

   n=0;
   while (!feof(fpin)) {
      fgets(line, 300, fpin);
      if(line[0]!='#' && line[0]!='\n')
	{
	  n=n+1;
	}
  }

  return n-1;
}


result mind(int nDat, epoch *data, char band)
{
  int i, j, jb, k, nSteps, nStepsFlux, nsf3, nsf6, cycle;
  double ayMin, ayMax, ayAvr, step=-1., probUncert, trueFlux, stepFlux=-1., miMax=-1., miMin=-1., probFlux, expSum, tmpProb, tmpUncert, tmpFlux, norm;
  double miBest, likelihoodMax, sf3, aversf3, modsf3, sf6, aversf6, modsf6;
  double *ax, *ay, *az, *mInd, *likelihood, **prob;
  char line[300];
  result res;
  FILE *fpIn;

  ax=(double *)malloc(nDat*sizeof(double));
  ay=(double *)malloc(nDat*sizeof(double));
  az=(double *)malloc(nDat*sizeof(double));
  
  ayMax=0.;
  ayMin=DBL_MAX;
  ayAvr=0.;
  j=0;
  for(i=0; i<nDat; i++)
    {
      if(band=='C' && data[i].bC>0.)
	{
	  ax[j]=data[i].mjd;
	  ay[j]=data[i].bC;
	  az[j]=data[i].errC;
	  if(ay[j]>ayMax)
	    ayMax=ay[j];
	  if(ay[j]<ayMin)
	    ayMin=ay[j];
	  ayAvr=ayAvr+ay[j];
	  
	  j=j+1;
	}
      if(band=='X' && data[i].bX>0.)
	{
	  ax[j]=data[i].mjd;
	  ay[j]=data[i].bX;
	  az[j]=data[i].errX;
	  //printf("%f   %f   %f\n", ax[j], ay[j], az[j]);
	  if(ay[j]>ayMax)
	    ayMax=ay[j];
	  if(ay[j]<ayMin)
	    ayMin=ay[j];
	  ayAvr=ayAvr+ay[j];
	  
	  j=j+1;
	}
      if(band=='K' && data[i].bK>0.)
	{
	  ax[j]=data[i].mjd;
	  ay[j]=data[i].bK;
	  az[j]=data[i].errK;
	  if(ay[j]>ayMax)
	    ayMax=ay[j];
	  if(ay[j]<ayMin)
	    ayMin=ay[j];
	  ayAvr=ayAvr+ay[j];
	  
	  j=j+1;
	}
      if(band=='V' && data[i].bV>0.)
	{
	  ax[j]=data[i].mjd;
	  ay[j]=data[i].bV;
	  az[j]=data[i].errV;
	  if(ay[j]>ayMax)
	    ayMax=ay[j];
	  if(ay[j]<ayMin)
	    ayMin=ay[j];
	  ayAvr=ayAvr+ay[j];
	  
	  j=j+1;
	}
    }
  nDat=j;
  ayAvr=ayAvr/(double)nDat;

  if(nDat>10)
    {
      //printf("%f   %f   %f\n", ayMin, ayMax, ayAvr);

      for(cycle=0; cycle<2; cycle++)
	{
	  likelihoodMax=0.;
	  if(cycle==0)
	    {
	      if(ayMin<DBL_MIN)
		{
		  for(i=0; i<nDat; i++)
		    ay[i]=ay[i]-ayMin+0.1;
		  ayMax=ayMax-ayMin+0.1;
		  ayMin=0.1;
		}
	      if(miMax<0.)
		miMax=(ayMax-ayMin)/ayAvr;
	      if(miMin<0.)
		miMin=0.;
	      if(step<0.)
		step=(miMax-miMin)/100.;
	      if(stepFlux<0.)
		stepFlux=(ayMax-ayMin)/100.;
	      
	      nSteps=(int)((miMax-miMin)/step+1);
	      nStepsFlux=(int)((ayMax-ayMin)/stepFlux+1);
	      mInd=(double *)malloc(nSteps*sizeof(double));
	      likelihood=(double *)malloc(nSteps*sizeof(double));
	      prob=(double **)malloc(nSteps*nStepsFlux*sizeof(double));
	      for(i=0; i<nSteps; i++)
		{
		  prob[i]=(double *)malloc(nStepsFlux*sizeof(double));
		}      
	      //printf("# nSteps: %d   nStepsFlux: %d   miMin: %f   miMax: %f\n", nSteps, nStepsFlux, miMin, miMax);
	    }
	  else
	    {
	      miMin=miBest-2*step;
	      if(miMin<0.)
		miMin=0.;
	      miMax=miBest+2*step;
	      step=(miMax-miMin)/100.;
	      //printf("# nSteps: %d   nStepsFlux: %d   miMin: %f   miMax: %f\n", nSteps, nStepsFlux, miMin, miMax);
	    }
	  
	  norm=0.;
	  for(i=0; i<nSteps; i++)
	    {
	      mInd[i]=miMin+i*step;
	      for(k=0; k<nStepsFlux; k++)
		{
		  trueFlux=ayMin+stepFlux*k;
		  probUncert=1.;
		  expSum=0.;
		  prob[i][k]=log10(trueFlux*stepFlux);
		  //printf("%e   %e   %e\n", mInd[i], trueFlux, prob[i][k]);
		  for(j=0; j<nDat; j++)
		    {
		      tmpUncert=1./(sqrt(2*M_PI*(pow(mInd[i]*trueFlux, 2.)+pow(az[j], 2.))));
		      tmpFlux=pow(ay[j]-trueFlux, 2.)/(pow(mInd[i]*trueFlux, 2.)+pow(az[j], 2.));
		      
		      prob[i][k]=prob[i][k]+log10(tmpUncert*exp(-0.5*tmpFlux));
		    }
		  if(norm<prob[i][k])
		    {
		      norm=prob[i][k];
		      //printf("%e   %e   %e   %e\n", mInd[i], trueFlux, prob[i][k], norm);
		    }
		}
	    }
	  for(i=0; i<nSteps; i++)
	    {
	      likelihood[i]=0.;
	      for(k=0; k<nStepsFlux; k++)
		{
		  prob[i][k]=exp(prob[i][k]-norm);
		  likelihood[i]=likelihood[i]+prob[i][k];
		}
	    }
	  miBest=0.;
	  likelihoodMax=0.;
	  for(i=0; i<nSteps; i++)
	    {
	      //printf("%e   %e\n", mInd[i], likelihood[i]);
	      if(likelihood[i]>likelihoodMax)
		{
		  likelihoodMax=likelihood[i];
		  miBest=mInd[i];
		}
	    }
	  //printf("MI best: %e   %e\n", miBest, likelihoodMax);
	}

      /* Mod SF at three hours */
      sf3=0.;
      nsf3=0;
      aversf3=0.;
      sf6=0.;
      nsf6=0;
      aversf6=0.;

      for(j=0; j<nDat; j++)
	{
	  for(jb=0; jb<nDat; jb++)
	    {
	      ///*
	      //if(fabs(ax[j]-ax[jb])>2.9*365.25 && fabs(ax[j]-ax[jb])<3.1*365.25)
	      if(fabs(ax[j]-ax[jb])>1.4*365.25 && fabs(ax[j]-ax[jb])<1.6*365.25)
		{
		  sf3=sf3+pow(ay[j]-ay[jb], 2.);
		  aversf3=aversf3+ay[j]+ay[jb];
		  nsf3=nsf3+1;
		}
	      //if(fabs(ax[j]-ax[jb])>5.8*365.25 && fabs(ax[j]-ax[jb])<6.2*365.25)
	      if(fabs(ax[j]-ax[jb])>2.9*365.25 && fabs(ax[j]-ax[jb])<3.1*365.25)
		{
		  sf6=sf6+pow(ay[j]-ay[jb], 2.);
		  aversf6=aversf6+ay[j]+ay[jb];
		  nsf6=nsf6+1;
		}
	      //*/
	      /*
	      if(fabs(ax[j]-ax[jb])>1.4*365.25 && fabs(ax[j]-ax[jb])<1.6*365.25)
		{
		  if(fabs(ay[j]-ay[jb])>sqrt(pow(az[j], 2.)+pow(az[jb], 2.)))
		    {
		      sf3=sf3+pow(fabs(ay[j]-ay[jb])-sqrt(pow(az[j], 2.)+pow(az[jb], 2.)), 2.);
		    }
		  aversf3=aversf3+ay[j]+ay[jb];
		  nsf3=nsf3+1;
		}
	      if(fabs(ax[j]-ax[jb])>2.9*365.25 && fabs(ax[j]-ax[jb])<3.1*365.25)
		{
		  if(fabs(ay[j]-ay[jb])>sqrt(pow(az[j], 2.)+pow(az[jb], 2.)))
		    {
		      sf6=sf6+pow(fabs(ay[j]-ay[jb])-sqrt(pow(az[j], 2.)+pow(az[jb], 2.)), 2.);
		    }
		  aversf6=aversf6+ay[j]+ay[jb];
		  nsf6=nsf6+1;
		}
	      */
	    }
	}
      sf3=sf3/(double)nsf3;
      if(nsf3>0)
	{
	  aversf3=aversf3/(2*(double)nsf3);
	  modsf3=sqrt(sf3/2.)/aversf3;
	}
      sf6=sf6/(double)nsf6;
      if(nsf6>0)
	{
	  aversf6=aversf6/(2*(double)nsf6);
	  modsf6=sqrt(sf6/2.)/aversf6;
	}
      /* End Mod SF at three hours */
      
    }
  res.ndat=nDat;
  res.aversf3=aversf3;
  res.aversf6=aversf6;
  if(nDat>10)
    {
      res.mind=miBest;
    }
  else
    {
      res.mind=-0.0999;     
    }
  if(nDat>10 && nsf3>0)
    res.modsf3=modsf3;
  else
    res.modsf3=-0.0999;
  if(nDat>10 && nsf6>0)
    res.modsf6=modsf6;
  else
    res.modsf6=-0.0999;
 
  return res;
}
