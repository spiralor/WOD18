#include <stdlib.h>
#include <stdio.h>

/***********************************************************

 DEFINED CONSTANTS

 MAXTAX: MAXIMUM NUMBER OF TAXA DELIMITORS IN A TAXA SET
 MAXSEC: MAXIMUM NUMBER OF SECONDARY HEADERS
 MAXBIO MAXIMUM NUMBER OF BIOLOGICAL HEADERS
 MAXPARM: MAXIMUM NUMBER OF MEASURED AND/OR CALCULATED VARIABLES
 MAXSECP: MAXIMUM NUMBER OF VARIABLE SPECIFIC SECONDARY HEADERS

************************************************************/

#define maxtax 30
#define maxsec 100
#define maxbio 50
#define maxparm 100
#define maxpsec 25 * maxparm

/**********************************************************

 *FP - WOD DATA FILE TO BE OPENED FOR READING
 *FPOUT - FILE TO BE OPENED FOR WRITING OUTPUT

***********************************************************/

FILE *fp,*fpout;

/**********************************************************

 HEADER INFORMATION FOR EACH CAST

 CC - COUNTRY CODE OF COUNTRY RESPONSIBLE FOR TAKING DATA AT CAST
 ICRUISE - NODC/WOD ASSIGNED CRUISE NUMBER (UNIQUE WITH COUNTRY CODE)
 OCAST - WOD UNIQUE CAST NUMBER
 YEAR, MONTH, DAY, HOUR - TIME AND DATA OF CAST
 LATITUDE, LONGITUDE - POSITION OF CAST
 LEVELS - NUMBER OF DEPTH LEVELS MEASURED AT CAST
 ISOOR - SET TO ZERO FOR OBSERVED LEVEL DATA, ONE FOR STANDARD LEVEL
 IP2 - VARIABLE CODES FOR VARIABLES MEASURED AT CAST
 IPERROR - ERROR FLAGS FOR EACH VARIABLE MEASURED AT CAST 

 HTOTFIG - TOTAL NUMBER OF FIGURES FOR 1: HOUR, 2: LATITUDE, 3: LONGITUDE
 HSIGFIG - NUMBER SIGNIFICANT FIGURES FOR 1: HOUR, 2: LATITUDE, 3: LONGITUDE
 HTOTFIG - NUMBER OF FIGURES TO THE RIGHT OF THE DECIMAL FOR
           1: HOUR, 2: LATITUDE, 3: LONGITUDE

**************************************************************/

char cc[2];
int icruise=0, ostation=0, year=0, month=0, day=0;
int hour,longitude,latitude;
int levels,isoor,nparm,ip2[maxparm],iperror[maxparm];
int htotfig[3],hsigfig[3],hrightfig[3];

/*************************************************************

 CHARACTER AND PRIMARY INVESTIGATOR DATA

 ORIGCFIG - NUMBER OF CHARACTERS IN ORIGINATORS CRUISE CODE
            (ZERO IF NOT GIVEN)
 ORIGC - ORIGINATORS CRUISE CODE
 ORIGSFIG - NUMBER OF CHARACTERS IN ORIGINATORS STATION CODE
            (ZERO IF NOT GIVEN)
 ORIGS - ORIGINATORS STATION CODE

 NPI - NUMBER OF PRIMARY INVESTIGATORS
 IPIP - VARIABLE FOR PRIMARY INVESTIGATOR
 IPI - PRIMARY INVESTIGATORY

*****************************************************************/

int origcfig,origsfig;
char origc[30],origs[30];
int ipip[maxparm],ipi[maxparm],npi;

/*************************************************************

 SECONDARY HEADER INFORMATION

 NSEC - NUMBER OF SECONDARY HEADERS
 STOTFIG - TOTAL FIGURES FOR EACH SECONDARY HEADER
 SSIGFIG - NUMBER OF SIGNIFICANT FIGURES FOR EACH SECONDARY HEADER
 SRIGHTFIG - NUMBER OF FIGURES RIGHT OF DECIMAL FOR EACH SECONDARY HEADER 
 SECCODE - SECONDARY HEADER CODE
 SECVAL - SECONDARY HEADER VALUE

**************************************************************/

int nsec;
int stotfig[maxsec],ssigfig[maxsec],srightfig[maxsec];
int seccode[maxsec],secval[maxsec];

/*************************************************************

 VARIABLE SPECIFIC SECONDARY HEADER INFORMATION

 NPSEC - NUMBER OF VARIABLE SPECIFIC SECONDARY HEADERS
 PSTOTFIG - TOTAL FIGURES FOR EACH VARIABLE SPECIFIC SECONDARY HEADER
 PSSIGFIG - NUMBER OF SIGNIFICANT FIGURES FOR EACH VARIABLE
            SPECIFIC SECONDARY HEADER
 PSRIGHTFIG - NUMBER OF FIGURES RIGHT OF DECIMAL FOR EACH VARIABLE
              SPECIFIC SECONDARY HEADER 
 PSECPARM - VALUE OF VARIABLE SPECIFIC SECONDARY HEADER
 PSECCODE - VARIABLE SPECIFIC SECONDARY HEADER CODE
 PSECVAL - VARIABLE SPECIFIC SECONDARY HEADER VALUE

**************************************************************/

int npsec;
int pstotfig[maxpsec],pssigfig[maxpsec],psrightfig[maxpsec];
int psecparm[maxpsec],pseccode[maxpsec],psecval[maxpsec];

/*************************************************************

 BIOLOGICAL HEADER INFORMATION

 NBIO - NUMBER OF BIOLOGICAL HEADERS
 BTOTFIG - TOTAL FIGURES FOR EACH BIOLOGICAL HEADER
 BSIGFIG - NUMBER OF SIGNIFICANT FIGURES FOR EACH BIOLOGICAL HEADER
 BRIGHTFIG - NUMBER OF FIGURES RIGHT OF DECIMAL FOR EACH BIOLOGICAL HEADER 
 BIOCODE - BIOLOGICAL HEADER CODE
 BIOVAL - BIOLOGICAL HEADER VALUE

**************************************************************/

int nbio;
int btotfig[maxbio],bsigfig[maxbio],brightfig[maxbio];
int biocode[maxbio],bioval[maxbio];

/*************************************************************

 TAXA SET INFORMATION

 NTSETS - NUMBER OF TAXA SETS AT CAST
 NTLOC - NUMBER OF DELIMITORS FOR EACH TAXA SET
 NTCODE - TAXA CODE FOR EACH TAXA VARIABLE
 NTVAL  - VALUE FOR EACH TAXA CODE
 NTERR - ERROR FLAG FOR EACH TAXA VALUE
 NTOERR - ORIGINATORS FLAG FOR EACH TAXA VALUE
 NTTOTFIG - TOTAL FIGURES FOR EACH TAXA VALUE
 NTSIGFIG - NUMBER OF SIGNIFICANT FIGURES FOR EACH TAXA VALUE
 NTRIGHTFIG - NUMBER OF FIGURES RIGHT OF THE DECIMAL FOR EACH TAXA VALUE

***************************************************************/

int ntsets;
int *ntloc,*ntcode,*ntval,*nterr,*ntoerr,*nttotfig,*ntsigfig,*ntrightfig;

/*************************************************************

 DEPTH INFORMATION

 DEPTH  - DEPTH AT EACH DEPTH LEVEL
 ZERR - DEPTH ERROR FLAG
 ZOERR - DEPTH ORIGINATORS FLAG
 ZTOTFIG - TOTAL FIGURES FOR EACH DEPTH VALUE
 ZSIGFIG - NUMBER OF SIGNIFICANT FIGURES FOR EACH DEPTH
 ZRIGHTFIG - NUMBER OF FIGURES RIGHT OF THE DECIMAL FOR EACH DEPTH

***************************************************************/

int *depth,*zerr,*zoerr,*ztotfig,*zsigfig,*zrightfig;

/*************************************************************

 MEASURED VARIABLE INFORMATION

 DVAL  - DATA VALUE AT EACH DEPTH LEVEL
 DERR - ERROR FLAG FOR EACH VARIABLE AT EACH LEVEL
 DOERR - ORIGINATORS FLAG FOR EACH VARIABLE AT EACH LEVEL
 DTOTFIG - TOTAL FIGURES FOR EACH DATA VALUE
 DSIGFIG - NUMBER OF SIGNIFICANT FIGURES FOR EACH DATA VALUE
 DRIGHTFIG - NUMBER OF FIGURES RIGHT OF THE DECIMAL FOR EACH DATA VALUE

***************************************************************/

int *dataval,*derr,*doerr,*dtotfig,*dsigfig,*drightfig;

/*************************************************************

 SIZE DELIMITORS

 ISIZE - PRESENT ARRAY SIZE NEEDED TO FIT ALL MEASURED VARIABLES
 ZSIZE - PRESENT ARRAY SIZE NEEDED TO FIT ALL DEPTHS
 NTSETSMAX - MAXIMUM NUMBER OF TAXA SETS
 ISIZEMAX - MAXIMUM ARRAY SIZE YET ENCOUNTERED FOR MEASURED VARIABLES
 ZSIZEMAX - MAXIMUM ARRAY SIZE YET ENCOUNTERED FOR DEPTHS

/*************************************************************/

int isize,zsize;
int ntsetsmax=0, isizemax=0,zsizemax=0;
int iVERSflag;

/*************************************************************

 INTERNAL ARRAYS

 TENP - POWERS OF TEN
 SDEPTH - STANDARD DEPTH LEVELS

MODIFIED SEP. 25, 2008 AS PER EMAIL FROM TERRY DEVEAU

*************************************************************/

 double tenp[] = { 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000,
                  100000000 };

 int sdepth[] = {  0., 10., 20., 30., 50., 75., 100., 125., 150.,
           200., 250., 300., 400., 500., 600., 700., 800., 900.,
           1000., 1100., 1200., 1300., 1400., 1500., 1750., 2000.,
           2500., 3000., 3500., 4000., 4500., 5000., 5500., 6000.,
           6500., 7000., 7500., 8000., 8500., 9000. };

main()

{

 char filename[80];
 int i=0, j, k, s, iend=0;
 int ncast=0;

 printf(" Enter input file name\n");
 scanf("%s",filename);

 if ((fp = fopen(filename,"rb+\0")) == NULL)
  printf("UNABLE TO OPEN FILE\n");

 else {

/********************************************************

 INITIALIZE DYNAMIC ARRAYS

*********************************************************/

  spacer(1);

/*   GET USER INFORMATION (NUMBER OF CASTS, OUTPUT FILE NAME) */

  printf(" Enter output file name\n");
  scanf("%s",filename);

  if ((fpout = fopen(filename,"w\0")) == NULL) {
   printf("UNABLE TO OPEN FILE\n");
  }

  printf(" ENTER NUMBER OF CASTS TO VIEW");
  printf (" (0 FOR ALL CASTS IN FILE)\n");
  if ( (s = scanf("%d",&ncast)) == 0 ) ncast=0;

  if ( ncast == 0 ) ncast=100000000;

  while ( !feof(fp) && (i++) < ncast ) {

/********************************************************

 READ IN CAST

**********************************************************/

   if ( ( iend = oclread() ) == -1 ) printf(" END OF FILE REACHED\n");

   else {

/********************************************************

 ENTER STANDARD LEVEL DEPTHS, IN CASE THIS IS STANDARD LEVEL DATA
 ONLY FOR FORMAT VERSIONS OLDER THAN WOD13

********************************************************/

    if (  iVERSflag != 2 )
     for ( j = 0; j < 40; j++ ) *(depth+j)= *(sdepth+j);

    printstation(i);
 
   }

  }

 }

}
oclread()

{

 int i,j;

 char wodform;
 int totfig=0, sigfig=0, rightfig=0;
 int nbytet,ntypec,nbytec,nbytes,nbyteb;
 int ninfc,ntoff, doff;
 int missing=-9999;
 int npinfs=0,npinfe=0,npinf;
 int iend=0;

/**********************************************************

 READ IN WOD FORMAT CODE: 'C' IS FOR WOD13 FORMAT.
 'B' FOR WOD05/WOD09 FORMAT
 'A' FOR WOD01 FORMAT.

***********************************************************/

 totfig= 1;
 if ( (iend = extractc(0,&totfig,&wodform)) == -1 ) return iend;

 if ( wodform == 'C' ) iVERSflag=2;

/**********************************************************

 READ IN NUMBER OF BYTES IN ASCII CAST 

***********************************************************/

 if ( ( iend = extracti(0,&totfig,&sigfig,&rightfig,&nbytet,missing))
          == -1 ) return iend;

/**********************************************************

 READ IN WOD CAST NUMBER

***********************************************************/

 if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&ostation,missing))
       == -1 ) return iend;

/*********************************************************

 READ IN NODC COUNTRY CODE

**********************************************************/

 totfig= 2;
 if ( (iend = extractc(0,&totfig,cc)) == -1 ) return iend;

/**********************************************************

 READ IN CRUISE NUMBER

***********************************************************/

 if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&icruise,
             missing)) == -1) return iend;

/**********************************************************

 READ IN YEAR, MONTH, DAY, TIME

***********************************************************/

 totfig=4;
 if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&year,missing))
       == -1) return iend;
 totfig=2;
 if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&month,missing))
       == -1) return iend;
 totfig=2;
 if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&day,missing))
       == -1) return iend;

 if ( (iend = extracti(1,htotfig,hsigfig,hrightfig,
              &hour, 9999)) == -1) return iend;

/**********************************************************

 READ IN LATITUDE AND LONGITUDE

***********************************************************/

 if ( (iend = extracti(1,(htotfig+1),(hsigfig+1),(hrightfig+1),
              &latitude, -9999)) == -1) return iend;

 if ( (iend = extracti(1,(htotfig+2),(hsigfig+2),(hrightfig+2),
              &longitude, -99999)) == -1) return iend;

/**********************************************************

 READ IN NUMBER OF LEVELS

***********************************************************/

 if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&levels,
             missing)) == -1) return iend;

/**********************************************************

 READ IN OBSERVED (0) OR STANDARD (1) LEVELS

***********************************************************/

 totfig=1;
 if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&isoor,missing))
       == -1) return iend;

/**********************************************************

 READ IN NUMBER OF VARIABLES AT THIS CAST (NOT INCLUDING
 BIOLOGY)

***********************************************************/

 totfig=2;
 if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&nparm,missing))
       == -1) return iend;

/**********************************************************

 READ IN EACH VARIABLE CODE AND PROFILE ERROR CODE

***********************************************************/

 for ( i =0; i< nparm; i++ ) {

  if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(ip2+i),
             missing)) == -1) return iend;
  totfig=1;
  if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(iperror+i),
             missing)) == -1) return iend;

/*******************************************************************

 READ NUMBER OF VARIABLE SPECIFIC SECOND HEADER VARIABLES

********************************************************************/

  if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&npinf,missing))
       == -1 ) return iend;
  npinfe += npinf;

/*******************************************************************

 READ IN EACH VARIABLE SPECIFIC SECOND HEADER VARIABLE

********************************************************************/

  for ( j = npinfs; j < npinfe; j++ ) {

   *(psecparm+j) = *(ip2+i);
   if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(pseccode+j),
             missing)) == -1) return iend;
   if ( (iend = extracti(1,(pstotfig+j),(pssigfig+j),(psrightfig+j),
                        (psecval+j), missing)) == -1) return iend;

  }

  npinfs += npinf;

 }

 npsec = npinfe;

/***************************************************************

 READ IN NUMBER OF BYTES IN CHARACTER AND PRIMARY INVESTIGATOR FIELDS

****************************************************************/

 if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&nbytec,
             missing)) == -1) return iend;

/**********************************************************

 READ IN NUMBER OF INFORMATION TYPES (MAX 3: FOR CRUISE CODE,
 STATION CODE AND PI INFORMATION)

***********************************************************/

 origcfig= 0;
 origsfig= 0;
 npi= 0;

 if ( nbytec > 0 ) {

  totfig=1;
  if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&ninfc,missing))
       == -1) return iend;

/**********************************************************

 READ IN TYPE OF INFORMATION : 1= CRUISE CODE,
 2=STATION CODE AND 3=PI INFORMATION)

***********************************************************/

  for ( i= 0; i < ninfc; i++) {

   totfig=1;
   if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&ntypec,missing))
       == -1) return iend;

/***********************************************************

 READ IN ORIGINATORS CRUISE CODE

************************************************************/

   if ( ntypec == 1) {

    totfig=2;
    if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&origcfig,
     missing)) == -1) return iend;
    if ( (iend = extractc(0,&origcfig,origc)) == -1 ) return iend;

   }

/***********************************************************

 READ IN ORIGINATORS STATION CODE

************************************************************/

   else if ( ntypec == 2 ) {

    totfig=2;
    if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&origsfig,
     missing)) == -1) return iend;
    if ( (iend = extractc(0,&origsfig,origs)) == -1 ) return iend;

   }

/***********************************************************

 READ IN PRIMARY INVESTIGATOR INFORMATION

************************************************************/

   else if ( ntypec == 3 ) {

/**********************************************************

 READ IN NUMBER OF PRIMARY INVESTIGATORS FOR THIS CAST

***********************************************************/

    totfig=2;
    if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,&npi,missing))
         == -1) return iend;

/**********************************************************

 READ IN EACH VARIABLE CODE AND PRIMARY INVESTIGATOR CODE
 FOR THAT VARIABLE

***********************************************************/

    for ( j =0; j< npi; j++ ) {

     if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(ipip+j),
              missing)) == -1) return iend;
     if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(ipi+j),
              missing)) == -1) return iend;

    }

   }

  }

 }

/***************************************************************

 READ IN NUMBER OF BYTES IN SECONDARY HEADER FIELDS

****************************************************************/

 if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&nbytes,
             missing)) == -1) return iend;

 nsec = 0;
 if ( nbytes > 0 ) {

/**************************************************************

 READ IN NUMBER OF SECONDARY HEADER VARIABLES

***************************************************************/

  if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&nsec,
              missing)) == -1) return iend;

/**********************************************************

 READ IN EACH SECONDARY HEADER VARIABLE CODE AND VALUE

***********************************************************/

  for ( i =0; i< nsec; i++ ) {

   if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(seccode+i),
             missing)) == -1) return iend;
   if ( (iend = extracti(1,(stotfig+i),(ssigfig+i),(srightfig+i),
                        (secval+i), missing)) == -1) return iend;

  }
 
 }

/***************************************************************

 READ IN NUMBER OF BYTES IN BIOLOGY HEADER AND TAXA SET FIELDS

****************************************************************/

 if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&nbyteb,
             missing)) == -1) return iend;

 nbio= 0;
 ntsets= 0;
 if ( nbyteb > 0 ) {

/**************************************************************

 READ IN NUMBER OF BIOLOGY HEADER VARIABLES

***************************************************************/

  if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&nbio,
              missing)) == -1) return iend;

/**********************************************************

 READ IN EACH BIOLOGY HEADER VARIABLE CODE AND VALUE

***********************************************************/

  for ( i =0; i< nbio; i++ ) {

   if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(biocode+i),
             missing)) == -1) return iend;
   if ( (iend = extracti(1,(btotfig+i),(bsigfig+i),(brightfig+i),
                        (bioval+i), missing)) == -1) return iend;

  }
  
/***************************************************************

 READ IN NUMBER OF TAXA SETS

****************************************************************/

  if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,&ntsets,
             missing)) == -1) return iend;

/***************************************************************

 IF NUMBER OF TAXA SETS IS GREATER THAN THE PREVIOUS MAXIMUM,
 REALLOCATE SPACE ACCORDINGLY

****************************************************************/

  if ( ntsets > ntsetsmax ) {

   ntsetsmax= ntsets;
   spacer(2);

  }

  if ( ntsets > 0 ) {

/**************************************************************

 READ IN NUMBER OF ENTRIES FOR THIS TAXA SET, CALCULATE OFFSET
 FOR READING IN DATA

***************************************************************/
  
   for ( j = 0; j < ntsets; j++) {

    if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(ntloc+j),
               missing)) == -1) return iend;

    ntoff= maxtax * j;

/**********************************************************

 READ IN EACH TAXA SET VARIABLE CODE AND VALUE,
 ERROR FLAG AND ORIGINATORS FLAG

***********************************************************/

    for ( i =0; i< *(ntloc+j); i++ ) {

     if ( (iend = extracti(0,&totfig,&sigfig,&rightfig,(ntcode+ntoff+i),
              missing)) == -1) return iend;
     if ( (iend = extracti(1,(nttotfig+ntoff+i),(ntsigfig+ntoff+i),
                 (ntrightfig+ntoff+i), (ntval+ntoff+i), missing))
                  == -1) return iend;
     totfig=1;
     if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(nterr+ntoff+i),
              missing)) == -1) return iend;
     totfig=1;
     if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(ntoerr+ntoff+i),
              missing)) == -1) return iend;

    }

   }
  
  }

 }

/***************************************************************

 IF NUMBER OF DEPTHS IS GREATER THAN THE PREVIOUS MAXIMUM,
  REALLOCATE SPACE ACCORDINGLY

****************************************************************/

 zsize= levels;
/*if ( isoor == 0 && zsize > zsizemax ) { */
 if ( zsize > zsizemax ) {

  zsizemax=zsize;
  spacer(3);

 }

/***************************************************************

 IF NUMBER OF DEPTHS MULTIPLIED BY NUMBER OF VARIABLES IS
 GREATER THAN THE PREVIOUS MAXIMUM, REALLOCATE SPACE ACCORDINGLY

****************************************************************/

 isize= nparm * levels;
 if ( isize > isizemax ) {

  isizemax=isize;
  spacer(4);

 }

/**********************************************************

 READ IN EACH DEPTH VALUE, ERROR FLAG, AND ORIGINATORS FLAG

***********************************************************/

 for ( j = 0; j < levels; j++ ) { 

   if ( isoor == 0 || wodform == 'C' ) {   

   if ( (iend = extracti(1,(ztotfig+j),(zsigfig+j),
               (zrightfig+j), (depth+j), missing))
               == -1) return iend;
   totfig=1;
   if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(zerr+j),
            missing)) == -1) return iend;
   totfig=1;
   if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(zoerr+j),
            missing)) == -1) return iend;

  }   

/**********************************************************

 READ IN EACH DATA VALUE

***********************************************************/

  for ( i =0; i< nparm; i++ ) {

   doff= i * levels;


   if ( (iend = extracti(1,(dtotfig+doff+j),(dsigfig+doff+j),
               (drightfig+doff+j), (dataval+doff+j), missing))
                == -1) return iend;

   if ( *(dtotfig+doff+j) > 0 ) {

    totfig=1;
    if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(derr+doff+j),
             missing)) == -1) return iend;
    totfig=1;
    if ( (iend = extracti(2,&totfig,&sigfig,&rightfig,(doerr+doff+j),
             missing)) == -1) return iend;
   }

   else {

    *(derr+doff+j)=0;
    *(doerr+doff+j)=0;
    *(drightfig+doff+j)=2;

   }


  }

 }

/***********************************************************

 READ TO END OF CAST

************************************************************/

 while ( ( i = fgetc(fp)) != '\n' && !feof(fp) );
 return iend;

}

/***************************************************************

              FUNCTION PRINTSTATION

 PRINTSTATION PRINTS CAST TO THE SCREEN

***************************************************************/

printstation(int i)

{

 int j, k, offs;
 float xlon, xlat, xhour, x;


/***************************************************************

 PRINT OUT HEADER

****************************************************************/

 fprintf(fpout," \n-----------------------------------------------------\n");
 fprintf(fpout," OUTPUT FROM ASCII FILE: CAST #%d\n",i);
 fprintf(fpout,"-----------------------------------------------------\n\n");

 xhour= (hour/ tenp[ *(hrightfig) ]);
 xlat= (latitude/ tenp[ *(hrightfig+1) ]);
 xlon= (longitude/ tenp[ *(hrightfig+2) ]);

 fprintf(fpout," CC  Cruise Latitude Longitude YYYY MM DD");
 fprintf(fpout,"  Time Cast    #Levels\n");
 fprintf(fpout," %2s  %6d  %7.3f  %8.3f %4d %2d %2d %5.2f %8d %4d\n\n\n",
   cc,icruise,xlat,xlon,year,month,day,xhour,ostation,levels);

/**************************************************************

 PRINT OUT CHARACTER DATA AND PRIMARY INVESTIGATORS DATA

***************************************************************/

 *(origc+origcfig)= '\0';
 *(origs+origsfig)= '\0';

 if ( origcfig > 0 ) fprintf(fpout," ORIGINATORS CRUISE CODE: %s\n",origc);
 if ( origsfig > 0 ) fprintf(fpout," ORIGINATORS STATION CODE: %s\n",origs);
 if ( origcfig > 0 || origsfig > 0 ) fprintf(fpout,"\n");

 for ( j= 0; j < npi; j++ ) {

  fprintf(fpout," PRIMARY INVESTIGATOR: %4d ... for variable #%2d\n",
   *(ipi+j),*(ipip+j));

 }
 
 if ( npi > 0 ) fprintf(fpout,"\n");

/**************************************************************

 PRINT OUT MEASURED VARIABLES

***************************************************************/

 if ( levels > 0 ) {

  fprintf(fpout,"      z          ");
  for ( j = 0; j < nparm; j++ ) fprintf(fpout,"  %2d            ",*(ip2+j));
  fprintf(fpout,"\n\n");

  for ( k = 0; k < levels; k++ ) {

   x = (*(depth+k)/ tenp[*(zrightfig+k)]);
   fprintf(fpout," %7.1f (%d) [%d%d]",x, *(zsigfig+k), *(zerr+k),*(zoerr+k));

   for ( j = 0; j < nparm; j++ ) {

    offs= j * levels;

    x = (*(dataval+offs+k)/ tenp[*(drightfig+offs+k)]);
    fprintf(fpout,"%7.3f (%d) [%d%d]",x, *(dsigfig+offs+k),
            *(derr+offs+k),*(doerr+offs+k));

   }

   fprintf(fpout,"\n");

  }

/************************************************************


 PRINT OUT PROFILE ERROR FLAGS

*************************************************************/

  fprintf(fpout, "\n ERR:           "); 
  for ( j = 0; j < nparm; j++ ) fprintf(fpout,"              %d ",
   *(iperror+j));
  fprintf(fpout,"\n\n");

 }

/************************************************************

 PRINT OUT SECONDARY HEADER INFORMATION

*************************************************************/

 for ( j = 0; j < nsec; j++ ) {

  fprintf(fpout,"   Second header #%3d  ", *(seccode+j));

  if ( *(srightfig+j) > 0 ) {

   x= (*(secval+j)/tenp[*(srightfig+j)]);
   fprintf(fpout,"%9.3f (%d)\n",x,*(ssigfig+j));

  }

  else fprintf(fpout,"%9d (%d)\n",*(secval+j),*(ssigfig+j));

 }
 fprintf(fpout,"\n");

/************************************************************

 PRINT OUT VARIABLE SPECIFIC SECONDARY HEADER INFORMATION

*************************************************************/

 for ( j = 0; j < npsec; j++ ) {

  fprintf(fpout,"   Variable #%3d Second header #%3d  ", *(psecparm+j),
   *(pseccode+j));

  if ( *(psrightfig+j) > 0 ) {

   x= (*(psecval+j)/tenp[*(psrightfig+j)]);
   fprintf(fpout,"%9.3f (%d)\n",x,*(pssigfig+j));

  }

  else fprintf(fpout,"%9d (%d)\n",*(psecval+j),*(pssigfig+j));

 }
 fprintf(fpout,"\n");

/************************************************************

 PRINT OUT BIOLOGY HEADER INFORMATION

*************************************************************/

 for ( j = 0; j < nbio; j++ ) {

  fprintf(fpout,"   Biology header #%3d  ", *(biocode+j));

  if ( *(brightfig+j) > 0 ) {

   x= (*(bioval+j)/tenp[*(brightfig+j)]);
   fprintf(fpout,"%9.3f (%d)\n",x,*(bsigfig+j));

  }

  else fprintf(fpout,"%9d (%d)\n",*(bioval+j),*(bsigfig+j));

 }

 fprintf(fpout,"\n");

/************************************************************

 PRINT OUT TAXA INFORMATION

*************************************************************/

 for ( j=0; j < ntsets; j++ ) {

  offs= j*maxtax; 

  for ( k=0; k < *(ntloc+j); k++ ) {

   if ( *(ntcode + offs + k) == 1 ) {

    fprintf(fpout,"\n Taxa-set %3d: Taxonomic Code[1]# %10d (%d) [%d%d]\n",
    j+1,*(ntval+ offs + k),*(ntsigfig+ offs+ k),*(nterr+ offs + k),
    *(ntoerr+ offs + k));

   }

   else {

    if ( *(ntrightfig+ offs + k) > 0 ) {

     x= (*(ntval+ offs + k)/tenp[*(ntrightfig+ offs + k)]);
     fprintf(fpout,"        Code #%2d        %8.3f (%d) [%d%d]\n",
      *(ntcode+ offs + k), x, *(ntsigfig+ offs + k),*(nterr+offs+k),
      *(ntoerr+offs+k));

    }

    else { 

     fprintf(fpout,"        Code #%2d        %8d (%d) [%d%d]\n",
      *(ntcode+ offs + k), *(ntval+ offs+k), *(ntsigfig+ offs + k),
      *(nterr+offs+k),*(ntoerr+offs+k));
    }

   }

  }

 }

}

/************************************************************

 SPACER.C SETS UP ORIGINAL SPACING FOR ALL DYNAMIC ARRAYS

*************************************************************/

spacer(
 
 int intime      /* SET TO ONE TO INITIALIZE ALL DYNAMIC ARRAYS,
                    SET TO TWO TO REDIMENSION TAXA ARRAYS,
                    SET TO THREE TO REDIMENSION DEPTH, 
                    SET TO FOUR TO REDIMENSION MEASURED VARIABLE ARRAYS
                 */

      )

{

 if ( intime == 1 ) { 

/***************************************************************

 ALLOCATE SPACE FOR TAXA SET DATA

****************************************************************/

  ntsets=1;
  ntsetsmax=1;
  if ( (ntcode =calloc( ntsets * maxtax, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntloc =calloc( ntsets * maxtax, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (nttotfig =calloc( ntsets * maxtax, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntsigfig =calloc( ntsets * maxtax, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntrightfig =calloc( ntsets * maxtax, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntval =calloc( ntsets * maxtax, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (nterr =calloc( ntsets * maxtax, sizeof(int)) )
        == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntoerr =calloc( ntsets * maxtax, sizeof(int)) )
        == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
 
/***************************************************************

 ALLOCATE SPACE FOR DEPTH AND MEASURED VARIABLES

****************************************************************/

  zsize= 137;
  zsizemax= 137;
  if ( (ztotfig =calloc(zsize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (zsigfig =calloc(zsize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (zrightfig =calloc(zsize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (depth =calloc( zsize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (zerr =calloc( zsize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (zoerr =calloc( zsize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
 
  isize=137;
  isizemax= 137;
  if ( (dtotfig =calloc(isize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA\n");
  if ( (dsigfig =calloc(isize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA\n");
  if ( (drightfig =calloc(isize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA\n");
  if ( (dataval =calloc( isize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA\n");
  if ( (derr =calloc( isize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA\n");
  if ( (doerr =calloc( isize, sizeof(int)) )
       == NULL )
    printf( " NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA\n");
 
 }

/***************************************************************

 REALLOCATE SPACE FOR TAXA DATA

****************************************************************/

 else if ( intime == 2 ) {

   ntsetsmax = ntsets;
   if ( (ntcode =realloc( ntcode, ntsets * maxtax * sizeof(int)) )
        == NULL )
     printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntloc =realloc( ntloc, ntsets * maxtax * sizeof(int)) )
        == NULL )
     printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (nttotfig =realloc( nttotfig, ntsets * maxtax * sizeof(int)) )
        == NULL )
     printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntsigfig =realloc( ntsigfig, ntsets * maxtax * sizeof(int)) )
        == NULL )
     printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntrightfig =realloc( ntrightfig, ntsets * maxtax * sizeof(int)))
        == NULL )
     printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntval =realloc( ntval, ntsets * maxtax * sizeof(int)) )
        == NULL )
     printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (nterr =realloc( nterr, ntsets * maxtax * sizeof(int)) )
         == NULL )
     printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);
  if ( (ntoerr =realloc( ntoerr, ntsets * maxtax * sizeof(int)) )
         == NULL )
     printf( " NOT ENOUGH SPACE IN MEMORY FOR %d TAXA SETS\n", ntsets);

 }

/***********************************************************

 REALLOCATE SPACE FOR DEPTH

************************************************************/

 else if ( intime == 3 ) {

  if ( (ztotfig =realloc( ztotfig, zsize * sizeof(int)) )
       == NULL )
    printf( " #1 NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (zsigfig =realloc( zsigfig, zsize * sizeof(int)) )
       == NULL )
    printf( " #2 NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (zrightfig =realloc( zrightfig, zsize * sizeof(int)) )
       == NULL )
    printf( " #3 NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (depth =realloc( depth, zsize * sizeof(int)) )
       == NULL )
    printf( " #4 NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (zerr =realloc( zerr, zsize * sizeof(int)) )
       == NULL )
    printf( " #5 NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");
  if ( (zoerr =realloc( zoerr, zsize * sizeof(int)) )
       == NULL )
    printf( " #6 NOT ENOUGH SPACE IN MEMORY FOR DEPTH DATA\n");

 }

/************************************************************

 REALLOCATE SPACE FOR MEASURED VARIABLES

*************************************************************/

 else if ( intime == 4 ) {

  if ( (dtotfig =realloc( dtotfig, isize * sizeof(int)) )
       == NULL )
    printf(" NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA");
  if ( (dsigfig =realloc( dsigfig, isize * sizeof(int)) )
       == NULL )
    printf(" NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA");
  if ( (drightfig =realloc( drightfig, isize * sizeof(int)) )
       == NULL )
    printf(" NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA");
  if ( (dataval =realloc( dataval, isize * sizeof(int)) )
       == NULL )
    printf(" NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA");
  if ( (derr =realloc( derr, isize * sizeof(int)) )
       == NULL )
    printf(" NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA");
  if ( (doerr =realloc( doerr, isize * sizeof(int)) )
       == NULL )
    printf(" NOT ENOUGH SPACE IN MEMORY FOR MEASURED DATA");
 
 }

}
/**********************************************************

                 FUNCTION EXTRACTI

 EXTRACT1 EXTRACTS AN INTEGER VALUE FROM FILE.  THERE ARE THREE
 TYPES OF EXTRACTIONS.

 THE FIRST IS ACCOMPLISHED BY EXTRACTING THE NUMBER OF BYTES IN
 THE INTEGER VALUE AND THEN USING THIS INFORMATION
 TO EXTRACT THE INTEGER VALUE ITSELF.  THIS IS FOR INTEGERS USED
 INTERNALLY TO READ THE DATA.

 THE SECOND TYPE EXTRACTS NUMBER OF SIGNIFICANT FIGURES, NUMBER
 OF TOTAL FIGURES, AND NUMBER OF FIGURES TO THE RIGHT OF THE
 DECIMAL BEFORE USING NUMBER OF TOTAL FIGURES TO EXTRACT THE ACTUAL
 VALUE.  IF THE NUMBER OF SIGNIGICANT FIGURES IS NOT AN INTEGER,
 BUT A NEGATIVE SIGN (-), THE DATA VALUE IS A MISSING VALUE.

 THE THIRD TYPE EXTRACTS A VALUE USING AN INPUT NUMBER OF
 FIGURES (TOTFIG)

***********************************************************/

extracti(

 int type,                    /* TYPE OF EXTRACTION:
                                0=INTERNAL INTEGER VALUE
                                1=OUTPUT DATA VALUE 
                                2= OUTPUT DATA VALUE FROM
                                   GIVEN NUMBER OF FIGURES */

 int *totfig,                  /* NUMBER OF FIGURES IN THE
                                 ASCII REPRESENTATION OF THE VALUE 
                                 BEING EXTRACTED */

 int *sigfig,                  /* NUMBER OF SIGNIFICANT FIGURES IN
                                 THE VALUE BEING EXTRACTED */

 int *rightfig,                /* NUMBER OF PLACES TO THE RIGHT OF
                                 THE DECIMAL IN THE VALUE BEING
                                 EXTRACTED */

 int *value,                   /* ACTUAL VALUE BEING EXTRACTED, IN
                                 INTEGER FORM */

 int missing                  /* MISSING VALUE MARKER */

           )

{

 int sign,j,i;

/********************************************************

 SKIP IF THIS IS END OF LINE CHARACTER

*********************************************************/

 if ( type != 2 ) i=nocrfgetc();

/********************************************************

 IF THIS IS THE END OF FILE, SET EOF TO ONE TO NOTIFY
 MAIN PROGRAM

*********************************************************/

 if ( feof(fp) ) return -1;

 else {

/*******************************************************

 IF THIS IS A MISSING VALUE (I='-'), SET VALUE ACCORDINGLY.

********************************************************/

  if ( i == '-' ) {

   *value = missing;
   *totfig = 0;
   *sigfig = 0;
   *rightfig = 0;
   return 0;

  }

  else {

/*********************************************************

 ELSE READ IN AND/OR SET NUMBER OF SIGNIFICANT FIGURES,
 NUMBER OF TOTAL FIGURES, AND NUMBER OF FIGURES RIGHT
 OF THE DECIMAL IF A TYPE 1, AND NUMBER OF TOTAL FIGURES
 IF TYPE 0.

**********************************************************/

   if ( type == 1 ) {

    *sigfig = i - '0';

    i = nocrfgetc();
    *totfig= i - '0';

    i = nocrfgetc();
    *rightfig= i - '0';

   }

   else if ( type == 0 ) *totfig= i - '0'; 

/**********************************************************

 READ IN VALUE, INCLUDING SIGN

***********************************************************/

   *value= 0;
   for ( j = 1; j <= *totfig; j++ ) {

    i = nocrfgetc();

    if ( j > 1 ) *value= 10 * *value + ( i - '0' );
    else {
     sign = (i == '-') ? -1 : 1;
     if ( sign == 1 && i != ' ') *value = (i - '0');
    }

   }

   *value *= sign;

  }

 }

return 0;

}
 
/*************************************************************

                 FUNCTION EXTRACTC

 EXTRACTC EXTRACTS CHARACTER DATA FROM GIVEN FILE

**************************************************************/

extractc(

 int type,              /* SET TO ZERO IF TOTFIG SUPPLIED,
                           SET TO ONE IF TOTFIG READ IN */

 int *totfig,            /* NUMBER OF BYTES IN CHARACTER DATA */

 char *cdata            /* ARRAY OF CHARACTER DATA */

        )

{

 int i,j;

 if ( type !=0 ) i = nocrfgetc();
 
 if ( feof(fp) ) return -1;

 else {

  if ( type !=0 ) *totfig = i - '0';

  for ( j = 0; j < *totfig; j++ ) {

   i = nocrfgetc();
   *(cdata+j) = i;

  }

 }

 if ( feof(fp) ) return -1;
 else return 0;

}

/***********************************************

          FUNCTION NOCRFGETC

 NOCRFGETC (NO CARRIAGE RETURN FGETC) READS WOD FORMAT
 SKIPPING END OF LINE CHARACTERS IN A WAY WHICH WILL 
 WORK IN PC AND UNIX ENVIRONMENT.  THE PC END OF
 LINE CHARACTER (^M), WHICH IS PRESENT IN FILES ON
 THE CDs, IS FOLLOWED BY THE \n CHARACTER ON SOME
 UNIX PLATFORMS

 FP - FILE IDENTIFIER (GLOBAL VARIABLE)

 RETURNS INPUT CHARACTER OR -1 FOR END OF FILE

************************************************/

nocrfgetc()

{

 int i;

 while ( !feof(fp) && isprint ( (i=fgetc(fp)) ) == 0 );

 if ( feof(fp) ) i = -1;

 return i;

} 
