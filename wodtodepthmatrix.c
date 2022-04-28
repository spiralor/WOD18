#include <stdlib.h>
#include <stdio.h>

/***********************************************************

 DEFINED CONSTANTS

 MAXTAX: MAXIMUM NUMBER OF TAXA DELIMITORS IN A TAXA SET
 MAXSEC: MAXIMUM NUMBER OF SECONDARY HEADERS
 MAXBIO MAXIMUM NUMBER OF BIOLOGICAL HEADERS
 MAXPARM: MAXIMUM NUMBER OF MEASURED AND/OR CALCULATED VARIABLES
 MAXSECP: MAXIMUM NUMBER OF VARIABLE SPECIFIC SECONDARY HEADERS
 IDIM: LONGITUDE DIMENSIONS
 JDIM: LATITUDE DIMENSIONS
 KDIM: STANDARD DEPTH LEVELS

************************************************************/

#define maxtax 30
#define maxsec 100
#define maxbio 50
#define maxparm 100
#define maxpsec 25 * maxparm
#define idim 360
#define jdim 180
#define kdim 102
#define kdimax 137
#define maxchoice 8

/**********************************************************

 *FP - WOD DATA FILE TO BE OPENED FOR READING
 *FPOUT - FILE TO BE OPENED FOR WRITING OUTPUT
 *FPLIST - FILE WITH LIST OF WOA FILES TO READ
 *FPWOA - WOA FILES TO READ IN

***********************************************************/

FILE *fp,*fpout,*fplist, *fpwoa;

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

 NVARS - MAXIMUM NUMBER OF VARIABLES
 DVAL  - DATA VALUE AT EACH DEPTH LEVEL
 DERR - ERROR FLAG FOR EACH VARIABLE AT EACH LEVEL
 DOERR - ORIGINATORS FLAG FOR EACH VARIABLE AT EACH LEVEL
 DTOTFIG - TOTAL FIGURES FOR EACH DATA VALUE
 DSIGFIG - NUMBER OF SIGNIFICANT FIGURES FOR EACH DATA VALUE
 DRIGHTFIG - NUMBER OF FIGURES RIGHT OF THE DECIMAL FOR EACH DATA VALUE

***************************************************************/

int nvars=43;
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

/*************************************************************

 USER SET VARIABLES

 XCHOICE - CHOICE OF DEPTH SLICES
 XINTERV1,XINTERV2 - SHALLOWEST/DEEPEST DEPTHS REQUESTED

*************************************************************/

char xchoice[2];
float xinterv1,xinterv2;


/*************************************************************

 INTERNAL ARRAYS

 TENP - POWERS OF TEN
 SDEPTH - STANDARD DEPTH LEVELS

*************************************************************/

 double tenp[] = { 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000,
                  100000000, 1000000000 };

 int sdepth[] = { 0, 5 , 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70,
                  75, 80, 85, 90, 95, 100, 125, 150, 175, 200, 225, 250,
                  275, 300, 325, 350, 375, 400,425, 450, 475, 500, 550,
                  600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100,
                  1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600,
                  1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2100,
                  2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100,
                  3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 
                  4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900,
                  5000, 5100, 5200, 5300, 5400, 5500, 5600, 5700, 5800, 5900,
                  6000, 6100, 6200, 6300, 6400, 6500, 6600, 6700, 6800, 6900,
                  7000, 7100, 7200, 7300, 7400, 7500, 7600, 7700, 7800, 7900,
                  8000, 8100, 8200, 8300, 8400, 8500, 8600, 8700, 8800, 8900,
                  9000 };

 char *namevar[] = { "Temp","Sal","Oxy","Phos","dum5","Sil","dum7",
                     "NO3","pH","dum10","Chl","dum12","dum13","dum14",
                     "dum15","dum16","Alk","dum18","dum19","pCO2",
                     "DIC","dum22","dum23","BAC","dum25","dum26",
                     "dum27","dum28","dum29","dum30","dum31","dum32",
                     "Trit","He","dHE3","dC14","dC13","Arg","Neo","CFC11",
                     "CFC12","CFC113","O18" };


main()

{

 char filename[80];
 int i=0, j, k, s, jchoice, iend=0;
 int ncast=0,dchoice;

 printf(" Enter input file name\n");
 scanf("%s",filename);

 if ((fp = fopen(filename,"rb+\0")) == NULL)
  printf("UNABLE TO OPEN FILE\n");

 else {

  printf( "Which variable would you like to see:\n");
  printf( " 0 - All\n");
  printf( " 1 - Temperature\n");
  printf( " 2 - Salinity\n");
  printf( " 3 - Oxygen\n");
  printf( " 4 - Phosphate\n");
  printf( " 6 - Silicate\n");
  printf( " 8 - Nitrate\n");
  printf( "11 - Chlorophyll\n");
  printf( "17 - Alkalinity\n");
  printf( "20 - pCO2\n");
  printf( "21 - tCO2\n");
  printf( "24 - BAC\n");
  printf( "33 - Tritium\n");
  printf( "34 - Helium\n");
  printf( "35 - deltaHE3\n");
  printf( "36 - deltaC14\n");
  printf( "37 - deltaC13\n");
  printf( "38 - Argon\n");
  printf( "39 - Neon\n");
  printf( "40 - CFC11\n");
  printf( "41 - CFC12\n");
  printf( "42 - CFC113\n");
  printf( "43 - Oxy18\n");
  scanf("%d",&jchoice);

  printf( "Which depths would you like output:\n");
  printf( " A - All\n");
  printf( " S - Surface [Shallowest measurement only]\n");
  printf( " B - Bottom [Deepest measurement only]\n");
  printf( " I - Interval [measurements only between min/max depth]\n");
  scanf("%s",xchoice); 

  if ( *xchoice == 'I' ) {
   printf("Enter shallowest depth\n");
   scanf("%f",&xinterv1);
   printf("Enter deepest depth\n");
   scanf("%f",&xinterv2);
  }

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

  printheader(jchoice);

  while ( !feof(fp) && iend != -1 && (i++) < ncast ) {

/********************************************************

 READ IN CAST

**********************************************************/

   if ( ( iend = oclread() ) == -1 ) printf(" END OF FILE REACHED\n");

   else printstation(i,jchoice);

  }

  i = fclose(fp); 
  i = fclose(fpout); 
  printf("iii %d\n",i);

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

 READ IN WOD FORMAT CODE: 'C' FOR WOD13 FORMAT
 'B' FOR WOD05/WOD09 FORMAT
 'A' FOR WOD01 FORMAT

***********************************************************/

 totfig= 1;
 if ( (iend = extractc(0,&totfig,&wodform)) == -1 ) return iend;

 if ( wodform != 'C' ) {
  printf("file is not in WOD13 format\n");
  printf("csvfromwod cannot translate data\n");
  return -1;
 }

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

             FUNCTION PRINTHEAD

 PRINTHEAD PRINTS OUT FILE HEADER

***************************************************************/

printheader( int jchoice )

{

 int k;

 fprintf(fpout,"ISO_country,Cruise_ID,Latitude,Longitude,Year,Month,Day,");
 fprintf(fpout,"Time,WOD_unique,depth(m),qc_flag");
 if ( jchoice != 0 ) 
   fprintf(fpout,",%s,qc_flag\n",namevar[jchoice-1]);
 else
  {
   for ( k = 0; k < nvars; k++ ) {
    fprintf(fpout,",%s,qc_flag",namevar[k]);
   }
    fprintf(fpout,"\n");
  } 

}
 
/***************************************************************

              FUNCTION PRINTSTATION

 PRINTSTATION PRINTS CAST TO THE SCREEN

***************************************************************/

printstation(int i, int jchoice)

{

 int ii, ij, j, k, k0, offs,hastmp, yright,yall,zright;
 int lat,lon,itseas, ixa, inbp;
 int iwritten=0, k2,loopvars,j0;
 int xright[maxparm],xall[maxparm];
 float x[maxparm];
 int ichoice[maxchoice],ilevelwrite;
 int kblock= idim * jdim, iblock= kdim * kblock;
 float xlon, xlat, xhour, y, z, xa, spaced=1.+1.E-07;

 ichoice[0]=1;
 ichoice[1]=2;
 ichoice[2]=3;
 ichoice[3]=4;
 ichoice[4]=6;
 ichoice[5]=8;
 ichoice[6]=9;
 ichoice[7]=17;


/***************************************************************

 PRINT OUT HEADER

****************************************************************/

 xhour= (hour/ tenp[ *(hrightfig) ]);
 xlat= (latitude/ tenp[ *(hrightfig+1) ]);
 xlon= (longitude/ tenp[ *(hrightfig+2) ]);

 lat= ((xlat + 90.)/spaced)+1.;

 if ( xlon < 0. ) lon=((xlon+360.)/spaced)+1.;
 else lon=(xlon/spaced)+1.;
 itseas=((month-1)/3)+13;

 if ( levels > 0 ) {

  if ( jchoice > 0 ) {
   hastmp=-1;
   loopvars=0;
   for ( j = 0; j < nparm; j++ ) {
    if ( *(ip2+ j) == jchoice ) hastmp=j;
   }

  }
  else {

   hastmp=0;
   loopvars=nvars-1;

  }
  if ( hastmp > -1 ) {

   ilevelwrite=0;

   for ( k0 = 0; k0 < levels; k0++ ) {

    k=k0;
    if ( *xchoice == 'B' ) k= (levels - k0 - 1 );

    iwritten=1;
    if ( jchoice > 0 ) {
     iwritten=0;
     offs=hastmp * levels;
     if ( *(dsigfig+offs+k) > 0 ) iwritten=1;
    }

    zright = *(zrightfig+k);
    z = (*(depth+k)/ tenp[zright]);
    if ( zright > 6 ) zright=6;

    if ( *xchoice == 'B' ||
         *xchoice == 'S' ) {
     if ( ilevelwrite != 0 ) iwritten=0; 
    }
    else if ( *xchoice == 'I' ) {
     if ( z < xinterv1 ||
          z > xinterv2 ) iwritten=0;
    }
 
    if ( iwritten == 1 ) {
     ilevelwrite=1;

     fprintf(fpout,"%2s,%d,%.3f,%.3f,%4d,%2d,%2d,",
      cc,icruise,xlat,xlon,year,month,day);
     if ( xhour >= 0.0 && xhour <= 24.0 ) 
      fprintf(fpout,"%.2f,",xhour);
     else
      fprintf(fpout,",");
     fprintf(fpout,"%d,",ostation);
 
     fprintf(fpout,"%.*f",*(zrightfig+k),z);
     fprintf(fpout,",%d",*(zerr+k));

     for ( j = 0; j <= loopvars; j++ ) {

      if ( jchoice == 0 ) {
       hastmp=-1;
       for ( j0 = 0; j0 < nparm; j0++ ) {
        if ( *(ip2+ j0) == j+1 ) hastmp=j0;
       }
      }

      if ( hastmp > -1 ) {

       offs= hastmp * levels;
      
       yright=0;
       yall=0;

       yright = *(drightfig+offs+k);
       yall = *(dtotfig+offs+k);
       y = (*(dataval+offs+k)/ tenp[yright]);
       if ( yright > 6 ) yright=6;

       fprintf(fpout,",%.*f",yright,y);
       fprintf(fpout,",%d",*(derr+offs+k));

      }

      else 
       fprintf(fpout,",,");
     
     }

     fprintf(fpout,"\n");

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

  zsize= kdimax;
  zsizemax= kdimax;
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
 
  isize=kdimax;
  isizemax= kdimax;
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

extractf(

 int *totfig,                  /* NUMBER OF FIGURES IN THE
                                 ASCII REPRESENTATION OF THE VALUE 
                                 BEING EXTRACTED */

 float *value                    /* ACTUAL VALUE BEING EXTRACTED, IN
                                 FLOATING POINT FORM */

           )

{

 int sign,j,i;
 float value2;

/**********************************************************

 READ IN VALUE, INCLUDING SIGN

***********************************************************/

 *value= 0.0;
 value2=0.0;
 for ( j = 1; j <= *totfig; j++ ) {

  i = nocrfgetcwoa();
  if ( i == -1 ) {
   *value=-99.999;
   return -1;
  }

  if ( j > 1 && j <= 3 && i != ' ' )
    *value= 10 * *value + ( i - '0' );
  else if ( j > 4 && j < 8) value2= 10 * value2 + ( i - '0' );
  else if ( j == 1 ) {
   sign = (i == '-') ? -1 : 1;
   if ( sign == 1 && i != ' ') *value = (i - '0');
  }

 }

 *value += value2/1000.;
 *value *= sign;

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

nocrfgetcwoa()

{

 int i;

 while ( !feof(fpwoa) && isprint ( (i=fgetc(fpwoa)) ) == 0 );

 if ( feof(fpwoa) ) i = -1;

 return i;

} 
