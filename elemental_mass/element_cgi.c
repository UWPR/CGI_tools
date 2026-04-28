#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "AminoAcidMasses.h"

#define MAX_SEQUENCE 100000
#define MAX_CHARGE   10
#define PROTON_MASS  1.00727646688


void EXTRACT_QUERY_STRING(int *iC,
      int *iH,
      int *iN,
      int *iO,
      int *iP,
      int *iS);

extern char* read_POST(void);
extern void getword(char *word, char *line, char stop);
extern void plustospace(char *str);
extern void unescape_url(char *url);

int main(int argc, char **argv)
{
   int iMassType;

   FILE  *fp;

   double pdMassAA[128];

   int iC = 0;
   int iH = 0;
   int iN = 0;
   int iO = 0;
   int iP = 0;
   int iS = 0;
 
   printf("Content-type: text/html\n\n");


   // header
   if ( (fp=fopen("/net/pr/vol3/www/html/__header.php","r"))!=NULL)
   {
      char szBuf[1000];

      while (fgets(szBuf, 1000, fp))
      {
         // strip out php commands
         if (strstr(szBuf, "<?php"))
         {
            int i;

            for (i=0; i<(int)strlen(szBuf); i++)
            {
               if (!strncmp(szBuf+i, "<?php", 5))
               {
                  while (strncmp(szBuf+i, "?>", 2))
                     i++;
                  i += 1;
               }
               else
                  printf("%c", szBuf[i]);
            }

         }
         else
         {
            printf("%s", szBuf);
         }
      }
      fclose(fp);
   }
   printf("\n");

   printf("    <div id=\"page\" class=\"container\">\n");
   printf("       <section>\n");
   printf("          <header class=\"major\">\n");
   printf("             <h2>Element mass</h2>\n");
   printf("          </header>\n");


   fflush(stdout);


   iMassType = 1;  // 1=mono, 0=avg

   EXTRACT_QUERY_STRING(&iC, &iH, &iN, &iO, &iP, &iS);

   {
      const char *szScriptName = getenv("SCRIPT_NAME");
      printf("         <form action=\"%s\" name=\"fragmentForm\" method=\"post\">\n", szScriptName ? szScriptName : "");
   }

   printf("<br>C <input type=\"text\" name=\"C\" size=\"2\" value=\"%d\">\n", iC);
   printf("H <input type=\"text\" name=\"H\" size=\"2\" value=\"%d\">\n", iH);
   printf("N <input type=\"text\" name=\"N\" size=\"2\" value=\"%d\">\n", iN);
   printf("O <input type=\"text\" name=\"O\" size=\"2\" value=\"%d\">\n", iO);
   printf("P <input type=\"text\" name=\"P\" size=\"2\" value=\"%d\">\n", iP);
   printf("S <input type=\"text\" name=\"S\" size=\"2\" value=\"%d\">\n", iS);

   printf("\n&nbsp; &nbsp; <input type=\"submit\" value=\"Go\"><br>\n");

   printf("         </form>\n\n");
   printf("         <div id=\"results\">\n");

   if (1)
   {
      double dMass;

      INITIALIZE_MASS(pdMassAA, iMassType);

      dMass = iC * pdMassAA['c']
         +  iH * pdMassAA['h']
         +  iN * pdMassAA['n']
         +  iO * pdMassAA['o']
         +  iP * pdMassAA['p']
         +  iS * pdMassAA['s'];

      printf("<br><tt>monoisotopic neutral mass:  %0.6f", dMass);
      printf("<br><br>1+: %0.6f\n", (dMass + 1*PROTON_MASS)/1);
      printf("<br>2+: %0.6f\n", (dMass + 2*PROTON_MASS)/2);
      printf("<br>3+: %0.6f\n", (dMass + 3*PROTON_MASS)/3);
      printf("<br>4+: %0.6f\n", (dMass + 4*PROTON_MASS)/4);
      printf("<br>5+: %0.6f\n", (dMass + 5*PROTON_MASS)/5);

      printf("<br><br>");
      if (iC > 0)
         printf("C%d:  %f<br>\n", iC, iC * pdMassAA['c']);
      if (iH > 0)
         printf("H%d:  %f<br>\n", iH, iH * pdMassAA['h']);
      if (iN > 0)
         printf("N%d:  %f<br>\n", iN, iN * pdMassAA['n']);
      if (iO > 0)
         printf("O%d:  %f<br>\n", iO, iO * pdMassAA['o']);
      if (iP > 0)
         printf("P%d:  %f<br>\n", iP, iP * pdMassAA['p']);
      if (iS > 0)
         printf("S%d:  %f<br>\n", iS, iS * pdMassAA['s']);
   }
   printf("</div>\n");
   printf("</div>\n");
   printf("</div>\n");

   // footer
   if ( (fp=fopen("/net/pr/vol3/www/html/__footer.php","r"))!=NULL)
   {
      char szBuf[1000];

      while (fgets(szBuf, 1000, fp))
      {
         // strip out php commands
         if (strstr(szBuf, "<?php"))
         {
            int i;

            for (i=0; i<(int)strlen(szBuf); i++)
            {
               if (!strncmp(szBuf+i, "<?php", 5))
               {
                  while (strncmp(szBuf+i, "?>", 2))
                     i++;
                  i += 1;
               }
               else
                  printf("%c", szBuf[i]);
            }

         }
         else
         {
            printf("%s", szBuf);
         }
      }
      fclose(fp);
   }

   exit(EXIT_SUCCESS);
}


void EXTRACT_QUERY_STRING(int *iC,
      int *iH,
      int *iN,
      int *iO,
      int *iP,
      int *iS)
{
   char *pStr = getenv("REQUEST_METHOD");

   if (pStr==NULL)
   {
      printf(" Error - this is a CGI program that cannot be\n");
      printf(" run from the command line.\n\n");
      exit(EXIT_FAILURE);
   }
   else  // get or post
   {
      int  i;
      int  iContentLength=0;
      char *szQuery,
           *szWord;

      if (!strcmp(pStr, "GET"))
      {
         szQuery = getenv("QUERY_STRING");
         if (szQuery == NULL)
            szQuery = (char*)"";
         iContentLength = strlen(szQuery);
      }
      else
      {
         szQuery = read_POST();
         if (szQuery == NULL)
            szQuery = (char*)"";
         iContentLength = strlen(szQuery);
      }

      if (strlen(szQuery)>0)
      {
         /* +1 so getword can always append a null terminator within the buffer */
         if ((szWord = (char*)malloc(iContentLength + 1))==NULL)
         {
            printf("<P>Error, cannot malloc szWord (size %d).\n", iContentLength);
            printf("</BODY>\n</HTML>\n");
            exit(EXIT_FAILURE);
         }

         szQuery[iContentLength]='\0';

         for (i=0; szQuery[0] != '\0'; i++)
         {
            getword(szWord, szQuery, '=');
            plustospace(szWord);
            unescape_url(szWord);

            if (!strcmp(szWord, "C") )
            {
               int iVal;
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               sscanf(szWord, "%d", &iVal);
               *iC = iVal;
            }
            else if (!strcmp(szWord, "H") )
            {
               int iVal;
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               sscanf(szWord, "%d", &iVal);
               *iH = iVal;
            }
            else if (!strcmp(szWord, "N") )
            {
               int iVal;
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               sscanf(szWord, "%d", &iVal);
               *iN = iVal;
            }
            else if (!strcmp(szWord, "O") )
            {
               int iVal;
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               sscanf(szWord, "%d", &iVal);
               *iO = iVal;
            }
            else if (!strcmp(szWord, "P") )
            {
               int iVal;
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               sscanf(szWord, "%d", &iVal);
               *iP = iVal;
            }
            else if (!strcmp(szWord, "S") )
            {
               int iVal;
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               sscanf(szWord, "%d", &iVal);
               *iS = iVal;
            }
            else
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
            }
         }
      }
   }
}
