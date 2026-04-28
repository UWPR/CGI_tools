#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* HTML-encode a string to stdout, escaping chars that are special in HTML attributes/content. */
static void print_html_encoded(const char *s)
{
   for (; *s; s++)
   {
      switch (*s)
      {
         case '&':  fputs("&amp;",  stdout); break;
         case '<':  fputs("&lt;",   stdout); break;
         case '>':  fputs("&gt;",   stdout); break;
         case '"':  fputs("&quot;", stdout); break;
         case '\'': fputs("&#39;",  stdout); break;
         default:   fputc(*s, stdout);       break;
      }
   }
}

/* Enzyme cut/no-cut strings may only contain amino-acid letters and '-'. */
static int is_valid_enzyme_str(const char *s)
{
   for (; *s; s++)
      if (!isalpha((unsigned char)*s) && *s != '-')
         return 0;
   return 1;
}

#define RULESFILE "/net/pr/vol3/www/cgi-bin/digest.cgi.rules"
#define DIGESTBINARY "/net/pr/vol1/ProteomicsResource/bin/digestdb_for_cgi"

struct ProteaseStruct
{
   char szName[96];
   char szCut[24];
   char szNoCut[24];
   int  iSense;    /* 0=n-term, 1=c-term */
} *pProtease, *pTmp;


void EXTRACT_QUERY_STRING(int *iWhichInput,
      int *iEnzymeNum,
      int *iVal,
      int *iCustomEnzymeSense,
      char *szCustomEnzymeCut,
      char *szCustomEnzymeNoCut,
      char *szInputSequence,
      double *dMinMass,
      double *dMaxMass,
      int *iNumAllowedMissed,
      int *iMassType);
struct ProteaseStruct* INITIALIZE_PROTEASE(int *iNumEnzymes,
      int *iNumAllocated);

extern char* read_POST(void);
extern void getword(char *word, char *line, char stop);
extern void plustospace(char *str);
extern void unescape_url(char *url);

int main(int argc, char **argv)
{
   int i;
   int iVal;
   int iNumEnzymes;
   int iNumAllocated;
   int iMassType;
   int iWhichInput;        /* 0=drop down form, 1=custom */
   double dMinMass;
   double dMaxMass;
   int iNumAllowedMissed;
   int iEnzymeNum;
   char szInputSequence[100000];
   FILE  *fp;
   int iCustomEnzymeSense;
   char szCustomEnzymeCut[48];
   char szCustomEnzymeNoCut[48];
   int iEnzymeSense;       /* these hold actual enzyme rules passed to digestdb */
   char szEnzymeCut[48];   /* these hold actual enzyme rules passed to digestdb */
   char szEnzymeNoCut[48]; /* these hold actual enzyme rules passed to digestdb */
 
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
         else if (strstr(szBuf, "<body>"))
         {
            printf("   <style>\n");
            printf("      table, th, td {width: 1%%; font-family:\"Courier New\", Courier, monospace; font-size:90%%; border: 1px solid #ADD8E6; padding-left: 0.4rem; padding-right: 0.4rem; text-align: center}\n");
            printf("      th {background-color: #F0FFFF}\n");
            printf("      a:link {text-decoration: none}\n");
            printf("   </style>\n");
            printf("%s", szBuf);
         }
         else
         {
            if (strstr(szBuf, "</head>"))
               printf("     <script type=\"text/javascript\" src=\"/js/sorttable.js\"></script>\n");

            printf("%s", szBuf);
         }
      }
      fclose(fp);
   }
   printf("\n");

   printf("    <div id=\"page\" class=\"container\">\n");
   printf("       <section>\n");
   printf("          <header class=\"major\">\n");
   printf("             <h1>Protein digestion</h1>\n");
   printf("          </header>\n");


   printf("  <script language=\"javascript\">\n");
   printf("     function pasteExample() {\n");
   printf("       document.getElementById('protein').value='MASFRLFLLCLAGLVFVSEAGSVGAGEPKCPLMVKVLDAV \
RGSPAANVGVKVFKKAADETWEPFASGKTSESGELHGLTT \
EDKFVEGLYKVELDTKSYWKSLGISPFHEFAEVVFTANDS \
GPRHYTIAALLSPYSYSTTALVSSPKA'\n");
   printf("     }\n");
   printf("\n");
   printf("     function clearExample() {\n");
   printf("       document.getElementById('protein').value=''\n");
   printf("     }\n");
   printf("\n");
   printf("     function whichInput() {\n");
   printf("        document.getElementById('whichinput').value=0;\n");
   printf("        document.getElementById('whichinput').checked=true;\n");
   printf("     }\n");
   printf("     function whichInput2() {\n");
   printf("        document.getElementById('whichinput2').value=1;\n");
   printf("        document.getElementById('whichinput2').checked=true;\n");
   printf("     }\n");
   printf("  </script>\n\n");

   iVal=0;
   strcpy(szInputSequence, "");

   dMinMass=50.0;
   dMaxMass=8000.0;
   iNumAllowedMissed=1;
   iMassType=1;  /*1=mono, 0=avg*/
   iEnzymeNum = 0;
   iWhichInput = 0;
   strcpy(szCustomEnzymeCut, "K");
   strcpy(szCustomEnzymeNoCut, "-");
   iCustomEnzymeSense=1;
   EXTRACT_QUERY_STRING(&iWhichInput, &iEnzymeNum, &iVal, &iCustomEnzymeSense,
         szCustomEnzymeCut, szCustomEnzymeNoCut, szInputSequence,
         &dMinMass, &dMaxMass, &iNumAllowedMissed, &iMassType);

   if (strlen(szCustomEnzymeCut)==0)
      strcpy(szCustomEnzymeCut, "-");
   if (strlen(szCustomEnzymeNoCut)==0)
      strcpy(szCustomEnzymeNoCut, "-");

   iNumAllocated = 0;
   iNumEnzymes = 0;
   pProtease = INITIALIZE_PROTEASE(&iNumEnzymes, &iNumAllocated);

   {
      const char *szScriptName = getenv("SCRIPT_NAME");
      printf("<form action=\"%s\" name=\"digestForm\" method=\"post\">", szScriptName ? szScriptName : "");
   }

   printf("<div id=\"left40entry\"><b>Select enzyme OR create own rule:</b><br>\n");
 
   printf("<input type=\"radio\" aria-label=\"select-enzyme\" name=\"whichinput\" id=\"whichinput\" value=\"0\"%s>\n", (iWhichInput==0?" checked":""));
   printf("<select name=\"enzyme\" style=\"width: 420px\" onfocus=\"whichInput()\">\n");
   iEnzymeSense=0;
   for (i=0; i<iNumEnzymes; i++)
   {
      char szSelected[50];
      char szDontCleave[512];

      szSelected[0]='\0';
      if (i==iEnzymeNum)
      {
         strcpy(szSelected, " selected");
         if (iWhichInput==0)
         {
            strcpy(szEnzymeCut, pProtease[i].szCut);
            strcpy(szEnzymeNoCut, pProtease[i].szNoCut);
            iEnzymeSense = pProtease[i].iSense;
         }
      }

      if (strstr(pProtease[i].szNoCut, "-"))
         szDontCleave[0] = '\0';

      else
         sprintf(szDontCleave, ", don't cleave \"%s\"", pProtease[i].szNoCut);

      printf("   <option value=\"%d\"%s>%s, %s, cleave \"%s\"%s\n",
            i,
            szSelected,
            pProtease[i].szName,
            pProtease[i].iSense==0 ? "C-term" : "N-term",
            pProtease[i].szCut,
            szDontCleave);
   }

   printf("</select>\n");

   if (iWhichInput==1)
   {
      strcpy(szEnzymeCut, szCustomEnzymeCut);
      strcpy(szEnzymeNoCut, szCustomEnzymeNoCut);
      iEnzymeSense = iCustomEnzymeSense;
   }

   printf("<br><input type=\"radio\" aria-label=\"select-custom-enzyme\" name=\"whichinput\" id=\"whichinput2\" value=\"1\"%s>\n", (iWhichInput==1?" checked":""));
   printf("cut:<input type=\"text\" aria-label=\"enzymecut\" name=\"enzymecut\" size=\"2\" value=\"");
   print_html_encoded(szCustomEnzymeCut);
   printf("\" onfocus=\"whichInput2()\"> ");
   printf("nocut:<input type=\"text\" aria-label=\"enzymenocut\" name=\"enzymenocut\" size=\"1\" value=\"");
   print_html_encoded(szCustomEnzymeNoCut);
   printf("\"onfocus=\"whichInput2()\">\n");
/*
   printf("sense:<input type=\"radio\" name=\"enzymesense\" onfocus=\"whichInput2()\" value=\"0\"%s>n", (iCustomEnzymeSense==0?" checked":""));
   printf("<input type=\"radio\" name=\"enzymesense\" onfocus=\"whichInput2()\" value=\"1\"%s>c\n", (iCustomEnzymeSense==1?" checked":""));
*/

   printf("<p><p><b>Output mass range:</b>\n");
   printf("<input type=\"text\" aria-label=\"minmass\" name=\"minmass\" size=\"5\" value=\"%0.2f\">-", dMinMass);
   printf("<input type=\"text\" aria-label=\"maxmass\" name=\"maxmass\" size=\"5\" value=\"%0.2f\">\n", dMaxMass);

   printf("<br><b>Allowed missed cleavages:</b>\n");
   printf("<select name=\"allowedmissed\" style=\"width: 40px\">\n");
   for (i=0; i<=5; i++)
      printf("   <option value=\"%d\"%s>%d\n", i, (iNumAllowedMissed==i?" selected":""), i);
   printf("</select>\n");

   printf("<br><input type=\"radio\" aria-label=\"monomass\" name=\"masstype\" value=\"1\"%s>mono ",
         (iMassType==1?" checked":""));
   printf("<input type=\"radio\" aria-label=\"avgmass\" name=\"masstype\" value=\"0\"%s>avg masses<BR>",
         (iMassType==0?" checked":""));
   printf("</div>\n");


   printf("<div id=\"right60entry\">\n");
   printf("<b>Paste in protein sequence here (sequence itself only please):</b><br>\n");
   printf("<textarea name=\"sequence\" rows=\"3\" cols=\"80\" wrap=\"virtual\" id=\"protein\">");
   print_html_encoded(szInputSequence);
   printf("</textarea>\n");
   printf("<br><a onclick=\"pasteExample();\"><u><font size=\"-2\">Click here to paste in a sample sequence</font></u></a>\n");
   printf(" &nbsp; &nbsp; &nbsp; <a onclick=\"clearExample();\"><u><font size=\"-2\">Clear sequence</font></u></a>\n");
   printf("<br><input type=\"submit\" aria-label=\"performdigest\" value=\"Digest!\"><br>\n");

   printf("</div>\n\n");
 
   printf("</form>\n\n");

   printf("<div id=\"results\">\n");
   if (strlen(szInputSequence)>0)
   {
      char szCommand[1024];
      FILE *fp;
      char szTmpFile[512];
      char szBuf[5000];

      for (i=0; i<(int)strlen(szInputSequence); i++)
         szInputSequence[i]=toupper(szInputSequence[i]);

      sprintf(szTmpFile, "/tmp/digestdb.fasta");
//    sprintf(szTmpFile, "/net/pr/vol3/software/digestdb.fasta");
      if ( (fp=fopen(szTmpFile, "w"))==NULL)
      {
         printf(" Error writing tmp db\n");
         exit(EXIT_FAILURE);
      }
      fprintf(fp, ">seq\n%s", szInputSequence);
      fclose(fp);

//    sprintf(szCommand, "/usr/bin/rm -f %s.out", szTmpFile);
//    system(szCommand);
      if (!is_valid_enzyme_str(szEnzymeCut) || !is_valid_enzyme_str(szEnzymeNoCut))
      {
         printf(" Error - invalid characters in enzyme specification.\n");
         exit(EXIT_FAILURE);
      }

      sprintf(szCommand, "%s -r%s -n%s -l%0.6f -h%0.6f -m%d -t%d -d%d %s | /usr/bin/sort -k7,7n -k8,8n",
            DIGESTBINARY, szEnzymeCut, szEnzymeNoCut, dMinMass, dMaxMass, iNumAllowedMissed, iMassType, iEnzymeSense, szTmpFile);

      //printf("<p>DEBUG: running command<br><tt>%s</tt><p>\n", szCommand);

//    system(szCommand);
//    strcat(szTmpFile, ".out");

//      printf("<p>DEBUG: display contents of file %s:<p>\n", szTmpFile);
      if ((fp = popen(szCommand, "r")) == NULL)
      {
         printf(" Error popen %s\n", szCommand);
         exit(EXIT_FAILURE);
      }

/*
      while (fgets(szBuf, 5000, fp))
         printf("<p>szBuf %s  len %d\n", szBuf, strlen(szBuf));
*/


      printf("<br><center><p style='font-family:\"Courier New\", Courier, monospace; font-size:90%%;'>Click on column headers to sort results. &#9830; is peptide fragmentation link.</p>\n");

      printf("<table class=\"sortable\" cellpadding=\"4\">\n");
      printf("<thead>");
      printf("<tr>");
      printf("<th>start</th>");
      printf("<th>end</th>");
      printf("<th>len</th>");
      printf("<th>neutral_mass</th>");
      printf("<th>prevAA</th>");
      printf("<th>peptide</th>");
      printf("<th>nextAA</th>");
      printf("<th>pI</th>");
      printf("</tr>");
      printf("</thead>\n");
      printf("<tbody>\n");

      while (fgets(szBuf, 5000, fp))
      {
         int iLen;
         int iStart;
         int iEnd;
         char szTmp[50];
         double dMass;
         double dPI;
         char cPre;
         char cPost;
         char szPep[5000];

         sscanf(szBuf, "%d %s %lf %c %s %c %d %d %lf\n",
               &iLen, szTmp, &dMass, &cPre, szPep, &cPost, &iStart, &iEnd, &dPI);
         printf("<tr>");
         printf("<td>%d</td>", iStart+1);
         printf("<td>%d</td>", iEnd+1);
         printf("<td>%d</td>", iLen);
         printf("<td style='text-align: right'>%0.6lf</td>", dMass - 1.00727646688);
         printf("<td>%c</td>", cPre);

         printf("<td style='text-align: left'>%s", szPep);
         printf("<a href=\"/cgi-bin/fragment.cgi?sequence=%s\">&#9830;</a>", szPep);
         printf("</td>");

         printf("<td>%c</td>", cPost);
         printf("<td>%0.2lf</td>", dPI);
         printf("</tr>\n");

//       printf("fget: %s<br>\n", szBuf);
      }

      printf("</tbody>\n");
      printf("</table></center>\n");
      pclose(fp);
   }

   printf("</div>\n");
   printf("</div>\n");
   printf("</div>\n");

   if ( (fp=fopen("/net/pr/vol3/www/html/__footer.php","r"))!=NULL)
   {
      char szBuf[1000];

      while (fgets(szBuf, 1000, fp))
      {
         printf("%s", szBuf);
      }
      fclose(fp);
   }

   exit(EXIT_SUCCESS);
} /*main*/


void EXTRACT_QUERY_STRING(int *iWhichInput,
      int *iEnzymeNum,
      int *iVal,
      int *iCustomEnzymeSense,
      char *szCustomEnzymeCut,
      char *szCustomEnzymeNoCut,
      char *szInputSequence,
      double *dMinMass,
      double *dMaxMass,
      int *iNumAllowedMissed,
      int *iMassType)
{
   char *pStr = getenv("REQUEST_METHOD");

   if (pStr==NULL)
   {
      printf(" Error - this is a CGI program that cannot be\n");
      printf(" run from the command line.\n\n");
      exit(EXIT_FAILURE);
   }
   else  /* get or post */
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
         if ((szWord = (char *)malloc(iContentLength + 1))==NULL)
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

            if (!strcmp(szWord, "i"))
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               sscanf(szWord, "%d", iVal);
            }
            else if (!strcmp(szWord, "sequence") )
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               strncpy(szInputSequence, szWord, 99999);
               szInputSequence[99999] = '\0';
            }
            else if (!strcmp(szWord, "whichinput") )
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               *iWhichInput= atoi(szWord);
            }
            else if (!strcmp(szWord, "enzyme") )
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               *iEnzymeNum = atoi(szWord);
            }
            else if (!strcmp(szWord, "enzymesense") )
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               sscanf(szWord, "%d", iCustomEnzymeSense);
            }
            else if (!strcmp(szWord, "enzymecut") )
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               strncpy(szCustomEnzymeCut, szWord, 47);
               szCustomEnzymeCut[47] = '\0';
            }
            else if (!strcmp(szWord, "enzymenocut") )
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               strncpy(szCustomEnzymeNoCut, szWord, 47);
               szCustomEnzymeNoCut[47] = '\0';
            }
            else if (!strcmp(szWord, "minmass") )
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               sscanf(szWord, "%lf", dMinMass);
            }
            else if (!strcmp(szWord, "maxmass") )
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               sscanf(szWord, "%lf", dMaxMass);
            }
            else if (!strcmp(szWord, "masstype") )
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               sscanf(szWord, "%d", iMassType);
               if (*iMassType < 0 || *iMassType>1)
                  *iMassType = 1;
            }
            else if (!strcmp(szWord, "allowedmissed") )
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               sscanf(szWord, "%d", iNumAllowedMissed);
               if (*iNumAllowedMissed <0)
                  *iNumAllowedMissed=0;
               if (*iNumAllowedMissed >5)
                  *iNumAllowedMissed=5;
            }
            else
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
            }
         } /*for*/
      }
   }
} /*EXTRACT_QUERY_STRING*/


/*
 * Protease definitions stored in cgi-bin/digest.cgi.rules
 */
struct ProteaseStruct* INITIALIZE_PROTEASE(int *iNumEnzymes,
      int *iNumAllocated) 
{
   FILE *fp;
   char szBuf[5000];
   struct ProteaseStruct *pTmp;

   *iNumAllocated = 15;
   pTmp = (struct ProteaseStruct *)malloc(sizeof(struct ProteaseStruct) * (*iNumAllocated));

   if (pTmp==NULL)
   {
      printf(" Error malloc pTmp\n");
      exit(EXIT_FAILURE);
   }

   if ( (fp=fopen(RULESFILE, "r"))==NULL)
   {
      printf(" Error opening enzyme config file %s\n\n", RULESFILE);
      exit(EXIT_FAILURE);
   }

   while (fgets(szBuf, 5000, fp))
   {
      int iSense;
      char szCut[24];
      char szNoCut[24];
      char szName[96];

      if (szBuf[0]!='#')
      {
         sscanf(szBuf, "%95s %d %23s %23s", szName, &iSense, szCut, szNoCut);
   
         if (*iNumEnzymes == *iNumAllocated)
         {
            struct ProteaseStruct *pTmp2;
            *iNumAllocated += 10;
            pTmp2 = (struct ProteaseStruct *)realloc(pTmp, (*iNumAllocated)*sizeof(struct ProteaseStruct));
            if (pTmp2 == NULL)
            {
               printf(" Error realloc\n\n");
               exit(EXIT_FAILURE);
            }
            pTmp = pTmp2;
         }
   
         strcpy((pTmp + *iNumEnzymes)->szName, szName);
         strcpy((pTmp + *iNumEnzymes)->szCut, szCut);
         strcpy((pTmp + *iNumEnzymes)->szNoCut, szNoCut);
         (pTmp + *iNumEnzymes)->iSense = iSense;
   
         (*iNumEnzymes)++;
      }
   }

   fclose(fp);

   return pTmp;

} /*INITIALIZE*/
