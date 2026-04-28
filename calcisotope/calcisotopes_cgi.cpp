//
// USAGE:  ./calcisotopes "K[325.13]K[170.11]FTENPKAG"
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "isotopes.h"

/* HTML-encode a string to stdout, escaping chars special in HTML attributes/content. */
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

#define PROTON_MASS      1.00727646688
#define MAX_SEQUENCE     512
#define SAMPLEPEPTIDE    "DIGYSESTEDQAMEDIK"

int piCompC[MAX_SEQUENCE];
int piCompH[MAX_SEQUENCE];
int piCompN[MAX_SEQUENCE];
int piCompO[MAX_SEQUENCE];
int piCompS[MAX_SEQUENCE];

void INIT_COMP(int *piCompC,
      int *piCompH,
      int *piCompN,
      int *piCompO,
      int *piCompS);
void EXTRACT_QUERY_STRING(char *szInputSequence,
                          int *iCharge);

extern char* read_POST(void);
extern void getword(char *word, char *line, char stop);
extern void plustospace(char *str);
extern void unescape_url(char *url);


int main(int argc, char **argv)
{
   FILE *fp;

   int iC = 0;
   int iH = 0;
   int iN = 0;
   int iO = 0;
   int iS = 0;

   int iCharge = 1;

   // sets element counts for each residue
   INIT_COMP(piCompC, piCompH, piCompN, piCompO, piCompS);

   char szComp[128];
   char szInput[1024];

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
   printf("             <h2>Peptide Isotope Calculator</h2>\n");
   printf("          </header>\n");


   printf("  <script language=\"javascript\">\n");
   printf("     function pasteExample() {\n");
   printf("       document.getElementById('peptide').value='%s'\n", SAMPLEPEPTIDE);
   printf("     }\n");
   printf("\n");
   printf("     function clearExample() {\n");
   printf("       document.getElementById('peptide').value=''\n");
   printf("     }\n");
   printf("\n");
   printf("     function whichInput() {\n");
   printf("        document.getElementById('whichinput').value=0;\n");
   printf("        document.getElementById('whichinput').checked=true;\n");
   printf("     }\n");
   printf("\n");
   printf("     function whichInput2() {\n");
   printf("        document.getElementById('whichinput2').value=1;\n");
   printf("        document.getElementById('whichinput2').checked=true;\n");
   printf("     }\n");
   printf("  </script>\n\n");

   fflush(stdout);


   szInput[0] = 0;

   EXTRACT_QUERY_STRING(szInput, &iCharge);

   {
      const char *szScriptName = getenv("SCRIPT_NAME");
      printf("         <form action=\"%s\" name=\"fragmentForm\" method=\"post\">\n", szScriptName ? szScriptName : "");
   }
   printf("\n");
   printf("         <div id=\"left40entry\">");
   printf("         This program calculates the isotope distribution of peptides using");
   printf("         the <a href=\"http://www.kombyonyx.com/isotopes/\">Isotope Distribution Calculator</a> code by James Redman at Cardiff Unversity.\n");
   printf("         It is described as an algorithm modeled after <a href=\"https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/jms.4498\">J. A. Yergey's procedure published in 1983</a>.\n");
   printf("         </div>\n");
   printf("         <p>\n");
   printf("         <div id=\"left40entry\"><b>Enter sequence here:</b>\n");
   printf("            <br> <a onclick=\"clearExample();\"><u><font size=\"-2\">Clear sequence</font></u></a>\n");
   printf("            &nbsp; &nbsp; &nbsp;<a onclick=\"pasteExample();\"><u><font size=\"-2\">Click here to paste in a sample peptide</font></u></a>\n");
   printf("            <br>peptide: <input type=\"text\" name=\"sequence\" id=\"peptide\" size=\"40\" value=\"");
   print_html_encoded(szInput);
   printf("\">\n");
   printf("            <br>charge: <input type=\"text\" name=\"charge\" id=\"charge\" size=\"2\" value=\"%d\"> ... precursor charge state is used only in the calculation of m/z values\n", iCharge);
   printf("            <p><br>&nbsp; <input type=\"submit\" value=\"Calculate\">\n");
   printf("         </div>\n");
   printf("\n");
   printf("         </form>\n\n");

   printf("         <div id=\"results\">\n");

   if (strlen(szInput)>0 && isalpha(szInput[0]))
   {
      // peptide specified on command line
      char cRes;
      int iLen = strlen(szInput);

      vector<double> vLabel;
      vector<double> vAbun;

      iH = 3;
      iO = 1;

      for (int i=0; i<iLen; i++)
      {
         cRes = szInput[i];
         int idx = (int)(unsigned char)cRes;  /* prevent negative index from signed char */

         iC += piCompC[idx];
         iH += piCompH[idx];
         iN += piCompN[idx];
         iO += piCompO[idx];
         iS += piCompS[idx];
      }

      szComp[0] = '\0';
      if (iC>0)
         sprintf(szComp+strlen(szComp), "C%d ", iC);
      if (iH>0)
         sprintf(szComp+strlen(szComp), "H%d ", iH);
      if (iN>0)
         sprintf(szComp+strlen(szComp), "N%d ", iN);
      if (iO>0)
         sprintf(szComp+strlen(szComp), "O%d ", iO);
      if (iS>0)
         sprintf(szComp+strlen(szComp), "S%d ", iS);

      int errnr, prec = 8, z;
      bool  masstocharge = false;
      double thr = 0.0001, res = 0.01;
      IsoCalc mycalc;

      masstocharge = true;     //choose whether to output mass or mass-to-charge ratio
      z = iCharge;

      errnr = mycalc.ReadAtomTable("/net/pr/vol1/ProteomicsResource/bin/isotopestable.txt");
      if (errnr)
      {
         cout << "Error: cannot read isotope table: /net/pr/vol1/ProteomicsResource/bin/isotopetable.txt\n";
         exit(1);
      }

      errnr = mycalc.SetComposition(szComp);
      if (errnr)
      {
         cout << "Error: cannot parse formula\n";
         exit(1);
      }
      mycalc.SetCharge(z);
      mycalc.SetThr(thr);
      mycalc.SetMassToCharge(masstocharge);
      mycalc.SetDegen(res);
      cout.precision(prec);     // set the precision
      errnr = mycalc.Calculate();
      if (errnr)
      {
         cout << "Error: cannot calculate distributions\n";
         exit(1);
      }                         // this should never happen - but just in case

      mycalc.Normalize();       // normalize max intensity to 100

      int npeaks=0;
      double mass, abun;
      mycalc.GetNPeaks(npeaks);

      printf("<p style='font-family: monospace; font-size:80%%'>\n");
      printf("sequence: ");
      print_html_encoded(szInput);
      printf("\n");
      printf("<br>composition:  %s\n", szComp);

      printf("<table style='font-family: monospace; font-size:80%%; width:50%%'>\n");
      printf("<tr><td>peak</td><td>m/z</td><td>relative abundance</td></tr>\n");

      double dMaxAbundance = 0;
      for (int i = 0; i < npeaks; i++)
      {
         mycalc.Peak(i, mass, abun);
         if (abun > dMaxAbundance)
            dMaxAbundance = abun;
      }
      for (int i = 0; i < npeaks; i++)
      {
         mycalc.Peak(i, mass, abun);

         printf("<tr>\n");
         if (i==0)
            printf("<td>Mono</td><td>%0.5f</td><td>%0.2f</td><br>", mass, 100.0 * abun / dMaxAbundance);
         else
            printf("<td>M+%d</td><td>%0.5f</td><td>%0.2f</td>", i,  mass, 100.0 * abun / dMaxAbundance);

         vLabel.push_back(mass);
         vAbun.push_back(100.0 * abun / dMaxAbundance);

         printf("</tr>\n");
      }
      printf("</table>\n\n");


      // https://www.chartjs.org/docs/latest/getting-started/
      printf("   <div class=\"chart-container\" style=\"width:600px\"><canvas id=\"myChart\"></canvas></div>\n");
      printf("<script src=\"https://cdn.jsdelivr.net/npm/chart.js\"></script>\
<script>\n\
  const ctx = document.getElementById('myChart');\n\
  new Chart(ctx, {\n\
    type: 'bar',\n\
    data: {\n\
      labels: [");

      for (auto it=vLabel.begin(); it!=vLabel.end(); ++it)
      {
         if (it != vLabel.begin())
            printf(", ");
         printf("'%0.4lf'", (*it));
      }

      printf("],\n\
      datasets: [{\n\
        label: 'relative isotope distribution',\n\
        data: [");

      for (auto it=vAbun.begin(); it!=vAbun.end(); ++it)
      {
         if (it != vAbun.begin())
            printf(", ");
         printf("'%0.4lf'", (*it));
      }

      printf("],\n        borderWidth: 1\n\
      }]\n\
    },\n\
    options: {\n\
      scales: {\n\
        y: { beginAtZero: true, grid: { display:false }  },\n\
        x: { display: true, text: 'm/z',  grid: { display:false }},\n\
      },\n\
      barThickness: 5,\n\
      plugins: { legend: { display: false } }\n\
    }\n\
  });\n\
</script>\n");

   }

   if ( (fp=fopen("/net/pr/vol3/www/html/__footer.php","r"))!=NULL)
   {
      char szBuf[1000];

      while (fgets(szBuf, 1000, fp))
      {
         // strip out php commands
         if (strstr(szBuf, "<?php"))
         {
            int i;
            int iLen = strlen(szBuf);

            for (i=0; i<iLen; i++)
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


void INIT_COMP(int *piCompC,
      int *piCompH,
      int *piCompN,
      int *piCompO,
      int *piCompS)
{
   int i;

   for (i=0; i<128; i++)
   {  
      piCompC[i]=0;
      piCompH[i]=0;
      piCompN[i]=0;
      piCompO[i]=0;
      piCompS[i]=0;
   }

   piCompC['G'] = 2  ;
   piCompC['A'] = 3  ;
   piCompC['S'] = 3  ;
   piCompC['P'] = 5  ;
   piCompC['V'] = 5  ;
   piCompC['T'] = 4  ;
   piCompC['C'] = 3  ;
   piCompC['L'] = 6  ;
   piCompC['I'] = 6  ;
   piCompC['N'] = 4  ;
   piCompC['D'] = 4  ;
   piCompC['Q'] = 5  ;
   piCompC['K'] = 6  ;
   piCompC['E'] = 5  ;
   piCompC['M'] = 5  ;
   piCompC['H'] = 6  ;
   piCompC['F'] = 9  ;
   piCompC['R'] = 6  ;
   piCompC['Y'] = 9  ;
   piCompC['W'] = 11 ;
   piCompC['O'] = 5  ;

   piCompH['G'] = 3  ;
   piCompH['A'] = 5  ;
   piCompH['S'] = 5  ;
   piCompH['P'] = 7  ;
   piCompH['V'] = 9  ;
   piCompH['T'] = 7  ;
   piCompH['C'] = 5  ;
   piCompH['L'] = 11 ;
   piCompH['I'] = 11 ;
   piCompH['N'] = 6  ;
   piCompH['D'] = 5  ;
   piCompH['Q'] = 8  ;
   piCompH['K'] = 12 ;
   piCompH['E'] = 7  ;
   piCompH['M'] = 9  ;
   piCompH['H'] = 7  ;
   piCompH['F'] = 9  ;
   piCompH['R'] = 12 ;
   piCompH['Y'] = 9  ;
   piCompH['W'] = 10 ;
   piCompH['O'] = 12 ;

   piCompN['G'] = 1 ;
   piCompN['A'] = 1 ;
   piCompN['S'] = 1 ;
   piCompN['P'] = 1 ;
   piCompN['V'] = 1 ;
   piCompN['T'] = 1 ;
   piCompN['C'] = 1 ;
   piCompN['L'] = 1 ;
   piCompN['I'] = 1 ;
   piCompN['N'] = 2 ;
   piCompN['D'] = 1 ;
   piCompN['Q'] = 2 ;
   piCompN['K'] = 2 ;
   piCompN['E'] = 1 ;
   piCompN['M'] = 1 ;
   piCompN['H'] = 3 ;
   piCompN['F'] = 1 ;
   piCompN['R'] = 4 ;
   piCompN['Y'] = 1 ;
   piCompN['W'] = 2 ;
   piCompN['O'] = 2 ;

   piCompO['G'] = 1 ;
   piCompO['A'] = 1 ;
   piCompO['S'] = 2 ;
   piCompO['P'] = 1 ;
   piCompO['V'] = 1 ;
   piCompO['T'] = 2 ;
   piCompO['C'] = 1 ;
   piCompO['L'] = 1 ;
   piCompO['I'] = 1 ;
   piCompO['N'] = 2 ;
   piCompO['D'] = 3 ;
   piCompO['Q'] = 2 ;
   piCompO['K'] = 1 ;
   piCompO['E'] = 3 ;
   piCompO['M'] = 1 ;
   piCompO['H'] = 1 ;
   piCompO['F'] = 1 ;
   piCompO['R'] = 1 ;
   piCompO['Y'] = 2 ;
   piCompO['W'] = 1 ;
   piCompO['O'] = 2 ;

   piCompS['C'] = 1;
   piCompS['M'] = 1;

}


void EXTRACT_QUERY_STRING(char *szInputSequence,
                          int *iCharge)
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

            if (!strcmp(szWord, "sequence") )
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               if (strlen(szWord)>MAX_SEQUENCE)
               {
                  printf(" Error - input string is greater than %d characters.\n", MAX_SEQUENCE);
                  break;
               }
               strcpy(szInputSequence, szWord);
            }
            else if (!strcmp(szWord, "charge") )
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               sscanf(szWord, "%d", iCharge);
            }
            else
            {
               getword(szWord, szQuery, '&');
            }
         }
      }
   }
} // EXTRACT_QUERY_STRING

