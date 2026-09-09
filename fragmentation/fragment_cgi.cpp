#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <vector>
#include <algorithm>

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

#include "AminoAcidMasses.h"
#include "cgi-page.h"

#define SAMPLEPEPTIDE "DIGSESTEDQAMEDIK"
#define MAX_SEQUENCE 100000
#define MAX_CHARGE   10


void EXTRACT_QUERY_STRING(char *szInputSequence,
      char *szUserMods,
      int *iMassType,
      int *iCharge,
      int *iIonSeries,
      int *bModified,
      int *bRunCalculator);
double COMPUTE_PI(char *seq,
      unsigned long seq_length,
      int charge_increment);
void INIT_COMP(int *piCompC,
      int *piCompH,
      int *piCompN,
      int *piCompO,
      int *piCompS,
      int *piCompSe);

extern char* read_POST(void);
extern void getword(char *word, char *line, char stop);
extern void plustospace(char *str);
extern void unescape_url(char *url);

int main(int argc, char **argv)
{
   int i;
   int iLenPeptide;
   int iMassType;
   int iCharge;
   int iIonSeries;
   int bModified;
   int bRunCalculator;

   char szInputSequence[MAX_SEQUENCE];
   char szUserMods[MAX_SEQUENCE];


   double pdMassAA[128];
   double pdPositionMod[MAX_SEQUENCE];

   int piCompC[MAX_SEQUENCE];
   int piCompH[MAX_SEQUENCE];
   int piCompN[MAX_SEQUENCE];
   int piCompO[MAX_SEQUENCE];
   int piCompS[MAX_SEQUENCE];
   int piCompSe[MAX_SEQUENCE];

   int iC;
   int iH;
   int iN;
   int iO;
   int iS;
   int iSe;
 
   printf("Content-type: text/html\n\n");

//       printf("\n#myStyledTable td { padding-left: 0.4rem; padding-right: 0.4rem; text-align: right; border: 1px solid; border-color: #ADD8E6; }\n");

   // header
   PRINT_PAGE_HEADER("Peptide fragmentation", NULL);
   printf("\n");

   printf("    <div id=\"page\" class=\"container\">\n");
   printf("       <section>\n");
   printf("          <header class=\"major\">\n");
   printf("             <h2>Peptide fragmentation</h2>\n");
   printf("             <p class=\"lede\">Theoretical fragment ion masses for a peptide sequence, with optional fixed or user defined modifications.</p>\n");
   printf("          </header>\n");


   printf("  <script>\n");
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


// strcpy(szInputSequence, SAMPLEPEPTIDE);
   szInputSequence[0]=0;
   szUserMods[0]=0;
   iMassType = 1;  // 1=mono, 0=avg
   iCharge = 1;
   iIonSeries = 18;
   bModified = 0;
   bRunCalculator = 1;
   memset(pdPositionMod, 0, sizeof(pdPositionMod));

   EXTRACT_QUERY_STRING(szInputSequence, szUserMods, &iMassType, &iCharge, &iIonSeries, &bModified, &bRunCalculator);

   {
      const char *szScriptName = getenv("SCRIPT_NAME");
      printf("         <form action=\"%s\" name=\"fragmentForm\" method=\"post\">\n", szScriptName ? szScriptName : "");
   }

   printf("\n");
   printf("         <div class=\"form-grid\">\n");

   // ---- left panel: sequence
   printf("         <div class=\"panel\">\n");
   printf("            <label class=\"field-label\" for=\"peptide\">Peptide sequence</label>\n");
   printf("            <div class=\"actions\" style=\"margin-top: 0\">\n");
   printf("               <input type=\"text\" name=\"sequence\" id=\"peptide\" class=\"mono\" size=\"32\" value=\"");
   print_html_encoded(szInputSequence);
   printf("\">\n");
   printf("               <input type=\"submit\" value=\"Fragment!\">\n");
   printf("            </div>\n");
   printf("            <div class=\"link-row\">\n");
   printf("               <button type=\"button\" class=\"link-btn\" onclick=\"pasteExample();\">Paste a sample peptide</button>\n");
   printf("               <button type=\"button\" class=\"link-btn\" onclick=\"clearExample();\">Clear sequence</button>\n");
   printf("            </div>\n");
   printf("            <details class=\"more\">\n");
   printf("               <summary>Notes on custom residues</summary>\n");
   printf("               <p class=\"hint\">Letters B, J, X and Z have zero mass and can be used as custom residues by specifying a modification mass for them. U is selenocysteine and O is pyrrolysine. ");
   printf("All lower case letters (except c, e, h, n, o, p, s, which carry elemental masses) also have zero mass and can be used as custom residues.</p>\n");
   printf("               <p class=\"hint\">For example, 2-hydroxyproline is C5H9NO3 with a free amino acid mass of 131.058243, calculated with the <a href=\"https://proteomicsresource.washington.edu/cgi-bin/element.cgi\">Elemental Mass calculator</a>. ");
   printf("This calculator needs residue masses, so subtract water: 131.058243 - 18.01056 = 113.047683. Use lower case <code>v</code> for 2-hydroxyproline in the sequence and specify <code>113.047683@v</code> as a user defined modification.</p>\n");
   printf("               <p class=\"hint\">If there is an extension that would be helpful in this calculator, <a href=\"mailto:engj@uw.edu\">let me know</a>.</p>\n");
   printf("            </details>\n");
   printf("         </div>\n");

   // ---- right panel: parameters
   printf("         <div class=\"panel\">\n");
   printf("            <fieldset class=\"options\">\n");
   printf("               <legend>Mass type</legend>\n");
   printf("               <div class=\"choice-row\">\n");
   printf("                  <label><input type=\"radio\" name=\"masstype\" value=\"1\"%s>monoisotopic</label>\n", (iMassType==1?" checked":""));
   printf("                  <label><input type=\"radio\" name=\"masstype\" value=\"0\"%s>average</label>\n", (iMassType==0?" checked":""));
   printf("               </div>\n");
   printf("            </fieldset>\n");

   printf("            <fieldset class=\"options\">\n");
   printf("               <legend>Product ion charge</legend>\n");
   printf("               <div class=\"choice-row\">\n");
   for (i=1; i<=MAX_CHARGE; i++)
      printf("                  <label><input type=\"radio\" name=\"chargestate\" value=\"%d\"%s>%d+</label>\n", i, (iCharge==i?" checked":""), i);
   printf("               </div>\n");
   printf("            </fieldset>\n");

   printf("            <fieldset class=\"options\">\n");
   printf("               <legend>Fragment ion types</legend>\n");
   printf("               <div class=\"choice-row\">\n");
   printf("                  <label><input type=\"checkbox\" name=\"ionseries\" value=\"1\"%s>a</label>\n", (iIonSeries&1?" checked":""));
   printf("                  <label><input type=\"checkbox\" name=\"ionseries\" value=\"2\"%s>b</label>\n", (iIonSeries&2?" checked":""));
   printf("                  <label><input type=\"checkbox\" name=\"ionseries\" value=\"4\"%s>c</label>\n", (iIonSeries&4?" checked":""));
   printf("                  <label><input type=\"checkbox\" name=\"ionseries\" value=\"8\"%s>x</label>\n", (iIonSeries&8?" checked":""));
   printf("                  <label><input type=\"checkbox\" name=\"ionseries\" value=\"16\"%s>y</label>\n", (iIonSeries&16?" checked":""));
   printf("                  <label><input type=\"checkbox\" name=\"ionseries\" value=\"32\"%s>z</label>\n", (iIonSeries&32?" checked":""));
   printf("                  <label><input type=\"checkbox\" name=\"ionseries\" value=\"64\"%s>z&#149;</label>\n", (iIonSeries&64?" checked":""));
   printf("               </div>\n");
   printf("            </fieldset>\n");

   printf("            <fieldset class=\"options\">\n");
   printf("               <legend>Modifications</legend>\n");
   printf("               <div class=\"choice-col\">\n");
   printf("                  <label><input type=\"radio\" name=\"modified\" value=\"0\"%s>None</label>\n", (bModified==0?" checked":""));
   printf("                  <label><input type=\"radio\" name=\"modified\" value=\"3\"%s>Carbamidomethyl C (57.021464 mono, 57.0513 avg)</label>\n", (bModified==3?" checked":""));
   printf("                  <label><input type=\"radio\" name=\"modified\" value=\"4\"%s>User defined: <input type=\"text\" name=\"usermods\" class=\"mono\" value=\"", (bModified==4?" checked":""));
   print_html_encoded(szUserMods);
   printf("\"></label>\n");
   printf("                  <p class=\"hint\">Enter <code>mass@AA</code> or <code>mass@position</code>, space separated, e.g. <code>15.995@M 57.0215@3</code>. ");
   printf("Use <code>[</code> for the N-terminus and <code>]</code> for the C-terminus, e.g. <code>16.0@[</code>.</p>\n");
   printf("               </div>\n");
   printf("            </fieldset>\n");
   printf("         </div>\n");

   printf("         </div>\n");  // form-grid
   printf("         </form>\n\n");

   printf("         <div id=\"results\">\n");

   if (bRunCalculator)
   {
      INITIALIZE_MASS(pdMassAA, iMassType);
      pdMassAA['['] = 0.0;  // initialize n-term mass
      pdMassAA[']'] = 0.0;  // initialize c-term mass

      if (bModified==1)
      {
         pdMassAA['K'] += 8.014199;
         pdMassAA['R'] += 10.008269;
      }
      else if (bModified==2)
      {
         if (iMassType==1) // mono
            pdMassAA['R'] += 3.98814;
         else
            pdMassAA['R'] += 3.9737;
      }
      else if (bModified==3)
      {
         if (iMassType==1) // mono
            pdMassAA['C'] += 57.021464;
         else
            pdMassAA['C'] += 57.0513;
      }
      else if (bModified==4)
      {
         char *tok;

         tok = strtok(szUserMods, " ");
         while (tok != NULL)
         {
            char *pStr;
            double dMass;
            char szResidue[100];

            pStr=strchr(tok, '@');
            if (!pStr) { tok = strtok(NULL, " "); continue; }
            *pStr = ' ';

            sscanf(tok, "%lf %99s", &dMass, szResidue);

            // first check if residue or number is entered
            if (strspn(szResidue, "0123456789")==strlen(szResidue))
            {
               int iPos=0;

               sscanf(szResidue, "%d", &iPos);
               if (iPos>=1 && iPos<=MAX_SEQUENCE)
                  pdPositionMod[iPos-1] = dMass;
            }
            else
            {
               pdMassAA[(unsigned char)szResidue[0]] += dMass;
            }

            tok = strtok(NULL, " ");
         }
      }

      if (1)
      {
         int i;
         int j;

         char szTmp[MAX_SEQUENCE];

         j=0;
         iLenPeptide = strlen(szInputSequence);
         for (i = 0; i < iLenPeptide; i++)
         {
            if (isalpha(szInputSequence[i]))
            {
               szTmp[j++] = szInputSequence[i];
               szTmp[j] = '\0';
            }
         }

         strcpy(szInputSequence, szTmp);
      }

      iLenPeptide = strlen(szInputSequence);

      if (iLenPeptide>0)
      {
         int i;
//       double dPI = COMPUTE_PI(szInputSequence, iLenPeptide, 0);
         double dAion = 0.0,
                dBion = 0.0,
                dCion = 0.0,
                dXion = 0.0,
                dYion = 0.0,
                dZion = 0.0,
                dZdotion = 0.0,
                dNterm = pdMassAA['h'] + pdMassAA['['],
                dCterm = pdMassAA['o'] + pdMassAA['h'] + pdMassAA[']'],
                dPepMass,
                dProton = PROTON_MASS,
                dCO = pdMassAA['c'] + pdMassAA['o'],
                dH2 = pdMassAA['h'] + pdMassAA['h'],
                dNH3 = pdMassAA['n'] + pdMassAA['h'] + pdMassAA['h'] + pdMassAA['h'];
         char szCharge[24];
         std::vector<double> vIons;

         if (iCharge==1)
            strcpy(szCharge, "<sup>+</sup>");
         else if (iCharge==2)
            strcpy(szCharge, "<sup>++</sup>");
         else   // charge state 3
         {
            sprintf(szCharge, "<sup>%d+</sup>", iCharge);
         }
   
         dPepMass = dNterm + dCterm + dProton;
         for (i=0; i<iLenPeptide; i++)
         {
            dPepMass += pdMassAA[(int)(szInputSequence[i])];
            dPepMass += pdPositionMod[i];
         }
   
         dBion = dNterm - pdMassAA['h'] + dProton;
         dYion = dPepMass;

         printf("<div class=\"table-wrap\">\n");
         printf("<table class=\"results\">\n");
         printf("<thead>");
         //printf("<tr align=\"center\">");
         printf("<tr>");
         if (iIonSeries&1)
            printf("<th class=\"center\">a%s</th>", szCharge);
         if (iIonSeries&2)
            printf("<th class=\"center\">b%s</th>", szCharge);
         if (iIonSeries&4)
            printf("<th class=\"center\">c%s</th>", szCharge);
         printf("<th class=\"aa\">#</th>");
         printf("<th class=\"aa\">AA</th>");
         printf("<th class=\"aa\">#</th>");
         if (iIonSeries&8)
            printf("<th class=\"center\">x%s</th>", szCharge);
         if (iIonSeries&16)
            printf("<th class=\"center\">y%s</th>", szCharge);
         if (iIonSeries&32)
            printf("<th class=\"center\">z%s</th>", szCharge);
         if (iIonSeries&64)
            printf("<th class=\"center\">z&#149;%s</th>", szCharge);
         printf("</tr>");
         printf("</thead>\n");
         printf("<tbody>\n");
   
         for (i=0; i<iLenPeptide; i++)
         {
            if (i<iLenPeptide-1)
            {
               dBion += pdMassAA[(int)(szInputSequence[i])];
               dBion += pdPositionMod[i];
               dAion = dBion - dCO;
               dCion = dBion + dNH3;
            }
            if (i>0)
            {
               dYion -= pdMassAA[(int)(szInputSequence[i-1])];
               dYion -= pdPositionMod[i-1];
               dXion = dYion + dCO - dH2;
               dZion = dYion - dNH3;
               dZdotion = dZion + pdMassAA['h'];
            }
            if (i==0)
               dYion -= pdMassAA['['];
   
            printf("<tr>");
   
            // print A-ions
            if (iIonSeries&1)
            {
               if (i<iLenPeptide-1)
               {
                  printf("<td class=\"num\">%0.6f</td>", (dAion + (iCharge-1)*dProton)/iCharge);
                  vIons.push_back((dAion + (iCharge-1)*dProton)/iCharge);
               }
               else
                  printf("<td></td>");
            }
   
            // print B-ions
            if (iIonSeries&2)
            {
               if (i<iLenPeptide-1)
               {
                  printf("<td class=\"num\">%0.6f</td>", (dBion + (iCharge-1)*dProton)/iCharge);
                  vIons.push_back((dBion + (iCharge-1)*dProton)/iCharge);
               }
               else
                  printf("<td></td>");
            }
   
            // print C-ions
            if (iIonSeries&4)
            {
               if (i<iLenPeptide-1)
               {
                  printf("<td class=\"num\">%0.6f</td>", (dCion + (iCharge-1)*dProton)/iCharge);
                  vIons.push_back((dCion + (iCharge-1)*dProton)/iCharge);
               }
               else
                  printf("<td></td>");
            }
   
            printf("<td class=\"aa\">%d</td>", i+1);
            printf("<td class=\"aa\">%c</td>", szInputSequence[i]);
            printf("<td class=\"aa\">%d</td>", iLenPeptide - i);
      
            // print X-ions
            if (iIonSeries&8)
            {
               if (i>0)
               {
                  printf("<td class=\"num\">%0.6f</td>", (dXion + (iCharge-1)*dProton)/iCharge);
                  vIons.push_back((dXion + (iCharge-1)*dProton)/iCharge);
               }
               else
                  printf("<td></td>");
            }
   
            // print Y-ions
            if (iIonSeries&16)
            {
               if (i>0)
               {
                  printf("<td class=\"num\">%0.6f</td>", (dYion + (iCharge-1)*dProton)/iCharge);
                  vIons.push_back((dYion + (iCharge-1)*dProton)/iCharge);
               }
               else
                  printf("<td></td>");
            }
   
            // print Z-ions
            if (iIonSeries&32)
            {
               if (i>0)
               {
                  printf("<td class=\"num\">%0.6f</td>", (dZion + (iCharge-1)*dProton)/iCharge);
                  vIons.push_back((dZion + (iCharge-1)*dProton)/iCharge);
               }
               else
                  printf("<td></td>");
            }
   
            // print Z-dot-ions
            if (iIonSeries&64)
            {
               if (i>0)
                  printf("<td class=\"num\">%0.6f</td>", (dZdotion + (iCharge-1)*dProton)/iCharge);
               else
                  printf("<td></td>");
            }
   
            printf("</tr>\n");
         }
   
         printf("</tbody>\n");
         printf("</table>\n");
         printf("</div>\n");

         printf("<pre class=\"readout\">");

         iC=0;
         iH=0;
         iN=0;
         iO=0;
         iS=0;
         iSe=0;

         INIT_COMP(piCompC, piCompH, piCompN, piCompO, piCompS, piCompSe);

         for (i=0; i<iLenPeptide; i++)
         {
            iC += piCompC[(int)szInputSequence[i]];
            iH += piCompH[(int)szInputSequence[i]];
            iN += piCompN[(int)szInputSequence[i]];
            iO += piCompO[(int)szInputSequence[i]];
            iS += piCompS[(int)szInputSequence[i]];
            iSe += piCompSe[(int)szInputSequence[i]];
         }
         printf("composition: ");
         if (iC > 0)
            printf("C<sub>%d</sub>", iC);
         if (iH > 0)
            printf("H<sub>%d</sub>", iH);
         if (iN > 0)
            printf("N<sub>%d</sub>", iN);
         if (iO > 0)
            printf("O<sub>%d</sub>", iO);
         if (iS > 0)
            printf("S<sub>%d</sub>", iS);
         if (iSe > 0)
            printf("Se<sub>%d</sub>", iSe);
         printf("\n");


         printf("neutral mass: %12.6f\n", dPepMass - dProton);
         printf("   [M +  H+]: %12.6f\n", dPepMass);
         printf("   [M + 2H+]: %12.6f\n", dPepMass + dProton);
         printf("   [M + 3H+]: %12.6f\n", dPepMass + 2*dProton);
         printf("   [M + 4H+]: %12.6f\n", dPepMass + 3*dProton);
         for (i=1; i<=MAX_CHARGE; i++)
         {
            if (i >= 10)
               printf("     +%2d m/z: %12.6f\n", i, (dPepMass + (i-1)*dProton)/i);
            else
               printf("      +%d m/z: %12.6f\n", i, (dPepMass + (i-1)*dProton)/i);
         }
         printf("</pre>\n");

         sort(vIons.begin(), vIons.end());

         printf("<details class=\"more\"><summary>Peaks in MS2 format</summary>\n");
         printf("<pre>");
         printf("S 1 1 %0.4lf\n", dPepMass);
         printf("Z 1 %0.4f\n", dPepMass);
         for (auto x : vIons)
         {
            printf("%lf 100\n", x);
         }
         printf("</pre></details>\n");
      }
   }
   printf("</div>\n");  // results
   printf("       </section>\n");
   printf("    </div>\n");  // page

   // footer
   PRINT_PAGE_FOOTER();

   exit(EXIT_SUCCESS);
}


void EXTRACT_QUERY_STRING(char *szInputSequence,
      char *szUserMods,
      int *iMassType,
      int *iCharge,
      int *iIonSeries,
      int *bModified,
      int *bRunCalculator)
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

         if (strstr(szQuery, "chargestate"))
            *iCharge = 0;
         if (strstr(szQuery, "ionseries"))
            *iIonSeries= 0;

         szQuery[iContentLength]='\0';

         for (i=0; szQuery[0] != '\0'; i++)
         {
            getword(szWord, szQuery, '=');
            plustospace(szWord);
            unescape_url(szWord);

            if (!strcmp(szWord, "sequence") )
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               if (strlen(szWord)>=MAX_SEQUENCE)
               {
                  printf(" Error - input string is greater than %d characters.\n", MAX_SEQUENCE);
                  *bRunCalculator = 0;
                  break;
               }
               strcpy(szInputSequence, szWord);
            }
            else if (!strcmp(szWord, "usermods") )
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               if (strlen(szWord) >= MAX_SEQUENCE)
               {
                  printf(" Error - usermods string is greater than %d characters.\n", MAX_SEQUENCE);
                  *bRunCalculator = 0;
                  break;
               }
               strcpy(szUserMods, szWord);
            }
            else if (!strcmp(szWord, "masstype") )
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               sscanf(szWord, "%d", iMassType);
               if (*iMassType < 0 || *iMassType>1)
                  *iMassType = 1;
            }
            else if (!strcmp(szWord, "chargestate") )
            {
               int iVal;
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               sscanf(szWord, "%d", &iVal);
               *iCharge = iVal;
            }
            else if (!strcmp(szWord, "ionseries") )
            {
               int iVal;
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               sscanf(szWord, "%d", &iVal);
               *iIonSeries += iVal;
            }
            else if (!strcmp(szWord, "modified") )
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
               sscanf(szWord, "%d", bModified);
            }
            else
            {
               getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
            }
         }
      }
      else
      {
         *bRunCalculator=0;
      }

      if (strlen(szInputSequence) == 0)
         *bRunCalculator = 0;
   }
} // EXTRACT_QUERY_STRING


void INIT_COMP(int *piCompC,
      int *piCompH,
      int *piCompN,
      int *piCompO,
      int *piCompS,
      int *piCompSe)
{
   int i;

   for (i=0; i<128; i++)
   {
      piCompC[i]=0;
      piCompH[i]=0;
      piCompN[i]=0;
      piCompO[i]=0;
      piCompS[i]=0;
      piCompSe[i]=0;
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
   piCompC['O'] = 12 ;
   piCompC['U'] = 3  ;

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
   piCompH['O'] = 19 ;
   piCompH['U'] = 5  ;

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
   piCompN['O'] = 3 ;
   piCompN['U'] = 1 ;

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
   piCompO['U'] = 1 ;

   piCompS['C'] = 1;
   piCompS['M'] = 1;

   piCompSe['U'] = 1;
}
