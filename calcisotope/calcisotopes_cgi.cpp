//
// USAGE:  ./calcisotopes "K[325.13]K[170.11]FTENPKAG"
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "isotopes.h"
#include "cgi-page.h"

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
   // results layout: summary + table on the left, bar chart on the right
   PRINT_PAGE_HEADER("Peptide Isotope Calculator",
      "   <style>\n"
      "      .results-grid { display: flex; flex-wrap: wrap; align-items: flex-start; gap: 1.25rem 2.5rem; }\n"
      "      .results-grid .results-data { flex: 0 1 auto; }\n"
      "      .results-grid .chart-container { flex: 1 1 360px; max-width: 640px; margin-top: 0; }\n"
      "      .chart-controls { display: flex; align-items: center; gap: .6rem; margin-top: .75rem; font-size: .88rem; color: var(--muted); }\n"
      "      .chart-controls label { font-weight: 600; }\n"
      "   </style>\n");
   printf("\n");

   printf("    <div id=\"page\" class=\"container\">\n");
   printf("       <section>\n");
   printf("          <header class=\"major\">\n");
   printf("             <h2>Peptide Isotope Calculator</h2>\n");
   printf("             <p class=\"lede\">Calculates the isotope distribution of a peptide using the <a href=\"http://www.kombyonyx.com/isotopes/\">Isotope Distribution Calculator</a> code by James Redman at Cardiff University, an algorithm modeled after <a href=\"https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/jms.4498\">J. A. Yergey's procedure published in 1983</a>.</p>\n");
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


   szInput[0] = 0;

   EXTRACT_QUERY_STRING(szInput, &iCharge);

   {
      const char *szScriptName = getenv("SCRIPT_NAME");
      printf("         <form action=\"%s\" name=\"fragmentForm\" method=\"post\">\n", szScriptName ? szScriptName : "");
   }
   printf("         <div class=\"panel\">\n");
   printf("            <div class=\"field-row\">\n");
   printf("            <div class=\"field\">\n");
   printf("               <label class=\"field-label\" for=\"peptide\">Peptide sequence</label>\n");
   printf("               <input type=\"text\" name=\"sequence\" id=\"peptide\" class=\"mono\" size=\"40\" value=\"");
   print_html_encoded(szInput);
   printf("\">\n");
   printf("            </div>\n");
   printf("            <div class=\"field\">\n");
   printf("               <label class=\"field-label\" for=\"charge\">Precursor charge</label>\n");
   printf("               <select name=\"charge\" id=\"charge\">\n");
   for (int z = 1; z <= 20; z++)
      printf("                  <option value=\"%d\"%s>%d</option>\n", z, (z == iCharge ? " selected" : ""), z);
   printf("               </select>\n");
   printf("            </div>\n");
   printf("            <input type=\"submit\" value=\"Calculate\">\n");
   printf("            </div>\n");  // field-row
   printf("            <div class=\"link-row\">\n");
   printf("               <button type=\"button\" class=\"link-btn\" onclick=\"pasteExample();\">Paste a sample peptide</button>\n");
   printf("               <button type=\"button\" class=\"link-btn\" onclick=\"clearExample();\">Clear sequence</button>\n");
   printf("            </div>\n");
   printf("         </div>\n");
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

      printf("<div class=\"results-grid\">\n");
      printf("<div class=\"results-data\">\n");
      printf("<dl class=\"summary\">\n");
      printf("<dt>sequence</dt><dd>");
      print_html_encoded(szInput);
      printf("</dd>\n");
      printf("<dt>composition</dt><dd>%s</dd>\n", szComp);
      printf("</dl>\n");

      printf("<div class=\"table-wrap\">\n");
      printf("<table class=\"results\">\n");
      printf("<thead><tr><th>peak</th><th class=\"num\">m/z</th><th class=\"num\">relative abundance</th></tr></thead>\n");
      printf("<tbody>\n");

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
            printf("<td>Mono</td><td class=\"num\">%0.5f</td><td class=\"num\">%0.2f</td>", mass, 100.0 * abun / dMaxAbundance);
         else
            printf("<td>M+%d</td><td class=\"num\">%0.5f</td><td class=\"num\">%0.2f</td>", i,  mass, 100.0 * abun / dMaxAbundance);

         vLabel.push_back(mass);
         vAbun.push_back(100.0 * abun / dMaxAbundance);

         printf("</tr>\n");
      }
      printf("</tbody>\n</table>\n</div>\n");
      printf("</div>\n\n");  // results-data


      // Profile-mode rendering of the centroid distribution: each isotope peak
      // becomes a Gaussian whose width follows the chosen resolving power
      // (FWHM = m/z / R), the Gaussians are summed on a fine m/z grid, and the
      // trace is drawn as a filled line with Chart.js.
      // https://www.chartjs.org/docs/latest/getting-started/
      printf("   <div class=\"chart-container\">\n");
      printf("      <canvas id=\"myChart\"></canvas>\n");
      printf("      <div class=\"chart-controls\">\n");
      printf("         <label for=\"resolution\">Resolving power</label>\n");
      printf("         <select id=\"resolution\" onchange=\"drawSpectrum()\">\n");
      printf("            <option value=\"1000\">1,000</option>\n");
      printf("            <option value=\"5000\">5,000</option>\n");
      printf("            <option value=\"10000\" selected>10,000</option>\n");
      printf("            <option value=\"30000\">30,000</option>\n");
      printf("            <option value=\"100000\">100,000</option>\n");
      printf("         </select>\n");
      printf("      </div>\n");
      printf("   </div>\n");
      printf("</div>\n");  // results-grid
      printf("<script src=\"https://cdn.jsdelivr.net/npm/chart.js\"></script>\n");
      printf("<script>\n");

      printf("  var peakMz = [");
      for (auto it=vLabel.begin(); it!=vLabel.end(); ++it)
      {
         if (it != vLabel.begin())
            printf(", ");
         printf("%0.5lf", (*it));
      }
      printf("];\n");

      printf("  var peakAbun = [");
      for (auto it=vAbun.begin(); it!=vAbun.end(); ++it)
      {
         if (it != vAbun.begin())
            printf(", ");
         printf("%0.4lf", (*it));
      }
      printf("];\n");

      printf("\
  var spectrumChart = null;\n\
  function profileTrace(resolvingPower) {\n\
    var n = peakMz.length;\n\
    var lo = peakMz[0], hi = peakMz[n - 1];\n\
    var sigmaMax = 0, sigmaMin = Infinity, i, j;\n\
    for (i = 0; i < n; i++) {\n\
      var s = (peakMz[i] / resolvingPower) / 2.3548;   /* FWHM -> sigma */\n\
      if (s > sigmaMax) sigmaMax = s;\n\
      if (s < sigmaMin) sigmaMin = s;\n\
    }\n\
    var spacing = (n > 1) ? (peakMz[1] - peakMz[0]) : 1.0;\n\
    var pad = Math.max(0.6 * spacing, 4 * sigmaMax);\n\
    var xmin = lo - pad, xmax = hi + pad;\n\
    /* sample finely enough that even narrow peaks are drawn smoothly */\n\
    var steps = Math.min(20000, Math.max(1500, Math.ceil((xmax - xmin) / (sigmaMin / 4))));\n\
    var dx = (xmax - xmin) / steps;\n\
    var pts = [], ymax = 0;\n\
    for (i = 0; i <= steps; i++) {\n\
      var x = xmin + i * dx, y = 0;\n\
      for (j = 0; j < n; j++) {\n\
        var sig = (peakMz[j] / resolvingPower) / 2.3548;\n\
        var d = (x - peakMz[j]) / sig;\n\
        if (d > -6 && d < 6) y += peakAbun[j] * Math.exp(-0.5 * d * d);\n\
      }\n\
      if (y > ymax) ymax = y;\n\
      pts.push({ x: x, y: y });\n\
    }\n\
    if (ymax > 0) for (i = 0; i < pts.length; i++) pts[i].y = 100 * pts[i].y / ymax;\n\
    return pts;\n\
  }\n\
  function drawSpectrum() {\n\
    var R = parseFloat(document.getElementById('resolution').value);\n\
    var pts = profileTrace(R);\n\
    if (spectrumChart) {\n\
      spectrumChart.data.datasets[0].data = pts;\n\
      spectrumChart.options.scales.x.min = pts[0].x;\n\
      spectrumChart.options.scales.x.max = pts[pts.length - 1].x;\n\
      spectrumChart.update();\n\
      return;\n\
    }\n\
    var ctx = document.getElementById('myChart');\n\
    spectrumChart = new Chart(ctx, {\n\
      type: 'line',\n\
      data: {\n\
        datasets: [{\n\
          label: 'relative intensity',\n\
          data: pts,\n\
          borderColor: '#4b2e83',\n\
          backgroundColor: 'rgba(75,46,131,0.12)',\n\
          borderWidth: 1.5,\n\
          fill: 'origin',\n\
          pointRadius: 0,\n\
          tension: 0\n\
        }]\n\
      },\n\
      options: {\n\
        animation: false,\n\
        aspectRatio: 1.6,\n\
        interaction: { mode: 'nearest', axis: 'x', intersect: false },\n\
        scales: {\n\
          x: { type: 'linear', min: pts[0].x, max: pts[pts.length - 1].x,\n\
               title: { display: true, text: 'm/z' },\n\
               grid: { display: false },\n\
               ticks: { maxTicksLimit: 8, callback: function (v) { return Number(v).toFixed(2); } } },\n\
          y: { beginAtZero: true, max: 105,\n\
               title: { display: true, text: 'relative intensity' },\n\
               grid: { display: false },\n\
               ticks: { stepSize: 25 } }\n\
        },\n\
        plugins: {\n\
          legend: { display: false },\n\
          tooltip: { callbacks: {\n\
            title: function (items) { return 'm/z ' + items[0].parsed.x.toFixed(4); },\n\
            label: function (item) { return item.parsed.y.toFixed(1); } } }\n\
        }\n\
      }\n\
    });\n\
  }\n\
  drawSpectrum();\n\
</script>\n");

   }
   printf("</div>\n");  // results
   printf("       </section>\n");
   printf("    </div>\n");  // page

   // footer
   PRINT_PAGE_FOOTER();

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
               if (strlen(szWord)>=MAX_SEQUENCE)
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
               if (*iCharge < 1)
                  *iCharge = 1;
               else if (*iCharge > 20)
                  *iCharge = 20;
            }
            else
            {
               getword(szWord, szQuery, '&');
            }
         }
      }
   }
} // EXTRACT_QUERY_STRING

