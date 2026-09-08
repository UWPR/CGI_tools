/*
 * cgi-page.h
 *
 * Self-contained HTML page header and footer for the UWPR CGI tools.
 *
 * The tools used to read /net/pr/vol3/www/html/__header.php and
 * __footer.php at run time, which tied them to the web site's file
 * layout and to its external CSS/JS assets.  The markup and styling
 * below are compiled into each CGI so the executables are stand-alone:
 * they can run on any host without any supporting files.
 *
 * Usage:
 *    PRINT_PAGE_HEADER("Tool title", szExtraHead);   // after Content-type
 *    ... page body ...
 *    PRINT_PAGE_FOOTER();
 *
 * szExtraHead may be NULL or a string of additional markup (e.g. <style>
 * or <script> blocks) that is inserted verbatim just before </head>.
 *
 * This file is plain C so it can be included from both .c and .cpp sources.
 */

#ifndef CGI_PAGE_H
#define CGI_PAGE_H

#include <stdio.h>

static const char *g_szPageCSS =
"   <style>\n"
"      html { box-sizing: border-box; }\n"
"      *, *:before, *:after { box-sizing: inherit; }\n"
"      body { margin: 0; background: #e9edf1; color: #333; font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; font-size: 14px; line-height: 1.5; }\n"
"      a { color: #2f6f9f; }\n"
"      a[onclick] { cursor: pointer; }\n"
"      .wrapper { max-width: 1200px; margin: 0 auto; padding: 0 1em; }\n"
"      #header { display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 0.5em; padding: 1em 0; }\n"
"      #header #logo a { font-size: 1.4em; font-weight: 700; color: #333; text-decoration: none; letter-spacing: 1px; text-transform: uppercase; }\n"
"      #header nav ul { list-style: none; margin: 0; padding: 0; display: flex; gap: 1.5em; }\n"
"      #header nav a { text-decoration: none; color: #555; text-transform: uppercase; font-size: 0.85em; letter-spacing: 1px; }\n"
"      #header nav a:hover { color: #2f6f9f; }\n"
"      #page { position: relative; background: #fff; padding: 3em 4em; border-radius: 6px; margin-bottom: 2em; }\n"
"      header.major { padding-bottom: 2em; text-transform: uppercase; }\n"
"      header.major h1, header.major h2 { margin: 0; padding: 0; line-height: 1em; font-weight: 800; font-size: 2.4em; }\n"
"      #left40entry, #right60entry { margin: 0.5em 0 1em 0; }\n"
"      input[type=text], select, textarea { font: inherit; padding: 0.2em 0.4em; border: 1px solid #bbb; border-radius: 3px; }\n"
"      input[type=submit] { font: inherit; padding: 0.35em 1em; background: #2f6f9f; color: #fff; border: none; border-radius: 3px; cursor: pointer; }\n"
"      input[type=submit]:hover { background: #245a82; }\n"
"      footer { padding: 1em 0 3em 0; text-align: center; font-size: 0.85em; color: #666; }\n"
"      footer a { color: #666; }\n"
"      @media (max-width: 700px) { #page { padding: 1.5em 1em; } header.major h1, header.major h2 { font-size: 1.8em; } }\n"
"   </style>\n";

static void PRINT_PAGE_HEADER(const char *szTitle, const char *szExtraHead)
{
   printf("<!DOCTYPE html>\n");
   printf("<html lang=\"en\">\n");
   printf("<head>\n");
   printf("   <meta charset=\"utf-8\">\n");
   printf("   <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n");
   printf("   <meta name=\"description\" content=\"University of Washington's Proteomics Resource\">\n");
   printf("   <meta name=\"keywords\" content=\"proteomics peptide protein mass spectrometry\">\n");
   if (szTitle && szTitle[0])
      printf("   <title>%s - UW Proteomics Resource</title>\n", szTitle);
   else
      printf("   <title>UW Proteomics Resource</title>\n");
   fputs(g_szPageCSS, stdout);
   if (szExtraHead)
      fputs(szExtraHead, stdout);
   printf("</head>\n");
   printf("<body>\n");
   printf("   <div class=\"wrapper\">\n");
   printf("   <main>\n");
   printf("      <div id=\"header\">\n");
   printf("         <div id=\"logo\"><a href=\"https://proteomicsresource.washington.edu/\">UW Proteomics Resource</a></div>\n");
   printf("         <nav id=\"nav\">\n");
   printf("            <ul>\n");
   printf("               <li><a href=\"https://proteomicsresource.washington.edu/\">home</a></li>\n");
   printf("            </ul>\n");
   printf("         </nav>\n");
   printf("      </div>\n");
   printf("\n");
}

static void PRINT_PAGE_FOOTER(void)
{
   printf("   </main>\n");
   printf("   <footer>\n");
   printf("      <div id=\"copyright\">\n");
   printf("         <p id=\"legal\">&copy; UWPR 2007. All Rights Reserved.\n");
   printf("            &nbsp; <a href=\"https://www.washington.edu/online/terms/\">Terms</a>\n");
   printf("            &nbsp; <a href=\"https://www.washington.edu/online/privacy/\">Privacy</a></p>\n");
   printf("      </div>\n");
   printf("   </footer>\n");
   printf("   </div>\n");
   printf("</body>\n");
   printf("</html>\n");
}

#endif /* CGI_PAGE_H */
