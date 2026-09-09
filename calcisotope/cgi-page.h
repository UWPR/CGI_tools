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
 * Shared CSS classes available to the tools:
 *    #page                 white card that holds the tool
 *    header.major h1/h2    tool title;  p.lede  short description under it
 *    .form-grid            responsive two-column layout for form panels
 *    .panel                a group of related fields
 *    .field / .field-label / .hint / .actions / .link-btn / .btn
 *    fieldset.options / legend / .choice-row   radio and checkbox groups
 *    table.results         data table (no outer border, zebra rows, hover)
 *       td.num right-aligned numbers, td.aa / th.aa highlighted residue cells
 *    .readout              monospace <pre> block for mass summaries
 *    .readout-hero         large single headline value
 *    details.more          collapsible secondary content
 *    .summary / .kv        key/value rows above a results table
 *    .note                 small muted note
 *
 * This file is plain C so it can be included from both .c and .cpp sources.
 */

#ifndef CGI_PAGE_H
#define CGI_PAGE_H

#include <stdio.h>

static const char *g_szPageCSS =
"   <style>\n"
"      :root {\n"
"         --bg: #f3f4f7; --surface: #ffffff; --surface-2: #f7f6fa;\n"
"         --text: #1f2430; --muted: #5f6672; --border: #e2e5eb; --border-strong: #c7ccd6;\n"
"         --accent: #4b2e83; --accent-dark: #39206a; --accent-soft: #ece7f5; --gold: #b7a57a;\n"
"         --mono: 'SFMono-Regular', Menlo, Consolas, 'Liberation Mono', monospace;\n"
"         --sans: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;\n"
"         --radius: 10px; --shadow: 0 1px 2px rgba(20,24,40,.05), 0 8px 24px rgba(20,24,40,.06);\n"
"      }\n"
"      *, *::before, *::after { box-sizing: border-box; }\n"
"      html { -webkit-text-size-adjust: 100%; }\n"
"      body { margin: 0; background: var(--bg); color: var(--text); font-family: var(--sans); font-size: 15px; line-height: 1.55; }\n"
"      a { color: var(--accent); text-decoration: none; }\n"
"      a:hover { text-decoration: underline; }\n"
"      .wrapper { max-width: 1160px; margin: 0 auto; padding: 0 1.25rem 2rem; }\n"
"      /* top bar */\n"
"      #header { display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: .75rem; padding: 1.1rem 0 1.25rem; }\n"
"      #header #logo a { display: inline-flex; align-items: center; gap: .6rem; font-weight: 700; font-size: 1.05rem; letter-spacing: .01em; color: var(--text); }\n"
"      #header #logo a::before { content: ''; width: 10px; height: 26px; border-radius: 3px; background: linear-gradient(180deg, var(--accent) 55%, var(--gold) 55%); }\n"
"      #header nav ul { list-style: none; margin: 0; padding: 0; display: flex; gap: .25rem; }\n"
"      #header nav a { display: inline-block; padding: .35rem .75rem; border-radius: 999px; color: var(--muted); font-size: .85rem; font-weight: 600; letter-spacing: .04em; text-transform: uppercase; }\n"
"      #header nav a:hover { background: var(--accent-soft); color: var(--accent); text-decoration: none; }\n"
"      /* tool card */\n"
"      #page { background: var(--surface); border: 1px solid var(--border); border-radius: var(--radius); box-shadow: var(--shadow); padding: 2.25rem 2.5rem 2.5rem; }\n"
"      header.major { margin-bottom: 1.75rem; }\n"
"      header.major h1, header.major h2 { margin: 0; font-size: 1.75rem; line-height: 1.2; font-weight: 700; letter-spacing: -.01em; }\n"
"      header.major h1::after, header.major h2::after { content: ''; display: block; width: 48px; height: 4px; margin-top: .6rem; border-radius: 2px; background: var(--accent); }\n"
"      .lede { margin: .9rem 0 0; color: var(--muted); }\n"
"      /* forms */\n"
"      form { margin: 0; }\n"
"      .form-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 1.25rem 2rem; }\n"
"      .panel { background: var(--surface-2); border: 1px solid var(--border); border-radius: var(--radius); padding: 1.25rem 1.4rem; }\n"
"      .field { margin: 0 0 1rem; }\n"
"      .field:last-child { margin-bottom: 0; }\n"
"      .field-label, legend { display: block; font-weight: 600; font-size: .9rem; margin-bottom: .4rem; }\n"
"      .hint { color: var(--muted); font-size: .84rem; line-height: 1.5; margin: .35rem 0 0; }\n"
"      .hint code, .inline-code { font-family: var(--mono); font-size: .9em; background: var(--accent-soft); padding: .05em .35em; border-radius: 4px; }\n"
"      input[type=text], input[type=number], select, textarea { font: inherit; color: var(--text); background: var(--surface); border: 1px solid var(--border-strong); border-radius: 7px; padding: .45rem .65rem; max-width: 100%; }\n"
"      input[type=text].mono, textarea.mono { font-family: var(--mono); font-size: .92rem; }\n"
"      input[type=text].short { width: 4.5rem; text-align: center; }\n"
"      input[type=text].wide, textarea.wide, select.wide { width: 100%; }\n"
"      textarea { resize: vertical; line-height: 1.45; }\n"
"      input:focus, select:focus, textarea:focus, button:focus-visible { outline: 3px solid rgba(75,46,131,.28); outline-offset: 1px; border-color: var(--accent); }\n"
"      input[type=radio], input[type=checkbox] { accent-color: var(--accent); width: 1rem; height: 1rem; margin: 0 .35rem 0 0; vertical-align: -.15em; }\n"
"      fieldset.options { border: 0; padding: 0; margin: 0 0 1rem; min-width: 0; }\n"
"      .choice-row { display: flex; flex-wrap: wrap; gap: .35rem 1rem; }\n"
"      .choice-row label, .choice-col label { display: inline-flex; align-items: center; gap: .1rem; font-size: .92rem; cursor: pointer; }\n"
"      .choice-col label { display: flex; margin: .3rem 0; }\n"
"      .choice-col label input[type=text] { margin-left: .5rem; flex: 1; }\n"
"      .actions { display: flex; flex-wrap: wrap; align-items: center; gap: .75rem; margin-top: 1rem; }\n"
"      .btn, input[type=submit] { font: inherit; font-weight: 600; color: #fff; background: var(--accent); border: 1px solid var(--accent); border-radius: 7px; padding: .5rem 1.25rem; cursor: pointer; transition: background .15s; }\n"
"      .btn:hover, input[type=submit]:hover { background: var(--accent-dark); }\n"
"      .link-btn { font: inherit; font-size: .85rem; color: var(--accent); background: none; border: 0; padding: 0; cursor: pointer; text-decoration: underline; text-underline-offset: 2px; }\n"
"      .link-btn:hover { color: var(--accent-dark); }\n"
"      .link-row { display: flex; flex-wrap: wrap; gap: 1rem; margin-top: .4rem; }\n"
"      .elements { display: flex; flex-wrap: wrap; align-items: flex-end; gap: .75rem 1rem; }\n"
"      .elements input[type=submit] { margin-left: .5rem; }\n"
"      .field-row { display: flex; flex-wrap: wrap; align-items: flex-end; gap: 1rem 2rem; margin: 0 0 1rem; }\n"
"      .field-row .field { margin: 0; }\n"
"      .field-row .field-label, .field-row .field label { display: block; }\n"
"      .elements label { display: flex; flex-direction: column; align-items: center; gap: .3rem; font-weight: 600; font-size: .9rem; }\n"
"      /* results */\n"
"      #results { margin-top: 2rem; }\n"
"      #results:empty { margin-top: 0; }\n"
"      .summary { display: grid; grid-template-columns: max-content 1fr; gap: .25rem 1.25rem; margin: 0 0 1.25rem; font-size: .92rem; }\n"
"      .summary dt { color: var(--muted); font-weight: 600; }\n"
"      .summary dd { margin: 0; font-family: var(--mono); word-break: break-all; }\n"
"      .note { color: var(--muted); font-size: .84rem; margin: 0 0 .75rem; }\n"
"      .table-wrap { overflow-x: auto; }\n"
"      table.results { border-collapse: collapse; font-family: var(--mono); font-size: .86rem; margin: 0; }\n"
"      table.results th, table.results td { padding: .4rem .8rem; border-bottom: 1px solid var(--border); text-align: left; white-space: nowrap; }\n"
"      table.results th { background: var(--accent-soft); color: var(--accent-dark); font-family: var(--sans); font-size: .8rem; font-weight: 600; letter-spacing: .03em; border-bottom: 1px solid var(--border-strong); }\n"
"      table.results th:first-child { border-top-left-radius: 7px; }\n"
"      table.results th:last-child { border-top-right-radius: 7px; }\n"
"      table.results tbody tr:nth-child(even) { background: var(--surface-2); }\n"
"      table.results tbody tr:hover { background: #f1edf8; }\n"
"      table.results .num { text-align: right; font-variant-numeric: tabular-nums; }\n"
"      table.results .center { text-align: center; }\n"
"      table.results .aa { text-align: center; background: #fbf8ee; font-weight: 600; }\n"
"      table.results tbody tr:hover .aa { background: #f5efd9; }\n"
"      table.results th.aa { background: var(--gold); color: #fff; }\n"
"      table.sortable th { cursor: pointer; user-select: none; }\n"
"      table.sortable th:hover { background: #e2d9f1; }\n"
"      .readout { font-family: var(--mono); font-size: .9rem; line-height: 1.6; background: var(--surface-2); border: 1px solid var(--border); border-radius: var(--radius); padding: 1rem 1.25rem; margin: 1.25rem 0 0; overflow-x: auto; }\n"
"      .readout pre { margin: 0; font: inherit; }\n"
"      .readout-hero { margin: 0 0 1rem; }\n"
"      .readout-hero span { display: block; color: var(--muted); font-size: .85rem; font-weight: 600; }\n"
"      .readout-hero strong { display: block; font-family: var(--mono); font-size: 1.7rem; font-weight: 600; letter-spacing: -.01em; }\n"
"      details.more { margin-top: 1rem; }\n"
"      details.more summary { cursor: pointer; color: var(--accent); font-size: .88rem; font-weight: 600; }\n"
"      details.more pre { font-family: var(--mono); font-size: .78rem; color: var(--muted); background: var(--surface-2); border: 1px solid var(--border); border-radius: var(--radius); padding: .9rem 1.1rem; margin: .6rem 0 0; overflow-x: auto; }\n"
"      .chart-container { max-width: 640px; margin-top: 1.5rem; }\n"
"      /* footer */\n"
"      footer { padding: 1.75rem 0 0; text-align: center; font-size: .82rem; color: var(--muted); }\n"
"      footer a { color: var(--muted); text-decoration: underline; text-underline-offset: 2px; }\n"
"      @media (max-width: 720px) {\n"
"         #page { padding: 1.5rem 1.1rem 1.75rem; }\n"
"         header.major h1, header.major h2 { font-size: 1.45rem; }\n"
"         #header { padding-bottom: 1rem; }\n"
"      }\n"
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
