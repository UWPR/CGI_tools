/*  motif_search.cpp  —  Protein Motif Search Tool (C++ CGI)
 *
 *  Compile:
 *      g++ -std=c++17 -O2 -o motif_search.cgi motif_search.cpp
 *      (needs cgi-page.h, the page header/stylesheet shared with the other
 *      UWPR CGI tools, in the same directory)
 *
 *  Deploy:
 *      cp motif_search.cgi /usr/lib/cgi-bin/
 *      chmod 755 /usr/lib/cgi-bin/motif_search.cgi
 *
 *  Apache — enable CGI execution for the directory that holds the binary:
 *      Options +ExecCGI
 *      AddHandler cgi-script .cgi
 *
 *  No third-party libraries required; uses C++17 standard library only.
 *
 *  Security measures:
 *    - POST body read with hard byte cap (MAX_UPLOAD_BYTES).
 *    - Motif validated against [A-Za-z]+ whitelist.
 *    - FASTA validated line-by-line: binary-byte check + IUPAC character
 *      whitelist. Only one sequence held in memory at a time.
 *    - No shell calls; no temporary files written to disk.
 *    - All user-supplied text passed through htmlEscape() before output.
 *    - Multipart boundary validated/used only as a plain string search token.
 */

#include <algorithm>
#include <cctype>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include "cgi-page.h"   // shared UWPR page chrome and stylesheet

/* ──────────────────────────────────────────────────────────────
   CONFIGURATION
   ────────────────────────────────────────────────────────────── */
static const size_t MAX_UPLOAD_BYTES    = 50ULL * 1024ULL * 1024ULL;  // 50 MB
static const size_t MAX_SEQUENCES       = 500000;
static const size_t MAX_SEQ_BYTES       = 10ULL * 1024ULL * 1024ULL; // 10 MB per sequence
static const size_t MAX_RESULT_ENTRIES  = 2000000; // max hidden result_* fields
static const size_t MAX_BOUNDARY_LEN    = 256;     // RFC 2046 §5.1.1 max 70 chars; we allow 256
static const char*  APP_VERSION         = "1.2";

/* ──────────────────────────────────────────────────────────────
   DATA STRUCTURES
   ────────────────────────────────────────────────────────────── */
struct Hit {
    int         pos;    // 1-based residue position
    std::string match;
};

struct ProteinResult {
    std::string      desc;
    std::vector<Hit> hits;
};

struct ParsedForm {
    std::string motif;
    std::string fastaText;      // pasted text (may be empty)
    std::string fastaFilename;  // original filename from upload
    std::string fastaFileData;  // raw bytes of uploaded file
    bool        hasFile = false;

    // Pre-computed hit data carried in the download re-POST.
    // Parallel arrays: resultDesc[i], resultPos[i], resultMatch[i] = one hit.
    std::vector<std::string> resultDesc;
    std::vector<std::string> resultPos;
    std::vector<std::string> resultMatch;
    std::string              resultSource;  // original filename / "pasted input"
    std::string              resultTotal;   // total proteins scanned (as string)
    bool                     downloadGene = false; // gene-only download requested
};

/* ──────────────────────────────────────────────────────────────
   STRING UTILITIES
   ────────────────────────────────────────────────────────────── */

// Escape for safe HTML output
std::string htmlEscape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 16);
    for (unsigned char c : s) {
        switch (c) {
            case '&':  out += "&amp;";  break;
            case '<':  out += "&lt;";   break;
            case '>':  out += "&gt;";   break;
            case '"':  out += "&quot;"; break;
            case '\'': out += "&#39;";  break;
            default:   out += static_cast<char>(c); break;
        }
    }
    return out;
}

// Strip leading/trailing ASCII whitespace including \r
std::string trim(const std::string& s) {
    const std::string ws = " \t\r\n";
    size_t b = s.find_first_not_of(ws);
    if (b == std::string::npos) return "";
    size_t e = s.find_last_not_of(ws);
    return s.substr(b, e - b + 1);
}

// ASCII upper-case copy
std::string toUpper(std::string s) {
    for (char& c : s)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

// URL percent-decode (for application/x-www-form-urlencoded)
std::string urlDecode(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '+') {
            out += ' ';
        } else if (s[i] == '%' && i + 2 < s.size() &&
                   std::isxdigit(static_cast<unsigned char>(s[i+1])) &&
                   std::isxdigit(static_cast<unsigned char>(s[i+2]))) {
            out += static_cast<char>(std::stoi(s.substr(i + 1, 2), nullptr, 16));
            i += 2;
        } else {
            out += s[i];
        }
    }
    return out;
}

// UTC timestamp string
std::string currentTimestamp() {
    std::time_t t = std::time(nullptr);
    char buf[64];
    struct tm tmBuf{};
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S UTC", gmtime_r(&t, &tmBuf));
    return buf;
}

/* ──────────────────────────────────────────────────────────────
   POST BODY READER
   ────────────────────────────────────────────────────────────── */

// Read stdin up to MAX_UPLOAD_BYTES.  Returns false if body exceeded cap.
bool readPostBody(std::string& body) {
    const char* clEnv = std::getenv("CONTENT_LENGTH");
    size_t contentLength = 0;
    if (clEnv) {
        try { contentLength = static_cast<size_t>(std::stoull(clEnv)); }
        catch (...) { contentLength = 0; }
    }
    if (contentLength > MAX_UPLOAD_BYTES) return false;

    size_t toRead = (contentLength > 0) ? contentLength : MAX_UPLOAD_BYTES + 1;
    body.resize(toRead);
    std::cin.read(&body[0], static_cast<std::streamsize>(toRead));
    body.resize(static_cast<size_t>(std::cin.gcount()));
    return body.size() <= MAX_UPLOAD_BYTES;
}

/* ──────────────────────────────────────────────────────────────
   MULTIPART FORM-DATA PARSER
   ────────────────────────────────────────────────────────────── */

// Extract boundary token from Content-Type header value
std::string extractBoundary(const std::string& ct) {
    const std::string key = "boundary=";
    size_t pos = ct.find(key);
    if (pos == std::string::npos) return "";
    std::string b = ct.substr(pos + key.size());
    if (!b.empty() && b.front() == '"') b = b.substr(1);
    size_t end = b.find_first_of("\";\r\n \t");
    if (end != std::string::npos) b = b.substr(0, end);
    // FIX: cap boundary length to prevent O(n) string searches with huge boundary
    if (b.size() > MAX_BOUNDARY_LEN) return "";
    return b;
}

// Pull name= and filename= from a Content-Disposition header block
void parseMimeHeaders(const std::string& hdr,
                      std::string& name, std::string& filename) {
    name.clear(); filename.clear();
    std::regex nameRe("name=\"([^\"]*)\"",     std::regex::icase);
    std::regex fileRe("filename=\"([^\"]*)\"", std::regex::icase);
    std::smatch m;
    if (std::regex_search(hdr, m, nameRe)) name     = m[1].str();
    if (std::regex_search(hdr, m, fileRe)) filename = m[1].str();
}

// Parse multipart/form-data.  Extracts motif, fasta_text, fasta_file, download.
ParsedForm parseMultipart(const std::string& body,
                          const std::string& boundary,
                          bool& hasDownload) {
    ParsedForm form;
    hasDownload = false;

    const std::string CRLF      = "\r\n";
    const std::string delim     = "--" + boundary;
    const std::string delimEnd  = "--" + boundary + "--";
    const std::string hdrSep    = "\r\n\r\n";

    size_t pos = body.find(delim);
    if (pos == std::string::npos) return form;

    while (true) {
        // Advance past the boundary line + CRLF
        size_t lineEnd = body.find(CRLF, pos);
        if (lineEnd == std::string::npos) break;
        pos = lineEnd + 2;

        // Stop at end boundary
        size_t safePos = (pos >= 2 + delimEnd.size()) ? pos - 2 - delimEnd.size() : 0;
        if (body.compare(safePos, delimEnd.size(), delimEnd) == 0) break;
        // Also stop if the next thing IS the end boundary
        if (body.compare(pos, delimEnd.size(), delimEnd) == 0) break;

        // Locate end of MIME headers
        size_t hdrEnd = body.find(hdrSep, pos);
        if (hdrEnd == std::string::npos) break;
        std::string headers = body.substr(pos, hdrEnd - pos);
        pos = hdrEnd + 4;

        // Locate start of next boundary (data ends just before \r\n--)
        size_t nextBound = body.find(CRLF + delim, pos);
        if (nextBound == std::string::npos) break;
        std::string data = body.substr(pos, nextBound - pos);
        pos = nextBound + 2 + delim.size();

        std::string fieldName, fieldFile;
        parseMimeHeaders(headers, fieldName, fieldFile);

        if      (fieldName == "motif")         form.motif     = trim(data);
        else if (fieldName == "fasta_text")    form.fastaText = data;
        else if (fieldName == "download")      hasDownload    = true;
        else if (fieldName == "result_source") form.resultSource = data;
        else if (fieldName == "result_total")  form.resultTotal  = data;
        else if (fieldName == "result_desc" &&
                 form.resultDesc.size() < MAX_RESULT_ENTRIES)
                                               form.resultDesc.push_back(data);
        else if (fieldName == "result_pos" &&
                 form.resultPos.size() < MAX_RESULT_ENTRIES)
                                               form.resultPos.push_back(data);
        else if (fieldName == "result_match" &&
                 form.resultMatch.size() < MAX_RESULT_ENTRIES)
                                               form.resultMatch.push_back(data);
        else if (fieldName == "download_gene") form.downloadGene = true;
        else if (fieldName == "fasta_file" && !fieldFile.empty() && !data.empty()) {
            form.hasFile       = true;
            form.fastaFilename = fieldFile;
            form.fastaFileData = std::move(data);
        }

        // Check for end boundary immediately following
        if (pos + 2 <= body.size() && body[pos] == '-' && body[pos+1] == '-') break;
    }
    return form;
}

// Parse application/x-www-form-urlencoded (used by the download re-POST)
ParsedForm parseUrlEncoded(const std::string& body, bool& hasDownload) {
    ParsedForm form;
    hasDownload = false;
    std::istringstream ss(body);
    std::string pair;
    while (std::getline(ss, pair, '&')) {
        size_t eq = pair.find('=');
        if (eq == std::string::npos) continue;
        std::string k = urlDecode(pair.substr(0, eq));
        std::string v = urlDecode(pair.substr(eq + 1));
        if      (k == "motif")         form.motif        = trim(v);
        else if (k == "fasta_text")    form.fastaText    = v;
        else if (k == "download")      hasDownload       = true;
        else if (k == "result_source") form.resultSource = v;
        else if (k == "result_total")  form.resultTotal  = v;
        else if (k == "result_desc" &&
                 form.resultDesc.size() < MAX_RESULT_ENTRIES)
                                       form.resultDesc.push_back(v);
        else if (k == "result_pos" &&
                 form.resultPos.size() < MAX_RESULT_ENTRIES)
                                       form.resultPos.push_back(v);
        else if (k == "result_match" &&
                 form.resultMatch.size() < MAX_RESULT_ENTRIES)
                                       form.resultMatch.push_back(v);
        else if (k == "download_gene")  form.downloadGene = true;
    }
    return form;
}

/* ──────────────────────────────────────────────────────────────
   VALIDATION
   ────────────────────────────────────────────────────────────── */

bool validateMotif(const std::string& raw, std::string& out, std::string& err) {
    out = toUpper(trim(raw));
    if (out.empty())   { err = "Motif cannot be empty."; return false; }
    if (out.size()>100){ err = "Motif is too long (max 100 characters)."; return false; }
    for (char c : out) {
        if (!std::isalpha(static_cast<unsigned char>(c))) {
            err = "Motif may only contain letters A\u2013Z (use X as wildcard).";
            return false;
        }
    }
    return true;
}

// IUPAC amino acid characters (standard 20 + ambiguous + stop * + gap -)
bool isValidAA(unsigned char c) {
    static const std::string valid =
        "ACDEFGHIKLMNPQRSTVWYBZXUOacdefghiklmnpqrstvwybzxuo*-";
    return valid.find(static_cast<char>(c)) != std::string::npos;
}

// Per-line binary/control byte check
bool hasBinaryBytes(const std::string& line) {
    for (unsigned char c : line)
        if (c<=0x08 || c==0x0B || c==0x0C || (c>=0x0E&&c<=0x1F) || c==0x7F)
            return true;
    return false;
}

/* ──────────────────────────────────────────────────────────────
   MOTIF → REGEX PATTERN
   ────────────────────────────────────────────────────────────── */

// X/x → [A-Za-z]   all other letters → literal (regex-escaped for safety)
std::string motifToPattern(const std::string& motif) {
    std::string pat;
    static const std::string meta = "^$.|?*+()[]{}\\";
    for (char c : motif) {
        if (c == 'X' || c == 'x') {
            pat += "[A-Za-z]";
        } else {
            if (meta.find(c) != std::string::npos) pat += '\\';
            pat += c;
        }
    }
    return pat;
}

/* ──────────────────────────────────────────────────────────────
   FASTA STREAMING SEARCH
   ────────────────────────────────────────────────────────────── */

/*
 * Read source line-by-line.  For each complete sequence:
 *   1. Validate characters.
 *   2. Search with pattern; collect hits.
 *   3. Immediately discard the sequence string before reading the next one.
 *
 * Returns false on validation error, setting `error`.
 */
bool streamSearchFasta(std::istream&               source,
                       const std::regex&            pattern,
                       size_t                       motifLen,
                       std::vector<ProteinResult>&  results,
                       size_t&                      total,
                       std::string&                 error) {
    std::string desc;
    std::string seq;
    bool        seenHeader = false;
    size_t      lineNum    = 0;

    // Search the buffered sequence and push any hits into results.
    // Clears seq immediately afterward to free memory.
    auto flushSeq = [&]() -> bool {
        if (desc.empty()) return true;
        if (seq.empty()) {
            error = "Sequence for '" + desc + "' is empty.";
            return false;
        }
        total++;
        if (total > MAX_SEQUENCES) {
            error = "File contains too many sequences (max " +
                    std::to_string(MAX_SEQUENCES) + ").";
            return false;
        }

        std::string upper = toUpper(seq);
        seq.clear();
        seq.shrink_to_fit();   // release memory now

        ProteinResult pr;
        pr.desc = desc;
        auto it  = std::sregex_iterator(upper.begin(), upper.end(), pattern);
        auto end = std::sregex_iterator();
        for (; it != end; ++it) {
            Hit h;
            h.pos   = static_cast<int>(it->position(0)) + 1;
            h.match = upper.substr(
                static_cast<size_t>(it->position(0)), motifLen);
            pr.hits.push_back(std::move(h));
        }
        if (!pr.hits.empty()) results.push_back(std::move(pr));
        return true;
    };

    std::string raw;
    while (std::getline(source, raw)) {
        lineNum++;
        if (!raw.empty() && raw.back() == '\r') raw.pop_back();  // strip \r
        if (raw.empty()) continue;

        if (hasBinaryBytes(raw)) {
            error = "Line " + std::to_string(lineNum) +
                    ": binary or non-printable content \u2014 not a FASTA file.";
            return false;
        }

        if (raw[0] == '>') {
            if (!flushSeq()) return false;
            seenHeader = true;
            std::string d = trim(raw.substr(1));
            desc = d.empty() ? ("Sequence_" + std::to_string(total + 1)) : d;
        } else {
            if (!seenHeader) {
                error = "Line " + std::to_string(lineNum) +
                        ": sequence data before any FASTA header.";
                return false;
            }
            for (size_t i = 0; i < raw.size(); ++i) {
                if (!isValidAA(static_cast<unsigned char>(raw[i]))) {
                    error = "Invalid character in sequence near line " +
                            std::to_string(lineNum) + ": \"" +
                            raw.substr(0, std::min(raw.size(), size_t(40))) + "\"";
                    return false;
                }
            }
            seq += raw;
            // FIX: cap single-sequence size to prevent memory exhaustion
            if (seq.size() > MAX_SEQ_BYTES) {
                error = "Sequence for '" + desc + "' exceeds the per-sequence "
                        "size limit (" + std::to_string(MAX_SEQ_BYTES / (1024*1024))
                        + " MB).";
                return false;
            }
        }
    }

    if (!flushSeq()) return false;

    if (total == 0) {
        error = "No valid FASTA records found in the input.";
        return false;
    }
    return true;
}

/* ──────────────────────────────────────────────────────────────
   HTTP OUTPUT HELPERS
   ────────────────────────────────────────────────────────────── */

void printHeaders(const std::string& ct = "text/html; charset=utf-8") {
    std::cout << "Content-Type: " << ct << "\r\n"
              // FIX: add security headers on all HTML responses
              << "X-Content-Type-Options: nosniff\r\n"
              << "X-Frame-Options: DENY\r\n"
              << "Referrer-Policy: strict-origin-when-cross-origin\r\n"
              << "Content-Security-Policy: default-src 'self'; "
                 "style-src 'self' 'unsafe-inline'; "
                 "script-src 'unsafe-inline'; "
                 "frame-ancestors 'none'\r\n"
              << "\r\n";
}

void printDownloadHeaders(const std::string& filename) {
    std::cout << "Content-Type: text/plain; charset=utf-8\r\n"
              << "Content-Disposition: attachment; filename=\""
              << filename << "\"\r\n"
              // FIX: add security headers on download responses
              << "X-Content-Type-Options: nosniff\r\n"
              << "X-Frame-Options: DENY\r\n"
              << "Content-Security-Policy: default-src 'none'\r\n"
              << "\r\n";
}

/* ──────────────────────────────────────────────────────────────
   HTML PAGE SECTIONS

   The page chrome (top bar, card, footer) and the shared stylesheet
   come from cgi-page.h, the same header used by the other UWPR CGI
   tools.  Only the few motif-specific rules live here.
   ────────────────────────────────────────────────────────────── */

static const char* g_szMotifCSS =
"   <style>\n"
"      /* motif search additions to the shared stylesheet */\n"
"      .error-box { background: #fdf3f3; border: 1px solid #e3b4b4; border-left: 4px solid #b23a3a; border-radius: var(--radius); padding: .9rem 1.1rem; margin: 0 0 1.5rem; color: #7a1010; font-size: .92rem; }\n"
"      .error-box p { margin: 0 0 .3rem; font-weight: 600; }\n"
"      .error-box ul { margin: 0; padding-left: 1.3em; }\n"
"      input[type=text].motif { width: 100%; max-width: 20rem; font-size: 1.15rem; letter-spacing: .12em; text-transform: uppercase; }\n"
"      input[type=file] { font: inherit; font-size: .9rem; color: var(--text); max-width: 100%; }\n"
"      input[type=file]::file-selector-button { font: inherit; font-size: .85rem; font-weight: 600; color: var(--accent); background: var(--surface); border: 1px solid var(--border-strong); border-radius: 7px; padding: .35rem .8rem; margin-right: .75rem; cursor: pointer; }\n"
"      input[type=file]::file-selector-button:hover { background: var(--accent-soft); }\n"
"      .or-divider { display: flex; align-items: center; gap: .75rem; margin: 1rem 0; color: var(--muted); font-size: .78rem; font-weight: 600; letter-spacing: .08em; text-transform: uppercase; }\n"
"      .or-divider::before, .or-divider::after { content: ''; flex: 1; height: 1px; background: var(--border-strong); }\n"
"      .btn.secondary { color: var(--accent); background: var(--surface); border-color: var(--accent); }\n"
"      .btn.secondary:hover { background: var(--accent-soft); }\n"
"      .actions form { display: inline; }\n"
"      table.hits { width: 100%; }\n"
"      table.hits th.pos, table.hits td.pos { width: 7rem; }\n"
"      table.hits tbody th { font-family: var(--mono); font-size: .84rem; font-weight: 600; color: var(--accent-dark); background: var(--surface-2); border-top: 1px solid var(--border-strong); border-bottom: 1px solid var(--border); white-space: normal; word-break: break-all; letter-spacing: 0; }\n"
"      table.hits tbody tr:hover th { background: var(--surface-2); }\n"
"      table.hits td.match { letter-spacing: .08em; }\n"
"      /* wildcard residue: colour plus dotted underline so it is not colour alone */\n"
"      .wc { background: #fbf8ee; color: #6b5a2a; border-radius: 2px; padding: 0 .1em; font-weight: 600; text-decoration: underline dotted; }\n"
"      .empty { background: var(--surface-2); border: 1px dashed var(--border-strong); border-radius: var(--radius); padding: 1.5rem; text-align: center; color: var(--muted); }\n"
"   </style>\n";

void printPageTop(const char* version) {
    std::string extra = std::string("   <!-- motif_search v") + version + " -->\n" + g_szMotifCSS;
    // cgi-page.h writes with printf; std::cout and stdio are kept in sync
    // (sync_with_stdio is left at its default), so the two may be mixed.
    std::cout.flush();
    PRINT_PAGE_HEADER("Motif Hunter", extra.c_str());
    fflush(stdout);

    std::cout <<
        "    <div id=\"page\" class=\"container\">\n"
        "       <section>\n"
        "          <header class=\"major\">\n"
        "             <h1>Motif hunter</h1>\n"
        "             <p class=\"lede\">Scan protein sequences in FASTA format for a short sequence motif. "
        "Use X as a single-residue wildcard. Each match is reported with its position, "
        "and the results can be downloaded as a text report.</p>\n"
        "          </header>\n\n";
}

void printErrors(const std::vector<std::string>& errors) {
    if (errors.empty()) return;
    std::cout <<
        "<div class=\"error-box\" role=\"alert\" "
        "aria-live=\"assertive\" aria-atomic=\"true\">\n"
        "  <p>Please correct the following:</p>\n"
        "  <ul>\n";
    for (const auto& e : errors)
        std::cout << "    <li>" << htmlEscape(e) << "</li>\n";
    std::cout << "  </ul>\n</div>\n";
}

void printForm(const std::string& motifVal, const std::string& fastaText) {
    std::cout <<
        "<form method=\"POST\" enctype=\"multipart/form-data\" "
        "id=\"searchForm\" aria-label=\"Protein motif search\">\n"
        "<div class=\"form-grid\">\n\n"

        // ── Motif panel ─────────────────────────────────────────────
        "<div class=\"panel\">\n"
        "<fieldset class=\"options\">\n"
        "<legend>Motif</legend>\n"
        "<input type=\"text\" id=\"motif\" name=\"motif\" class=\"mono motif\"\n"
        "       placeholder=\"FLXLFX\"\n"
        "       value=\"" << htmlEscape(motifVal) << "\"\n"
        "       autocomplete=\"off\" spellcheck=\"false\" maxlength=\"100\"\n"
        "       aria-describedby=\"motif-hint\" aria-required=\"true\" required>\n"
        "<p class=\"hint\" id=\"motif-hint\">Letters A&#8211;Z only. X matches any single amino acid. "
        "<br>Examples: <code>FLXLFX</code>, <code>GXXGXG</code>, <code>RGD</code>, <code>WXXXW</code></p>\n"
        "</fieldset>\n"
        "<div class=\"actions\"><input type=\"submit\" name=\"search\" value=\"Search\"></div>\n"
        "</div>\n\n"  // panel

        // ── Sequences panel ─────────────────────────────────────────
        "<div class=\"panel\">\n"
        "<label class=\"field-label\" for=\"fasta_text\">Protein sequences (FASTA)</label>\n"
        "<p class=\"hint\" style=\"margin: 0 0 .5rem\">Each sequence must start with a "
        "<code>&gt;</code> header line.</p>\n"
        "<textarea id=\"fasta_text\" name=\"fasta_text\" class=\"wide mono\" rows=\"10\"\n"
        "          aria-describedby=\"fasta-paste-hint\" spellcheck=\"false\"\n"
        "          placeholder=\"&gt;sp|P12345|MYO_HUMAN Myosin heavy chain"
        "&#10;MSSTKIHLEQHVKEIDISQ...\">"
        << htmlEscape(fastaText) <<
        "</textarea>\n"
        "<p class=\"hint\" id=\"fasta-paste-hint\">Paste sequences above, or upload a file below.</p>\n"
        "<div class=\"or-divider\" aria-hidden=\"true\">or</div>\n"
        "<div class=\"field\">\n"
        "<label class=\"field-label\" for=\"fasta_file\">Upload a FASTA file</label>\n"
        "<input type=\"file\" id=\"fasta_file\" name=\"fasta_file\"\n"
        "       accept=\".fasta,.fa,.faa,.txt\" aria-describedby=\"file-hint\">\n"
        "<p class=\"hint\" id=\"file-hint\">Accepted: .fasta .fa .faa .txt &nbsp;&middot;&nbsp; "
        "Max 50&nbsp;MB. <br>An uploaded file takes precedence over pasted text.</p>\n"
        "</div>\n"
        "</div>\n\n"  // panel

        "</div>\n"    // form-grid
        "</form>\n\n";
}

// Render one matched sequence string with wildcard positions highlighted
std::string renderMatch(const std::string& match, const std::string& motif) {
    std::string out;
    for (size_t i = 0; i < match.size() && i < motif.size(); ++i) {
        std::string aa = htmlEscape(std::string(1, match[i]));
        if (motif[i] == 'X') {
            out += "<span class=\"wc\" aria-label=\"wildcard residue "
                +  aa + "\">" + aa + "</span>";
        } else {
            out += aa;
        }
    }
    return out;
}

void printResults(const std::vector<ProteinResult>& results,
                  const std::string& motif,
                  size_t totalProteins, size_t totalHits,
                  const std::string& fastaSource) {

    std::cout <<
        "<dl class=\"summary\" aria-label=\"Search summary\">\n"
        "  <dt>Motif</dt><dd class=\"stat-val\">" << htmlEscape(motif) << "</dd>\n"
        "  <dt>Source</dt><dd class=\"stat-val\">" << htmlEscape(fastaSource) << "</dd>\n"
        "  <dt>Proteins scanned</dt><dd class=\"stat-val\">" << totalProteins << "</dd>\n"
        "  <dt>Proteins with hits</dt><dd class=\"stat-val\">" << results.size() << "</dd>\n"
        "  <dt>Matches found</dt><dd class=\"stat-val\">" << totalHits << "</dd>\n"
        "</dl>\n\n";

    if (results.empty()) {
        std::cout <<
            "<div class=\"empty\" role=\"status\">\n"
            "  No sequences matched the motif <strong>"
            << htmlEscape(motif) << "</strong>.\n"
            "  Check the motif pattern or try a different FASTA input.\n"
            "</div>\n";
        return;
    }

    // ── Download buttons ────────────────────────────────────────────
    // Shared hidden fields are repeated in each form because browsers only
    // submit fields belonging to the form that was submitted.
    // Each form is self-contained so either button works independently.
    auto emitHiddenFields = [&]() {
        std::cout
            << "    <input type=\"hidden\" name=\"motif\" value=\""
            << htmlEscape(motif) << "\">\n"
            << "    <input type=\"hidden\" name=\"result_source\" value=\""
            << htmlEscape(fastaSource) << "\">\n"
            << "    <input type=\"hidden\" name=\"result_total\" value=\""
            << totalProteins << "\">\n";
        for (const auto& r : results) {
            for (const auto& h : r.hits) {
                std::cout
                    << "    <input type=\"hidden\" name=\"result_desc\" value=\""
                    << htmlEscape(r.desc) << "\">\n"
                    << "    <input type=\"hidden\" name=\"result_pos\" value=\""
                    << h.pos << "\">\n"
                    << "    <input type=\"hidden\" name=\"result_match\" value=\""
                    << htmlEscape(h.match) << "\">\n";
            }
        }
    };

    std::cout << "<div class=\"actions\" style=\"margin: 0 0 1.25rem\">\n";

    std::cout <<
        "  <form method=\"POST\" enctype=\"multipart/form-data\" "
        "aria-label=\"Download full results\">\n"
        "    <input type=\"hidden\" name=\"download\" value=\"1\">\n";
    emitHiddenFields();
    std::cout <<
        "    <button type=\"submit\" class=\"btn secondary\">"
        "Download results (.txt)</button>\n"
        "  </form>\n";

    std::cout <<
        "  <form method=\"POST\" enctype=\"multipart/form-data\" "
        "aria-label=\"Download protein and gene names only\">\n"
        "    <input type=\"hidden\" name=\"download_gene\" value=\"1\">\n";
    emitHiddenFields();
    std::cout <<
        "    <button type=\"submit\" class=\"btn secondary\">"
        "Download protein &amp; gene list (.txt)</button>\n"
        "  </form>\n";

    std::cout << "</div>\n\n";  // actions

    // ── Hit table: one tbody per protein, description as a group row ──
    std::cout <<
        "<p class=\"note\">Highlighted residues are at wildcard (X) positions of the motif.</p>\n"
        "<div class=\"table-wrap\">\n"
        "<table class=\"results hits\">\n"
        "<thead><tr>"
        "<th class=\"num pos\" scope=\"col\">position</th>"
        "<th scope=\"col\">matched sequence</th>"
        "</tr></thead>\n";

    for (const auto& r : results) {
        std::cout <<
            "<tbody class=\"hit-card\">\n"
            "<tr><th colspan=\"2\" scope=\"rowgroup\">" << htmlEscape(r.desc) << "</th></tr>\n";
        for (const auto& h : r.hits) {
            std::cout <<
                "<tr>"
                "<td class=\"num pos\">" << h.pos << "</td>"
                "<td class=\"match\">" << renderMatch(h.match, motif) << "</td>"
                "</tr>\n";
        }
        std::cout << "</tbody>\n";
    }

    std::cout <<
        "</table>\n"
        "</div>\n";     // table-wrap
}

void printPageBottom() {
    std::cout <<
        "       </section>\n"
        "    </div>\n"  // page
        "<script>\n"
        "/* Auto-uppercase motif field */\n"
        "var motifInput = document.getElementById('motif');\n"
        "motifInput.addEventListener('input', function () {\n"
        "  var pos = motifInput.selectionStart;\n"
        "  motifInput.value = motifInput.value.toUpperCase();\n"
        "  motifInput.setSelectionRange(pos, pos);\n"
        "});\n"
        "/* Move focus to results after submit */\n"
        "var rs = document.getElementById('results');\n"
        "if (rs && rs.children.length) { rs.setAttribute('tabindex','-1'); rs.focus(); }\n"
        "</script>\n";
    std::cout.flush();
    PRINT_PAGE_FOOTER();
    fflush(stdout);
}


/* ──────────────────────────────────────────────────────────────
   PROTEIN / GENE NAME EXTRACTION
   ────────────────────────────────────────────────────────────── */

/*
 * Parse a FASTA description line into an accession token and a gene name.
 *
 * Accession rules:
 *   If the first word (whitespace-delimited) contains '|' characters,
 *   split on '|' and take the third field (index 2).
 *   e.g. "sp|P01893|HLAH_HUMAN ..." → "HLAH_HUMAN"
 *   Otherwise the first word itself is the accession.
 *   e.g. "MY_PROTEIN description ..." → "MY_PROTEIN"
 *
 * Gene name rules:
 *   Search for " GN=<token> " anywhere in the description (case-sensitive,
 *   standard UniProt format). The token ends at the next space or end-of-string.
 *   Returns "" if no GN= field is present.
 */
std::pair<std::string, std::string> parseProteinGene(const std::string& desc) {
    // ── accession ──────────────────────────────────────────────
    std::string accession;
    size_t spacePos = desc.find_first_of(" \t");
    std::string firstWord = (spacePos == std::string::npos)
                            ? desc : desc.substr(0, spacePos);

    size_t pipe1 = firstWord.find('|');
    if (pipe1 != std::string::npos) {
        size_t pipe2 = firstWord.find('|', pipe1 + 1);
        if (pipe2 != std::string::npos) {
            // Third field exists
            accession = firstWord.substr(pipe2 + 1);
        } else {
            // Only one pipe: take second field
            accession = firstWord.substr(pipe1 + 1);
        }
    } else {
        accession = firstWord;
    }

    // ── gene name ──────────────────────────────────────────────
    std::string gene;
    const std::string tag = " GN=";
    size_t gnPos = desc.find(tag);
    if (gnPos != std::string::npos) {
        size_t start = gnPos + tag.size();
        size_t end   = desc.find(' ', start);
        gene = (end == std::string::npos)
               ? desc.substr(start)
               : desc.substr(start, end - start);
    }

    return {accession, gene};
}

/*
 * Stream a gene-only download report: one line per unique protein with hits,
 * formatted as "<accession>\t<gene>" (tab-separated; gene may be empty).
 * Proteins are de-duplicated — each appears only once regardless of hit count.
 */
void sendDownloadGeneOnly(const std::vector<ProteinResult>& results,
                          const std::string& motif,
                          const std::string& source) {
    std::string safe;
    for (char c : motif)
        safe += std::isalnum(static_cast<unsigned char>(c)) ? c : '_';

    printDownloadHeaders("motif_" + safe + "_gene_results.txt");

    std::cout
        << "# Protein Gene Report\n"
        << "# Generated : " << currentTimestamp() << "\n"
        << "# Source    : " << source << "\n"
        << "# Motif     : " << motif << "  (X = any amino acid)\n"
        << "# Columns   : Accession<TAB>Gene (GN= field; blank if absent)\n"
        << "# Proteins with hits: " << results.size() << "\n"
        << std::string(72, '-') << "\n";

    for (const auto& r : results) {
        auto [accession, gene] = parseProteinGene(r.desc);
        std::cout << accession << "\t" << gene << "\n";
    }
}

/* ──────────────────────────────────────────────────────────────
   DOWNLOAD MODE
   ────────────────────────────────────────────────────────────── */

void sendDownload(const std::vector<ProteinResult>& results,
                  const std::string& motif,
                  size_t totalProteins,
                  const std::string& source) {
    std::string safe;
    for (char c : motif)
        safe += std::isalnum(static_cast<unsigned char>(c)) ? c : '_';

    printDownloadHeaders("motif_" + safe + "_results.txt");

    size_t totalHits = 0;
    for (const auto& r : results) totalHits += r.hits.size();

    std::cout
        << "# Protein Motif Search Report\n"
        << "# Generated : " << currentTimestamp() << "\n"
        << "# Source    : " << source << "\n"
        << "# Motif     : " << motif << "  (X = any amino acid)\n"
        << "# Proteins scanned  : " << totalProteins << "\n"
        << "# Proteins with hits: " << results.size() << "\n"
        << "# Total matches     : " << totalHits << "\n"
        << std::string(72, '-') << "\n\n";

    for (const auto& r : results) {
        std::cout << ">" << r.desc << "\n";
        for (const auto& h : r.hits)
            std::cout << "  Position " << std::setw(6) << h.pos
                      << " : " << h.match << "\n";
        std::cout << "\n";
    }
}

/* ──────────────────────────────────────────────────────────────
   MAIN
   ────────────────────────────────────────────────────────────── */
int main() {

    const char* requestMethod  = std::getenv("REQUEST_METHOD");
    const char* contentTypeEnv = std::getenv("CONTENT_TYPE");
    std::string contentType    = contentTypeEnv ? contentTypeEnv : "";

    /* ── GET: serve blank form ── */
    if (!requestMethod || std::string(requestMethod) != "POST") {
        printHeaders();
        printPageTop(APP_VERSION);
        printForm("", "");
        std::cout << "<div id=\"results\" aria-live=\"polite\" aria-label=\"Search results\"></div>\n\n";
        printPageBottom();
        return 0;
    }

    /* ── POST ── */
    std::string body;
    bool sizeOk = readPostBody(body);

    std::vector<std::string> errors;
    std::string motifDisplay;
    std::vector<ProteinResult> results;
    size_t totalProteins = 0;
    size_t totalHits     = 0;
    std::string fastaSource;
    bool hasDownload = false;

    if (!sizeOk) {
        errors.push_back("Uploaded data exceeds the 50 MB size limit.");
    }

    // Parse form
    ParsedForm form;
    std::string boundary = extractBoundary(contentType);
    if (!boundary.empty()) {
        form = parseMultipart(body, boundary, hasDownload);
    } else {
        form = parseUrlEncoded(body, hasDownload);
    }
    body.clear();
    body.shrink_to_fit();   // free POST body — we have what we need

    // 1. Validate motif
    std::string motifErr;
    if (errors.empty()) {
        if (!validateMotif(form.motif, motifDisplay, motifErr))
            errors.push_back(motifErr);
    }

    // 2. Determine FASTA source.
    //    Skip entirely for download re-POSTs that carry pre-computed hit fields —
    //    those posts never include FASTA data and don't need it.
    bool hasInput = false;
    std::string fastaData;
    const bool isDownloadRepost = (hasDownload || form.downloadGene) && !form.resultDesc.empty();

    if (errors.empty() && !isDownloadRepost) {
        if (form.hasFile && !form.fastaFileData.empty()) {
            fastaData   = std::move(form.fastaFileData);
            fastaSource = form.fastaFilename;
            hasInput    = true;
        } else if (!trim(form.fastaText).empty()) {
            fastaData   = form.fastaText;
            fastaSource = "pasted input";
            hasInput    = true;
        } else {
            errors.push_back(
                "Please paste a FASTA sequence or upload a FASTA file.");
        }
    }

    // 3. Stream-search (skipped for download re-POSTs)
    std::string streamError;
    if (errors.empty() && hasInput && !isDownloadRepost) {
        std::istringstream stream(fastaData);
        fastaData.clear();
        fastaData.shrink_to_fit();

        std::string patStr = motifToPattern(motifDisplay);
        std::regex  pattern;
        try {
            pattern = std::regex(patStr, std::regex_constants::icase);
        } catch (const std::regex_error& re) {
            errors.push_back(std::string("Internal regex error: ") + re.what());
        }

        if (errors.empty()) {
            bool ok = streamSearchFasta(stream, pattern,
                                        motifDisplay.size(),
                                        results, totalProteins, streamError);
            if (!ok)
                errors.push_back("FASTA validation failed: " + streamError);
            else
                for (const auto& r : results) totalHits += r.hits.size();
        }
    }

    // 4. Download mode
    //    Fast path: if the form contains pre-computed result hidden fields
    //    (result_desc/pos/match), reconstruct ProteinResult objects from them
    //    and stream the report — no FASTA re-parsing needed.
    const bool wantDownload     = hasDownload || form.downloadGene;
    if (wantDownload && errors.empty()) {
        // Check whether we have pre-computed hit data from hidden fields
        if (!form.resultDesc.empty() &&
            form.resultDesc.size() == form.resultPos.size() &&
            form.resultDesc.size() == form.resultMatch.size()) {

            // Reconstruct ProteinResult vector from parallel arrays
            std::vector<ProteinResult> dlResults;
            for (size_t i = 0; i < form.resultDesc.size(); ++i) {
                Hit h;
                try {
                    long long v = std::stoll(form.resultPos[i]);
                    // FIX: clamp to a sane range; reject obviously bogus values
                    if (v < 0 || v > 100000000LL) v = 0;
                    h.pos = static_cast<int>(v);
                }
                catch (...) { h.pos = 0; }
                h.match = form.resultMatch[i];

                // Group hits under the same desc into one ProteinResult
                if (!dlResults.empty() && dlResults.back().desc == form.resultDesc[i]) {
                    dlResults.back().hits.push_back(std::move(h));
                } else {
                    ProteinResult pr;
                    pr.desc = form.resultDesc[i];
                    pr.hits.push_back(std::move(h));
                    dlResults.push_back(std::move(pr));
                }
            }

            size_t dlTotal = 0;
            try { dlTotal = std::stoull(form.resultTotal); } catch (...) {}
            std::string dlSource = form.resultSource.empty()
                                   ? "unknown" : form.resultSource;

            if (form.downloadGene) {
                sendDownloadGeneOnly(dlResults, motifDisplay, dlSource);
            } else {
                sendDownload(dlResults, motifDisplay, dlTotal, dlSource);
            }
            return 0;
        }
        // If no pre-computed hits but download was requested, fall through
        // to the normal search path (handles the pasted-text re-POST case).
    }


    // 5. Render HTML response
    printHeaders();
    printPageTop(APP_VERSION);
    printErrors(errors);
    printForm(motifDisplay, form.fastaText);

    // The results container is always present (empty on a blank form),
    // matching the other tools; #results:empty collapses its top margin.
    std::cout << "<div id=\"results\" aria-live=\"polite\" aria-label=\"Search results\">\n";
    if (errors.empty() && totalProteins > 0)
        printResults(results, motifDisplay, totalProteins, totalHits, fastaSource);
    std::cout << "</div>\n\n";  // results

    printPageBottom();
    return 0;
}
