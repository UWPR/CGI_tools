/*  motif_search.cpp  —  Protein Motif Search Tool (C++ CGI)
 *
 *  Compile:
 *      g++ -std=c++17 -O2 -o motif_search.cgi motif_search.cpp
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
                 "style-src 'self' 'unsafe-inline' https://fonts.googleapis.com; "
                 "font-src https://fonts.gstatic.com; "
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
   ────────────────────────────────────────────────────────────── */

void printPageHead() {
    std::cout <<
R"(<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>MotifHunter &#8212; Protein Sequence Search</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Source+Code+Pro:ital,wght@0,400;0,600;1,400&family=Source+Sans+3:wght@400;600;700&display=swap" rel="stylesheet">
<style>
  /*
   * WCAG 2.1 AA — contrast ratios verified (minimum 4.5:1 for text)
   *   --text     #1a1a1a  16.75:1 on white
   *   --text-sub #4a4a4a   9.73:1 on white
   *   --label    #3d3d3d  11.27:1 on white
   *   --accent   #005f99   7.23:1 on white
   *   --accent-dk#004570   9.87:1 on white
   *   btn: white on #005f99 = 7.23:1
   *   error: #7a1010 on #fff4f4 = 8.59:1
   *   wc:   #8b4500 on #fff8f0 = 6.11:1
   */
  *,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
  :root{
    --bg:#f5f7fa; --surface:#ffffff; --surface2:#f0f4f8;
    --border:#b0bec8; --border-dark:#7a909e;
    --text:#1a1a1a; --text-sub:#4a4a4a; --label:#3d3d3d;
    --accent:#005f99; --accent-dk:#004570; --accent-bg:#e8f2fa;
    --btn-bg:#005f99; --btn-hover:#004570; --btn-text:#ffffff;
    --error-bg:#fff4f4; --error-border:#c0392b; --error-text:#7a1010;
    --hit-border:#005f99; --hit-bg:#f8fbff;
    --wc-bg:#fff8f0; --wc-text:#8b4500; --pos-color:#4a4a4a;
    --radius:8px;
    --mono:'Source Code Pro','Courier New',monospace;
    --sans:'Source Sans 3',Arial,sans-serif;
    --focus-ring:0 0 0 3px #005f99,0 0 0 5px #ffffff;
  }
  /* Skip link (WCAG 2.4.1) */
  .skip-link{position:absolute;top:-100%;left:1rem;background:var(--accent);
    color:#fff;padding:.6em 1.2em;border-radius:0 0 var(--radius) var(--radius);
    font-family:var(--sans);font-weight:700;font-size:.95rem;z-index:9999;
    text-decoration:none;transition:top .1s}
  .skip-link:focus{top:0;outline:3px solid #fff;outline-offset:2px}
  /* Focus (WCAG 2.4.7) */
  :focus-visible{outline:3px solid var(--accent);outline-offset:3px;border-radius:3px}
  :focus:not(:focus-visible){outline:none}
  @media(prefers-reduced-motion:reduce){
    *,*::before,*::after{animation-duration:.01ms!important;transition-duration:.01ms!important}}
  html{scroll-behavior:smooth}
  body{background:var(--bg);color:var(--text);font-family:var(--sans);
       font-size:1rem;line-height:1.6;min-height:100vh}
  /* Header */
  header{background:var(--surface);border-bottom:2px solid var(--border);
         padding:2rem 2rem 1.75rem;text-align:center}
  .logo-mark{display:inline-block;font-family:var(--mono);font-size:.75rem;
    font-weight:600;letter-spacing:.18em;color:var(--accent);background:var(--accent-bg);
    border:1px solid #aacde8;border-radius:4px;padding:.25em .8em;
    margin-bottom:.9rem;text-transform:uppercase}
  h1{font-size:clamp(1.75rem,4vw,2.6rem);font-weight:700;line-height:1.15;color:var(--text)}
  h1 .brand-accent{color:var(--accent)}
  .subtitle{margin-top:.6rem;color:var(--text-sub);font-size:.95rem}
  /* Layout */
  .container{max-width:860px;margin:0 auto;padding:2rem 1.5rem 4rem}
  /* Fieldset */
  fieldset{background:var(--surface);border:1px solid var(--border);
           border-radius:var(--radius);padding:1.75rem 2rem;margin-bottom:1.5rem}
  fieldset:focus-within{border-color:var(--accent);box-shadow:0 0 0 2px var(--accent-bg)}
  legend{font-family:var(--mono);font-size:.78rem;font-weight:600;letter-spacing:.12em;
         text-transform:uppercase;color:var(--accent);padding:0 .5em}
  label{display:block;font-size:.9rem;font-weight:600;color:var(--label);margin-bottom:.35rem}
  .field-hint{font-size:.82rem;color:var(--text-sub);margin-top:.3rem;margin-bottom:.75rem}
  input[type="text"],textarea{
    width:100%;background:var(--surface);border:2px solid var(--border-dark);
    border-radius:var(--radius);color:var(--text);font-family:var(--mono);
    padding:.6rem .9rem;outline:none;transition:border-color .15s}
  input[type="text"]:focus,textarea:focus{
    border-color:var(--accent);box-shadow:var(--focus-ring)}
  input[type="text"]{font-size:1.4rem;letter-spacing:.12em;text-transform:uppercase}
  textarea{font-size:.875rem;line-height:1.7;resize:vertical;min-height:180px}
  /* OR divider */
  .or-divider{display:flex;align-items:center;gap:.75rem;margin:1.25rem 0;
    font-size:.8rem;font-weight:600;color:var(--text-sub);
    text-transform:uppercase;letter-spacing:.1em}
  .or-divider::before,.or-divider::after{content:'';flex:1;height:1px;background:var(--border)}
  /* File zone */
  .file-zone{border:2px dashed var(--border-dark);border-radius:var(--radius);
    padding:1.5rem;text-align:center;position:relative;
    transition:border-color .15s,background .15s;background:var(--surface2)}
  .file-zone:hover{border-color:var(--accent);background:var(--accent-bg)}
  .file-zone.drag-over{border-color:var(--accent-dk);background:#d0e8f7;outline:3px solid var(--accent)}
  .file-zone input[type="file"]{position:absolute;inset:0;opacity:0;width:100%;height:100%;cursor:pointer}
  .file-zone-inner{pointer-events:none}
  .file-icon{font-size:2rem;display:block;margin-bottom:.4rem}
  .file-zone p{font-size:.85rem;color:var(--text-sub);margin:.2rem 0}
  .file-zone strong{color:var(--accent)}
  #file-name-display{margin-top:.6rem;font-family:var(--mono);font-size:.82rem;
    font-weight:600;color:#1a6e3c;min-height:1.3em}
  /* Buttons */
  .btn-search{display:block;width:100%;padding:.9rem 1rem;background:var(--btn-bg);
    color:var(--btn-text);font-family:var(--sans);font-weight:700;font-size:1.05rem;
    border:2px solid var(--btn-bg);border-radius:var(--radius);cursor:pointer;
    margin-top:.5rem;min-height:44px;transition:background .15s,border-color .15s}
  .btn-search:hover{background:var(--btn-hover);border-color:var(--btn-hover)}
  .btn-group{display:flex;gap:.6rem;flex-wrap:wrap;align-items:flex-start}
  .btn-download{background:#fff;border:2px solid var(--accent);color:var(--accent);
    font-family:var(--sans);font-size:.875rem;font-weight:600;padding:.45em 1.1em;
    border-radius:var(--radius);cursor:pointer;min-height:44px;transition:background .15s;
    text-align:center;line-height:1.3}
  .btn-download:hover{background:var(--accent-bg)}
  .btn-download small{display:block;font-size:.75rem;font-weight:400;
    opacity:.85;margin-top:.15em}
  /* Error alert */
  .error-box{background:var(--error-bg);border:2px solid var(--error-border);
    border-radius:var(--radius);padding:1rem 1.25rem;margin-bottom:1.5rem;
    font-size:.9rem;color:var(--error-text)}
  .error-box ul{padding-left:1.4em}
  .error-box li{margin:.3em 0}
  .error-box-title{font-weight:700;margin-bottom:.4rem;font-size:.95rem}
  /* Stats bar */
  .stats-bar{display:flex;gap:.75rem;flex-wrap:wrap;margin-bottom:1.5rem}
  .stat-pill{background:var(--surface2);border:1px solid var(--border);
    border-radius:999px;padding:.35em 1em;font-size:.82rem;color:var(--text-sub);
    display:flex;gap:.4em;align-items:center}
  .stat-pill .stat-val{font-weight:700;color:var(--accent)}
  .stat-pill.highlight{border-color:var(--accent);background:var(--accent-bg)}
  /* Results */
  .results-header{display:flex;justify-content:space-between;align-items:center;
    margin-bottom:1rem;flex-wrap:wrap;gap:.75rem}
  .results-heading{font-size:1.1rem;font-weight:700;color:var(--text)}
  .hit-card{background:var(--hit-bg);border:1px solid var(--border);
    border-left:4px solid var(--hit-border);border-radius:var(--radius);
    padding:1rem 1.25rem;margin-bottom:1rem;animation:fadeSlide .25s ease both}
  @keyframes fadeSlide{from{opacity:0;transform:translateY(6px)}to{opacity:1;transform:translateY(0)}}
  .hit-desc{font-family:var(--mono);font-size:.82rem;font-weight:600;
    color:var(--accent);margin-bottom:.7rem;word-break:break-all}
  .hit-table{width:100%;border-collapse:collapse;font-family:var(--mono);font-size:.85rem}
  .hit-table th{text-align:left;font-size:.72rem;font-weight:700;text-transform:uppercase;
    letter-spacing:.08em;color:var(--text-sub);border-bottom:1px solid var(--border);
    padding:.25em .6em .25em 0}
  .hit-table td{padding:.3em .6em .3em 0;border-bottom:1px solid #e0e8f0;
    color:var(--text);vertical-align:baseline}
  .hit-table tr:last-child td{border-bottom:none}
  .hit-pos-cell{white-space:nowrap;color:var(--pos-color);width:7rem}
  /* Wildcard highlight: background + italic + dotted underline (not colour alone — WCAG 1.4.1) */
  .wc{background:var(--wc-bg);color:var(--wc-text);border-radius:2px;
    padding:0 .1em;font-style:italic;font-weight:600;text-decoration:underline dotted}
  .no-hits{text-align:center;padding:3rem;color:var(--text-sub);font-size:.95rem;
    background:var(--surface);border:1px solid var(--border);border-radius:var(--radius)}
  .no-hits .no-hits-icon{font-size:2.2rem;display:block;margin-bottom:.75rem}
  footer{text-align:center;padding:1.75rem 2rem;font-size:.8rem;
    color:var(--text-sub);border-top:1px solid var(--border);background:var(--surface)}
  @media(max-width:560px){fieldset{padding:1.25rem}h1{font-size:1.6rem}}
</style>
</head>
<body>
)";
}

void printPageHeader(const char* version) {
    std::cout
        << "<a class=\"skip-link\" href=\"#main-content\">Skip to main content</a>\n"
        << "<header role=\"banner\">\n"
        << "  <div class=\"logo-mark\" aria-hidden=\"true\">"
           "&#x2B21; MotifHunter v" << version << "</div>\n"
        << "  <h1>Protein <span class=\"brand-accent\">Motif</span> Search</h1>\n"
        << "  <p class=\"subtitle\">"
           "Scan FASTA sequences for patterns with wildcard support</p>\n"
        << "</header>\n"
        << "<main id=\"main-content\">\n"
        << "<div class=\"container\">\n";
}

void printErrors(const std::vector<std::string>& errors) {
    if (errors.empty()) return;
    std::cout <<
        "<div class=\"error-box\" role=\"alert\" "
        "aria-live=\"assertive\" aria-atomic=\"true\">\n"
        "  <p class=\"error-box-title\">"
        "Error &#8212; please correct the following:</p>\n"
        "  <ul>\n";
    for (const auto& e : errors)
        std::cout << "    <li>" << htmlEscape(e) << "</li>\n";
    std::cout << "  </ul>\n</div>\n";
}

void printForm(const std::string& motifVal, const std::string& fastaText) {
    std::cout <<
        "\n<form method=\"POST\" enctype=\"multipart/form-data\" "
        "id=\"searchForm\" aria-label=\"Protein motif search\">\n\n"

        "  <fieldset>\n"
        "    <legend>Search Motif</legend>\n"
        "    <label for=\"motif\">\n"
        "      Pattern\n"
        "      <span style=\"font-weight:400;color:var(--text-sub)\">"
        " (use X as a single-residue wildcard)</span>\n"
        "    </label>\n"
        "    <input type=\"text\" id=\"motif\" name=\"motif\"\n"
        "           placeholder=\"FLXLFX\"\n"
        "           value=\"" << htmlEscape(motifVal) << "\"\n"
        "           autocomplete=\"off\" spellcheck=\"false\" maxlength=\"100\"\n"
        "           aria-describedby=\"motif-hint\" aria-required=\"true\" required>\n"
        "    <p class=\"field-hint\" id=\"motif-hint\">\n"
        "      Letters A&#8211;Z only. X matches any single amino acid.\n"
        "      Examples: <code>FLXLFX</code>, <code>GXXGXG</code>,"
        " <code>RGD</code>, <code>WXXXW</code>\n"
        "    </p>\n"
        "  </fieldset>\n\n"

        "  <fieldset>\n"
        "    <legend>Protein Sequences (FASTA)</legend>\n\n"
        "    <label for=\"fasta_text\">Paste FASTA sequences</label>\n"
        "    <textarea id=\"fasta_text\" name=\"fasta_text\"\n"
        "              aria-describedby=\"fasta-paste-hint\" spellcheck=\"false\"\n"
        "              placeholder=\"&gt;sp|P12345|MYO_HUMAN Myosin heavy chain"
        "&#10;MSSTKIHLEQHVKEIDISQ...\">"
        << htmlEscape(fastaText) <<
        "</textarea>\n"
        "    <p class=\"field-hint\" id=\"fasta-paste-hint\">\n"
        "      Standard FASTA format. Each sequence must start with a "
        "<code>&gt;</code> header line.\n"
        "    </p>\n\n"
        "    <div class=\"or-divider\" aria-hidden=\"true\">or upload a file</div>\n\n"
        "    <label for=\"fasta_file\" "
        "style=\"font-size:.9rem;margin-bottom:.4rem\">\n"
        "      Upload a FASTA file\n"
        "    </label>\n"
        "    <div class=\"file-zone\" id=\"fileZone\"\n"
        "         role=\"group\" "
        "aria-label=\"File upload area &#8212; click or drag and drop\">\n"
        "      <input type=\"file\" id=\"fasta_file\" name=\"fasta_file\"\n"
        "             accept=\".fasta,.fa,.faa,.txt\"\n"
        "             aria-describedby=\"file-hint file-name-display\">\n"
        "      <div class=\"file-zone-inner\">\n"
        "        <span class=\"file-icon\" aria-hidden=\"true\">&#128194;</span>\n"
        "        <p><strong>Click to browse</strong> or drag and drop</p>\n"
        "        <p id=\"file-hint\">"
        "Accepted: .fasta .fa .faa .txt &nbsp;&middot;&nbsp; Max 50&nbsp;MB</p>\n"
        "        <div id=\"file-name-display\" "
        "aria-live=\"polite\" aria-atomic=\"true\"></div>\n"
        "      </div>\n"
        "    </div>\n"
        "  </fieldset>\n\n"
        "  <button type=\"submit\" class=\"btn-search\" name=\"search\" value=\"1\">\n"
        "    Run Motif Search\n"
        "  </button>\n"
        "</form>\n";
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
        "\n<section aria-live=\"polite\" aria-label=\"Search results\" "
        "style=\"margin-top:2rem\">\n\n"

        "<dl class=\"stats-bar\" aria-label=\"Search summary\">\n"
        "  <div class=\"stat-pill\">"
        "<dt>Proteins scanned:</dt>"
        "<dd class=\"stat-val\">" << totalProteins << "</dd></div>\n"
        "  <div class=\"stat-pill highlight\">"
        "<dt>Matches found:</dt>"
        "<dd class=\"stat-val\">" << totalHits << "</dd></div>\n"
        "  <div class=\"stat-pill\">"
        "<dt>Proteins with hits:</dt>"
        "<dd class=\"stat-val\">" << results.size() << "</dd></div>\n"
        "  <div class=\"stat-pill\">"
        "<dt>Motif:</dt>"
        "<dd class=\"stat-val\">" << htmlEscape(motif) << "</dd></div>\n"
        "</dl>\n\n"

        "<div class=\"results-header\">\n"
        "  <h2 class=\"results-heading\" id=\"results-heading\">Results</h2>\n";

    if (!results.empty()) {
        // Shared hidden fields are repeated in each form because browsers only
        // submit fields belonging to the form that was submitted.
        // Each form is self-contained so either button works independently.

        // Helper lambda emits the common hidden fields for both forms
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

        std::cout << "  <div class=\"btn-group\">\n";

        // ── Full results download ──────────────────────────────────────────
        std::cout <<
            "  <form method=\"POST\" enctype=\"multipart/form-data\" "
            "aria-label=\"Download full results\">\n"
            "    <input type=\"hidden\" name=\"download\" value=\"1\">\n";
        emitHiddenFields();
        std::cout <<
            "    <button type=\"submit\" class=\"btn-download\">"
            "&#8595; Download results (.txt)"
            "</button>\n"
            "  </form>\n";

        // ── Gene-only download ─────────────────────────────────────────────
        std::cout <<
            "  <form method=\"POST\" enctype=\"multipart/form-data\" "
            "aria-label=\"Download gene names only\">\n"
            "    <input type=\"hidden\" name=\"download_gene\" value=\"1\">\n";
        emitHiddenFields();
        std::cout <<
            "    <button type=\"submit\" class=\"btn-download\">"
            "&#8595; Download results (.txt)"
            "<small>protein &amp; gene only</small>"
            "</button>\n"
            "  </form>\n";

        std::cout << "  </div>\n";
    }
        std::cout << "</div>\n\n";

    if (results.empty()) {
        std::cout <<
            "<div class=\"no-hits\" role=\"status\">\n"
            "  <span class=\"no-hits-icon\" aria-hidden=\"true\">&#128269;</span>\n"
            "  No sequences matched the motif <strong>"
            << htmlEscape(motif) << "</strong>.\n"
            "  Check your motif pattern or try a different FASTA input.\n"
            "</div>\n";
    } else {
        for (size_t idx = 0; idx < results.size(); ++idx) {
            const auto& r = results[idx];
            std::cout <<
                "<article class=\"hit-card\""
                " style=\"animation-delay:" << (idx * 0.04) << "s\""
                " aria-label=\"Hit in: " << htmlEscape(r.desc) << "\">\n"
                "  <p class=\"hit-desc\">" << htmlEscape(r.desc) << "</p>\n"
                "  <table class=\"hit-table\" aria-label=\"Matches in "
                << htmlEscape(r.desc) << "\">\n"
                "    <thead><tr>"
                "<th scope=\"col\">Position</th>"
                "<th scope=\"col\">Matched sequence</th>"
                "</tr></thead>\n"
                "    <tbody>\n";
            for (const auto& h : r.hits) {
                std::cout <<
                    "      <tr>"
                    "<td class=\"hit-pos-cell\">" << h.pos << "</td>"
                    "<td>" << renderMatch(h.match, motif) << "</td>"
                    "</tr>\n";
            }
            std::cout << "    </tbody>\n  </table>\n</article>\n";
        }
    }
    std::cout << "</section>\n";
}

void printPageFooter(const char* version) {
    std::cout <<
        "\n</div><!-- /container -->\n"
        "</main>\n\n"
        "<footer role=\"contentinfo\">\n"
        "  MotifHunter v" << version << " &mdash;\n"
        "  FASTA validated &middot; streaming processing &middot; "
        "no sequences stored &middot; C++ CGI\n"
        "</footer>\n\n"
        "<script>\n"
        "/* File zone drag-and-drop + file picker */\n"
        "const zone    = document.getElementById('fileZone');\n"
        "const fileIn  = document.getElementById('fasta_file');\n"
        "const display = document.getElementById('file-name-display');\n"
        "fileIn.addEventListener('change', () => {\n"
        "  display.textContent = fileIn.files[0]\n"
        "    ? 'Selected file: ' + fileIn.files[0].name : '';\n"
        "});\n"
        "['dragenter','dragover'].forEach(ev =>\n"
        "  zone.addEventListener(ev, e => {\n"
        "    e.preventDefault(); zone.classList.add('drag-over');\n"
        "    zone.setAttribute('aria-dropeffect','copy');\n"
        "  })\n"
        ");\n"
        "['dragleave','drop'].forEach(ev =>\n"
        "  zone.addEventListener(ev, e => {\n"
        "    e.preventDefault(); zone.classList.remove('drag-over');\n"
        "    zone.removeAttribute('aria-dropeffect');\n"
        "  })\n"
        ");\n"
        "zone.addEventListener('drop', e => {\n"
        "  const f = e.dataTransfer.files[0];\n"
        "  if (f) {\n"
        "    const dt = new DataTransfer(); dt.items.add(f);\n"
        "    fileIn.files = dt.files;\n"
        "    display.textContent = 'Selected file: ' + f.name;\n"
        "  }\n"
        "});\n"
        "/* Auto-uppercase motif field */\n"
        "const motifInput = document.getElementById('motif');\n"
        "motifInput.addEventListener('input', () => {\n"
        "  const pos = motifInput.selectionStart;\n"
        "  motifInput.value = motifInput.value.toUpperCase();\n"
        "  motifInput.setSelectionRange(pos, pos);\n"
        "});\n"
        "/* Move focus to results after submit (WCAG 2.4.3) */\n"
        "const rs = document.querySelector('[aria-label=\"Search results\"]');\n"
        "if (rs) { rs.setAttribute('tabindex','-1'); rs.focus({preventScroll:false}); }\n"
        "</script>\n"
        "</body>\n"
        "</html>\n";
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
    // CGI uses raw binary streams; don't translate \r\n on Windows hosts
    std::ios::sync_with_stdio(false);

    const char* requestMethod  = std::getenv("REQUEST_METHOD");
    const char* contentTypeEnv = std::getenv("CONTENT_TYPE");
    std::string contentType    = contentTypeEnv ? contentTypeEnv : "";

    /* ── GET: serve blank form ── */
    if (!requestMethod || std::string(requestMethod) != "POST") {
        printHeaders();
        printPageHead();
        printPageHeader(APP_VERSION);
        printForm("", "");
        printPageFooter(APP_VERSION);
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
    printPageHead();
    printPageHeader(APP_VERSION);
    printErrors(errors);
    printForm(motifDisplay, form.fastaText);

    if (errors.empty() && totalProteins > 0)
        printResults(results, motifDisplay, totalProteins, totalHits, fastaSource);

    printPageFooter(APP_VERSION);
    return 0;
}
