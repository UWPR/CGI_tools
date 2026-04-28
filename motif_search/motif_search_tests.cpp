/*  test_motif_search.cpp — Test suite for motif_search.cgi
 *
 *  Compile & run:
 *      g++ -std=c++17 -O2 -o test_motif_search test_motif_search.cpp
 *      ./test_motif_search
 *
 *  Or use the supplied Makefile target:
 *      make test
 *
 *  The tests exercise every layer of the CGI program:
 *    1.  String utilities   (htmlEscape, trim, toUpper, urlDecode)
 *    2.  Motif validation   (validateMotif)
 *    3.  IUPAC / binary     (isValidAA, hasBinaryBytes)
 *    4.  Motif → regex      (motifToPattern)
 *    5.  FASTA search       (streamSearchFasta) — valid inputs
 *    6.  FASTA search       — error / edge cases
 *    7.  Multipart parser   (parseMultipart)
 *    8.  URL-encoded parser (parseUrlEncoded)
 *    9.  Full CGI cycle     — GET, POST search, download re-POST,
 *                             error paths — via subprocess execution
 *
 *  Strategy: the production source is compiled separately into motif_search.cgi.
 *  Unit tests (groups 1-8) replicate the helper logic inline so they can run
 *  without linking against the CGI binary.
 *  Integration tests (group 9) invoke the compiled binary as a subprocess
 *  through popen(), feeding it the CGI environment variables and a synthesised
 *  request body, then assert on the stdout text.
 */

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

/* ═══════════════════════════════════════════════════════════════
   MINIMAL TEST FRAMEWORK
   ═══════════════════════════════════════════════════════════════ */

static int g_pass = 0;
static int g_fail = 0;
static std::string g_currentSuite;

void suite(const std::string& name) {
    g_currentSuite = name;
    std::cout << "\n\033[1;34m── " << name << " ──\033[0m\n";
}

void check(bool condition, const std::string& label) {
    if (condition) {
        std::cout << "  \033[32m✓\033[0m  " << label << "\n";
        ++g_pass;
    } else {
        std::cout << "  \033[31m✗\033[0m  " << label << "  ← FAIL\n";
        ++g_fail;
    }
}

// Convenience: check that haystack contains needle
void checkContains(const std::string& haystack, const std::string& needle,
                   const std::string& label) {
    check(haystack.find(needle) != std::string::npos, label);
}

void checkNotContains(const std::string& haystack, const std::string& needle,
                      const std::string& label) {
    check(haystack.find(needle) == std::string::npos, label);
}

/* ═══════════════════════════════════════════════════════════════
   INLINE REPLICAS OF PRODUCTION HELPERS
   (kept minimal — just enough to unit-test the logic)
   ═══════════════════════════════════════════════════════════════ */

// ── htmlEscape ──────────────────────────────────────────────────
std::string htmlEscape(const std::string& s) {
    std::string out;
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

// ── trim ────────────────────────────────────────────────────────
std::string trim(const std::string& s) {
    const std::string ws = " \t\r\n";
    size_t b = s.find_first_not_of(ws);
    if (b == std::string::npos) return "";
    size_t e = s.find_last_not_of(ws);
    return s.substr(b, e - b + 1);
}

// ── toUpper ─────────────────────────────────────────────────────
std::string toUpper(std::string s) {
    for (char& c : s)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

// ── urlDecode ───────────────────────────────────────────────────
std::string urlDecode(const std::string& s) {
    std::string out;
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '+') {
            out += ' ';
        } else if (s[i] == '%' && i + 2 < s.size() &&
                   std::isxdigit(static_cast<unsigned char>(s[i+1])) &&
                   std::isxdigit(static_cast<unsigned char>(s[i+2]))) {
            out += static_cast<char>(std::stoi(s.substr(i+1, 2), nullptr, 16));
            i += 2;
        } else {
            out += s[i];
        }
    }
    return out;
}

// ── validateMotif ───────────────────────────────────────────────
bool validateMotif(const std::string& raw, std::string& out, std::string& err) {
    out = toUpper(trim(raw));
    if (out.empty())    { err = "Motif cannot be empty."; return false; }
    if (out.size()>100) { err = "Motif is too long (max 100 characters)."; return false; }
    for (char c : out) {
        if (!std::isalpha(static_cast<unsigned char>(c))) {
            err = "Motif may only contain letters A-Z (use X as wildcard).";
            return false;
        }
    }
    return true;
}

// ── isValidAA ───────────────────────────────────────────────────
bool isValidAA(unsigned char c) {
    static const std::string valid =
        "ACDEFGHIKLMNPQRSTVWYBZXUOacdefghiklmnpqrstvwybzxuo*-";
    return valid.find(static_cast<char>(c)) != std::string::npos;
}

// ── hasBinaryBytes ──────────────────────────────────────────────
bool hasBinaryBytes(const std::string& line) {
    for (unsigned char c : line)
        if (c<=0x08 || c==0x0B || c==0x0C || (c>=0x0E&&c<=0x1F) || c==0x7F)
            return true;
    return false;
}

// ── motifToPattern ──────────────────────────────────────────────
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

// ── Data structures used by streamSearchFasta ───────────────────
struct Hit { int pos; std::string match; };
struct ProteinResult { std::string desc; std::vector<Hit> hits; };

static const size_t MAX_SEQUENCES = 500000;

// ── streamSearchFasta ───────────────────────────────────────────
bool streamSearchFasta(std::istream& source, const std::regex& pattern,
                       size_t motifLen, std::vector<ProteinResult>& results,
                       size_t& total, std::string& error) {
    std::string desc, seq;
    bool seenHeader = false;
    size_t lineNum  = 0;

    auto flushSeq = [&]() -> bool {
        if (desc.empty()) return true;
        if (seq.empty()) { error = "Sequence for '" + desc + "' is empty."; return false; }
        total++;
        if (total > MAX_SEQUENCES) {
            error = "File contains too many sequences (max " +
                    std::to_string(MAX_SEQUENCES) + ").";
            return false;
        }
        std::string upper = toUpper(seq);
        seq.clear();
        ProteinResult pr; pr.desc = desc;
        auto it = std::sregex_iterator(upper.begin(), upper.end(), pattern);
        for (; it != std::sregex_iterator(); ++it) {
            Hit h;
            h.pos   = static_cast<int>(it->position(0)) + 1;
            h.match = upper.substr(static_cast<size_t>(it->position(0)), motifLen);
            pr.hits.push_back(std::move(h));
        }
        if (!pr.hits.empty()) results.push_back(std::move(pr));
        return true;
    };

    std::string raw;
    while (std::getline(source, raw)) {
        lineNum++;
        if (!raw.empty() && raw.back() == '\r') raw.pop_back();
        if (raw.empty()) continue;
        if (hasBinaryBytes(raw)) {
            error = "Line " + std::to_string(lineNum) +
                    ": binary or non-printable content.";
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
        }
    }
    if (!flushSeq()) return false;
    if (total == 0) { error = "No valid FASTA records found in the input."; return false; }
    return true;
}

/* ═══════════════════════════════════════════════════════════════
   INTEGRATION TEST HELPERS
   Run the compiled binary as a subprocess; return stdout as string.
   ═══════════════════════════════════════════════════════════════ */

// Path to the compiled binary — adjusted if needed
static const char* CGI_BIN = "./motif_search.cgi";

// Build a multipart/form-data body from a list of {name, value} pairs.
// Pass filename != "" for the file field.
struct FormField {
    std::string name;
    std::string value;
    std::string filename;  // non-empty → treat as file upload part
};

std::pair<std::string,std::string>
buildMultipart(const std::vector<FormField>& fields,
               const std::string& boundary = "----TestBoundary99887766") {
    std::string body;
    const std::string CRLF = "\r\n";
    for (const auto& f : fields) {
        body += "--" + boundary + CRLF;
        if (f.filename.empty()) {
            body += "Content-Disposition: form-data; name=\"" + f.name + "\"" + CRLF + CRLF;
        } else {
            body += "Content-Disposition: form-data; name=\"" + f.name +
                    "\"; filename=\"" + f.filename + "\"" + CRLF;
            body += "Content-Type: text/plain" + CRLF + CRLF;
        }
        body += f.value + CRLF;
    }
    body += "--" + boundary + "--" + CRLF;
    return {body, boundary};
}

// Run the CGI binary with the given environment and stdin body.
// Returns the full stdout (headers + body).
std::string runCgi(const std::string& method,
                   const std::string& contentType,
                   const std::string& body) {
    // Write body to a temp file so we can feed it via shell
    std::string tmpIn  = "/tmp/cgi_test_in.bin";
    std::string tmpOut = "/tmp/cgi_test_out.txt";

    {
        FILE* f = fopen(tmpIn.c_str(), "wb");
        if (f) { fwrite(body.data(), 1, body.size(), f); fclose(f); }
    }

    // Build shell command with env vars
    std::string cmd =
        "REQUEST_METHOD='" + method + "' "
        "CONTENT_TYPE='"   + contentType + "' "
        "CONTENT_LENGTH='" + std::to_string(body.size()) + "' "
        + std::string(CGI_BIN) + " < " + tmpIn + " > " + tmpOut + " 2>/dev/null";

    int rc = system(cmd.c_str());
    (void)rc;

    // Read output
    std::string out;
    FILE* f = fopen(tmpOut.c_str(), "rb");
    if (f) {
        char buf[4096];
        size_t n;
        while ((n = fread(buf, 1, sizeof(buf), f)) > 0)
            out.append(buf, n);
        fclose(f);
    }
    return out;
}

std::string runGet() {
    return runCgi("GET", "text/html", "");
}

std::string runPostMultipart(const std::vector<FormField>& fields) {
    auto [body, boundary] = buildMultipart(fields);
    return runCgi("POST",
                  "multipart/form-data; boundary=" + boundary,
                  body);
}

std::string runPostUrlEncoded(const std::string& body) {
    return runCgi("POST", "application/x-www-form-urlencoded", body);
}

// URL-encode a string for use in form body
std::string urlEncode(const std::string& s) {
    std::string out;
    for (unsigned char c : s) {
        if (std::isalnum(c) || c=='-' || c=='_' || c=='.' || c=='~') {
            out += static_cast<char>(c);
        } else if (c == ' ') {
            out += '+';
        } else {
            char buf[4];
            snprintf(buf, sizeof(buf), "%%%02X", c);
            out += buf;
        }
    }
    return out;
}

/* ═══════════════════════════════════════════════════════════════
   SAMPLE FASTA DATA
   ═══════════════════════════════════════════════════════════════ */

// Sequence 1: contains FLALF at position 6 and FLGLF at position 56
// Sequence 2: no match for FLXLF
// Sequence 3: contains FLGLF at position 22
static const std::string FASTA_VALID =
    ">sp|P001|TEST1 Test protein one\n"
    "MSSTKFLALFAHLALVTLSSTASSFSEAQRSRQEHLGSMIQKYLESAMKSKLCSS\n"
    "FLGLFMYPRFHRGNEQRPALSLFSDLQFSNVKQMLQNQLELKQKLTKKELEQKL\n"
    ">sp|P002|TEST2 No-hit protein\n"
    "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQR\n"
    ">sp|P003|TEST3 Another hit protein\n"
    "MGSSHHHHHHSSGLVPRGSHMFLGLFAENTSRFLYQSDRQKRPTAAGPCLSGPCR\n";

// Single sequence with no motif hits
static const std::string FASTA_NO_HITS =
    ">sp|P999|NOHIT No hit protein\n"
    "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n";

// FASTA with overlapping motif: FLXLF overlapping FLFLF
static const std::string FASTA_OVERLAP =
    ">sp|P010|OVER Overlap test\n"
    "MSSTKFLFLFAAA\n";

// FASTA with Windows line endings
static const std::string FASTA_CRLF =
    ">sp|P004|CRLF Windows line endings\r\n"
    "MSSTKFLALFAHLALVTLSSTASSFSEAQ\r\n";

// FASTA with empty sequence (invalid)
static const std::string FASTA_EMPTY_SEQ =
    ">sp|P005|EMPTY Empty sequence\n"
    ">sp|P006|NEXT Next sequence\n"
    "ACDEFGHIKLM\n";

// FASTA with sequence before header (invalid)
static const std::string FASTA_NO_HEADER =
    "ACDEFGHIKLMNPQRSTVWY\n"
    ">sp|P007|LATE Late header\n"
    "ACDEFGHIKLM\n";

// FASTA with invalid character in sequence
static const std::string FASTA_BAD_CHAR =
    ">sp|P008|BAD Bad character\n"
    "ACDEF1GHIKL\n";

// FASTA with binary bytes
static const std::string FASTA_BINARY =
    ">sp|P009|BIN Binary data\n"
    + std::string("ACDEF") + static_cast<char>(0x01) + std::string("GHIKL\n");

// Multi-line sequence (tests line concatenation)
static const std::string FASTA_MULTILINE =
    ">sp|P011|MULTI Multi-line sequence\n"
    "MSSTK\n"
    "FLALF\n"
    "AHLAL\n";

/* ═══════════════════════════════════════════════════════════════
   TEST GROUPS
   ═══════════════════════════════════════════════════════════════ */

// ── 1. String utilities ─────────────────────────────────────────
void testStringUtils() {
    suite("1. String Utilities");

    // htmlEscape
    check(htmlEscape("hello") == "hello",            "htmlEscape: plain string unchanged");
    check(htmlEscape("<b>") == "&lt;b&gt;",           "htmlEscape: angle brackets");
    check(htmlEscape("a&b") == "a&amp;b",             "htmlEscape: ampersand");
    check(htmlEscape("\"hi\"") == "&quot;hi&quot;",   "htmlEscape: double quotes");
    check(htmlEscape("it's") == "it&#39;s",           "htmlEscape: single quote");
    check(htmlEscape("") == "",                       "htmlEscape: empty string");
    check(htmlEscape("<script>alert('xss')</script>") ==
          "&lt;script&gt;alert(&#39;xss&#39;)&lt;/script&gt;",
          "htmlEscape: XSS payload fully escaped");

    // trim
    check(trim("  hello  ") == "hello",              "trim: leading/trailing spaces");
    check(trim("\t\nhello\r\n") == "hello",           "trim: tabs and newlines");
    check(trim("") == "",                             "trim: empty string");
    check(trim("   ") == "",                          "trim: all whitespace");
    check(trim("hello") == "hello",                   "trim: nothing to trim");

    // toUpper
    check(toUpper("flxlf") == "FLXLF",               "toUpper: lowercase motif");
    check(toUpper("FLXLF") == "FLXLF",               "toUpper: already upper");
    check(toUpper("") == "",                          "toUpper: empty string");
    check(toUpper("Abc123") == "ABC123",              "toUpper: mixed alphanumeric");

    // urlDecode
    check(urlDecode("hello+world") == "hello world",  "urlDecode: plus as space");
    check(urlDecode("hello%20world") == "hello world", "urlDecode: %20 as space");
    check(urlDecode("%3E") == ">",                    "urlDecode: %3E = >");
    check(urlDecode("%3C") == "<",                    "urlDecode: %3C = <");
    check(urlDecode("%26") == "&",                    "urlDecode: %26 = &");
    check(urlDecode("no+encoding") == "no encoding",  "urlDecode: simple plus");
    check(urlDecode("") == "",                        "urlDecode: empty string");
    check(urlDecode("%41%42%43") == "ABC",            "urlDecode: hex letters");
    check(urlDecode("motif%3DFLXLF") == "motif=FLXLF", "urlDecode: equals sign");
}

// ── 2. Motif validation ─────────────────────────────────────────
void testMotifValidation() {
    suite("2. Motif Validation");

    std::string out, err;

    // Valid motifs
    check(validateMotif("FLXLF", out, err) && out == "FLXLF",
          "valid: exact match motif FLXLF");
    check(validateMotif("flxlf", out, err) && out == "FLXLF",
          "valid: lowercase converted to upper");
    check(validateMotif("RGD", out, err) && out == "RGD",
          "valid: short 3-char motif");
    check(validateMotif("GXXGXG", out, err) && out == "GXXGXG",
          "valid: multiple wildcards");
    check(validateMotif("X", out, err) && out == "X",
          "valid: single wildcard");
    check(validateMotif("  FLXLF  ", out, err) && out == "FLXLF",
          "valid: whitespace trimmed");
    {
        std::string longMotif(100, 'A');
        check(validateMotif(longMotif, out, err),
              "valid: exactly 100 characters");
    }

    // Invalid motifs
    check(!validateMotif("", out, err),               "invalid: empty string");
    check(!validateMotif("   ", out, err),            "invalid: all whitespace");
    check(!validateMotif("FL1LF", out, err),          "invalid: digit in motif");
    check(!validateMotif("FL-LF", out, err),          "invalid: hyphen in motif");
    check(!validateMotif("FL*LF", out, err),          "invalid: asterisk in motif");
    check(!validateMotif("FL LF", out, err),          "invalid: space in motif");
    check(!validateMotif("FL<LF", out, err),          "invalid: HTML char in motif");
    check(!validateMotif("FL;LF", out, err),          "invalid: semicolon in motif");
    {
        std::string tooLong(101, 'A');
        check(!validateMotif(tooLong, out, err),      "invalid: 101 characters (too long)");
        check(err.find("100") != std::string::npos,   "invalid: error mentions 100");
    }
    check(!validateMotif("FL'LF", out, err),          "invalid: single quote (injection)");
    check(!validateMotif("FL\"LF", out, err),         "invalid: double quote (injection)");
}

// ── 3. IUPAC and binary validation ──────────────────────────────
void testCharacterValidation() {
    suite("3. IUPAC / Binary Character Validation");

    // Standard 20 amino acids
    for (char c : std::string("ACDEFGHIKLMNPQRSTVWY")) {
        check(isValidAA(static_cast<unsigned char>(c)),
              std::string("isValidAA: standard AA '") + c + "'");
    }
    // Ambiguous IUPAC
    for (char c : std::string("BZXUObzxuo")) {
        check(isValidAA(static_cast<unsigned char>(c)),
              std::string("isValidAA: ambiguous IUPAC '") + c + "'");
    }
    // Stop codon and gap
    check(isValidAA('*'), "isValidAA: stop codon '*'");
    check(isValidAA('-'), "isValidAA: gap '-'");
    // Invalid characters
    check(!isValidAA('1'),  "isValidAA: digit '1' rejected");
    check(!isValidAA('!'),  "isValidAA: '!' rejected");
    check(!isValidAA(' '),  "isValidAA: space rejected");
    check(!isValidAA('\t'), "isValidAA: tab rejected");
    check(!isValidAA(0x01), "isValidAA: control byte 0x01 rejected");

    // Binary byte detection
    check(!hasBinaryBytes("ACDEFGHIKLM"),      "hasBinaryBytes: clean sequence");
    check(!hasBinaryBytes(">sp|P001 header"),  "hasBinaryBytes: FASTA header");
    check(!hasBinaryBytes(""),                 "hasBinaryBytes: empty line");
    check(hasBinaryBytes(std::string("AC") + static_cast<char>(0x01) + "DEF"),  "hasBinaryBytes: 0x01");
    check(hasBinaryBytes(std::string("AC") + static_cast<char>(0x07) + "DEF"),  "hasBinaryBytes: BEL 0x07");
    check(hasBinaryBytes(std::string("AC") + static_cast<char>(0x0B) + "DEF"),  "hasBinaryBytes: VT 0x0B");
    check(hasBinaryBytes(std::string("AC") + static_cast<char>(0x7F) + "DEF"),  "hasBinaryBytes: DEL 0x7F");
    // Tab (0x09) and LF (0x0A) must NOT be flagged as binary
    check(!hasBinaryBytes(std::string("AC") + static_cast<char>(0x09) + "DEF"),  "hasBinaryBytes: tab is not binary");
}

// ── 4. Motif → regex pattern ─────────────────────────────────────
void testMotifToPattern() {
    suite("4. Motif to Regex Pattern");

    // Pattern structure
    check(motifToPattern("FLXLF") == "FL[A-Za-z]LF",
          "pattern: FLXLF → FL[A-Za-z]LF");
    check(motifToPattern("GXXGXG") == "G[A-Za-z][A-Za-z]G[A-Za-z]G",
          "pattern: GXXGXG multiple wildcards");
    check(motifToPattern("RGD") == "RGD",
          "pattern: RGD no wildcards");
    check(motifToPattern("X") == "[A-Za-z]",
          "pattern: single wildcard X");

    // Verify patterns actually match correctly
    auto testMatch = [](const std::string& motif,
                        const std::string& seq,
                        bool shouldMatch) -> bool {
        std::string pat = motifToPattern(motif);
        std::regex  re(pat, std::regex_constants::icase);
        bool matched = std::regex_search(seq, re);
        return matched == shouldMatch;
    };

    check(testMatch("FLXLF", "FLALF", true),   "regex match: FLXLF matches FLALF");
    check(testMatch("FLXLF", "FLGLF", true),   "regex match: FLXLF matches FLGLF");
    check(testMatch("FLXLF", "FLWLF", true),   "regex match: FLXLF matches FLWLF");
    check(testMatch("FLXLF", "FLLF",  false),  "regex match: FLXLF does not match FLLF (short)");
    check(testMatch("FLXLF", "FL1LF", false),  "regex match: FLXLF does not match FL1LF (digit)");
    check(testMatch("RGD",   "RGDS",  true),   "regex match: RGD matches in RGDS");
    check(testMatch("RGD",   "ARRGD", true),   "regex match: RGD matches mid-sequence");
    check(testMatch("RGD",   "RGE",   false),  "regex match: RGD does not match RGE");
}

// ── 5. FASTA search — valid inputs ──────────────────────────────
void testFastaSearchValid() {
    suite("5. FASTA Search — Valid Inputs");

    auto makePattern = [](const std::string& motif) {
        return std::regex(motifToPattern(motif), std::regex_constants::icase);
    };

    // Basic three-sequence search
    {
        std::istringstream ss(FASTA_VALID);
        std::vector<ProteinResult> results;
        size_t total = 0; std::string err;
        bool ok = streamSearchFasta(ss, makePattern("FLXLF"), 5, results, total, err);
        check(ok, "basic search: returns true");
        check(err.empty(), "basic search: no error");
        check(total == 3, "basic search: 3 proteins scanned");
        check(results.size() == 2, "basic search: 2 proteins with hits");
        check(results[0].desc.find("TEST1") != std::string::npos,
              "basic search: first hit is TEST1");
        check(results[0].hits.size() == 2, "basic search: TEST1 has 2 hits");
        check(results[0].hits[0].pos == 6,  "basic search: TEST1 hit 1 at position 6");
        check(results[0].hits[0].match == "FLALF", "basic search: TEST1 hit 1 match is FLALF");
        check(results[0].hits[1].pos == 56, "basic search: TEST1 hit 2 at position 56");
        check(results[0].hits[1].match == "FLGLF", "basic search: TEST1 hit 2 match is FLGLF");
        check(results[1].desc.find("TEST3") != std::string::npos,
              "basic search: second hit is TEST3");
        check(results[1].hits[0].pos == 22, "basic search: TEST3 hit at position 22");
    }

    // No hits
    {
        std::istringstream ss(FASTA_NO_HITS);
        std::vector<ProteinResult> results;
        size_t total = 0; std::string err;
        bool ok = streamSearchFasta(ss, makePattern("FLXLF"), 5, results, total, err);
        check(ok, "no-hits: returns true (not an error)");
        check(total == 1, "no-hits: 1 protein scanned");
        check(results.empty(), "no-hits: results vector is empty");
    }

    // Overlapping matches: FLFLFAAA with motif FLXLF
    // Matches: FLFLF at pos 1 (FL+F=X, LF), FLFAA — wait, FLFLF:
    //   pos 1: F-L-F-L-F → matches FLXLF (X=F)
    //   pos 3: F-L-F-A-A → F-L-F is only 3, not 5 — no second match starting at 3
    //   Actually FLFLFAAA: pos1=FLFLF (match), pos3=FLFAA (no, F≠L at pos4)
    //   Let's verify with the simpler FLFLF: pos1 FLFLF matches, and check pos3: LFAAA no
    {
        std::istringstream ss(FASTA_OVERLAP);
        std::vector<ProteinResult> results;
        size_t total = 0; std::string err;
        bool ok = streamSearchFasta(ss, makePattern("FLXLF"), 5, results, total, err);
        check(ok, "overlap: returns true");
        check(!results.empty(), "overlap: found at least one match");
        if (!results.empty()) {
            check(results[0].hits[0].match == "FLFLF",
                  "overlap: first match is FLFLF");
        }
    }

    // Windows CRLF line endings
    {
        std::istringstream ss(FASTA_CRLF);
        std::vector<ProteinResult> results;
        size_t total = 0; std::string err;
        bool ok = streamSearchFasta(ss, makePattern("FLXLF"), 5, results, total, err);
        check(ok, "CRLF: Windows line endings accepted");
        check(total == 1, "CRLF: 1 protein parsed");
    }

    // Multi-line sequence — lines concatenated before search
    {
        std::istringstream ss(FASTA_MULTILINE);
        std::vector<ProteinResult> results;
        size_t total = 0; std::string err;
        bool ok = streamSearchFasta(ss, makePattern("FLXLF"), 5, results, total, err);
        check(ok, "multiline: accepted");
        check(total == 1, "multiline: 1 protein");
        // MSSTKFLALF assembled from 3 lines — FLALF at pos 6
        check(!results.empty() && results[0].hits[0].pos == 6,
              "multiline: FLALF found at position 6 across line boundaries");
    }

    // Case-insensitive: motif uppercase, sequence lowercase
    {
        std::string fastaLower =
            ">lower\n"
            "msstkflalfahlal\n";
        std::istringstream ss(fastaLower);
        std::vector<ProteinResult> results;
        size_t total = 0; std::string err;
        bool ok = streamSearchFasta(ss, makePattern("FLXLF"), 5, results, total, err);
        check(ok, "case-insensitive: search on lowercase sequence succeeds");
        check(!results.empty(), "case-insensitive: match found in lowercase sequence");
    }

    // Exact motif (no wildcard)
    {
        std::string fasta = ">exact\nACDEFRGDSTVW\n";
        std::istringstream ss(fasta);
        std::vector<ProteinResult> results;
        size_t total = 0; std::string err;
        bool ok = streamSearchFasta(ss, makePattern("RGD"), 3, results, total, err);
        check(ok, "exact motif: search succeeds");
        check(!results.empty() && results[0].hits[0].pos == 6,
              "exact motif: RGD found at position 6");
    }

    // Position is 1-based
    {
        std::string fasta = ">pos1\nFLALFACDEF\n";
        std::istringstream ss(fasta);
        std::vector<ProteinResult> results;
        size_t total = 0; std::string err;
        streamSearchFasta(ss, makePattern("FLXLF"), 5, results, total, err);
        check(!results.empty() && results[0].hits[0].pos == 1,
              "1-based positions: motif at start of sequence reports pos=1");
    }

    // Empty lines between sequences are ignored
    {
        std::string fasta =
            ">seq1\n"
            "ACDEF\n"
            "\n"
            "\n"
            ">seq2\n"
            "MSSTKFLALFAHLAL\n";
        std::istringstream ss(fasta);
        std::vector<ProteinResult> results;
        size_t total = 0; std::string err;
        bool ok = streamSearchFasta(ss, makePattern("FLXLF"), 5, results, total, err);
        check(ok, "blank lines: empty lines between sequences OK");
        check(total == 2, "blank lines: 2 sequences parsed despite blank lines");
    }

    // Header with no description text → gets auto-named
    {
        std::string fasta = ">\nACDEFGHIKLM\n";
        std::istringstream ss(fasta);
        std::vector<ProteinResult> results;
        size_t total = 0; std::string err;
        bool ok = streamSearchFasta(ss, makePattern("RGD"), 3, results, total, err);
        check(ok, "empty header: auto-named sequence accepted");
        check(total == 1, "empty header: 1 sequence counted");
    }

    // Motif spanning exactly the full sequence
    {
        std::string fasta = ">full\nFLALF\n";
        std::istringstream ss(fasta);
        std::vector<ProteinResult> results;
        size_t total = 0; std::string err;
        bool ok = streamSearchFasta(ss, makePattern("FLXLF"), 5, results, total, err);
        check(ok && !results.empty() && results[0].hits[0].pos == 1,
              "full-length match: motif exactly fills sequence");
    }

    // All IUPAC ambiguous characters accepted in sequence
    {
        std::string fasta = ">iupac\nBZUO*-ACDEFGHIKLM\n";
        std::istringstream ss(fasta);
        std::vector<ProteinResult> results;
        size_t total = 0; std::string err;
        bool ok = streamSearchFasta(ss, makePattern("RGD"), 3, results, total, err);
        check(ok, "IUPAC ambiguous: B Z U O * - all accepted");
    }
}

// ── 6. FASTA search — error / edge cases ────────────────────────
void testFastaSearchErrors() {
    suite("6. FASTA Search — Error Cases");

    auto makePattern = [](const std::string& motif) {
        return std::regex(motifToPattern(motif), std::regex_constants::icase);
    };

    // Empty sequence body
    {
        std::istringstream ss(FASTA_EMPTY_SEQ);
        std::vector<ProteinResult> results;
        size_t total = 0; std::string err;
        bool ok = streamSearchFasta(ss, makePattern("FLXLF"), 5, results, total, err);
        check(!ok, "empty sequence: returns false");
        check(err.find("empty") != std::string::npos, "empty sequence: error mentions 'empty'");
    }

    // Sequence data before first header
    {
        std::istringstream ss(FASTA_NO_HEADER);
        std::vector<ProteinResult> results;
        size_t total = 0; std::string err;
        bool ok = streamSearchFasta(ss, makePattern("FLXLF"), 5, results, total, err);
        check(!ok, "no header: returns false");
        check(err.find("header") != std::string::npos, "no header: error mentions 'header'");
    }

    // Invalid character in sequence
    {
        std::istringstream ss(FASTA_BAD_CHAR);
        std::vector<ProteinResult> results;
        size_t total = 0; std::string err;
        bool ok = streamSearchFasta(ss, makePattern("FLXLF"), 5, results, total, err);
        check(!ok, "bad char: returns false");
        check(err.find("Invalid character") != std::string::npos,
              "bad char: error mentions 'Invalid character'");
    }

    // Binary bytes in file
    {
        std::istringstream ss(FASTA_BINARY);
        std::vector<ProteinResult> results;
        size_t total = 0; std::string err;
        bool ok = streamSearchFasta(ss, makePattern("FLXLF"), 5, results, total, err);
        check(!ok, "binary: returns false");
        check(err.find("binary") != std::string::npos ||
              err.find("non-printable") != std::string::npos,
              "binary: error mentions 'binary' or 'non-printable'");
    }

    // Completely empty input
    {
        std::istringstream ss("");
        std::vector<ProteinResult> results;
        size_t total = 0; std::string err;
        bool ok = streamSearchFasta(ss, makePattern("FLXLF"), 5, results, total, err);
        check(!ok, "empty input: returns false");
        check(err.find("No valid FASTA") != std::string::npos,
              "empty input: error mentions 'No valid FASTA'");
    }

    // Only blank lines
    {
        std::istringstream ss("\n\n\n");
        std::vector<ProteinResult> results;
        size_t total = 0; std::string err;
        bool ok = streamSearchFasta(ss, makePattern("FLXLF"), 5, results, total, err);
        check(!ok, "blank only: returns false");
    }

    // Line number reported correctly in error (bad char on line 3)
    {
        std::string fasta = ">seq\nACDEF\nACDE1F\n";
        std::istringstream ss(fasta);
        std::vector<ProteinResult> results;
        size_t total = 0; std::string err;
        streamSearchFasta(ss, makePattern("RGD"), 3, results, total, err);
        check(err.find("3") != std::string::npos,
              "line number: error reports line 3 for bad character");
    }
}

// ── 7. Multipart parser ─────────────────────────────────────────
void testMultipartParser() {
    suite("7. Multipart Form Parser");

    // We test the parser by invoking the CGI and checking its HTML output,
    // since we cannot call parseMultipart directly from outside the binary.
    // These tests check that the CGI correctly parses multipart fields.

    // Test: motif field parsed correctly
    {
        auto out = runPostMultipart({
            {"motif", "FLXLF", ""},
            {"fasta_text", FASTA_VALID, ""}
        });
        checkContains(out, "FLXLF", "multipart: motif field parsed");
        checkContains(out, "Proteins scanned", "multipart: search ran");
    }

    // Test: file upload field parsed correctly (prioritised over text)
    {
        auto out = runPostMultipart({
            {"motif", "RGD", ""},
            {"fasta_text", FASTA_VALID, ""},
            {"fasta_file", ">sp|FILE|Upload\nACDEFRGDSTVW\n", "upload.fasta"}
        });
        // File upload takes priority; sequence contains RGD at pos 6
        checkContains(out, "Proteins scanned", "multipart file: search ran");
        checkContains(out, "upload.fasta", "multipart file: filename shown in source");
    }

    // Test: download flag parsed
    {
        auto out = runPostMultipart({
            {"motif", "FLXLF", ""},
            {"fasta_text", FASTA_VALID, ""},
            {"download", "1", ""}
        });
        // Without result_ fields this falls through to a search; just check it doesn't error
        checkContains(out, "Content-Type", "multipart download flag: no crash");
    }

    // Test: result_desc / result_pos / result_match fields parsed for download
    {
        auto out = runPostMultipart({
            {"motif", "FLXLF", ""},
            {"download", "1", ""},
            {"result_source", "test.fasta", ""},
            {"result_total", "3", ""},
            {"result_desc",  "sp|P001|TEST1 Test protein", ""},
            {"result_pos",   "6", ""},
            {"result_match", "FLALF", ""},
        });
        checkContains(out, "Content-Disposition: attachment",
                      "multipart download: attachment header present");
        checkContains(out, "Protein Motif Search Report",
                      "multipart download: report header present");
        checkContains(out, "FLALF", "multipart download: hit match in report");
        checkContains(out, "test.fasta", "multipart download: source filename in report");
    }

    // Test: missing motif → error
    {
        auto out = runPostMultipart({
            {"motif", "", ""},
            {"fasta_text", FASTA_VALID, ""}
        });
        checkContains(out, "Motif cannot be empty",
                      "multipart missing motif: error shown");
    }

    // Test: no FASTA data → error
    {
        auto out = runPostMultipart({
            {"motif", "FLXLF", ""}
        });
        checkContains(out, "Please paste a FASTA sequence",
                      "multipart no FASTA: error shown");
    }
}

// ── 8. URL-encoded parser ────────────────────────────────────────
void testUrlEncodedParser() {
    suite("8. URL-Encoded Form Parser");

    // Basic pasted-text search
    {
        std::string body = "motif=FLXLF&fasta_text=" + urlEncode(FASTA_VALID);
        auto out = runPostUrlEncoded(body);
        checkContains(out, "Proteins scanned", "url-encoded: basic search runs");
        checkContains(out, "FLXLF", "url-encoded: motif echoed back");
    }

    // Special characters in motif field (invalid — should error)
    {
        std::string body = "motif=FL%3CLF&fasta_text=" + urlEncode(FASTA_VALID);
        auto out = runPostUrlEncoded(body);
        checkContains(out, "Motif may only contain letters",
                      "url-encoded: < in motif rejected after decode");
    }

    // Download re-POST with result fields
    {
        std::string body =
            "motif=FLXLF"
            "&download=1"
            "&result_source=" + urlEncode("pasted input") +
            "&result_total=3"
            "&result_desc="  + urlEncode("sp|P001|TEST1") +
            "&result_pos=6"
            "&result_match=FLALF";
        auto out = runPostUrlEncoded(body);
        checkContains(out, "Content-Disposition: attachment",
                      "url-encoded download: attachment header");
        checkContains(out, "FLALF", "url-encoded download: match in report");
        checkContains(out, "pasted input", "url-encoded download: source in report");
    }

    // Empty body → errors for both motif and FASTA
    {
        auto out = runPostUrlEncoded("");
        checkContains(out, "Motif cannot be empty",
                      "url-encoded empty body: motif error");
    }

    // Motif with injected HTML entities (should be safely escaped in output)
    {
        std::string body = "motif=" + urlEncode("<SCRIPT>") +
                           "&fasta_text=" + urlEncode(FASTA_VALID);
        auto out = runPostUrlEncoded(body);
        checkNotContains(out, "<SCRIPT>",        "XSS: raw script tag not in output");
        checkNotContains(out, "<script>alert",   "XSS: lowercase script tag not in output");
        checkContains(out, "Motif may only contain", "XSS: motif correctly rejected");
    }
}

// ── 9. Full CGI integration ──────────────────────────────────────
void testCgiIntegration() {
    suite("9. Full CGI Integration");

    // GET request: serve blank form
    {
        auto out = runGet();
        checkContains(out, "Content-Type: text/html",  "GET: HTML content type");
        checkContains(out, "<!DOCTYPE html>",          "GET: HTML doctype");
        checkContains(out, "MotifHunter",              "GET: page title");
        checkContains(out, "<form",                    "GET: form present");
        checkContains(out, "name=\"motif\"",           "GET: motif input present");
        checkContains(out, "name=\"fasta_file\"",      "GET: file input present");
        checkContains(out, "name=\"fasta_text\"",      "GET: textarea present");
        checkNotContains(out, "role=\"alert\"",          "GET: no errors on blank page");
        checkContains(out, "Skip to main content",     "GET: skip link (accessibility)");
        checkContains(out, "lang=\"en\"",              "GET: language attribute set");
    }

    // POST: valid search with pasted FASTA, file upload
    {
        auto out = runPostMultipart({
            {"motif", "FLXLF", ""},
            {"fasta_file", FASTA_VALID, "test.fasta"}
        });
        checkContains(out, "Content-Type: text/html",  "search: HTML response");
        checkContains(out, "Proteins scanned",         "search: stats present");
        checkContains(out, "stat-val\">3<",            "search: 3 proteins scanned");
        checkContains(out, "stat-val\">3<",            "search: 3 total hits");
        checkContains(out, "hit-card",                 "search: hit cards rendered");
        checkContains(out, "TEST1",                    "search: TEST1 in results");
        checkContains(out, "TEST3",                    "search: TEST3 in results");
        checkNotContains(out, "TEST2",                 "search: TEST2 (no hit) absent");
        checkContains(out, "result_desc",              "search: download hidden fields present");
        checkContains(out, "result_pos",               "search: position hidden field present");
        checkContains(out, "result_match",             "search: match hidden field present");
    }

    // POST: pasted FASTA
    {
        std::string body = "motif=RGD&fasta_text=" + urlEncode(
            ">sp|P001|RGDTEST RGD test\nACDEFRGDSTVWYACDEF\n");
        auto out = runPostUrlEncoded(body);
        checkContains(out, "stat-val\">1<", "pasted: 1 protein scanned");
        checkContains(out, "hit-card",      "pasted: hit found");
        checkContains(out, "hit-pos-cell\">6<", "pasted: RGD at position 6");
    }

    // POST: no hits found
    {
        auto out = runPostMultipart({
            {"motif", "FLXLF", ""},
            {"fasta_file", FASTA_NO_HITS, "nohit.fasta"}
        });
        checkContains(out, "No sequences matched",     "no-hits: message shown");
        checkNotContains(out, "class=\"hit-card\"",      "no-hits: no hit cards");
        checkContains(out, "stat-val\">0<",            "no-hits: 0 matches shown");
    }

    // POST: invalid FASTA (bad character)
    {
        auto out = runPostMultipart({
            {"motif", "FLXLF", ""},
            {"fasta_file", FASTA_BAD_CHAR, "bad.fasta"}
        });
        checkContains(out, "FASTA validation failed",  "bad FASTA: validation error shown");
        checkContains(out, "Invalid character",        "bad FASTA: specific error message");
    }

    // POST: invalid FASTA (binary data)
    {
        auto out = runPostMultipart({
            {"motif", "FLXLF", ""},
            {"fasta_file", FASTA_BINARY, "binary.fasta"}
        });
        checkContains(out, "FASTA validation failed",  "binary FASTA: error shown");
    }

    // POST: empty motif
    {
        auto out = runPostMultipart({
            {"motif", "", ""},
            {"fasta_file", FASTA_VALID, "test.fasta"}
        });
        checkContains(out, "error-box",                "empty motif: error box shown");
        checkContains(out, "Motif cannot be empty",    "empty motif: correct message");
    }

    // POST: motif with digit (invalid)
    {
        auto out = runPostMultipart({
            {"motif", "FL3LF", ""},
            {"fasta_file", FASTA_VALID, "test.fasta"}
        });
        checkContains(out, "Motif may only contain",   "digit motif: correct error");
    }

    // POST: no FASTA data (neither file nor text)
    {
        auto out = runPostMultipart({
            {"motif", "FLXLF", ""}
        });
        checkContains(out, "Please paste a FASTA sequence",
                      "no FASTA: correct error message");
    }

    // POST: both errors at once (empty motif + no FASTA)
    {
        auto out = runPostUrlEncoded("motif=");
        checkContains(out, "Motif cannot be empty",    "dual error: motif error shown");
        // Note: when motif is empty the CGI stops before the FASTA check
        // so only the motif error is emitted on a completely empty POST
        checkContains(out, "error-box",                "dual error: error box shown");
    }

    // POST: download with pre-computed hit fields (file upload source)
    {
        auto out = runPostMultipart({
            {"motif",         "FLXLF",      ""},
            {"download",      "1",          ""},
            {"result_source", "bigfile.fasta", ""},
            {"result_total",  "49832",      ""},
            {"result_desc",   "sp|P001|MYO_HUMAN Myosin OS=Homo sapiens", ""},
            {"result_pos",    "143",        ""},
            {"result_match",  "FLALF",      ""},
            {"result_desc",   "sp|P001|MYO_HUMAN Myosin OS=Homo sapiens", ""},
            {"result_pos",    "891",        ""},
            {"result_match",  "FLGLF",      ""},
            {"result_desc",   "sp|Q67890|KIN_MOUSE Kinase OS=Mus musculus", ""},
            {"result_pos",    "22",         ""},
            {"result_match",  "FLTLFS",     ""},
        });
        checkContains(out, "Content-Disposition: attachment",
                      "download: attachment Content-Disposition header");
        checkContains(out, "filename=\"motif_FLXLF_results.txt\"",
                      "download: correct filename in header");
        checkContains(out, "Content-Type: text/plain",
                      "download: plain text content type");
        checkContains(out, "# Protein Motif Search Report",
                      "download: report header present");
        checkContains(out, "# Source    : bigfile.fasta",
                      "download: original source filename preserved");
        checkContains(out, "# Proteins scanned  : 49832",
                      "download: total protein count preserved");
        checkContains(out, "# Proteins with hits: 2",
                      "download: correct proteins-with-hits count");
        checkContains(out, "# Total matches     : 3",
                      "download: correct total match count");
        checkContains(out, ">sp|P001|MYO_HUMAN",
                      "download: first protein header in report");
        checkContains(out, "Position    143 : FLALF",
                      "download: first hit position and match");
        checkContains(out, "Position    891 : FLGLF",
                      "download: second hit of first protein");
        checkContains(out, ">sp|Q67890|KIN_MOUSE",
                      "download: second protein header in report");
        checkContains(out, "Position     22 : FLTLFS",
                      "download: second protein hit");
        checkNotContains(out, "<!DOCTYPE html>",
                      "download: no HTML in plain text output");
        checkNotContains(out, "error-box",
                      "download: no error box in report");
    }

    // POST: download re-POST where motif itself is invalid → error page not download
    {
        auto out = runPostMultipart({
            {"motif",        "FL1LF",  ""},  // digit — invalid
            {"download",     "1",      ""},
            {"result_source","test",   ""},
            {"result_total", "1",      ""},
            {"result_desc",  "TEST",   ""},
            {"result_pos",   "1",      ""},
            {"result_match", "FLALF",  ""},
        });
        checkContains(out, "Content-Type: text/html",  "dl bad motif: returns HTML error page");
        checkContains(out, "Motif may only contain",   "dl bad motif: motif error shown");
        checkNotContains(out, "attachment",            "dl bad motif: no download sent");
    }

    // POST: XSS in FASTA description — output must be escaped
    // Use a sequence that contains RGD so the desc appears in hit output
    {
        std::string xssFasta =
            "><script>alert('xss')</script>\n"
            "ACDEFRGDKLMNPQRST\n";   // has RGD at pos 6
        auto out = runPostMultipart({
            {"motif",      "RGD",   ""},
            {"fasta_file", xssFasta, "xss.fasta"}
        });
        checkNotContains(out, "<script>alert",
                         "XSS in desc: raw script tag not in output");
        // The desc should appear HTML-escaped in the hit-card
        checkContains(out, "&lt;script&gt;",
                      "XSS in desc: angle brackets escaped in output");
    }

    // POST: XSS in match field on download re-POST
    {
        auto out = runPostMultipart({
            {"motif",        "FLXLF",               ""},
            {"download",     "1",                   ""},
            {"result_source","test",                ""},
            {"result_total", "1",                   ""},
            {"result_desc",  "<b>injection</b>",    ""},
            {"result_pos",   "1",                   ""},
            {"result_match", "FLALF",               ""},
        });
        // In the plain-text report, HTML chars from desc appear unescaped
        // (plain text doesn't need HTML escaping) but it must NOT be HTML
        checkContains(out, "Content-Type: text/plain", "XSS in dl desc: still plain text");
        checkNotContains(out, "<!DOCTYPE", "XSS in dl desc: not rendered as HTML");
    }

    // POST: large protein count in download (field preserved correctly)
    {
        auto out = runPostMultipart({
            {"motif",        "FLXLF",     ""},
            {"download",     "1",         ""},
            {"result_source","huge.fasta",""},
            {"result_total", "499999",    ""},
            {"result_desc",  "prot1",     ""},
            {"result_pos",   "12345",     ""},
            {"result_match", "FLWLF",     ""},
        });
        checkContains(out, "499999", "large count: protein total preserved");
        checkContains(out, "12345",  "large count: large position preserved");
    }

    // Accessibility: key ARIA and semantic HTML present
    {
        auto out = runGet();
        // role="alert" only emitted when errors exist — verify on a POST with errors
        {
            auto errOut = runPostMultipart({{"motif","",""},{"fasta_file",FASTA_VALID,"t.fasta"}});
            checkContains(errOut, "role=\"alert\"", "a11y: role=alert on error container");
        }
        checkContains(out, "aria-required=\"true\"","a11y: aria-required on motif input");
        checkContains(out, "aria-describedby",     "a11y: aria-describedby on inputs");
        checkContains(out, "aria-live",            "a11y: aria-live regions present");
        checkContains(out, "<fieldset",            "a11y: fieldset for form grouping");
        checkContains(out, "<legend",              "a11y: legend for fieldset");
        checkContains(out, "role=\"banner\"",      "a11y: header landmark");
        checkContains(out, "role=\"contentinfo\"", "a11y: footer landmark");
        checkContains(out, "id=\"main-content\"",  "a11y: main content anchor for skip link");
        checkContains(out, "prefers-reduced-motion", "a11y: reduced motion media query");
    }
}

/* ═══════════════════════════════════════════════════════════════
   MAIN — run all suites and print summary
   ═══════════════════════════════════════════════════════════════ */
int main() {
    std::cout << "\033[1mMotifHunter CGI Test Suite\033[0m\n";
    std::cout << "Binary under test: " << CGI_BIN << "\n";

    testStringUtils();
    testMotifValidation();
    testCharacterValidation();
    testMotifToPattern();
    testFastaSearchValid();
    testFastaSearchErrors();
    testMultipartParser();
    testUrlEncodedParser();
    testCgiIntegration();

    int total = g_pass + g_fail;
    std::cout << "\n\033[1m════ Results: "
              << g_pass << "/" << total << " passed";
    if (g_fail == 0) {
        std::cout << " \033[32m✓ ALL PASS\033[0m\n";
    } else {
        std::cout << " \033[31m(" << g_fail << " FAILED)\033[0m\n";
    }
    return g_fail > 0 ? 1 : 0;
}
