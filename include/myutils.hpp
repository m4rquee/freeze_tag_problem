#ifndef MY_UTILS_DEFINE
#define MY_UTILS_DEFINE

#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <time.h>
#include <vector>

using namespace std;

extern double MY_EPS;    // used to verify if numbers are the same, within some error of MY_EPS
extern double MY_INF;    // all numbers above this threshold are treated as infinity
extern string PDF_READER;// Program to open a pdf file.
                         // You can change it for another pdf viewer,
                         // but the program must run in a terminal with the syntax:
                         // <pdfreader>  <pdffile>

//====================================================================================
//     * Dealing with PDF files and strings

int view_pdf_file(const string &filename);// opens a pdf file using the program defined at PDF_READER

inline bool is_space(char c) { return (c == ' ') || (c == '\t'); }

inline bool is_comment(string line, string prefix) {// if the given prefix appears after all blank characters of a line
    int i = 0;
    while (i < line.length() && is_space(line[i])) i++;// skip over blank characters

    if (i == line.length()) return true;                  // consider empty lines as comments
    if (i + prefix.length() > line.length()) return false;// cannot contain the prefix

    return !line.compare(i, prefix.length(), prefix);// has the prefix after all blank characters
}

inline bool is_suffix(string str, string suf) {
    return str.size() >= suf.size() && !str.compare(str.size() - suf.size(), suf.size(), suf);
}

inline bool is_prefix(string str, string pre) { return str.size() >= pre.size() && !str.compare(0, pre.size(), pre); }

inline int hex2int(string s) {
    int n = s.length(), d = 0;
    for (int i = 0; i < n; i++) {
        char c = tolower(s[i]);
        if (c >= '0' && c <= '9') d = d * 16 + (c - '0');
        else if (c >= 'a' && c <= 'f')
            d = d * 16 + (c - 'a' + 10);
        else
            throw invalid_argument("Character that is not a hexadecimal digit found.");
    }
    return d;
}

inline bool file_exists(const string &filename) {
    ifstream f(filename.c_str());
    bool ret = f.is_open();
    if (ret) f.close();
    return ret;
}

inline char char_tolower(char c) { return tolower(c); }

inline void lowercase(string &s) { transform(s.begin(), s.end(), s.begin(), char_tolower); }

static inline string &ltrim(string &str) {// left trim
    str.erase(str.begin(), find_if(str.begin(), str.end(), not1(ptr_fun<int, int>(isspace))));
    return str;
}

static inline string &rtrim(string &str) {// right trim
    str.erase(find_if(str.rbegin(), str.rend(), not1(ptr_fun<int, int>(isspace))).base(), str.end());
    return str;
}

static inline string &trim(string &str) {// left and right trim
    return ltrim(rtrim(str));
}

//====================================================================================
//     * Functions to test values

inline bool is_frac(double x) {// if x is fractional (within a certain small error)
    double f;
    f = ceil(x - MY_EPS) - x;
    if (f < MY_EPS || f > 1.0 - MY_EPS) return false;
    return true;// MY_EPS <= ceil(x)-x <= 1-MY_EPS
}

inline bool is_equal(double x, double y) {// if x and y are equal (within a certain small error)
    if (x > y + MY_EPS || x < y - MY_EPS) return false;
    return true;// y-eps <= x <= y+eps
}

// The next functions suppose that 0 <= x <= 1:
inline bool bin_is_one(double x) { return x > (double) (1.0 - MY_EPS); }
inline bool bin_is_zero(double x) { return x < MY_EPS; }

// ========================================================================
//     * Reading tables from text files

class StringTable {
public:
    StringTable(int nrows, ifstream &file);
    StringTable(int nrows, int ncols);
    int nrows, ncols;
    vector<string> header;
    vector<int> column_size;
    vector<vector<string>> lines;
    int column_index(const string &column_name);   // return the column index, or -1 if it does not exist
    string first(const string &column_name);       // return the first element of the column
    int first_int(const string &column_name);      // return the first element of the column as int
    double first_double(const string &column_name);// return the first element of the column as double
    bool read_column(const string &column_name, vector<string> &col);
    bool read_column(const string &column_name, vector<int> &col);
    bool read_column(const string &column_name, vector<double> &col);
    bool entry(int row, int col, string &entry);
    bool entry(int row, int col, double &entry);
    bool entry(int row, int col, int &entry);
    void print();
};

inline int StringTable::column_index(const string &column_name) {// return the column index, or -1 if it does not exist
    string cname = column_name;
    lowercase(cname);
    for (int col = 0; col < this->ncols; col++)
        if (cname == this->header[col]) return col;
    return -1;
}

inline bool StringTable::entry(int row, int col, string &entry) {
    if (row < 0 || row >= this->nrows) return false;
    if (col < 0 || col >= this->ncols) return false;
    entry = this->lines[row][col];
    return true;
}

inline bool StringTable::entry(int row, int col, double &entry) {
    if (row < 0 || row >= this->nrows) return false;
    if (col < 0 || col >= this->ncols) return false;
    entry = stod(this->lines[row][col]);
    return true;
}

inline bool StringTable::entry(int row, int col, int &entry) {
    if (row < 0 || row >= this->nrows) return false;
    if (col < 0 || col >= this->ncols) return false;
    entry = stoi(this->lines[row][col]);
    return true;
}

inline string get_path(string filename) {// the path (ending with "/") of a given filename
    int n = filename.size();
    if (n == 0) throw invalid_argument("Invalid filename.");
    n--;
    while (n >= 0 && filename[n] != '/') n--;
    if (n < 0) return "./";
    return filename.substr(0, n + 1);
}

#endif// MY_UTILS_DEFINE
