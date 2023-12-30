#include "myutils.hpp"
#include <iomanip>
#include <iostream>

double MY_EPS = 1E-1;
double MY_INF = 1E10;
string PDF_READER = "xdg-open";// the Linux will choose the default one

#define MAX_CMD_SIZE 1000// maximum number of characters of a command used by the system routine

//     * Utility functions
//     * Routines to read a table from text file

//====================================================================================
//     * Dealing with PDF files and strings

int view_pdf_file(const string &filename) {
    char cmd[MAX_CMD_SIZE];
    sprintf(cmd, "%s %s", PDF_READER.c_str(), filename.c_str());
    system(cmd);
    return 0;
}

// ========================================================================
//     * Reading tables from text files

// read a table containing a header (first line that is not a comment) and then,
// read nr rows, according to the header.
StringTable::StringTable(int nr, ifstream &file) {
    string word, h;
    this->nrows = nr;
    this->line.reserve(nr);
    line.resize(nr);

    this->ncols = 0;
    while (getline(file, h) && is_comment(h, "#"))
        ;
    replace(h.begin(), h.end(), '\t', ' ');
    istringstream token(h);

    // Read the header:
    while (getline(token, word, ' ')) {
        trim(word);
        if (word.empty()) continue;
        lowercase(word);
        this->header.push_back(word);
        this->column_size.push_back(word.length());
        this->ncols++;
    }
    // Read the rows:
    for (int i = 0; i < this->nrows; i++) {
        while (getline(file, h) && is_comment(h, "#"))
            ;
        replace(h.begin(), h.end(), '\t', ' ');
        istringstream token(h);

        // Read a single line:
        int c = 0;
        this->line[i].reserve(this->ncols);
        this->line[i].resize(this->ncols);
        while (getline(token, word, ' ')) {
            trim(word);
            if (word.empty()) continue;
            if (c == this->ncols) {
                stringstream buffer;
                buffer << "Number of cols, in line " << i + 1 << " is larger than defined in the header." << endl;
                throw runtime_error(buffer.str());
            }
            this->line[i][c] = word;
            if (word.length() > this->column_size[c]) this->column_size[c] = word.length();
            c++;
        }
        if (c < this->ncols) {
            stringstream buffer;
            buffer << "Number of cols, in line " << i + 1 << " is smaller than defined in the header (" << c << " < "
                   << this->ncols << ")" << endl;
            throw runtime_error(buffer.str());
        }
    }
}

StringTable::StringTable(int nrows, int ncols) {// empty startup
    string word, h;
    this->nrows = nrows;
    this->ncols = ncols;
    this->header.reserve(ncols);
    this->header.resize(ncols);
    this->line.reserve(nrows);
    line.resize(nrows);
    for (int i = 0; i < nrows; i++) {
        this->line[i].reserve(ncols);
        this->line[i].resize(ncols);
    }
}

string StringTable::first(const string &column_name) {
    int col = this->column_index(column_name);
    if (col != -1) return line[0][col];
    else
        throw invalid_argument("Invalid column name.");
}

int StringTable::first_int(const string &column_name) {
    int col = this->column_index(column_name);
    if (col != -1) return stoi(line[0][col]);
    else
        throw invalid_argument("Invalid column name.");
}

double StringTable::first_double(const string &column_name) {
    int col = this->column_index(column_name);
    if (col != -1) return stod(line[0][col]);
    else
        throw invalid_argument("Invalid column name.");
}

bool StringTable::read_column(const string &column_name, vector<string> &col) {
    int col_i = this->column_index(column_name);
    if (col_i == -1) return false;
    col.reserve(this->nrows);
    col.resize(this->nrows);
    for (int i = 0; i < this->nrows; i++) col[i] = this->line[i][col_i];
    return true;
}

bool StringTable::read_column(const string &column_name, vector<int> &col) {
    int col_i = this->column_index(column_name);
    if (col_i == -1) return false;
    col.reserve(this->nrows);
    col.resize(this->nrows);
    for (int i = 0; i < this->nrows; i++) col[i] = stoi(this->line[i][col_i]);
    return true;
}

bool StringTable::read_column(const string &column_name, vector<double> &col) {
    int col_i = this->column_index(column_name);
    if (col_i == -1) return false;
    col.reserve(this->nrows);
    col.resize(this->nrows);
    for (int i = 0; i < this->nrows; i++) col[i] = stod(this->line[i][col_i]);
    return true;
}

void StringTable::print() {// prints the whole table
    for (int j = 0; j < this->ncols; j++) {
        for (int k = 0; k < this->column_size[j]; k++) cout << "-";
        cout << "-+ ";
    }
    cout << endl;

    for (int j = 0; j < this->ncols; j++) cout << setw(this->column_size[j]) << header[j] << " | ";
    cout << endl;
    for (int j = 0; j < this->ncols; j++) {
        for (int k = 0; k < this->column_size[j]; k++) cout << "-";
        cout << " + ";
    }
    cout << endl;

    for (int i = 0; i < this->nrows; i++) {
        for (int j = 0; j < this->ncols; j++) cout << setw(this->column_size[j]) << this->line[i][j] << " | ";
        cout << endl;
    }

    for (int j = 0; j < this->ncols; j++) {
        for (int k = 0; k < this->column_size[j]; k++) cout << "-";
        cout << " + ";
    }
    cout << endl;
}
