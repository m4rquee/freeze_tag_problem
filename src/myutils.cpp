#include "myutils.hpp"

double MY_EPS = 1E-1;
double MY_INF = 1E10;
string PDF_READER = "xdg-open";// the Linux will choose the default one

#define MAX_CMD_SIZE 1000// maximum number of characters of a command used by the system routine

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

// Read a table containing a header (first line that is not a comment) and then,
// read n rows, according to the header:
StringTable::StringTable(int n, ifstream &file) {
    string word, line;
    this->nrows = n;
    this->lines.reserve(n);
    line.resize(n);

    this->ncols = 0;

    // Read the header:
    while (getline(file, line) && is_comment(line, "#")) {}
    replace(line.begin(), line.end(), '\t', ' ');
    istringstream tokenizer(line);
    while (getline(tokenizer, word, ' ')) {
        trim(word);
        if (word.empty()) continue;
        lowercase(word);
        this->header.push_back(word);
        this->column_size.push_back(word.length());
        this->ncols++;
    }

    // Read the rows:
    for (int i = 0; i < this->nrows; i++) {
        while (getline(file, line) && is_comment(line, "#")) {}
        replace(line.begin(), line.end(), '\t', ' ');
        tokenizer = istringstream(line);

        // Read a single line:
        int c = 0;
        this->lines[i].reserve(this->ncols);
        this->lines[i].resize(this->ncols);
        while (getline(tokenizer, word, ' ')) {
            trim(word);
            if (word.empty()) continue;
            if (c == this->ncols) {
                stringstream buffer;
                buffer << "Number of columns in line " << i + 1 << " is larger than defined in the header." << endl;
                throw runtime_error(buffer.str());
            }
            this->lines[i][c] = word;
            this->column_size[c] = max(this->column_size[c], (int) word.length());
            c++;
        }
        if (c < this->ncols) {
            stringstream buffer;
            buffer << "Number of columns in line " << i + 1 << " is smaller than defined in the header." << endl;
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
    this->lines.reserve(nrows);
    lines.resize(nrows);
    for (int i = 0; i < nrows; i++) {
        this->lines[i].reserve(ncols);
        this->lines[i].resize(ncols);
    }
}

string StringTable::first(const string &column_name) {
    int col = this->column_index(column_name);
    if (col != -1) return lines[0][col];
    else
        throw invalid_argument("Invalid column name.");
}

int StringTable::first_int(const string &column_name) {
    int col = this->column_index(column_name);
    if (col != -1) return stoi(lines[0][col]);
    else
        throw invalid_argument("Invalid column name.");
}

double StringTable::first_double(const string &column_name) {
    int col = this->column_index(column_name);
    if (col != -1) return stod(lines[0][col]);
    else
        throw invalid_argument("Invalid column name.");
}

bool StringTable::read_column(const string &column_name, vector<string> &col) {
    int col_i = this->column_index(column_name);
    if (col_i == -1) return false;
    col.reserve(this->nrows);
    col.resize(this->nrows);
    for (int i = 0; i < this->nrows; i++) col[i] = this->lines[i][col_i];
    return true;
}

bool StringTable::read_column(const string &column_name, vector<int> &col) {
    int col_i = this->column_index(column_name);
    if (col_i == -1) return false;
    col.reserve(this->nrows);
    col.resize(this->nrows);
    for (int i = 0; i < this->nrows; i++) col[i] = stoi(this->lines[i][col_i]);
    return true;
}

bool StringTable::read_column(const string &column_name, vector<double> &col) {
    int col_i = this->column_index(column_name);
    if (col_i == -1) return false;
    col.reserve(this->nrows);
    col.resize(this->nrows);
    for (int i = 0; i < this->nrows; i++) col[i] = stod(this->lines[i][col_i]);
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
        for (int j = 0; j < this->ncols; j++) cout << setw(this->column_size[j]) << this->lines[i][j] << " | ";
        cout << endl;
    }

    for (int j = 0; j < this->ncols; j++) {
        for (int k = 0; k < this->column_size[j]; k++) cout << "-";
        cout << " + ";
    }
    cout << endl;
}
