#ifndef MY_COLOR_DEFINE
#define MY_COLOR_DEFINE

#include "myutils.hpp"
#include <iomanip>
#include <map>

typedef struct {
    string color_name, rgb;
    int rgb_red, rgb_green, rgb_blue;
    int next_visually_distinct_color;
} color_element;

class ColorClass {
public:
    ColorClass();// use colors defined in the routine init
    string name(int color_code);
    int insert(const string &color_name, const string &rgb);
    string rgb(int color_code);
    bool exists(const string &color_name);
    int size() const;
    void init();
    void print();

    std::map<string, int> code;// dictionary that returns the code of a given color name
    vector<color_element> colors;
    vector<int> ith_visually_distinct_color;
    int first_visually_distinct_color{}, last_visually_distinct_color{}, n_visually_distinct_colors{};
};

extern ColorClass Colors;// like a static instance of ColorClass

inline int ColorCode(string color_name) { return Colors.code[color_name]; }
inline int ColorSize() { return Colors.colors.size(); }
inline string ColorRGB(string color_name) { return Colors.rgb(Colors.code[color_name]); }
inline string ColorRGB(int color_code) { return Colors.rgb(color_code); }
inline bool ColorExists(string color_name) { return Colors.exists(color_name); }
inline string color_name(int code) { return Colors.name(code); }
inline int ColorInsert(string color_name, string rgb) { return Colors.insert(color_name, rgb); }
inline void ColorPrint() { Colors.print(); }

int InsertVisualDistinctColor(const string &color_name);
int NextVisualDistinctColor(int color_code);
int FirstVisualDistinctColor();
int SizeVisualDistinctColor();
int ith_VisualDistinctColor(int i);
inline string ith_VisualDistinctColorName(int i) { return Colors.name(ith_VisualDistinctColor(i)); }

#endif// MY_COLOR_DEFINE
