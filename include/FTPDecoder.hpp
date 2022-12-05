#ifndef FTP_DECODER_DEFINE
#define FTP_DECODER_DEFINE

#include "ftp_utils.hpp"
#include <algorithm>
#include <list>
#include <vector>

class FTPDecoder {
    FTP_Instance &P;

public:
    explicit FTPDecoder(FTP_Instance &P);
    ~FTPDecoder();

    double decode(const std::vector<double> &chromosome) const;
    void decode(const std::vector<double> &chromosome, DNodeVector &Sol) const;
};

#endif// FTP_DECODER_DEFINE
