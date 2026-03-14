#pragma once

#include <vector>
#include <string>

struct Sobolparams
{
    int polynomial;
    int q;
    std::vector<int> minit;
};

std::vector<Sobolparams> extract_Sobolparams(const std::string &filename, int dimensions);

std::vector<std::vector<double>> Sobolpts(int n0, int npts, int d, const std::vector<Sobolparams> &S);