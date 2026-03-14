#include "sobol.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;
#define MSB 32

int rightmostbit(int n)
{
    int pos = 0;
    while (n & 1)
    {
        n >>= 1;
        pos++;
    }
    return pos;
}

vector<Sobolparams> extract_Sobolparams(const string &filename, int dimensions)
{
    ifstream file(filename);
    string line;
    getline(file, line);
    vector<Sobolparams> S;
    while (getline(file, line))
    {
        stringstream ss(line);
        int d, s, a;
        if (!(ss >> d >> s >> a))
            continue;
        if (d - 1 <= dimensions)
        {
            vector<int> mi_vector;
            int value;
            while (ss >> value)
            {
                mi_vector.push_back(value);
            }
            S.push_back({a, s, mi_vector});
        }
        else
        {
            return S;
        }
    }
    return {};
}

vector<vector<bool>> Sobolmat(bitset<MSB> c, int q, int r, vector<int> m)
{
    vector<bool> a(q - 1, 0);
    vector<vector<bool>> V(r, vector<bool>(r, 0));

    m.resize(r, 0);

    if (q)
    {
        for (int i = 0; i < q - 1; ++i)
        {
            a[i] = c[q - 2 - i];
        }
        for (int j = q - 1; j < r; ++j)
        {
            for (int i = 0; i < q - 1; ++i)
            {
                m[j] ^= (1 << (i + 1)) * a[i] * m[j - i - 1];
            }
            m[j] ^= ((1 << (q)) * m[j - q]) ^ m[j - q];
        }

        for (int j = 0; j < r; ++j)
        {
            int val = m[j];
            for (int i = 0; i <= j; ++i)
            {
                V[i][j] = ((val >> (j - i)) & 1);
            }
        }
    }
    else
    {
        for (int i = 0; i < r; ++i)
        {
            for (int j = 0; j < r; ++j)
            {
                if (i == j)
                {
                    V[i][j] = 1;
                }
            }
        }
    }

    return V;
}
vector<vector<double>> Sobolpts(int n0, int npts, int d, const vector<Sobolparams> &S)
{
    int nmax = npts + n0 - 1;
    int rmax = 1 + floor(log2(nmax));
    int r = (n0 > 1 ? 1 + floor(log2(n0 - 1)) : 1);
    vector<vector<double>> P(npts, vector<double>(d, 0));
    vector<vector<bool>> y(rmax, vector<bool>(d, 0));
    vector<vector<vector<bool>>> V(d, vector<vector<bool>>(rmax, vector<bool>(rmax, 0)));
    int qnext = (1 << r);

    bitset<MSB> g((n0 - 1) ^ ((n0 - 1) >> 1));

    for (int i = 0; i < d; ++i)
    {
        Sobolparams s = S[i];
        bitset<MSB> c(s.polynomial);
        V[i] = Sobolmat(c, s.q, rmax, s.minit);
    }
    for (int i = 0; i < d; ++i)
    {

        for (int m = 0; m < rmax; ++m)
        {
            for (int n = 0; n < rmax; ++n)
            {
                y[m][i] = y[m][i] ^ (V[i][m][n] & g[rmax - n - 1]);
            }
        }
    }
    int l;
    for (int k = n0; k < nmax; ++k)
    {
        if (k == qnext)
        {
            r += 1;
            l = 1;
            qnext *= 2;
        }
        else
        {
            l = r - rightmostbit(k - 1);
        }
        for (int i = 0; i < d; ++i)
        {
            for (int m = 0; m < r; ++m)
            {
                y[m][i] = y[m][i] ^ V[i][m][r - l];
            }
            for (int j = 0; j < r; ++j)
            {
                P[k - n0][i] += y[j][i] * pow(2.0, -(j + 1));
            }
        }
    }
    return P;
}
