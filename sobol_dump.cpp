#include "sobol.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

int main(int argc, char *argv[])
{
    if (argc < 5)
    {
        std::cerr << "Usage: ./sobol_dump <N> <d> <direction_file> <output_file>\n";
        return 1;
    }

    int N = std::stoi(argv[1]);
    int d = std::stoi(argv[2]);
    std::string dir_f = argv[3];
    std::string out_f = argv[4];

    std::vector<Sobolparams> params = extract_Sobolparams(dir_f, d);
    std::vector<std::vector<double>> pts = Sobolpts(1, N, d, params);

    std::ofstream out(out_f);
    if (!out)
    {
        std::cerr << "Could not open output file: " << out_f << "\n";
        return 1;
    }

    out << std::fixed << std::setprecision(10);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < d; ++j)
        {
            out << pts[i][j];
            if (j < d - 1)
                out << '\t';
        }
        out << '\n';
    }
    return 0;
}