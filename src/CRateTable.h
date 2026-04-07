#pragma once

#include <string>

void InitCRateTable(const std::string& filename,
                    double mD,
                    double CR,
                    double g,
                    double T);

double CRateInterp(double q2, double phi);
