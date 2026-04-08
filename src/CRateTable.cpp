#pragma once

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Interpolate.h"

namespace {
    constexpr size_t Nq = 180;
    constexpr size_t Nphi = 180;

    struct Row {
        double q;
        double phi;
        double C;
        double err;
    };

    Interpolate gCRateInterpolate(Nq, Nphi);
    bool gCRateReady = false;

    double gQ2Min = 0.0, gQ2Max = 0.0;
    double gPhiMin = 0.0, gPhiMax = 0.0;

    double WrapPhi(double phi) {
        const double twoPi = 2.0 * M_PI;
        while (phi < 0.0) phi += twoPi;
        while (phi >= twoPi) phi -= twoPi;
        return phi;
    }

    double ParseSpecialDouble(const std::string& s) {
        std::string t = s;
        for (char& c : t) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));

        if (t == "nan" || t == "+nan" || t == "-nan") {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (t == "inf" || t == "+inf" || t == "infinity" || t == "+infinity") {
            return std::numeric_limits<double>::infinity();
        }
        if (t == "-inf" || t == "-infinity") {
            return -std::numeric_limits<double>::infinity();
        }
        return std::stod(s);
    }
}

void InitCRateTable(const std::string& filename,
                    double mD,
                    double CR,
                    double g,
                    double T)
{
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::vector<Row> data;
    data.reserve(Nq * Nphi);

    std::string line;
    size_t lineNo = 0;

    while (std::getline(fin, line)) {
        ++lineNo;
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string sq, sphi, sC, serr;

        if (!(iss >> sq >> sphi >> sC >> serr)) {
            std::cerr << "Bad line format at line " << lineNo << ": " << line << std::endl;
            throw std::runtime_error("Failed to parse data line.");
        }

        Row row;
        row.q   = std::stod(sq);
        row.phi = std::stod(sphi);
        row.C   = ParseSpecialDouble(sC);
        row.err = ParseSpecialDouble(serr);

        data.push_back(row);
    }

    // std::cerr << "Read rows = " << data.size() << std::endl;
    // std::cerr << "Nq = " << Nq << ", Nphi = " << Nphi
    //          << ", expected = " << Nq * Nphi << std::endl;

    if (!data.empty()) {
        std::cerr << "First row: "
                  << data.front().q << " "
                  << data.front().phi << " "
                  << data.front().C << " "
                  << data.front().err << std::endl;

        std::cerr << "Last row: "
                  << data.back().q << " "
                  << data.back().phi << " "
                  << data.back().C << " "
                  << data.back().err << std::endl;
    }

    if (data.size() != Nq * Nphi) {
        throw std::runtime_error("InitCRateTable: file size does not match 180x180.");
    }

    // repair NaN/Inf ：using near point of phi values
    for (size_t i = 0; i < Nq; ++i) {
        for (size_t j = 0; j < Nphi; ++j) {
            size_t idx = i * Nphi + j;

            if (!std::isfinite(data[idx].C) || !std::isfinite(data[idx].err)) {
                size_t jl = (j == 0 ? Nphi - 1 : j - 1);
                size_t jr = (j + 1 == Nphi ? 0 : j + 1);

                size_t idxL = i * Nphi + jl;
                size_t idxR = i * Nphi + jr;

                bool leftOK  = std::isfinite(data[idxL].C) && std::isfinite(data[idxL].err);
                bool rightOK = std::isfinite(data[idxR].C) && std::isfinite(data[idxR].err);

                if (leftOK && rightOK) {
                    std::cerr << "Repair bad point at (i=" << i << ", j=" << j
                              << ") using neighbors." << std::endl;
                    data[idx].C   = 0.5 * (data[idxL].C   + data[idxR].C);
                    data[idx].err = 0.5 * (data[idxL].err + data[idxR].err);
                } else if (leftOK) {
                    std::cerr << "Repair bad point at (i=" << i << ", j=" << j
                              << ") using left neighbor only." << std::endl;
                    data[idx].C   = data[idxL].C;
                    data[idx].err = data[idxL].err;
                } else if (rightOK) {
                    std::cerr << "Repair bad point at (i=" << i << ", j=" << j
                              << ") using right neighbor only." << std::endl;
                    data[idx].C   = data[idxR].C;
                    data[idx].err = data[idxR].err;
                } else {
                    std::cerr << "Bad point cannot be repaired at line index " << idx << std::endl;
                    throw std::runtime_error("NaN/Inf point cannot be repaired from neighbors.");
                }
            }
        }
    }

    // q^2 
    for (size_t i = 0; i < Nq; ++i) {
        double q2 = data[i * Nphi].q * data[i * Nphi].q;
        gCRateInterpolate.setX(i, q2);

        if (i == 0) gQ2Min = q2;
        if (i == Nq - 1) gQ2Max = q2;
    }

    // phi 
    for (size_t j = 0; j < Nphi; ++j) {
        double phi = data[j].phi;
        gCRateInterpolate.setY(j, phi);

        if (j == 0) gPhiMin = phi;
        if (j == Nphi - 1) gPhiMax = phi;
    }

    // fill out 
    for (size_t i = 0; i < Nq; ++i) {
        for (size_t j = 0; j < Nphi; ++j) {
            size_t idx = i * Nphi + j;

            double C_tilde = (mD * mD) / (CR * g * g * T) * data[idx].C;
            gCRateInterpolate.setValues(i, j, C_tilde);
        }
    }

    gCRateInterpolate.init();
    gCRateReady = true;
}

double CRateInterp(double q2, double phi)
{
    if (!gCRateReady) {
        throw std::runtime_error("CRateInterp called before InitCRateTable.");
    }

    phi = WrapPhi(phi);
    phi = std::max(gPhiMin, std::min(phi, gPhiMax));
    // q2  = std::max(gQ2Min,  std::min(q2,  gQ2Max));
    if (q2 < gQ2Min) {
        return gCRateInterpolate(gQ2Min, phi);
    }
    else if (q2 > gQ2Max) {
        return gCRateInterpolate(gQ2Max, phi) * std::pow(gQ2Max / q2, 4);
    }

    return gCRateInterpolate(q2, phi);
}
