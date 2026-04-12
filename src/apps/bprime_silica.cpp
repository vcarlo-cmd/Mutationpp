/**
 * @file bprime_silica.cpp
 *
 * @brief Computes the B' table for silica (SiO2) ablation in air.
 *
 * Implements the surface mass balance for a pure SiO2 char material.
 * The tracked element is Si, and the char elemental composition is:
 *   y_Si = M_Si / M_SiO2 ≈ 0.4674
 *   y_O  = 2*M_O / M_SiO2 ≈ 0.5326
 *
 * Usage:
 *   bprime_silica -T T1:dT:T2 -P pressure_pa -b Bg -m mixture -bl BL_comp [-py pyro_comp]
 *
 * Example:
 *   bprime_silica -T 300:25:5000 -P 101325 -b 0 -m silica-air -bl air
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 */

#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

#include "mutation++.h"

#ifdef _GNU_SOURCE
#include <fenv.h>
#endif

using namespace std;
using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;

// ---------------------------------------------------------------------------
// Command-line option handling (same pattern as bprime.cpp)
// ---------------------------------------------------------------------------

typedef struct Options {
    double T1;
    double T2;
    double dT;
    double P1;
    double Bg;
    std::string mixture;
    std::string boundary_layer_comp;
    std::string pyrolysis_composition;
    bool pyrolysis_exist = false;
} Options;

bool optionExists(int argc, char** argv, const std::string& option) {
    return (std::find(argv, argv + argc, option) != argv + argc);
}

std::string getOption(int argc, char** argv, const std::string& option) {
    std::string value;
    char** ptr = std::find(argv, argv + argc, option);
    if (ptr == argv + argc || ptr + 1 == argv + argc)
        value = "";
    else
        value = *(ptr + 1);
    return value;
}

void printHelpMessage(const char* const name) {
    std::string tab("    ");
    cout.setf(std::ios::left, std::ios::adjustfield);
    cout << endl;
    cout << "Usage: " << name << " [OPTIONS]" << endl;
    cout << "Compute the B' table for SiO2 (silica) ablation using Mutation++." << endl;
    cout << endl;
    cout << tab << "-h, --help          prints this help message" << endl;
    cout << tab << "-T                  temperature range \"T1:dT:T2\" or T (default = 300:100:5000 K)" << endl;
    cout << tab << "-P                  pressure in Pa (default = 1 atm)" << endl;
    cout << tab << "-b                  pyrolysis non-dimensional mass blowing rate (default = 0)" << endl;
    cout << tab << "-m                  mixture name (default = silica-air)" << endl;
    cout << tab << "-bl                 boundary layer edge composition name" << endl;
    cout << tab << "-py                 pyrolysis composition name (optional)" << endl;
    cout << endl;
    cout << "Example:" << endl;
    cout << tab << name << " -T 300:25:5000 -P 101325 -b 0 -m silica-air -bl air" << endl;
    cout << endl;
    exit(0);
}

bool parseRange(const std::string& range, double& x1, double& x2, double& dx) {
    std::vector<std::string> tokens;
    String::tokenize(range, tokens, ":");
    if (!String::isNumeric(tokens)) return false;
    switch (tokens.size()) {
        case 1:
            x1 = atof(tokens[0].c_str());
            x2 = x1;
            dx = 1.0;
            break;
        case 3:
            x1 = atof(tokens[0].c_str());
            x2 = atof(tokens[2].c_str());
            dx = atof(tokens[1].c_str());
            break;
        default:
            return false;
    }
    if (dx == 0.0) { x2 = x1; dx = 1.0; }
    return true;
}

bool parseRange(const std::string& range, double& x1) {
    std::vector<std::string> tokens;
    String::tokenize(range, tokens, ":");
    if (!String::isNumeric(tokens)) return false;
    x1 = atof(tokens[0].c_str());
    return true;
}

Options parseOptions(int argc, char** argv) {
    Options opts;
    if (argc < 2) printHelpMessage(argv[0]);
    if (optionExists(argc, argv, "-h") || optionExists(argc, argv, "--help"))
        printHelpMessage(argv[0]);

    if (optionExists(argc, argv, "-T")) {
        if (!parseRange(getOption(argc, argv, "-T"), opts.T1, opts.T2, opts.dT)) {
            cout << "Bad format for temperature range!" << endl;
            printHelpMessage(argv[0]);
        }
    } else {
        opts.T1 = 300.0; opts.T2 = 5000.0; opts.dT = 100.0;
    }

    if (optionExists(argc, argv, "-P")) {
        if (!parseRange(getOption(argc, argv, "-P"), opts.P1)) {
            cout << "Bad format for pressure!" << endl;
            printHelpMessage(argv[0]);
        }
    } else {
        opts.P1 = ONEATM;
    }

    if (optionExists(argc, argv, "-b")) {
        if (!parseRange(getOption(argc, argv, "-b"), opts.Bg)) {
            cout << "Bad format for Bg!" << endl;
            printHelpMessage(argv[0]);
        }
    } else {
        opts.Bg = 0.0;
    }

    opts.mixture = optionExists(argc, argv, "-m")
        ? getOption(argc, argv, "-m") : "silica-air";

    if (optionExists(argc, argv, "-bl")) {
        opts.boundary_layer_comp = getOption(argc, argv, "-bl");
    } else {
        printHelpMessage(argv[0]);
    }

    if (optionExists(argc, argv, "-py")) {
        opts.pyrolysis_composition = getOption(argc, argv, "-py");
        opts.pyrolysis_exist = true;
    }

    return opts;
}

// ---------------------------------------------------------------------------
// Silica surface mass balance
// ---------------------------------------------------------------------------
//
// Elemental mass balance at the wall (per element k):
//   y_w,k * (Bc + Bg + 1) = Bc * y_c,k + Bg * y_g,k + y_e,k
//
// For SiO2 char :  y_c,Si = M_Si/M_SiO2,  y_c,O = 2*M_O/M_SiO2
// Tracking element: Si
//   Bc = (y_e,Si + Bg*y_g,Si - yw,Si*(1+Bg)) / (yw,Si - y_c,Si)
//
void surfaceMassBalanceSilica(
    Mixture& mix,
    const double* const Yke,   // elemental mass fractions, boundary layer edge
    const double* const Ykg,   // elemental mass fractions, pyrolysis gas
    const double T,
    const double P,
    const double Bg,
    double& Bc,
    double& hw,
    double* const Xw)          // species mole fractions at wall (output)
{
    const int ne = mix.nElements();
    const int ng = mix.nGas();
    const int ns = mix.nSpecies();

    // SiO2 elemental mass fractions (char composition)
    const double M_Si   = 28.086;
    const double M_O    = 15.999;
    const double M_SiO2 = M_Si + 2.0 * M_O;   // 60.084 g/mol

    const int iSi = mix.elementIndex("Si");
    const int iO  = mix.elementIndex("O");

    std::vector<double> Ychar(ne, 0.0);
    Ychar[iSi] = M_Si / M_SiO2;
    Ychar[iO]  = 2.0 * M_O / M_SiO2;

    // 1. Initialise wall elemental composition = BL edge + pyrolysis contribution
    std::vector<double> p_Xw(ne);
    double sum = 0.0;
    for (int i = 0; i < ne; ++i) {
        p_Xw[i] = Yke[i] + Bg * Ykg[i];
        sum += p_Xw[i];
    }

    // 2. Add large amount of SiO2 char to simulate infinite ablative surface
    double char_amount = std::max(100.0 * Bg, 200.0);
    for (int i = 0; i < ne; ++i)
        p_Xw[i] += char_amount * Ychar[i];
    sum += char_amount;

    for (int i = 0; i < ne; ++i)
        p_Xw[i] /= sum;

    // 3. Convert elemental mass fractions -> mole fractions, then equilibrate
    mix.convert<YE_TO_XE>(p_Xw.data(), p_Xw.data());
    mix.equilibriumComposition(T, P, p_Xw.data(), Xw, IN_PHASE);

    // 4. Gas-phase mean molecular weight and Si mass fraction at wall
    double mwg = 0.0, ywSi = 0.0;
    for (int j = 0; j < ng; ++j) {
        mwg  += mix.speciesMw(j) * Xw[j];
        ywSi += mix.elementMatrix()(j, iSi) * Xw[j];
    }
    ywSi *= mix.atomicMass(iSi) / mwg;

    // 5. B'c (silica-specific denominator: yw,Si - y_c,Si instead of yw,C - 1)
    double yck = Ychar[iSi];   // M_Si / M_SiO2 ≈ 0.4674
    Bc = (Yke[iSi] + Bg * Ykg[iSi] - ywSi * (1.0 + Bg)) / (ywSi - yck);
    Bc = std::max(Bc, 0.0);

    // 6. Wall enthalpy (J/kg)
    std::vector<double> h(ns);
    mix.speciesHOverRT(T, h.data());
    hw = 0.0;
    for (int i = 0; i < ng; ++i)
        hw += Xw[i] * h[i];
    hw *= RU * T / mwg;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[])
{
#ifdef _GNU_SOURCE
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

    Options opts = parseOptions(argc, argv);

    Mixture mix(opts.mixture);

    const int ne = mix.nElements();
    const int ns = mix.nSpecies();

    std::vector<double> Yke(ne, 0.0);
    std::vector<double> Ykg(ne, 0.0);
    std::vector<double> Xw(ns, 0.0);

    mix.getComposition(opts.boundary_layer_comp, Yke.data(), Composition::MASS);
    if (opts.pyrolysis_exist)
        mix.getComposition(opts.pyrolysis_composition, Ykg.data(), Composition::MASS);

    // Print header
    cout << setw(10) << "\"Tw[K]\"" << setw(15) << "\"B'c\"" << setw(15) << "\"hw[MJ/kg]\"";
    for (int i = 0; i < ns; ++i)
        cout << setw(25) << "\"" + mix.speciesName(i) + "\"";
    cout << endl;

    // Main loop over temperature
    double Bc, hw;
    for (double T = opts.T1; T < opts.T2 + 1.0e-6; T += opts.dT) {
        surfaceMassBalanceSilica(mix, Yke.data(), Ykg.data(),
                                 T, opts.P1, opts.Bg, Bc, hw, Xw.data());
        cout << setw(10) << T << setw(15) << Bc << setw(15) << hw / 1.0e6;
        for (int i = 0; i < ns; ++i)
            cout << setw(25) << Xw[i];
        cout << endl;
    }

    return 0;
}
