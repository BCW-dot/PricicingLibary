#pragma once

#include "stdafx.h"
#include "StockPriceModel.h"

std::vector<double> generatePoisson_path(double start, double end, double param);
std::vector<double> generatePoisson_Martingale(double start, double end, int n, double param);
std::vector<double> generate_Merton_paths(double start, double end, int n, double s0, double param, double mu_j, double sigma_j, double r, double sigma, std::vector<double>& s);

void testPoissonProcess();
