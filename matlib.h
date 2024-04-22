#pragma once

#include "stdafx.h"
#include "RealFunction.h"


/*  Create a linearly spaced vector */
std::vector<double> linspace( double from, double to, int numPoints );
/*  Compute the sum of a vector */
double sum( const std::vector<double>& v );
/*  Compute the mean of a vector */
double mean( const std::vector<double>& v );
/*  Compute the standard deviation of a vector */
double standardDeviation( const std::vector<double>& v, bool population=0 );
/*  Find the minimum of a vector */
double min( const std::vector<double>& v );
/*  Find the maximum of a vector */
double max( const std::vector<double>& v );
/*  Find the given percentile of a vector */
double prctile( const std::vector<double>& v, double percentage );
/*  Sort a vector */
std::vector<double> sort( const std::vector<double>&  v );

/*  Create uniformly distributed random numbers */
std::vector<double> randuniform( int n );
/*  Create normally distributed random numbers */
std::vector<double> randn( int n );
std::vector<double> my_randn( int n );
std::vector<double> my_randuniform( int n );
void my_rng();
/* open a plot/hist right away */
void open_plot(std::string file_name);
void open_hist(std::string file_name);
/* print a window of a vector */
void print_window(std::vector<double> v, int start, int end);
/*  Seeds the default random number generator */
void rng( const std::string& setting );

void print(std::vector<double> v);

/**
 *  Computes the cumulative
 *  distribution function of the
 *  normal distribution
 */
double normcdf( double x );

/* Computes the inverse of normcdf */
double norminv( double x ); 


/*  Create a line chart */
void plot( const std::string& fileName,
           const std::vector<double>& x,
           const std::vector<double>& y);

/*  Plot a histogram */
void hist( const std::string& fileName,
           const std::vector<double>& values,
           int numBuckets=10);

/*  Integrate using the rectangle rule */
double integral( RealFunction& f,
                 double a,
                 double b,
                 int nSteps );

double sum_parallel(std::vector<double>& v);
double sum_simple(std::vector<double> v);
void writeVectorsToCsv(const std::vector<double>& v1, const std::vector<double>& t, const std::vector<double>& v2, std::string& filename);
void fill_correlated_normals(std::vector<double>& n1, std::vector<double>& n2, double rho);





/**
 *  Test function
 */
void testMatlib();

