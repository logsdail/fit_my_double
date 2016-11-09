/*
 *  @Outputs.h
 *  fit_my_ecp
 *
 *  @brief Simple parent class to deal with testing of
 *  ECPs within the structure of the properties calculating environment
 *
 *  Created by Andrew Logsdail on 01/06/2012.
 *  Copyright 2012 University College London. All rights reserved.
 *
 */

#ifndef OUTPUTS_H
#define OUTPUTS_H

#include <iostream>
#include <cmath>
// Personal headers
#include "Structures.h"

class Outputs{
	
public:
	
	/*
	 Constructor
	 
	 No params
	 */
	Outputs(){}
	
	/*
	 Constructor
	 
	 @param[in/out] seed Initial seed for random number
	 */
	Outputs(int *seed)
	{
		idum = seed;
		converged = false;
	}
	
	/*
	 Deconstructor
	 
	 No params
	 */
	virtual ~Outputs(){}
	
	/*
	 Check if the calculation has converged
	 
	 @return bool Operator with true or false value
	 */
	virtual bool get_converged()
	{
		return converged;
	}
	
	/*
	 Prints details of method being used
	 
	 No params
	 */
	virtual void print_type()
	{
		std::cout << "Performing an analysis of outputs only" << std::endl;
	}
	
	/*
	 Set the initial parameters for the Powells optimisation
	 
	 @param[in] ss The step size to be initially used for A and zeta
	 @param[in] sr Step reduction rate for A and zeta
	 @param[in] ssm Convergence criteria, defined as the target step size to reach, for A and zeta
	 @param[in] min Minimum values for A and zeta
	 @param[in] max Maximum values for A and zeta
	 @param[in] msod Maximum steps in any one direction - by default this is disabled
	 */
	virtual void set_parameters(std::vector< std::vector<double> > ss,
				    std::vector< std::vector<double> > sr,
  				    std::vector< std::vector<double> > ssm,
				    std::vector<double> min,
				    std::vector<double> max,
				    int msod)
	{;}
	
	/*
	 Method to set the initial GA parameters.
	 Requirement if we are running a quick GA search.
	 
	 @param[in] ps The size of the GA population
	 @param[in] ms The size of the mutant population
	 @param[in] os The size of the offspring population
	 @param[in] cc Convergence criteria to terminate a GA search
	 @param[in] md Check for dynamic mutation
	 */
	virtual void set_ga_parameters(int ps,
								   int ms,
								   int os,
								   int cc,
								   bool md)
	{;}
	
	/*
	 Defines initial ECPs, and sets up class
	 
	 @param[in] v Vector of initial gaussian(s) read in from ecp.template
	 */
	virtual void set_starting_gaussians(gaussian v)
	{;}
	
	virtual std::vector<gaussian> get_ecps_to_test();
	
	/*
	 Takes in the details of the last search run, and saves data
	 
	 @param[in] v Vector of ECPs just tested
	 @param[in] n Index of the highest ranked ECP
	 */
	virtual void set_ecps_tested(std::vector<gaussian> v, int n)
	{;}
	
protected:

	// Random number pointer
	int *idum;
	// Boolean to check conerged state
	bool converged;
	// Standard vectors
	gaussian starting_gaussians;
	std::vector<gaussian> ecps_to_test;
	
	/*
	 Compare two ecps to see if they are the same
	 
	 @param[in] ecp1 First gaussian ECP
	 @param[in] ecp2 Second gaussian ECP
         @param[in] DELTA gap for comparison of values
	 @return bool True or false operator
	 */
	bool compare_ecps(gaussian ecp1, gaussian ecp2, double DELTA = 0.000001)
	{
		return compare_ecps(ecp1.values,ecp2.values,DELTA);
	}
	
	// Return a generic default gaussian
	gaussian get_default();
	
	// Compare two double arrays to see if they are the same
	bool compare_ecps(std::vector<gaussian_info> ecp1, std::vector<gaussian_info> ecp2, double DELTA = 0.000001);
	
};

#endif

