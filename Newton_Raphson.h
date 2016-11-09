/*
 *  @Newton_Raphson.h
 *  fit_my_ecp
 *
 *  @brief Implementation of Powell's minimisation method.
 *  Inherits some attributes from Outputs, which are rewritten here
 *
 *  Created by Andrew Logsdail on 16/08/2012.
 *  Copyright 2012 University College London. All rights reserved.
 *
 */

#ifndef NEWTON_RAPHSON_H
#define NEWTON_RAPHSON_H

#include "Utils.h"
#include "Outputs.h"
#include "lbfgs.h"

class Newton_Raphson : public Outputs {
	
public:
	
	/*
	 Constructor
	 
	 No params
	 */
	Newton_Raphson(){}

	Newton_Raphson(int *seed, bool flag = false);
	
	/*
	 Deconstructor
	 
	 No params
	 */
	~Newton_Raphson(){}
	
	virtual void set_parameters(std::vector< std::vector<double> > ss,
				    std::vector< std::vector<double> > sr,
				    std::vector< std::vector<double> > ssm,
  				    std::vector<double> min,
				    std::vector<double> max,
				    int msod);	
	
	/*
	 Defines initial ECPs, and sets up class
	 
	 @param[in] v Vector of initial gaussian(s) read in from ecp.template
	 */
	void set_starting_gaussians(gaussian v)
	{
		starting_gaussians = v; 
		ecps_to_test.push_back(v);
		previous_minimum = v; 
		calculate_ecps_to_test();
	}
	
	/*
	 Takes in the details of the last search run, and saves data
	 

	 @param[in] n Index of the highest ranked ECP
	 */
	void set_ecps_tested(std::vector<gaussian> v, int n)
	{
		ecps_tested = v; 
		number_one_ranked = n;
		check_converged();
	}
	
	
	/*
	 Method to print the search type being conducted, common with all other search methods
	 
	 No Param
	 */
	virtual void print_type()
	{
		if (lbfgs_flag)
		{
			std::cout << "Performing a LBFGS minimisation" << std::endl;
		}
		else
		{
			std::cout << "Performing a Newton-Raphson minimisation" << std::endl;
		}
	}
	
private:
	// Int	
	int number_one_ranked;
	// Vectors
	std::vector< std::vector<int> > search_vectors; // 2D vector. Search vectors (x) by factor (y)
	std::vector<gaussian> ecps_tested;
	// Floats
	std::vector< std::vector<double> > step_size;
        std::vector< std::vector<double> > step_size_min;
        std::vector<double> maximums;
        std::vector<double> minimums;
        // Boolean
        bool lbfgs_flag;
	// Others
        gaussian previous_minimum;
	
	virtual void calculate_ecps_to_test();
	
	virtual void check_converged();
	
	void check_duplicates();
	
	void resize_search_vectors(const unsigned int size);

        void PrintMinimiserOptions();

	scitbx::lbfgs::minimizer<double>* m_pMinimizer;
};

#endif

