/*
 *  @Linear.h
 *  fit_my_ecp
 *
 *  @brief Implementation of Linear Scan Method.
 *  Inherits some attributes from Outputs, which are rewritten here
 *
 *  Created by Andrew Logsdail on 06/11/2013.
 *  Copyright 2012 University College London. All rights reserved.
 *
 */

#ifndef LINEAR_H
#define LINEAR_H

#include "Utils.h"
#include "Outputs.h"

class Linear : public Outputs {
	
public:
	
	/*
	 Constructor
	 
	 No params
	 */
	Linear(){}
	
	Linear(int *seed);
	
	/*
	 Deconstructor
	 
	 No params
	 */
	~Linear(){}
	
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
		//ecps_to_test.push_back(v);
		//previous_minimum = v; 
		calculate_ecps_to_test();
	}
	
	/*
	 Takes in the details of the last search run, and saves data
	 
	 @param[in] v Vector of ECPs just tested
	 @param[in] n Index of the highest ranked ECP
	 */
	void set_ecps_tested(std::vector<gaussian> v, int n)
	{
		ecps_tested = v; 
		//number_one_ranked = n;
		check_converged();
	}
	
	
	/*
	 Method to print the search type being conducted, common with all other search methods
	 
	 No Param
	 */
	virtual void print_type()
	{
		std::cout << "Performing a Linear Scan" << std::endl;
	}
	
protected:
	
	//std::vector<int>::size_type vector_counter;
	//int max_steps_in_one_direction;
	int minimisation_count;
	//int number_one_ranked;
	// Vectors
	//std::vector< std::vector<int> > search_vectors; // 2D vector. Search vectors (x) by factor (y)
	//std::vector<int> vector_counts;
	//std::vector<int> powell_search_vectors;
	std::vector<gaussian> ecps_tested;
	// Floats
	std::vector< std::vector<double> > step_size;
	//std::vector< std::vector<double> > step_reduction;
	//std::vector< std::vector<double> > step_size_min;
	std::vector<double> minimums;
	std::vector<double> maximums;
	// Others
	//gaussian previous_minimum;
	
 	virtual void calculate_ecps_to_test();
	
	/*
	Set this to true after one loop - instant exit
 
	No param
	*/
	virtual void check_converged()
	{
        	converged = true;
	}
	
	void check_boundaries_duplicates();
	
private:
	
	void resize_search_vectors(gaussian g);
};

#endif

