/*
 *  @file Powells.cpp
 *  fit_my_ecp
 *
 *  Created by Andrew Logsdail on 01/06/2012.
 *  Copyright 2012 University College London. All rights reserved.
 *
 * EDITS:
 * 05/11/2013 : Fixed the resizing of search vectors for varying step sizes
 *
 */

#include "Powells.h"

using namespace std;

/*
 Constructor
 
 @param[in/out] seed Pointer to the seed
 */
Powells::Powells(int *seed) : Outputs(seed)
{	
	//idum = seed;
	vector_counter = 0;
	max_steps_in_one_direction = 0;
	minimisation_count = 0;
	number_one_ranked = 0;
	
	// Search area confines
	step_size.resize(2);
	step_reduction.resize(2);
	step_size_min.resize(2);
	minimums.resize(2,0.0);
	maximums.resize(2,0.0);
}

/*
 Set up vectors for creating Powell's search vector
 
 @param[in] Number of different directions possible
 */
void Powells::resize_search_vectors(gaussian g)
{
	const unsigned int size = g.values.size();

	// Intiate search vectors
	search_vectors.resize(size+1);
	
	// Set up search vectors
	for (unsigned int i = 0; i < size; i++)
	{
		search_vectors[i].resize(size,0);
		
		// Set one pathway as operative in initial vector
		search_vectors[i][i] = 1;
	}
	
	// Set up Powells vector
	search_vectors[size].resize(size,0);
	
	// Set up counter for search vectors
	vector_counts.resize(size+1,0);

        // Addition: We can resize the step size vectors here, if need be.
        for (int i = 0; i < 2; i++)
	{
		if (step_size[i].size() != size)
		{
			step_size[i].resize(size,step_size[i][0]);
		}

		if (step_reduction[i].size() != size)
                {
                        step_reduction[i].resize(size,step_reduction[i][0]);
                }

                if (step_size_min[i].size() != size)
                {
                        step_size_min[i].resize(size,step_size_min[i][0]);
                }
	}

        // Temporary values
        vector<double> new_step_size;
        vector<double> new_step_reduction;
	vector<double> new_step_size_min;
        int counter = 0;

        //cout << "RESIZING" << endl;

	//This is so cumbersome, but for now will do OK.
        // We can resize the step size vectors here, if need be.
        for (int i = 0; i < 2; i++)
        {
                new_step_size.clear();
                counter = 0;

                // Sizes aren't the same
                if (step_size[i].size() != size)
                {
                        // Check what type of Gaussian component we have
                        for (vector<gaussian_info>::size_type j = 0; j < size; j++)
                        {
                                // If not the same add padding
                                if (i != g.values[j].type)
                                {
                                        // step_size[j].insert(step_size[j].begin(),0.0);
                                        new_step_size.push_back(0.0);
                                }
                                else
                                {
                                        new_step_size.push_back(step_size[i][counter]);
                                        counter++;
                                }
                        }
                }

                // Did adding the padding fix everything?
                if (new_step_size.size() != size)
                {
                        // If not we just resize with a fixed value for everything
                        // Not pretty but saves trouble further on
                        step_size[i].resize(1);
                        step_size[i].resize(size,step_size[i][0]);
                }
                else
                {
                        step_size[i] = new_step_size;
                }

		counter = 0;
		new_step_reduction.clear();

		 // Sizes aren't the same
                if (step_reduction[i].size() != size)
                {
                        // Check what type of Gaussian component we have
                        for (vector<gaussian_info>::size_type j = 0; j < size; j++)
                        {
                                // If not the same add padding
                                if (i != g.values[j].type)
                                {
                                        new_step_reduction.push_back(0.0);
                                }
                                else
                                {
                                        new_step_reduction.push_back(step_reduction[i][counter]);
                                        counter++;
                                }
                        }
                }

                // Did adding the padding fix everything?
                if (new_step_reduction.size() != size)
                {
                        // If not we just resize with a fixed value for everything
                        // Not pretty but saves trouble further on
                        step_reduction[i].resize(1);
                        step_reduction[i].resize(size,step_reduction[i][0]);
                }
                else
                {
                        step_size[i] = new_step_reduction;
                }

                counter = 0;
                new_step_size_min.clear();

                 // Sizes aren't the same
                if (step_size_min[i].size() != size)
                {
                        // Check what type of Gaussian component we have
                        for (vector<gaussian_info>::size_type j = 0; j < size; j++)
                        {
                                // If not the same add padding
                                if (i != g.values[j].type)
                                {
                                        new_step_size_min.push_back(0.0);
                                }
                                else
                                {
                                        new_step_size_min.push_back(step_reduction[i][counter]);
                                        counter++;
                                }
                        }
                }

                // Did adding the padding fix everything?
                if (new_step_size_min.size() != size)
                {
                        // If not we just resize with a fixed value for everything
                        // Not pretty but saves trouble further on
                        step_size_min[i].resize(1);
                        step_size_min[i].resize(size,step_size_min[i][0]);
                }
                else
                {
                        step_size[i] = new_step_size_min;
                }
        }
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
void Powells::set_parameters(vector< vector<double> > ss,
			     vector< vector<double> > sr,
			     vector< vector<double> > ssm,
			     vector<double> min,
			     vector<double> max,
			     int msod)
{
	step_size = ss;
	step_reduction = sr;
	step_size_min = ssm;
	minimums = min;
	maximums = max;
	max_steps_in_one_direction = msod;
	
	// - step size
	if (step_size[0].size() == 0)
	{
		cout << "Using default starting step size for A: 1" << endl;
		step_size[0].push_back(1.0);
	}
	if (step_size[1].size() == 0)
	{
		cout << "Using default starting step size for Z: 1" << endl;
		step_size[1].push_back(1.0);
	}
	
	// - step size reduction
	if (step_reduction[0].size() == 0)
	{
		cout << "Using default step size reduction for A: 0.1 (i.e. step size is reduced by an order when we are at minima)" << endl;
		step_reduction[0].push_back(0.1);
	}
	if (step_reduction[1].size() == 0)
	{
		cout << "Using default step size reduction for Z: 0.1 (i.e. step size is reduced by an order when we are at minima)" << endl;
		step_reduction[1].push_back(0.1);
	}
	
	// - step size minimum
	if (step_size_min[0].size() == 0)
	{
		cout << "Using default step size minimum for convergence for A: 1.0" << endl;
		step_size_min[0].push_back(1.0);
	}
	if (step_size_min[1].size() == 0)
	{
		cout << "Using default step size minimum for convergence for Z: 1.0" << endl;
		step_size_min[1].push_back(1.0);
	}
	
	// - constraining values for A
	if (minimums[0] == 0.0)
	{
		cout << "Using default constraint of minimum value for A: 0.0" << endl;
		minimums[0] = 0.0;
	}
	if (maximums[0] == 0.0)
	{
		cout << "Using default constraint of maximum value for A: 1000" << endl;
		maximums[0] = 1000.0;
	}
	// - constraining values for Z
	if (minimums[1] == 0.0)
	{
		cout << "Using default constraint of minimum value for Z: 0.0" << endl;
		minimums[1] = 0.0;
	}
	if (maximums[1] == 0.0)
	{
		cout << "Using default constraint of maximum value for Z: 1000" << endl;
		maximums[1] = 1000.0;
	}
	
	// Maximum steps in one direction
	if (max_steps_in_one_direction == 0)
	{
		cout << "Using default max steps in one direction: OFF" << endl;
		max_steps_in_one_direction = 0;
	}
	// Defaults are set
}

/*
 Calculate ECPs to test from the initial data
 
 No param.
 */
void Powells::calculate_ecps_to_test()
{		
	//cout << ecps_to_test.size() << endl;	

	// We'll do this using a univariate method, coupled with Powell's if I get time
	gaussian g = ecps_to_test[0];
	
	// If not done already, initiate search vectors
	if (search_vectors.size() == 0)
	{
		resize_search_vectors(g);
	}
	
	// Dynamically scale all vectors
	bool values_allowed = false;
	while (!values_allowed)
	{
		// Count the number of zeros
		int zero_counter = 0;
		
		for (vector<gaussian_info>::size_type a = 0; a < g.values.size(); a++)
		{
			// g.values[a].value += search_vectors[vector_counter][g.values[a].type]*step_size[g.values[a].type][a];
			g.values[a].value += search_vectors[vector_counter][a]*step_size[g.values[a].type][a];
			
			if (g.values[a].value == 0)
			{
				zero_counter++;
			}
		}
		
		if (zero_counter == 0)
		{
			values_allowed = true;
		}
		else 
		{
			// Non-critical error. Using work around
			cout << "One of the ECP values is zero, which is unacceptable! Making another step" << endl;
		}
	}
	
	ecps_to_test.push_back(g);
	
	// Add variant with decreased value
	g = ecps_to_test[0];
	
	// Dynamically scale all vectors in opposite direction
	values_allowed = false;
	while (!values_allowed)
	{
		// Count the number of zeros
		int zero_counter = 0;
		
		for (vector<gaussian_info>::size_type a = 0; a < g.values.size(); a++)
		{			
			// g.values[a].value -= search_vectors[vector_counter][g.values[a].type]*step_size[g.values[a].type][a];
			g.values[a].value -= search_vectors[vector_counter][a]*step_size[g.values[a].type][a];
			
			if (g.values[a].value == 0)
			{
				zero_counter++;
			}
		}
		
		if (zero_counter == 0)
		{
			values_allowed = true;
		}
		else 
		{
			// Non-critical error. Using work around
			cout << "One of the ECP values is zero, which is unacceptable! Making another step." << endl;
		}
	}
	
	// Add new values
	ecps_to_test.push_back(g);	
	
	//cout << ecps_to_test.size() << endl;
	
	// Check new values for duplicates
	check_boundaries_duplicates();

	//cout << ecps_to_test.size() << endl;
}

/*
 Check if the powells method has converged - if not form new ECPs to search
 
 No param
 */
void Powells::check_converged()
{
	// Clear ECPs
	ecps_to_test.clear();
	
	// Firstly escape if we are just running a calculation
	if (ecps_tested[number_one_ranked].values.size() == 0)
	{
		minimisation_count = search_vectors.size();
	}
	// Define new variable previous_minimum and we'll compare to that each time through
	else if (compare_ecps(ecps_tested[number_one_ranked],previous_minimum))
	{
		// Pointer to value being tested
		vector_counter++;
		
		// Number of times we have been through the minima without change
		minimisation_count++;
		
		if (vector_counter >= search_vectors.size())
		{
			// If we've gone beyond powells go back to the start
			vector_counter = 0;

			// Sum of squares of search vectors to see if it is greater than 0
			int i = 0;
			
			for (vector<int>::size_type i = 0; i < search_vectors[search_vectors.size()-1].size(); i++)
			{
				i += (search_vectors[search_vectors.size()-1][i] * search_vectors[search_vectors.size()-1][i]);
			}
			
			// We will only save this vector if it is greater than 0
			// As otherwise it erases a potentially suitable search vector
			// And we could end up with stagnation of the vectors
			if (i > 0)
			{
				// Assign Powell's vector to most used vector previously, be it positive or negative
				int j = 0;
				for (vector<int>::size_type k = j + 1; k < vector_counts.size()-1; k++)
				{
					// Check if counter is bigger for the current counter
					if (vector_counts[k] > vector_counts[j])
					{
						// If so assign this value to i
						j = k;
					}
				}
				// We have the biggest value. Now we'll copy the vectors into this search vector
				search_vectors[j] = search_vectors[search_vectors.size()-1];
			}

			// Reset Powell's search vectors
			for (vector<int>::size_type i = 0; i < search_vectors[search_vectors.size()-1].size(); i++)
			{
				search_vectors[search_vectors.size()-1][i] = 0;
			}
			
			// Reset vector_counts
			for (vector<int>::size_type i = 0; i < vector_counts.size(); i++)
			{
				vector_counts[i] = 0;
			}
		}
	}
	else
	{
		// Increment Powell's search vector if not a Powell's search
		if (vector_counter < (search_vectors.size()-1))
		{
			// search vector could be negative
			int direction_factor = -1;
			int direction_sum = 0;
			
			// Check over all values for direction taken
			for (vector<gaussian_info>::size_type a = 0; a < ecps_tested[number_one_ranked].values.size(); a++)
			{
				//if (ecps_tested[number_one_ranked].values[a].value == 
				//	(previous_minimum.values[a].value + search_vectors[vector_counter][previous_minimum.values[a].type]*step_size[previous_minimum.values[a].type][a]))
				if (ecps_tested[number_one_ranked].values[a].value == 
					(previous_minimum.values[a].value + search_vectors[vector_counter][a]*step_size[previous_minimum.values[a].type][a]))
					
				{
					direction_sum++;
				}
			}
			// Compare total. If matches number of vectors we have moved in a positive direction
			int ecps_tested_one_size = ecps_tested[number_one_ranked].values.size();
			if (direction_sum == ecps_tested_one_size)
			{
				direction_factor = 1;
			}
			
			// Search vector direction has now been set correctly
			for (vector<int>::size_type i = 0; i < search_vectors[search_vectors.size()-1].size(); i++)
			{
				search_vectors[search_vectors.size()-1][i] += (direction_factor*search_vectors[vector_counter][i]);
			}
			
			// Let's put a break out here if we exceed the maximum steps in one direction
			// Default value of three? Will force an exploration of local search area as we go
			if (vector_counts[vector_counter] == max_steps_in_one_direction)
			{
				// Number of times we have been through the minima without change
				vector_counter++;
			}
		}
			
		// Add one to the count of how many times this search vector has been used
		vector_counts[vector_counter]++;
		// Restart minimisation counter
		minimisation_count = 0;
		// Save new minimum
		previous_minimum = ecps_tested[number_one_ranked];
	}
	
	// Get rank 1 and push back to ecps_to_test
	ecps_to_test.push_back(ecps_tested[number_one_ranked]);
	// ecps_to_test.push_back(previous_minimum);
	
	// To break the while loop, we need to decrease the step_size
	// This would depend on convergence at the current values
	
	int search_vectors_size = search_vectors.size();
	
	if (minimisation_count >= search_vectors_size) 
	{
		cout << "Reducing step sizes in search" << endl;
		for (vector<double>::size_type i = 0; i < step_size[0].size(); i++)
		{
			// Reduce the step size
			if (step_size[0][i] >= step_size_min[0][i])
			{
				cout << "Old Step for A [" << i << "] = " << step_size[0][i];
				step_size[0][i] *= step_reduction[0][i];
				cout << "; New Step for A [" << i << "] = " << step_size[0][i] << endl;
			}
			if (step_size[1][i] >= step_size_min[1][i])
			{
				cout << "Old Step for Z [" << i << "] = " << step_size[1][i];
				step_size[1][i] *= step_reduction[1][i];
				cout << "; New Step for Z [" << i << "] = " << step_size[1][i] << endl;
			}
		}
		// And reset minimisation_count
		minimisation_count = 0;
	}

	// Hope we are converged....
        converged = true;

        for (vector<double>::size_type i = 0; i < step_size[0].size(); i++)
        {
 	        // Check if we are above the convergence thresholds
                if ((step_size[0][i] >= step_size_min[0][i]) || (step_size[1][i] >= step_size_min[1][i]))
                {
			converged = false;	
                }
        }
	
//	if ((step_size[0] < step_size_min[0]) && (step_size[1] < step_size_min[1]))
//	{
//		//cout << "converged" << endl;
//		converged = true;
//	}
//	else

	// Sadly not. Repeat
        if (!converged)
	{
		//cout << "not converged" << endl;
		// converged = false;
		calculate_ecps_to_test(); 
	}
}

/*
 Check the boundary conditions to ensure everything is reasonable,
 and that we do not have any duplicates,
 or singularities at 0 which have no effect!
 
 No param
 */
void Powells::check_boundaries_duplicates()
{
	for (vector<gaussian>::size_type ecps_counter = 0; ecps_counter != ecps_to_test.size(); ecps_counter++)
	{
		//Check for duplicates
		for (vector<gaussian>::size_type ecps_counter_two = ecps_counter + 1; ecps_counter_two != ecps_to_test.size(); ecps_counter_two++)
		{
			if (compare_ecps(ecps_to_test[ecps_counter],ecps_to_test[ecps_counter_two]))
			{
				ecps_to_test.erase(ecps_to_test.begin()+ecps_counter_two);
				// Decrease a by one as we've just removed an item
				ecps_counter_two--;
			}
		}
		
		// Check for limits
		for (vector<gaussian_info>::size_type i = 0; i < ecps_to_test[ecps_counter].values.size(); i++)
		{
			if ((minimums[ecps_to_test[ecps_counter].values[i].type] > ecps_to_test[ecps_counter].values[i].value) ||
				(maximums[ecps_to_test[ecps_counter].values[i].type] < ecps_to_test[ecps_counter].values[i].value))
			{
				cout << "Erasing a value as it is outside the boundaries. Entry : " << ecps_to_test[ecps_counter].values[i].value << endl;;
				// Erase value
				ecps_to_test.erase(ecps_to_test.begin()+ecps_counter);
				ecps_counter--;
				i = ecps_to_test[ecps_counter].values.size();
			}
		}
	}
}
