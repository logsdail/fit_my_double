/*
 *  @file Linear.cpp
 *  fit_my_ecp
 *
 *  Created by Andrew Logsdail on 05/11/2013.
 *  Copyright 2012 University College London. All rights reserved.
 *
 */

#include "Linear.h"

using namespace std;

/*
 Constructor
 
 @param[in/out] seed Pointer to the seed
 */
Linear::Linear(int *seed) : Outputs(seed)
{	
	//idum = seed;
	
	// Search area confines
	step_size.resize(2);
	minimums.resize(2,0.0);
	maximums.resize(2,0.0);
}

/*
 Set up vectors for creating Powell's search vector
 
 @param[in] Number of different directions possible
 */
void Linear::resize_search_vectors(gaussian g)
{
	const unsigned int size = g.values.size();
	// Temporary values
	vector<double> new_step_size;
	int counter = 0;

	//cout << "RESIZING" << endl;

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
				//cout << g.values[j].value << " " << g.values[j].type << endl;
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
				//cout << i << " " << j << " " << new_step_size[j] << endl;
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
	}
}

/*
 Set the initial parameters for the Linear optimisation
 
 @param[in] ss The step size to be initially used for A and zeta
 @param[in] sr Step reduction rate for A and zeta
 @param[in] ssm Convergence criteria, defined as the target step size to reach, for A and zeta
 @param[in] min Minimum values for A and zeta
 @param[in] max Maximum values for A and zeta
 @param[in] msod Maximum steps in any one direction - by default this is disabled
 */
void Linear::set_parameters(vector< vector<double> > ss,
			     vector< vector<double> > sr,
			     vector< vector<double> > ssm,
			     vector<double> min,
			     vector<double> max,
			     int msod)
{
	step_size = ss;
	minimums = min;
	maximums = max;

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
	
	// Defaults are set
}

/*
 Calculate ECPs to test from the initial data
 
 No param.
 */
void Linear::calculate_ecps_to_test()
{		
	//cout << "Creating ECPS" << endl;
	//cout << ecps_to_test.size() << endl;	

	// We'll do this using a univariate method, coupled with Powell's if I get time
	gaussian g = starting_gaussians;
	
	// If not done already, initiate search vectors
	if (g.values.size() > 1)
	{
		resize_search_vectors(g);
	}

	//cout << step_size[0][0] << step_size[0][1] << step_size[1][0] << step_size[1][1] << endl;
	//cout << step_size[0].size() << endl;
	
        // Initialise with minima
        for (vector<gaussian_info>::size_type a = 0; a < g.values.size(); a++)
        {
		g.values[a].value = minimums[g.values[a].type];

		// Protect against an instance where the starting point is an invalid number
		if (g.values[a].value == 0.0)
		{
			g.values[a].value += step_size[g.values[a].type][a];
			//cout << g.values[a].value << " " << minimums[g.values[a].type] << " " << step_size[g.values[a].type][a] << endl;
		}
        }

        ecps_to_test.push_back(g);

        // Now lets do all the directions
	// Give ourselves a variable to store the ecps_to_test length
	// As this is altered in each loop

	unsigned int ecps_to_test_length = ecps_to_test.size();

	// cout << "G VALUES SIZE = " << g.values.size() << endl;

	// Add all possible combinations on to the previous collections of results
        for (vector<gaussian_info>::size_type a = 0; a < g.values.size(); a++)
	{
		ecps_to_test_length = ecps_to_test.size();
		// cout << "ECPS LENGTH = " << ecps_to_test_length << endl;

		for (vector<gaussian>::size_type i = 0; i < ecps_to_test_length; i++)
		{
			// cout << ecps_to_test_length << " " << i << endl;
			g = ecps_to_test[i];

			//cout << g.values[a].value << " " << maximums[g.values[a].type] << " " << step_size[g.values[a].type][a] << endl;

			while (g.values[a].value < maximums[g.values[a].type])
        		{
				// cout << g.values[a].value << " " << maximums[g.values[a].type]	<< endl;
	                	if (g.values[a].value != 0.0)
        	        	{
                	        	ecps_to_test.push_back(g);
                		}

	                	g.values[a].value += step_size[g.values[a].type][a];
        		}
		}
	}

	//cout << "Created ECPs" << endl;
	//cout << ecps_to_test.size() << endl;

	// Insert Starting Gaussian For Comparison
	ecps_to_test.insert(ecps_to_test.begin(),starting_gaussians);
	
	// Check new values for duplicates
	check_boundaries_duplicates();

	//cout << "Tested ECPs" << endl;
	//cout << ecps_to_test.size() << endl;
}

/*
 Check the boundary conditions to ensure everything is reasonable,
 and that we do not have any duplicates,
 
 No param
 */
void Linear::check_boundaries_duplicates()
{
	for (vector<gaussian>::size_type ecps_counter = 0; ecps_counter != ecps_to_test.size(); ecps_counter++)
	{
		//Check for duplicates
		for (vector<gaussian>::size_type ecps_counter_two = ecps_counter + 1; ecps_counter_two != ecps_to_test.size(); ecps_counter_two++)
		{
			if (compare_ecps(ecps_to_test[ecps_counter],ecps_to_test[ecps_counter_two]))
			{
				// cout << "ECP " << ecps_counter << " IS THE SAME AS " << ecps_counter_two << endl;
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
				cout << "Erasing a value as it is outside the boundaries. Entry : " << ecps_to_test[ecps_counter].values[i].value;
				cout << " of Gaussian " << i << endl;
				// Erase value
				ecps_to_test.erase(ecps_to_test.begin()+ecps_counter);
				ecps_counter--;
				i = ecps_to_test[ecps_counter].values.size();
			}
		}
	}
}
