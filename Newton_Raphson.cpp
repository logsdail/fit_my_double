/*
 *  @file Newton_Raphson.cpp
 *  fit_my_ecp
 *
 *  Created by Andrew Logsdail on 01/06/2012.
 *  Copyright 2012 University College London. All rights reserved.
 *
LOG:
19/10/2013: Edited L-BFGS to do max two calls to run()
            (Previously had a while loop which was dangerous)
 */

#include "Newton_Raphson.h"

using namespace std;

/*
 Constructor
 
 @param[in/out] seed Pointer to the seed
 @param[in] flag Determines if we are doing an LBFGS or Newton-Raphson Minimisation
 */
Newton_Raphson::Newton_Raphson(int *seed, bool flag) : Outputs(seed)
{	
	//idum = seed;
	
	// Search area confines
	step_size.resize(2);
        step_size_min.resize(2);

	// L-BFGS Flag
	lbfgs_flag = flag;
}

/*
 Set up vectors for creating search vectors
 
 @param[in] Number of different directions possible
 */
void Newton_Raphson::resize_search_vectors(const unsigned int size)
{
	// Intiate search vectors
	search_vectors.resize(size);
	
	// Set up search vectors
	for (unsigned int i = 0; i < size; i++)
	{
		search_vectors[i].resize(size,0);
		
		// Set one pathway as operative in initial vector
		search_vectors[i][i] = 1;
	}

        // Set up minimiser
        // We will try first with just initalising the minimiser using most defaults
        if (lbfgs_flag)
	{
		// size of system
		// number of steps to store
		// and max. function evaluations
	       	m_pMinimizer = new scitbx::lbfgs::minimizer<double>(size, 5, 20);

        	// Print options
        	this->PrintMinimiserOptions();
	}

        // Addition: We can resize the step size vectors here, if need be.
        for (int i = 0; i < 2; i++)
        {
                if (step_size[i].size() != size)
                {
                        step_size[i].resize(size,step_size[i][0]);
                }

                if (step_size_min[i].size() != size)
                {
                        step_size_min[i].resize(size,step_size_min[i][0]);
                }
        }

}

/*
 Set the initial parameters for the Newton_Raphson optimisation
 
 @param[in] ss The step size to be initially used for A and zeta
 @param[in] sr Step reduction rate for A and zeta
 @param[in] ssm Convergence criteria, defined as the target step size to reach, for A and zeta
 @param[in] min Minimum values for A and zeta
 @param[in] max Maximum values for A and zeta
 @param[in] msod Maximum steps in any one direction - by default this is disabled
 */
void Newton_Raphson::set_parameters(vector< vector<double> > ss,
				    vector< vector<double> > sr,
				    vector< vector<double> > ssm,
				    vector<double> min,
				    vector<double> max,
				    int msod)
{
	step_size = ss;
	//step_reduction = sr;
	step_size_min = ssm;
        minimums = min;
        maximums = max;
	
	// - step size
	if (step_size[0].size() == 0)
	{
		cout << "Using default starting step size for A: 0.001" << endl;
		step_size[0].push_back(0.001);
	}
	if (step_size[1].size() == 0)
	{
		cout << "Using default starting step size for Z: 0.001" << endl;
		step_size[1].push_back(0.001);
	}
	// - step size minimum
	if (step_size_min[0].size() == 0)
	{
		cout << "Using default step size minimum for convergence for A: 1.0" << endl;
		step_size_min[0].push_back(1.0);
	}
	if (step_size_min[1].size() == 0)
	{
		cout << "Using default step size minimum for convergence for Z: 0.001" << endl;
		step_size_min[1].push_back(0.001);
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
}

/*
 Calculate ECPs to test from the initial data
 
 No param.
 */
void Newton_Raphson::calculate_ecps_to_test()
{		
	// If not done already, initiate search vectors
	if (search_vectors.size() == 0)
	{
		resize_search_vectors(ecps_to_test[0].values.size());
	}

        // cout << ecps_to_test[0].values.size() << endl;
	
        for (vector<int>::size_type vector_counter = 0; vector_counter < search_vectors.size(); vector_counter++)
	{
		gaussian g_plus = ecps_to_test[0];
		gaussian g_minus = ecps_to_test[0];
 
        	for (vector<gaussian_info>::size_type a = 0; a < ecps_to_test[0].values.size(); a++)
		{
			// cout << a << endl;
			g_plus.values[a].value += search_vectors[vector_counter][a]*step_size[g_plus.values[a].type][a];
			g_minus.values[a].value -= search_vectors[vector_counter][a]*step_size[g_minus.values[a].type][a];
		}
	
		ecps_to_test.push_back(g_plus);
		ecps_to_test.push_back(g_minus);
        }
}

/*
 Check if the Newton_Raphson has converged - if not form new ECPs to search
 
 No param
 */
void Newton_Raphson::check_converged()
{
	// Clear ECPs
	ecps_to_test.clear();

        // Clear variables for Newton_Raphson
        // vector<double> x(search_vectors.size()+1);
        // vector<double> g(search_vectors.size()+1);
        // vector<double> g2(search_vectors.size()+1); 

        vector<double> x(search_vectors.size());
        vector<double> g(search_vectors.size());
        vector<double> g2(search_vectors.size());

        //We need to work out the gradient at previous minimum 
        for (vector<gaussian>::size_type i = 0; i < ecps_tested.size(); i++)
	{
		if (compare_ecps(previous_minimum,ecps_tested[i]))
		{
			previous_minimum = ecps_tested[i];
		}
	}	

	// Work out gradients
	for (vector< vector<int> >::size_type vector_counter = 0; vector_counter < search_vectors.size(); vector_counter++)
        {
                gaussian g_plus = previous_minimum;
                gaussian g_minus = previous_minimum;

                for (vector<gaussian_info>::size_type a = 0; a < previous_minimum.values.size(); a++)
                {
                        g_plus.values[a].value += search_vectors[vector_counter][a]*step_size[g_plus.values[a].type][a];
                        g_minus.values[a].value -= search_vectors[vector_counter][a]*step_size[g_minus.values[a].type][a];
                }

		for (vector<gaussian>::size_type i = 0; i < ecps_tested.size(); i++)
        	{
                	if (compare_ecps(g_plus,ecps_tested[i]))
                	{
                        	g_plus = ecps_tested[i];
                	}
			else if (compare_ecps(g_minus,ecps_tested[i]))
			{
				g_minus = ecps_tested[i];
			}
        	}

		x[vector_counter] = previous_minimum.values[vector_counter].value;
		g[vector_counter] = (g_plus.function - g_minus.function)/(g_plus.values[vector_counter].value - g_minus.values[vector_counter].value);

		//double gradient_after = (g_plus.function - previous_minimum.function)/(g_plus.values[vector_counter].value - previous_minimum.values[vector_counter].value); // gradient after		
	        //double gradient_before = (previous_minimum.function - g_minus.function)/(previous_minimum.values[vector_counter].value - g_minus.values[vector_counter].value); // gradient before
		//double distance_between = (((g_plus.values[vector_counter].value + previous_minimum.values[vector_counter].value)/2) - ((previous_minimum.values[vector_counter].value + g_minus.values[vector_counter].value)/2)); 

		//g2[vector_counter] = (gradient_after - gradient_before) / distance_between;

		cout << endl;
		cout << "value:   " << x[vector_counter] << endl;
		cout << "dy/dx:   " << g[vector_counter] << endl;

		// Alternative method. This equivalent to working d2y/dx2.
		double functions = g_plus.function + g_minus.function - (2*previous_minimum.function);
		double distance_between = previous_minimum.values[vector_counter].value - g_minus.values[vector_counter].value; 
		distance_between *= distance_between;
		g2[vector_counter] = functions / distance_between;

                cout << "d2y/dx2: " << g2[vector_counter] << endl;
                cout << endl;


		//cout << "Newton Raphson: intersect: " << x[vector_counter] - ( previous_minimum.function / g[vector_counter] ) << endl;
		//cout << "Newton-Raphson: min/max:   " << x[vector_counter] - ( g[vector_counter] / g2[vector_counter] ) << endl;
		//}
        }

        // Pointers to the first element of variable and gradient vectors
        double* ptrX = &(*(x.begin()));
        double* ptrG = &(*(g2.begin()));

        // Put in max_cycles criteria
        size_t max_cycles = 2000;

	// Then we need to run the lbfgs
	if (lbfgs_flag)
	{

	        // Run minimizer until it returns true or we run out of cycles
                // I got this so wrong! I need two calls to run(). If the first is true, we need the gradients so.
                // If it returns false we do the second call ro run(), then we get the gradients. So max two calls.
                if((!(m_pMinimizer->run(ptrX, previous_minimum.function, ptrG))) && (m_pMinimizer->nfun() < max_cycles))
        	{
                        m_pMinimizer->run(ptrX, previous_minimum.function, ptrG);
	        }

	}
	else
	{
	        for (vector< vector<int> >::size_type vector_counter = 0; vector_counter < search_vectors.size(); vector_counter++)
        	{
			// Update x values
			// cout << ( g[vector_counter] / g2[vector_counter] ) << endl;
                        if ( g[vector_counter] != 0 && g2[vector_counter] != 0 )
			{
                		x[vector_counter] -= ( g[vector_counter] / g2[vector_counter] );
			}
		}
	}

        gaussian new_minimum = previous_minimum;

	for (vector< vector<int> >::size_type vector_counter = 0; vector_counter < search_vectors.size(); vector_counter++)
        {
                new_minimum.values[vector_counter].value = x[vector_counter];
     		//cout << "New x value1: " <<  x[vector_counter] << endl;
        }

	double delta = 0;
//        vector<double>::iterator delta_A = min_element(step_size_min[0].begin(),step_size_min[0].end());
//        vector<double>::iterator delta_Z = min_element(step_size_min[1].begin(),step_size_min[1].end());
        double delta_A = *min_element(step_size_min[0].begin(),step_size_min[0].end());
        double delta_Z = *min_element(step_size_min[1].begin(),step_size_min[1].end());

        if (delta_A > delta_Z)
        {
                delta = delta_Z;
        }
        else
        {
                delta = delta_A;
        }

	// Also we want a convergence check we use the one that comes with this library. Epsilon should be soft-coded
        // Currently this is quite tight. Default is 1e-5
        // double epsilon = 1e-10; // Convergence criteria is ||g|| < epsilon * ||x||
        // scitbx::lbfgs::traditional_convergence_test<double> is_Converged(previous_minimum.values.size(),epsilon);
        // bool convergence_test = is_Converged(ptrX, ptrG);

        // Trying a different approach
        double epsilon = 1e-2;
        bool convergence_test = true;
        for (vector< vector<int> >::size_type vector_counter = 0; vector_counter < g2.size(); vector_counter++)
        {
                //cout << vector_counter << " " << g2[vector_counter] << " " << epsilon << " " << convergence_test << endl;
        	if ( abs(g2[vector_counter]) > epsilon )
		{
			convergence_test = false;
		}
		//cout << vector_counter << " " << g2[vector_counter] << " " << epsilon << " " << convergence_test << endl;
        }

        if (lbfgs_flag)
        {
        	convergence_test = ( convergence_test || m_pMinimizer->nfun() > max_cycles );
        }

        //cout << endl;
        for (vector< vector<int> >::size_type vector_counter = 0; vector_counter < search_vectors.size(); vector_counter++)
        {
               cout << "New value1: " <<  new_minimum.values[vector_counter].value << endl;
		cout << "Previous value1: " << previous_minimum.values[vector_counter].value << endl;
        }
        //cout << endl;
        //cout << compare_ecps(previous_minimum,new_minimum,delta) << " " << convergence_test << endl;

	//if (compare_ecps(previous_minimum,new_minimum,delta))
        if ( compare_ecps(previous_minimum,new_minimum,delta) || convergence_test ) 
        {
		if (lbfgs_flag)
		{
			if (m_pMinimizer->nfun() > max_cycles)
                	{
                        	cout << "Function evaluations has exceeded " << max_cycles << " so exiting" << endl;
                        	cout << "Check this as this may be due to oscillating minima" << endl;
                	} 
			else
			{
				cout << "Convergence obtained" << endl;
        	        	cout << "gnorm: " << m_pMinimizer->euclidean_norm(&(*(g2.begin()))) << endl;
				cout << "xnorm: " << m_pMinimizer->euclidean_norm(&(*(x.begin()))) << endl;
        		        cout << "nfunc: " << m_pMinimizer->nfun() << endl;
		                cout << "iter:  " << m_pMinimizer->iter() << endl;
			}
		}

		//cout << "converged" << endl;
		converged = true;

                previous_minimum = new_minimum;

		// Delete minimiser
		if (lbfgs_flag)
		{
        		delete m_pMinimizer;
		}
	}
	else
	{
		if (lbfgs_flag)
		{
			// Why is this double called? Test with this removed.
	                // m_pMinimizer->run(ptrX, previous_minimum.function, ptrG);
	
			gaussian new_minimum = previous_minimum;
 		       	// Print out values
        		 for (vector< vector<int> >::size_type vector_counter = 0; vector_counter < search_vectors.size(); vector_counter++)
	        	{
        	        	new_minimum.values[vector_counter].value = x[vector_counter];
                		// cout << "New value2: " <<  new_minimum.values[vector_counter].value << endl;
	                }
		}
        
        	// Then we need to assign the new values to previous_minimum
	       	for (vector< vector<int> >::size_type vector_counter = 0; vector_counter < search_vectors.size(); vector_counter++)
        	{
			cout << "Old value: " <<  previous_minimum.values[vector_counter].value << endl;
               		previous_minimum.values[vector_counter].value = x[vector_counter];
			cout << "New value: " <<  previous_minimum.values[vector_counter].value << endl;
	       	}

		ecps_to_test.push_back(previous_minimum);

		//cout << "not converged" << endl;
		// converged = false;
		calculate_ecps_to_test(); 
	}

	//these are both pointers to the first value of a vector
	//vector will handle itself so I can safely derefence them
	ptrG = NULL;
	ptrX = NULL;
}

void Newton_Raphson::PrintMinimiserOptions()
{
	//if(!m_pMinimizer)
	//{
	//	throw cError("No minimiser instantiated!");
	//}
	std::cout << "======= Minimiser Settings =======\n";
	std::cout << "n: " << m_pMinimizer->n() << "\n";
	std::cout << "m: " << m_pMinimizer->m() << "\n";
	std::cout << "xtol: " << m_pMinimizer->xtol() << "\n";
	std::cout << "gtol: " << m_pMinimizer->gtol() << "\n";
	std::cout << "stpmin: " << m_pMinimizer->stpmin() << "\n";
	std::cout << "stpmax: " << m_pMinimizer->stpmax() << "\n";
}

