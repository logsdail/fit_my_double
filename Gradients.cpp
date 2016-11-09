/*
 *  @file Gradients.cpp
 *  fit_my_ecp
 *
 *  Created by Andrew Logsdail on 03/06/2012.
 *  Copyright 2012 University College London. All rights reserved.
 *
 */

#include "Gradients.h"

using namespace std;

/*
 Method to take in the gradients, and separate out each region to find the details about each region
 such as min, max, average.
 
 @param[in] gradient_input Vector of lines n containing the gradients input data
 @param[in] total_centres Number of centres in the ChemShell model
 @param[in] centre_regions_total Array giving the total number of centres in each region
 @param[in] centre_regions Array of all centres, detailing their region.
 @return vector<regions_data> Vector of customised structure containing all the gradient information, by region
 */
vector<regions_data> digest_gradients(vector<string> gradient_input, int total_centres, 
				      vector<int> centre_regions_total, vector<int> centre_regions,
				      bool absolute_gradients)
{
	// Number of regions
	const int number_of_regions = 5;
	
	// Calculate useful information from output file, as used for convergence.
	// This should be pumped into an output file in a comprehendable setup.
	// First lets extract the gradients and calculate norms
	// We need some variables to store the gradients as we go
	double gradx = 0.0;
	double grady = 0.0;
	double gradz = 0.0;
	double gnorm = 0.0;
	//
	// Boolean check
	bool flag = false;
	// Setup gnorm counters
	vector<regions_data> regions(number_of_regions);
	// Initialise all counters
	for (int i = 0; i < number_of_regions; i++)
	{
		regions[i].gnorm = 0.0;
		regions[i].gnorm_max = 0.0;
		regions[i].gradx_max = 0.0;
		regions[i].grady_max = 0.0;
		regions[i].gradz_max = 0.0;
	}
	
	// Counter for current centre
	int centre_count = 0; // Actual centre minus one for array access
	// Then lets run a for loop over the data
	
	for (vector<string>:: size_type i = 0; i < gradient_input.size(); i++)
	{
		vector<string> tokens;
		// gradient_input_string = gradient_input[i];
		// string temp = gradient_input[i];
		// cout << i << " " << gradient_input.size() << endl;
		Tokenize(gradient_input[i],tokens,"= \n\t");
		// cout << "out" << endl;
		
		if ((!flag) && (tokens.size() > 2))
		{
			// block = dense_real_matrix = <number of atoms> is the line in the punch file which marks the atomic coordinates
			if (cmpStr(tokens[0],"block") && 
				cmpStr(tokens[1],"dense_real_matrix"))
			{
				flag = true;
			}
		}
		else if ((flag) && (centre_count < total_centres))
		{
			StringToNumber(gradient_input[i],gradx);
			StringToNumber(gradient_input[i+1],grady);
			StringToNumber(gradient_input[i+2],gradz);
			
			if (absolute_gradients) 
			{
				gnorm = sqrt((gradx*gradx) + (grady*grady) + (gradz*gradz));
			}
			else
			{
				gnorm = gradx + grady + gradz;
			}
			// We need to increment the line number so we are at the next variable
			// Just incremented by two as it'll be incremented again by the for loop
			i += 2;
			// We need to check what region the current atom is associated with
			regions[centre_regions[centre_count]-1].gnorm +=  gnorm;
			//cout << regions[centre_regions[centre_count]-1].gnorm << endl;
			
			// Keep a record of maxima for all regios
			if (gradx*gradx > regions[centre_regions[centre_count]-1].gradx_max*regions[centre_regions[centre_count]-1].gradx_max)
			{
				regions[centre_regions[centre_count]-1].gradx_max = gradx;
			}
			if (grady*grady >regions[centre_regions[centre_count]-1].grady_max*regions[centre_regions[centre_count]-1].grady_max)
			{
				regions[centre_regions[centre_count]-1].grady_max = grady;
			}
			if (gradz*gradz > regions[centre_regions[centre_count]-1].gradz_max*regions[centre_regions[centre_count]-1].gradz_max)
			{
				regions[centre_regions[centre_count]-1].gradz_max = gradz;
			} 
			// We don't need to square this one as it'll always be positive
			if (gnorm*gnorm > regions[centre_regions[centre_count]-1].gnorm_max*regions[centre_regions[centre_count]-1].gnorm_max)
			{
				regions[centre_regions[centre_count]-1].gnorm_max = gnorm;
			}
			// Increase centre count
			centre_count++;
		}
	}
	
	// If flag isn't true then we failed to find the gradient information
	// So this must be a system error. Mark as failed
	if (!flag)
	{	
		cout << "Failure: Converting Gradients" << endl;
		return regions;
	}
	else
	{
		cout << "Gradients: a.u. (eV/A)" << endl;
		// Now we've matched the total_regions vector with the regions vector in gaussian
		// We can run this as a loop, which is nicer to read
		for (int i = 0; i < number_of_regions; i++)
		{
			if (centre_regions_total[i] != 0)
			{
				regions[i].gnorm /= centre_regions_total[i];
				cout << "Region " << i+1 << " Gnorm Average: " << regions[i].gnorm;
				cout << "  \t (" << regions[i].gnorm*gradient_converter << ")" << endl;
			}
			else
			{
				regions[i].gnorm = -1;
			}
		}
		
		// Print the gnorm average over the whole cluster
		double gnorm_average = 0.0;
		int gnorm_regions = 0;
		for (int i = 0; i < number_of_regions; i++)
		{
			if (regions[i].gnorm != -1)
			{
				gnorm_average += regions[i].gnorm;
				gnorm_regions++;
			}
		}
		gnorm_average /= gnorm_regions;
		
		cout << "Gnorm Average:   \t" << gnorm_average << "  \t (" << gnorm_average*gradient_converter << ") over " << gnorm_regions << " regions" <<  endl;
		cout << endl;
		
		for (int i = 0; i < number_of_regions; i++)
		{
			if (regions[i].gnorm != -1)
			{
				cout << "Region " << i+1 << " Max Grad: x: " << regions[i].gradx_max;
				cout << "  \t y: " << regions[i].grady_max;
				cout << "  \t z: " << regions[i].gradz_max;
				cout << "  \t Gnorm: " << regions[i].gnorm_max << endl;
			}
		}
		cout << endl;
		
		return regions;
	}
}
