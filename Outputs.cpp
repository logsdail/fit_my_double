/*
 *  @file Outputs.cpp
 *  fit_my_ecp
 *
 *  Created by Andrew Logsdail on 01/06/2012.
 *  Copyright 2012 University College London. All rights reserved.
 *
 */

#include "Outputs.h"

using namespace std;

/*
 Returns a vector of type gaussian, containing ECPS to test
 
 @return vector<gaussian> ECPs to test
 */
vector<gaussian> Outputs::get_ecps_to_test()
{
	// We need to fool the system into performing just the analysis.
	// So we'll add a ghost gaussian variable to the starting_gaussian_ecps vector
	if (starting_gaussians.values.size() == 0)
	{
		ecps_to_test.push_back(get_default());
	}
		
	return ecps_to_test;
}
/*
 Prints out a default ECP, initialised for use
 
 @return gaussian ECP
 */
gaussian Outputs::get_default()
{
	gaussian g;
	g.function = 999999;
	g.failed = false;
	
	return g;
}

/*
 Compare two double arrays to see if they are the same
 
 @param[in] ecp1 Vector of type double with ECP values
 @param[in] ecp2 Vector of type double with ECP values
 @param[in] DELTA double with optional value of accuracy check
 @return bool Result of comparison check
 */
bool Outputs::compare_ecps(vector<gaussian_info> ecp1, vector<gaussian_info> ecp2, double DELTA /* = 0.000001*/)
{
	
	// Check sizes
	if (ecp1.size() != ecp2.size())
	{
		return false;
	}
	
	// Check if the values don't match
	for (vector<double>::size_type i = 0; i < ecp1.size(); i++)
	{
		// Compare to delta
		if (abs(ecp1[i].value - ecp2[i].value) > DELTA)
		{
			//cout << ecp1[i].value << "& " << ecp2[i].value << ", ";
			return false;
		}
		//cout << endl;
	}
	
	//cout << "BAD NEWS" << endl;
	// Return true as they are the same
	return true;
}

