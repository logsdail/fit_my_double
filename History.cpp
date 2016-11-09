/*
 *  @file History.cpp
 *  fit_my_ecp
 *
 *  Created by Andrew Logsdail on 02/06/2012.
 *  Copyright 2012 University College London. All rights reserved.
 *
 */

#include "History.h"

using namespace std;

/*
 Read in the data from an old file, and update history to reflect this
 
 @param[in] input Data read in from old log file
 */
// void History::insert_old_data(vector<string> input, Functions *func_calc)
void History::insert_old_data(vector<string> input, vector<string> regions_input)
{
         /** NOW THE ONLY WORKING HEADER TYPE
	 This will be the newest header:
	 Entry      |       Gnorm 1 (70)    Gnorm 2 (20)    Gnorm 3 (0.1)   An. Spread (1)         HOMO (0)       LUMO (0)  |       Function (Rank) 
	 Coupled with a regions.log of format:
	 Entry   |               Line            Type            Value
	 ================================================================================================
	 0       |               34              0               298.22733
		 |               34              1               715.1192
	 ================================================================================================
	 
	 **/

	if (input.size() > 0)
	{
		string buffer = "";
		vector< vector<gaussian_info> > regions_present;
		
		// Check if we have a regions file
		if (regions_input.size() > 0)
		{
			regions_present = decompose_regions(regions_input);
		}
		
		vector<string> tokens;
		Tokenize(buffer,tokens,"|= \n\t");

		/*********
		// Force recalculate if requested means this will not be checked
		// And all functions will be recalculated. Expensive?
	        // Let's compare the weights to our current values
	        vector<double> old_weights;
        	double temp_weight;
	        // Read in the weights
        	for (vector<string>::size_type i = 0; i < tokens.size(); i++)
	        {
        		//cout << tokens[i] << endl;
                	if (tokens[i][0] == '(' && tokens[i][tokens[i].length()-1] == ')')
	                {
        	        	buffer = tokens[i].substr(1,tokens[i].length()-2);
                	        StringToNumber(buffer,temp_weight);
	                        old_weights.push_back(temp_weight);
	                }
        	}
	        // Now let's check if we need to recalculate each loop
        	bool same_as_before = func_calc->check_weights(old_weights);
                // cout << same_as_before << endl; 

		// Check if we've included that targets in the output
		buffer = input[2];
		tokens.clear();
		Tokenize(buffer,tokens,"()|=* \n\t");
		bool targets_included = false;
                
		if (cmpStr(tokens[0],"Targets"))
		{
			//cout << same_as_before << endl;
			targets_included = true;
			// We need to compare the targets, otherwise tyhe function is wrong
			vector<double> old_targets;
			double temp_target;

			for (vector<string>::size_type i = 1; i < tokens.size(); i++)
			{
				StringToNumber(tokens[i],temp_target);
				old_targets.push_back(temp_target);
			}

			// Check the targets
			if (same_as_before)
			{
				same_as_before = func_calc->check_targets(old_targets,dataset);
				//cout << same_as_before << endl;
			}
		}

		// Force recalculation of Function.
		// This makes all of the above redundant, and should be restructured so that we don't check
		// unless we need to! But for a rainy day.
		if (force_recalc)
		{
			same_as_before = false;
		} 
		******/

		// Start from 1
		for (vector<string>::size_type i = 3; i < input.size(); i++)
		{
			//cout << i << " " << input.size() << endl;
			buffer = input[i];
			//cout << buffer << endl;
			tokens.clear();
			Tokenize(buffer,tokens,"()|= \n\t");
			// Check if we have any useful information
			if (tokens.size() > 1)
			{
				gaussian g;
				
				// Get gaussian_info from regions
                                // Assume always new file types
				StringToNumber(tokens[0],g.index);
					
				//cout << g.index << " " << regions_present.size() << endl;
				g.values = regions_present[g.index];
				
				// assign new index value
				// g.index = ecp_history.size();
				
				// Setup gnorm counters
				int number_of_regions = 5;
				g.regions.resize(number_of_regions);
                                g.orbital_spread.resize(1);
				// Read in
				StringToNumber(tokens[1],g.regions[0].gnorm);
				StringToNumber(tokens[2],g.regions[1].gnorm);
				StringToNumber(tokens[3],g.regions[2].gnorm);
				StringToNumber(tokens[4],g.orbital_spread[0].spread);
				g.orbital_spread[0].spread /= hartree_to_eV; // Scale back to AU
				StringToNumber(tokens[5],g.HOMO_value);
				g.HOMO_value /= hartree_to_eV; // Scale back to AU
                                // New to spin-polarised fitting
                                StringToNumber(tokens[6],g.LUMO_value);
                                g.LUMO_value /= hartree_to_eV; // Scale back to AU

				StringToNumber(tokens[7],g.function);
				// Tokens[8] should be the rank from the previous calculation
				// This needs improving but for now it works
				if (tokens.size() > 9)
				{
                                        // This only applies to a Gamess-UK calculation
					g.dma_spread.resize(1);
					StringToNumber(tokens[9],g.dma_spread[0].spread);
				}

				// Set rank to arbitrary value of 0
				g.rank = 0;
				g.failed = false;

                                //cout << "Getting Function Value" << endl;
				
				// If this is 0, this is a dud value
				if (g.function > 0) 
				{
					/******
					if (g.function != 888888)
					{
					        // Create a temporary vector to hold the DMA spreads
                                        	vector<double> temp(g.dma_spread.size());
                                        	for (vector<double>::size_type i = 0; i < temp.size(); i++)
                                        	{
                                                	temp[i] = g.dma_spread[i].spread;
                                        	}

						// Lets recalculate the function with new weightings IF NECESSARY
                                                if (!same_as_before)
						{
							g.function = func_calc->calculate_function(g.regions[0].gnorm,
												   g.regions[1].gnorm,
												   g.regions[2].gnorm,
												   g.orbital_spread[0].spread,
												   g.HOMO_value,
        	                                                                                   g.LUMO_value,
												   temp,
												   dataset);
						}
					}
					*******/
					// Add to history
					add_to_history(g);
				}
                                //cout << "Continue" << endl;
			}
			//cout << "NEXT LOOP" << endl;
		}
		//cout << "OUTSIDE LOOP" << endl;
	}
        //cout << "Exiting insert_old_data" << endl;
}

/*
 Recalculate functions at input for History contents 

 @param[in] *func_calc pointer to the function calculator
 @param[in] dataset counter for appropriate data in targets/weights
 @param[in] force_recalc boolean whether we should just recalc straight away
 */
void History::recalc_function(Functions *func_calc, int dataset, bool force_recalc)
{
	// Check the history is not empty
	if (size() > 0)
        {

	        // Check the functions are the same now as they were before
		if (!force_recalc)
		{
			// Create a temporary vector to hold the DMA spreads
                        vector<double> temp(ecp_history[0].dma_spread.size());
                        for (vector<double>::size_type i = 0; i < temp.size(); i++)
                        {
                        	temp[i] = ecp_history[0].dma_spread[i].spread;
                        }

			if (ecp_history[0].function != (func_calc->calculate_function(ecp_history[0].regions[0].gnorm,
                                                        ecp_history[0].regions[1].gnorm,ecp_history[0].regions[2].gnorm,
                                                        ecp_history[0].orbital_spread[0].spread,
                                                        ecp_history[0].HOMO_value,
                                                        ecp_history[0].LUMO_value,temp,dataset)))
			{
				// If the values don't match - recalc
				force_recalc = true;
			}
		}

		// Recalculate functions for all Gaussians if necessary
		if (force_recalc)
		{
			cout << "Recalculating all Functions in History file" << endl;

	                // Loop through the history file.
        	        for (vector<gaussian>::size_type i = 0; i < ecp_history.size(); i++)
                	{
			
				// Don't recalculate for duds
				if (ecp_history[i].function != 888888)
	                        {
	
        	        	      	// Create a temporary vector to hold the DMA spreads
                	                vector<double> temp(ecp_history[i].dma_spread.size());
                                	for (vector<double>::size_type j = 0; j < temp.size(); j++)
                                	{
                                		temp[j] = ecp_history[i].dma_spread[j].spread;
                                	}

	                                // Recalculate the function with new weightings IF NECESSARY
                        	       	ecp_history[i].function = func_calc->calculate_function(ecp_history[i].regions[0].gnorm,
                                                                             ecp_history[i].regions[1].gnorm,ecp_history[i].regions[2].gnorm,
                                               	                             ecp_history[i].orbital_spread[0].spread,
                                                       	                     ecp_history[i].HOMO_value,
                                                               	             ecp_history[i].LUMO_value,
                                                                             temp,dataset);
        	                }
			}
		}
	}
}

/*
 Add data from regions file into History

 @param[in] regions_input string containing all the information from regions.log
 @param[out] regions_output Decomposed information
 */
vector< vector<gaussian_info> > History::decompose_regions(vector<string> regions_input)
{
	const int max_regions = 100000;
	vector< vector<gaussian_info> > regions_output(max_regions);
	vector<gaussian_info> current_region;
	int index = 0;
	
	// Start from line number 2
	for (vector<string>::size_type i = 2; i < regions_input.size(); i++)
	{
		string buffer = regions_input[i];
		vector<string> tokens;
		
		if (buffer[0] != '=')
		{
			Tokenize(buffer,tokens,"()= \n\t");
			int offset = 0;
			
			// Check for index
			if (tokens[0][0] == '|')
			{
				offset = 1;
			}
			else if (tokens[1][0] == '|')
			{
				// Get old index value
				StringToNumber(tokens[0],index); 
			}
			
			// Assign values from regions information
			gaussian_info gaussian_new;
			StringToNumber(tokens[2-offset],gaussian_new.line_number);
			StringToNumber(tokens[3-offset],gaussian_new.type);
			StringToNumber(tokens[4-offset],gaussian_new.value);
			
			current_region.push_back(gaussian_new);
		}
		else
		{
			// Push back complete gaussian_info for this index
			// regions_output.push_back(current_region);
			regions_output[index] = current_region;
			current_region.clear();
		}
	}
	
	return regions_output;
}

/*
 Add gaussian to history contents
 
 @param[in] g Input gaussian
 */
void History::add_to_history(gaussian g)
{
	int does_exist = check_history(g);

        //cout << "Does exist? " << does_exist << " for " << g.values[0].value << endl;
	
	if (does_exist == -1)
	{
		// Create index if not set
		if (g.index == 0)
		{
			g.index = ecp_history.size();
		}
		// Push back on to file
		ecp_history.push_back(g);
	}
}

#define DELTA (0.000001)
/*
  Compare current gaussian against those in history
 
 @param[in] g Gaussian to check through the history for
 @return int Index of the gaussian if found
 */
int History::check_history(gaussian g)
{
	int does_exist = -1;
	
        //cout << size() << endl;

	if (size() > 0)
	{
		// Loop through the history file. 
		for (vector<gaussian>::size_type i = 0; i < ecp_history.size(); i++)
		{
			int k = 0;

 			//cout << i << " " << g.values.size() << " " << ecp_history.size() << endl;
			
			// Check if we are copying from history OK
			for (vector<double>::size_type j = 0; j < g.values.size(); j++)
			{
				//cout << g.values[j].value << endl;
                                //cout << ecp_history.size() << i << endl;
                                //cout << ecp_history[i].values.size() << endl;
                                //cout << ecp_history[i].values[j].value << endl;
                                //cout << abs(g.values[j].value-ecp_history[i].values[j].value) << endl;
                                //cout << DELTA << endl;
				if (abs(g.values[j].value-ecp_history[i].values[j].value) < DELTA)
				{
					//cout << "True" << endl;
					k++;
				}
			}
			
			int g_values_size = g.values.size();
			
			if (k == g_values_size)
			{
				does_exist = i;
			}
		}
	}

        //cout << does_exist << endl;
	
	return does_exist;
}
#undef DELTA
