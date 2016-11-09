/*
 *  @file Punch.cpp
 *  fit_my_ecp
 *
 *  Created by Andrew Logsdail on 03/06/2012.
 *  Copyright 2012 University of Birmingham. All rights reserved.
 *
 */

#include "Punch.h"

using namespace std;

/*
 Fireproof Constructor
 
 No params
 */
Punch::Punch()
{
	number_of_regions = 5;
	total_centres = 0;
	centre_regions_total.resize(number_of_regions,0);
	
	// Make a note of anions in region 1
	region_1_anions = 0;
}

/*
 Print method for the number of centres in each region
 
 No params
 */
void Punch::print_regions()
{
	cout << endl;
	cout << "Starting Quantity of Centres" << endl;
	// Initiate counters at 0
	for (int i = 0; i < number_of_regions; i++)
	{
		cout << "Region " << i+1 << " Centres: " << centre_regions_total[i] << endl;
	}
	cout << "Total Centres: " << total_centres << endl;
	cout << endl;
}

/*
 Print method for the minimum bond lengths between regions
 
 No params
 */
void Punch::print_bond_data()
{
	// Create Vector for minima
	vector< vector<double> > minima(number_of_regions);
	vector< vector<int> > minima_labels(number_of_regions);
	// Resize for the require regions
	const int initial_value = 99999;
	//
	for (int i = 0; i < number_of_regions; i++)
	{
		minima[i].resize(number_of_regions,initial_value);
		minima_labels[i].resize(number_of_regions);
	}

	// Calculate bond lengths
	for (vector<int>::size_type i = 0; i < centre_regions.size(); i++)
	{
		for (vector<int>::size_type j = i + 1; j < centre_regions.size(); j++)
		{
			double x_squared = centre_coords[i].x - centre_coords[j].x;
			double y_squared = centre_coords[i].y - centre_coords[j].y;
			double z_squared = centre_coords[i].z - centre_coords[j].z;
				
			x_squared *= x_squared;
			y_squared *= y_squared;
			z_squared *= z_squared;
			
			double new_value = sqrt(x_squared + y_squared + z_squared);
			
			if (minima[centre_regions[i]-1][centre_regions[j]-1] > new_value)
			{
				minima[centre_regions[i]-1][centre_regions[j]-1] = new_value;
				minima[centre_regions[j]-1][centre_regions[i]-1] = new_value;

				minima_labels[centre_regions[i]-1][centre_regions[j]-1] = i;
                                minima_labels[centre_regions[j]-1][centre_regions[i]-1] = j;
			}
		}
	}
	
	cout << endl;
	cout << "Minimum Bond Lengths:" << endl;
	// Initiate counters at 0
	for (int i = 0; i < number_of_regions-1; i++)
	{
		// Keep track of the number of print outs
		int print_out_counter = 0;
		//
		for (int j = i + 1; j < number_of_regions; j++)
		{
			if (minima[i][j] != initial_value)
			{
				print_out_counter++;
				// Check if they are the same types
				if ( i != j )
				{
					cout << "R" << i+1 << " - R" << j+1 << " ( " << centre_coords[minima_labels[i][j]].label << " - " << centre_coords[minima_labels[j][i]].label << " )";
				}
				else
				{
					cout << "R" << i+1 << " - R" << j+1;
				}
				cout << " : " << minima[i][j] << " a.u. (" << minima[i][j]*bohr_to_angstrom << " A)" << endl;
			}
		}
		// Tidy up printouts
		if (print_out_counter > 0)
		{
			cout << endl;
		}
	}
}

/*
 Method to search through each line in the input and decide if it is an atom,
 and if so what region it is in
 
 @param[in] p Punch template - this needs to be soft-coded as if we are fitting to a defect
 system then the numbers will be different and this will be re-analaysed
 */
void Punch::analyse(vector<string> p)
{
	int centres_start_point = 0;
	int centres_end_point = 0;
	int i;
	
	for (i = 0; i < size(); i++)
	{
		vector<string> tokens;
		Tokenize(p[i],tokens,"= \n\t");
		
		if ((total_centres == 0) && 
			(tokens.size() > 2))
		{
			// block = coordinates records = <number of atoms> is the line in the punch file which marks the atomic coordinates
			if (cmpStr(tokens[0],"block") && 
				cmpStr(tokens[1],"coordinates") && 
				cmpStr(tokens[2],"records"))
			{
				StringToNumber(tokens[3],total_centres);
				centre_regions.resize(total_centres,0);
				centre_coords.resize(total_centres);
				
				centres_start_point = i;
				centres_end_point = centres_start_point+total_centres;
				// cout << centres_end_point << " " << centres_start_point << endl;
			}
		} 
		// Identify if the atom of this line number is region 1-5
		else if (i <= centres_end_point)
		{
			int array_pointer = i-centres_start_point-1;
			// This'll give the region in one nice little move;
			StringToNumber(tokens[0].substr(tokens[0].length()-1),centre_regions[array_pointer]);

                        // Shrink tokens to just contain the label
                        //centre_coords[array_pointer].label = tokens[0].substr(0,tokens[0].length()-1);
                        tokens[0] = tokens[0].substr(0,tokens[0].length()-1);

			// Store data for geometry analysis
                        centre_coords[array_pointer].label = tokens[0];
			StringToNumber(tokens[1],centre_coords[array_pointer].x);
			StringToNumber(tokens[2],centre_coords[array_pointer].y);
			StringToNumber(tokens[3],centre_coords[array_pointer].z);
			
			// Add 1 to counter for this region
			centre_regions_total[centre_regions[array_pointer]-1]++;
			
			if (centre_regions[array_pointer] == 1)
			{
				// Make a note if this is an anionic species for determining spread of eigenvalues
                                bool r1_anion = true;
								
				for (size_t ichar = 0; ichar < tokens[0].length(); ichar++ )
                                {
					// Convert string to uppercase
					tokens[0][ichar] = toupper(tokens[0][ichar]);
			 	}

                                if (tokens[0].length() == anion_species.length())
				{
					for (size_t ichar = 0; ichar < tokens[0].length(); ichar++ )
                                	{
						//cout << anion_species << " " << tokens[0][ichar] << endl;
						if (anion_species[ichar] != tokens[0][ichar]) 
						{
							r1_anion = false;
						}
					}
                                }
				else
				{
					r1_anion = false;
				}

				if (r1_anion) 
				{
					region_1_anions++;
				}

				// Make a note of R1 species
                                bool r1_exists = false;
                                for (vector<string>::size_type i = 0; i < r1_species.size(); i++)
                                {
                                        // Reusing r1 anion to see if this species has been noted as present in R1
                                        r1_anion = true;

					//cout << tokens[0] << " " << tokens[0].length() << " " << r1_species[i] << " " << r1_species[i].length() << endl;

					// Get rid quickly of any mismatching sized strings
					if (tokens[0].length() == r1_species[i].length())
					{

                                		for (size_t ichar = 0; ichar < tokens[0].length(); ichar++ )
                                		{
                                 		       	// cout << r1_species[i] << " " << tokens[0] << endl;
                                        		if (tokens[0][ichar] != r1_species[i][ichar])
                                        		{
                                                		r1_anion = false;
                                        		}
                                		}

                        	                if (r1_anion)
                	                        {
							r1_exists = true;
	                                        }
					}
					//cout << r1_anion << endl;
                                }
                                if (!r1_exists)
                                {
					r1_species.push_back(tokens[0]);
                                }
			}
		}
		// We can then match this information to the gradients and energy at the end
	}
}

/*
 Compares two punch templates to see if a defect has been introduced
 
 @param[in] punch_output Updated punch file
 */
void Punch::compare(vector<string> punch_output)
{
	// This should only run once in the scenario where we have introduced a defect
	// Needs full testing but I'm confident
	if (punch_output.size() != punch_template.size())
	{
		cout << endl;
		cout << "Punch output and template do not have the same number of atoms" << endl;
		cout << "Recounting centres for analysis purposes, and saving these values" << endl;
		cout << endl;
		
		// We need to reset some counters
		// Initiate counters at 0
		total_centres = 0;
		centre_regions_total.resize(number_of_regions,0);
		region_1_anions = 0;
		
		// Clear data
		centre_regions.clear();
		centre_coords.clear();
		
		// Update the numbers to reflect the potentially defective system
		analyse(punch_output);
	}
}

