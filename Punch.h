/*
 *  @Punch.h
 *  fit_my_ecp
 *
 *  @brief Deals with reading through the punch file and deciphering it
 *
 *  Created by Andrew Logsdail on 03/06/2012.
 *  Copyright 2012 University College London. All rights reserved.
 *
 */

#ifndef PUNCH_H
#define PUNCH_H
#include <iostream>
#include <ctype.h>
// Personal headers
#include "Utils.h"
#include "Structures.h"

class Punch {

public:
	
	Punch();
	
	/*
	 Deconstructor
	 
	 No params
	 */
	~Punch(){}
	
	/*
	 Sets the contents of the punch template file, and calls analyse.
	 We need to save this for outputing in each run
	 
	 @param[in] s Vector containing the contents of the punch.template file
	 */
	void set_punch_template(std::vector<std::string> s, std::string a)
	{
		// Transform the anion_species to uppercase for comparisons
		anion_species = a;
                for (size_t ichar = 0; ichar < a.length(); ichar++ )
                {
                	// cout << anion_species << endl;
                        anion_species[ichar] = toupper(a[ichar]);
                }

		punch_template = s;
		analyse(s);
	}
	
	void compare(std::vector<std::string> punch_output);
	
	void print_regions();
	
	void print_bond_data();
	
	/*
	 Returns the total number of centres in the calculation
	 
	 @return int Total centres
	 */
	int get_total_centres()
	{
		return total_centres;
	}

        /*
         Return the anion type for region 1 (For parameterising)
         
         @return String anion type in region 1 
         */
        std::string get_region_1_anion_species()
        {
                return anion_species;
        }
	
	/*
	 Return the number of selected in region 1 (For parameterising)
	 
	 @return int Total Anions in Region 1
	 */
	int get_region_1_anions()
	{
		return region_1_anions;
	}
	
	/*
	 Return the regions of each centre
	 
	 @return vector<int> Vector array containing the region numbers, of size total_centres
	 */
	std::vector<int> get_centre_regions()
	{
		return centre_regions;
	}
	
	/*
	 Return the number of centres in each region
	 
	 @return vector<int> Vector array containing the quantities, in each region, of size 5
	 */
	std::vector<int> get_centre_regions_total()
	{
		return centre_regions_total;
	}
	
	/*
	 Return the template punch for runtime
	 
	 @return vector<string> String vector of the entire original template, without any distortions
	 */
	std::vector<std::string> get_punch_template()
	{
		return punch_template;
	}

	/*
	 Return an array containing the characters of the possible region 1 atoms

         @return vector<string> String vector with character list of atoms present in R1
         */
        std::vector<std::string> get_region_1_species()
	{
		return r1_species;
        }

	/*
	 Debug to print the atoms in region 1. Let's hide this away for now
	 */
        void print_region_1_species()
	{
		std::cout << "Atoms in Region 1: ";
	        for (std::vector<int>::size_type i = 0; i < r1_species.size(); i++)
        	{
			std::cout << r1_species[i] << " ";
		}
		std::cout << std::endl;
	}
	
private:
	
	struct xyz
	{
		std::string label;
		double x;
		double y;
		double z;
	};
	
	int number_of_regions;
	int total_centres;
	int region_1_anions;

        // anionic species
        std::string anion_species;
	
	std::vector<std::string> punch_template;
	
	std::vector<int> centre_regions;
	std::vector<xyz> centre_coords;
	std::vector<int> centre_regions_total;

        std::vector<std::string> r1_species;
	
	void analyse(std::vector<std::string>);
	
	/*
	 Returns the size of the punch template, in lines
	 
	 @return int Number of lines
	 */
	int size()
	{
		return punch_template.size();
	}
};

#endif
