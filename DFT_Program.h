/*
 *  @DFT_Program.h
 *  fit_my_ecp
 *
 *  @brief Deals with translating inputs/outputs for DFT_Programs
 *  Hopefully we can use this modular setup for other codes in the future
 *
 *  Created by Andrew Logsdail on 18/01/2013.
 *  Copyright 2012 University College London. All rights reserved.
 *
 */


#ifndef DFT_PROGRAM_H
#define DFT_PROGRAM_H
#include <iostream>
// Personal headers
#include "Utils.h"
#include "Structures.h"

class DFT_Program {
	
public:
	
        /*
         Constructor

	 Initialise variables
	*/
	DFT_Program()
	{
        	region_1_anion_offset = 0;
        	anion_species = "";
	}
	
	/*
	 Deconstructor
	 
	 No params
	 */
	virtual ~DFT_Program(){}

	/*
	Set region offset for anions
	
	@param[in] o Number of anions in region 1
	*/
        void set_anion_offset(int o)
	{
		region_1_anion_offset = o;
	}

        /*
        Set label for anions

        @param[in] a Label for anion species
        */
        void set_anion_species(std::string a)
        {
                anion_species = a;
        }
	
	/*
	 Saves the ecp template to file
	 
	 @param[in] s Input strings
	 */
	void set_ecp_template(std::vector<std::string> s)
	{
		ecp_template = s;
		calculate_starting_gaussian_ecps();
	}
	
	/*
	 Returns the ecp template
	 
	 @return vector<string> Data on the ecp template
	 */
	std::vector<std::string> get_ecp_template()
	{
		return ecp_template;
	}
	
        /*
         Stub for returning the ECP input

         To be over-written by inheriting class
         */
	virtual std::vector<std::string> get_ecp_template(gaussian g)
        {
  		return ecp_template;
        }
	
	/*
	 Return the starting gaussian ecps
	 
	 @return gaussian ECPs from input
	 */
	gaussian get_starting_gaussian_ecps()
	{
		// return starting_gaussian_ecps;
		
		return starting_gaussians;
	}

        /*
         Stub for analysing the contents of the QM output

         To be over-written by inheriting class
         */
        virtual gaussian digest_electronic(std::vector<std::string> qm_output, gaussian g, int r1_anion, std::vector<std::string> r1_species)
        {
		return g;
        }

        virtual std::string type()
	{
		return "NULL";
	}
	
protected:
	
	std::vector<std::string> ecp_template;
	
	gaussian starting_gaussians;
	
	int region_1_anion_offset;

        std::string anion_species;
	
	/*
	 Return the size of the ecp_template
	 
	 @return int Size of array
	 */
	int size()
	{
		return ecp_template.size();
	}

        /*
         Stub for calculation of starting gaussians

	 To be over-written by inheriting class
         */

        virtual void calculate_starting_gaussian_ecps(){}
	
};

#endif

