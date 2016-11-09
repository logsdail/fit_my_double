/*
 *  @Nwchem.h
 *  fit_my_ecp
 *
 *  @brief Deals with translating inputs/outputs with NWChem
 *  Hopefully we can use this modular setup for other codes in the future
 *
 *  Created by Andrew Logsdail on 03/06/2012.
 *  Copyright 2012 University College London. All rights reserved.
 *
 */


#ifndef NWCHEM_H
#define NWCHEM_H
#include <iostream>
#include <algorithm>
// Personal headers
#include "Utils.h"
#include "Structures.h"
#include "DFT_Program.h"

class Nwchem : public DFT_Program {
	
public:
	/*
         Constructor

         No params */
	Nwchem(){}
	
	/*
	 Deconstructor
	 
	 No params
	 */
	~Nwchem(){}

        std::vector<std::string> get_ecp_template(gaussian g);

	gaussian digest_electronic(std::vector<std::string> nwchem_output, gaussian g, int r1_anion, std::vector<std::string> r1_species);

        std::string type()
	{ return "NWCHEM"; }
	
private:
	
	void calculate_starting_gaussian_ecps();
	
	int digest_rhf(std::vector<std::string> tokens, 
				   std::vector<std::string> *nwchem_output, 
				   int *MO_line_number, gaussian *g, int i, int en);
	
	int digest_uhf(std::vector<std::string> tokens, 
				   std::vector<std::string> *nwchem_output, 
				   int *MO_line_number, gaussian *g, int i, int en, int spin_counter);

        int digest_dma(std::vector<std::string> tokens,
                                   std::vector<std::string> *nwchem_output,
                                   int *DMA_line_number, gaussian *g, int i);
	
};

#endif

