/*
 *  @file Nwchem.cpp
 *  fit_my_ecp
 *
 *  Created by Andrew Logsdail on 03/06/2012.
 *  Copyright 2012 University of Birmingham. All rights reserved.
 *
 */

#include "Nwchem.h"

using namespace std;

/*
 Method taking the input ECP file and separating out the input ECPs
 
 No params
 */
void Nwchem::calculate_starting_gaussian_ecps()
{
	// starting_gaussian_ecps.resize(0);
	// Maybe we can implement this to deal with contracted in the long term.
	// We will need this to deal with contracted Gaussians straight off!
	for (vector<string>::size_type i = 0; i < ecp_template.size(); i++)
	{
		// Break down the string
		vector<string> tokens;
		Tokenize(ecp_template[i],tokens,"!? \n\t");
		// Check for presence of exclamations
		size_t first = ecp_template[i].find_first_of('!');
		size_t last  = ecp_template[i].find_last_of('!');
		
		if ((ecp_template[i][0] == '!') ||
			((first != string::npos) && (last != string::npos) && (first != last)))
		{
			// We need to break this down. Firstly we'll do the Coefficient
			gaussian_info g;

			// Break down the string
			// vector<string> tokens;
			// Tokenize(ecp_template[i],tokens,"! \n\t");
			
			// Store line number
			g.line_number = i;
			
			// We need to make sure values is the right size
			StringToNumber(tokens[1],g.value);
			// This is a exponent 
			g.type = 1;
			
			starting_gaussians.values.push_back(g);
			
			// Now to do the coefficient
			// gaussian_info g;
			
			// Store line number
			g.line_number = i;
			
			// We need to make sure values is the right size
			StringToNumber(tokens[2],g.value);
			// This is an coefficient
			g.type = 0;
			
			starting_gaussians.values.push_back(g);
		}
		else if (first != string::npos)
		// An exclamation exists on this line
		{
			string temp = ecp_template[i].substr(first+1);
			
			// cout << temp << endl;
			
			// Break down the string
			vector<string> tokens_short;
			Tokenize(temp,tokens_short,"! \n\t");
			
			// cout << tokens_short[0] << endl;
			
			// Create gaussian
			gaussian_info g;
			
			// Store line number
			g.line_number = i;
			
			// We need to make sure values is the right size
			if (tokens_short.size() > 0)
			{
				StringToNumber(tokens_short[0],g.value);
			}
			else
			{
				cout << "Critical Error. We have a rogue exclamation mark in the ECP file, with no trailing value." << endl;
				// Exit as we won't manage to run
				exit(EXIT_FAILURE);	
			}
			
			// Check size. If size = 1, we have coefficient. If size = 2, we have an exponent
			if (tokens_short.size() == 1)
			{
				g.type = 0;
			}
			else
			{
				g.type = 1;
			}

			
			starting_gaussians.values.push_back(g);
		}
	}
}

/*
 Take an input ECP, add the data where appropriate into the ECP template file
 and then push this data back out in a vector<string> format
 
 @param[in] g Gaussian ECP
 @return vector<string> Contents of new ECP file
 */
vector<string> Nwchem::get_ecp_template(gaussian g)
{
	// Set the out data to ecp template and then edit
	vector<string> outData = ecp_template;;
	string sentence = "";
	
	for (vector<int>::size_type i = 0; i < g.values.size(); i++)
	{
		// Break down the string
		vector<string> tokens;
		Tokenize(outData[g.values[i].line_number],tokens,"! \n\t");

		// Originally, I designed this subroutines so that:
		// Coefficient = 0
		// Exponent = 1
		// As this is the ordering in the Gamess-UK ECP input

		// In NWChem these have been defined the other way round. But to make it all
		// transferable I have maintained the same coding.
		// However in NWCHEM the ecp template the order is exponent, coefficient.
		
		// Replace the value as needed
		NumberToString(g.values[i].value,tokens[abs(g.values[i].type-1)+1]);
		// This is plus one to offset the r factor at the start of each line
		
		// Now reconstruct string and submit back to outData
		ostringstream sentence;
		sentence.precision(8);
		for (vector<string>::size_type j = 0; j < tokens.size(); j++)
		{
			sentence << tokens[j] << "\t\t";
		}
		
		outData[g.values[i].line_number] = sentence.str();
	}
	 
	return outData;
}

/*
 Decodes the NWChem output
 
 @param[in] nwchem_output Vector containing string output
 @param[in] g Current gaussian being investigated
 @param[in] r1_anions Number of anions in R1 
 @param[in] r1_species Array of types for R1 species
 @return gaussian Updated information on the ECP
 */
gaussian Nwchem::digest_electronic(vector<string> nwchem_output, gaussian g, int r1_anions, vector<string> r1_species)
{
	// Zero all the counters we use
	int MO_line_number = 0;
        int DMA_line_number = 0;
	int number_electrons = 0;
        int alpha_electrons = 0;
        int beta_electrons = 0;
	int spin_counter = 0;
        int spin_polarised = 0;
	// Zero gaussian counters
	//g.anion_max = -1000.0;
	//g.anion_min = 1000.0;
	//g.anion_average = 0;
	//g.anion_spread = 0.0;
	g.HOMO_value = -100.0; // Arbitrarily low
	g.LUMO_value = 100.0; // Arbitrarily high
        // Zero anion 1s counters
        g.orbital_spread.resize(1);
        g.orbital_spread[0].max = -1000.0;
        g.orbital_spread[0].min =  1000.0;
        g.orbital_spread[0].average = 0.0;
        g.orbital_spread[0].spread  = 0.0;
        g.orbital_spread[0].quantity = r1_anions;
        g.orbital_spread[0].label    = anion_species;
        //Zero dma spread counters, and resize
        g.dma_spread.resize(r1_species.size());
	for (vector<min_max_spread>::size_type i = 0; i < g.dma_spread.size(); i++)
	{
		g.dma_spread[i].max = -1000.0;
	      	g.dma_spread[i].min =  1000.0;
	        g.dma_spread[i].average = 0.0;
	        g.dma_spread[i].spread  = 0.0;
        	g.dma_spread[i].quantity =  0;
	        g.dma_spread[i].label    = r1_species[i];
	}
	// We should move all of this g initialisation into a new class at some point, something like g_init()
	
	if (!g.failed)
	{
		for (vector<string>::size_type i = 0; i < nwchem_output.size(); i++)
		{
			vector<string> tokens;
			Tokenize(nwchem_output[i],tokens,"= \n\t");
			
			if ((tokens.size() > 3) && MO_line_number == 0)
			{
				// We need to find this line and get the number of electrons
				// Wavefunction type: closed shell. / spin polarized (one only)
				// Or such like as it appears before the MOs. Then round to nearest integer, as it is not explicitly defined
				// We only need to do this if number of electrons isn't defined, as it shouldn't change
				if (number_electrons == 0)
				{
					//if (cmpStr(tokens[0],"No.") &&
					//	cmpStr(tokens[1],"of") &&
					//	cmpStr(tokens[2],"electrons") &&
					//	cmpStr(tokens[3],":")) 
					if (cmpStr(tokens[0],"Wavefunction") &&
					    cmpStr(tokens[1],"type:"))
					{
                                                if (cmpStr(tokens[2],"spin") &&
                                                    cmpStr(tokens[3],"polarized."))
						{
							spin_polarised = 1;
							//cout << "Spin Polarised" << endl;
						}
						// Now lets gather all the other information in
						tokens.clear();
						Tokenize(nwchem_output[i+2],tokens,"= \n\t");

						// Save number of electrons. This seems to work fine
						StringToNumber(tokens[4],number_electrons);
						//cout << "Number of Electrons " << number_electrons << endl;

						// Lets force read the electrons in each channel here
						// First up Alpha electrons
						tokens.clear();
                                                Tokenize(nwchem_output[i+3],tokens,"= \n\t");
						StringToNumber(tokens[3],alpha_electrons);
						//cout << "Alpha Electrons " << alpha_electrons << endl;

						// And then Beta electrons
						tokens.clear();
						Tokenize(nwchem_output[i+4],tokens,"= \n\t");
                                                StringToNumber(tokens[3],beta_electrons);
						//cout << "Beta Electrons " << beta_electrons << endl;
					}
				}
				else if (spin_polarised == 0)
				{
					// Setting i to nwchem_output gives us an exit route
					// This is now disabled, and this process just returns i as is
					i = digest_rhf(tokens, &nwchem_output, &MO_line_number, 
								   &g, i, number_electrons);
				}
				else
				{
					// This is for unrestricted calculation
					if (spin_counter == 0)
					{
						spin_counter = digest_uhf(tokens, &nwchem_output, &MO_line_number, 
									   &g, i, alpha_electrons, spin_counter);

						//cout << spin_counter << " " << MO_line_number << endl;
	                                        if (spin_counter == 1)
        	                                {
                	                                MO_line_number = 0;
                        	                }
					}
					else
					{
						spin_counter = digest_uhf(tokens, &nwchem_output, &MO_line_number,
                                                                           &g, i, beta_electrons, spin_counter);
						
	                                        // If DMA is not needed this can be used
        	                                //else if (spin_counter == 2)
                	                        //{
                        	                //       i = nwchem_output.size();
                                	        //}
					}
				}
			} 
			// By checking the number of electrons we will only call this after
			// we have made a decent way through the output
			else if ((tokens.size() == 4 && number_electrons != 0))
			{
				// Now we are looking for distributed multipole analysis.
                                i = digest_dma(tokens, &nwchem_output, &DMA_line_number, &g, i);
			}
		}
		
		// MOs were not found and read in
		// We'll count this as failures
		if (MO_line_number == 0)
		{
			g.failed = true;
			cout << "Critical Failure: Getting MOs" << endl;
		}
		
		/* TO DO
		if (DMA_line_number == 0)
		{
			cout << "Non-critical Failure: Getting DMA" << endl;	
		}
		*/
	}
	
	// Print out all the electronic information
		
	if (!g.failed)
	{
		// So everything is calculated
		// Lets print our values and check them
		cout << "MO values: a.u. (eV)" << endl;
		cout << "HOMO Value: \t\t" << g.HOMO_value << "\t (" << g.HOMO_value*hartree_to_eV << ")" << endl;
		cout << "LUMO Value: \t\t" << g.LUMO_value << "\t (" << g.LUMO_value*hartree_to_eV << ")" << endl;
		cout << g.orbital_spread[0].label << " orbital average: \t" << g.orbital_spread[0].average << " \t (" ;
                cout << g.orbital_spread[0].average*hartree_to_eV << ") for " << r1_anions << " " << g.orbital_spread[0].label << " in region 1 " << endl;
		cout << g.orbital_spread[0].label << " orbital  spread: \t" << g.orbital_spread[0].spread  << " \t (" ;
                cout << g.orbital_spread[0].spread*hartree_to_eV << ")" << endl;
		cout << endl;

		/* TO DO
		if (DMA_line_number != 0)
		{
			cout << "**TESTING**" << endl;
			for (vector<min_max_spread>::size_type i = 0; i < g.dma_spread.size(); i++)
        		{
				cout << "DMA Spread for : " << g.dma_spread[i].label << " is " << g.dma_spread[i].spread;
				cout << "   \t( Min : " << g.dma_spread[i].min << ", Max : " << g.dma_spread[i].max << ", Quantity : " << g.dma_spread[i].quantity << " )" << endl;
			}
		}
		else
		{
			cout << "Non-critical Failure: DMA not present" << endl;
		}
		cout << endl;
		*/
	}
	else
		// Print out error about calculation failing
	{
		cout << endl;
		cout << "Bad News:" << endl;
		cout << "Looks like the calculation failed. No values have been saved" << endl;
		cout << endl;
		// Set g.function to an arbitary value
		g.function = 999999;
	}
	
	
	return g;
}

/*
 Method designed to take apart the restricted hartree-fock output from NWChem 
 
 @param[in] tokens Current line contents
 @param[in/out] nwchem_output Vector of output content
 @param[in/out] MO_line_number The position in the file where the MOs start
 @param[in/out] g Updating gaussin ECP
 @param[in] i Current index in the nwchem_output
 @param[in] en Number of electrons
 
 @return int Index for continued search
 */
int Nwchem::digest_rhf(vector<string> tokens, vector<string> *nwchem_output, int *MO_line_number, 
						  gaussian *g, int i, int en)
{
	/** This is the proceeding line from the outputs before we get the MOs

                       DFT Final Molecular Orbital Analysis
                       ------------------------------------

 Vector    1  Occ=2.000000D+00  E=-1.862267D+01
              MO Center=  2.8D+00,  2.8D+00,  8.4D+00, r^2= 5.7D-01
   Bfn.  Coefficient  Atom+Function         Bfn.  Coefficient  Atom+Function
  ----- ------------  ---------------      ----- ------------  ---------------
   344      0.636973   13 O  s              343      0.433767   13 O  s

 Vector    2  Occ=2.000000D+00  E=-1.862264D+01
              MO Center=  2.7D+00,  2.8D+00, -2.6D+00, r^2= 2.5D+00
   Bfn.  Coefficient  Atom+Function         Bfn.  Coefficient  Atom+Function
  ----- ------------  ---------------      ----- ------------  ---------------
    64      0.627347    3 O  s               63      0.427213    3 O  s


	 So we need to characterise this in the output **/
	
	// This will find the start of the MOs
	if (cmpStr(tokens[0],"DFT") &&
		cmpStr(tokens[1],"Final") &&
		cmpStr(tokens[2],"Molecular") &&
		cmpStr(tokens[3],"Orbital") &&
		cmpStr(tokens[4],"Analysis"))
	{
		*MO_line_number = i + 3; // Skip forward to where MOs start
  		double value = 0.0;
		vector<double> MO_energy;
		int last_MO_line_number = 0;
		
		// Different to Gamess-UK, these values aren't line after line.
		// So let's get all the values into an new vector, and read them out from there.
                for (vector<string>::size_type j = *MO_line_number; j < nwchem_output->size(); j++)
                {
			// cout << j << endl;
			tokens.clear();
			Tokenize(nwchem_output->at(j),tokens,"= \n\t");

                        // cout << nwchem_output->at(j) << " " << tokens.size() << endl;
			
			if (tokens.size() > 4)
			{
				if (cmpStr(tokens[0],"Vector"))
				{
                	                // cout << "in" << endl;
	                                last_MO_line_number = (int) j;
					// cout << tokens[5] << endl;
					replace(tokens[5].begin(),tokens[5].end(),'D','E');
					//cout << tokens[5] << endl;
					StringToNumber(tokens[5],value);
					MO_energy.push_back(value);
	
					// cout << last_MO_line_number << "  " << MO_energy[MO_energy.size()-1] << endl;		
				}
			}
                        
			// cout << j << " " << nwchem_output->size() << endl;
		}
		
		// Firstly lets work out the spread of Anion Orbitals 
		g->orbital_spread[0].min = MO_energy[region_1_anion_offset];
		
		// Find Anion maximum
                g->orbital_spread[0].max = MO_energy[region_1_anion_offset+(g->orbital_spread[0].quantity)-1];
                
		//cout << g->orbital_spread[0].min << " " << g->orbital_spread[0].max << " " << region_1_anion_offset << endl;
		
		g->orbital_spread[0].spread = g->orbital_spread[0].max - g->orbital_spread[0].min;
		
		// Lets work out the average of the Anion Orbital 
		for (int a = region_1_anion_offset; a < (g->orbital_spread[0].quantity)+region_1_anion_offset; a++)
		{
			//cout << tokens[2] << " " << value << " " << r1_anions << endl;
			g->orbital_spread[0].average += MO_energy[a];
		}
		g->orbital_spread[0].average /= g->orbital_spread[0].quantity;
		
		// Now lets set the HOMO and LUMO values
		g->HOMO_value = MO_energy[(en/2)-1];
		g->LUMO_value = MO_energy[(en/2)];
		// And we're done! Now to print these an see if they are right....
		// We'll do that outside the loop

		// Prior to DMA this was used. The use of DMA doubles the output_read process. Hey ho.	
		// This will break us out of the loop now we have the information we want
		// return nwchem_output->size();

		return last_MO_line_number;
	}
	
	// No luck finding results, return i
	return i;
}

/*
 Method designed to take apart the unrestricted hartree-fock output from NWChem 
 
 @param[in] tokens Current line contents
 @param[in/out] nwchem_output Vector of output content
 @param[in/out] MO_line_number The position in the file where the MOs start
 @param[in/out] g Updating gaussin ECP
 @param[in] i Current index in the nwchem_output
 @param[in] en Number of electrons in this channel
 @param[in] spin_counter Check if we are doing alpha or beta
 
 @return int Index for continued search
 */
int Nwchem::digest_uhf(vector<string> tokens, vector<string> *nwchem_output, int *MO_line_number, 
						  gaussian *g, int i, int en, int spin_counter)
{
	// We should be able to recalculate this from the input files.
	
	/** This is the proceeding line from the outputs before we get the MOs
                    DFT Final Alpha Molecular Orbital Analysis
                    ------------------------------------------

 Vector   19  Occ=1.000000D+00  E=-1.873084D+01
              MO Center=  2.8D+00,  2.8D+00,  2.8D+00, r^2= 1.5D-02
   Bfn.  Coefficient  Atom+Function         Bfn.  Coefficient  Atom+Function
  ----- ------------  ---------------      ----- ------------  ---------------
   120      0.639008    5 O  s              119      0.435300    5 O  s
	
	The header for the Beta orbitals is:
                     DFT Final Beta Molecular Orbital Analysis
                     -----------------------------------------

 Vector   19  Occ=1.000000D+00  E=-1.873084D+01
              MO Center=  2.8D+00,  2.8D+00,  2.8D+00, r^2= 1.5D-02
   Bfn.  Coefficient  Atom+Function         Bfn.  Coefficient  Atom+Function
  ----- ------------  ---------------      ----- ------------  ---------------
   120      0.639008    5 O  s              119      0.435300    5 O  s

	**/
	
	// This will find the start of the MOs
	if (cmpStr(tokens[0],"DFT") &&
		cmpStr(tokens[1],"Final") &&
		cmpStr(tokens[3],"Molecular") &&
		cmpStr(tokens[4],"Orbital") &&
		cmpStr(tokens[5],"Analysis"))
	{
                *MO_line_number = i + 3; // Skip forward to where MOs start
                double value = 0.0;
                vector<double> MO_energy;

                // Different to Gamess-UK, these values aren't line after line.
                // So let's get all the values into an new vector, and read them out from there.
                for (vector<string>::size_type j = *MO_line_number; j < nwchem_output->size(); j++)
                {
                        //cout << j << endl;
                        tokens.clear();
                        Tokenize(nwchem_output->at(j),tokens,"= \n\t");

                        if (cmpStr(tokens[0],"Vector"))
                        {
				//cout << tokens[5] << endl;
                                replace(tokens[5].begin(),tokens[5].end(),'D','E');
                                //cout << tokens[5] << endl;
                                StringToNumber(tokens[5],value);
                                MO_energy.push_back(value);

                                //cout << MO_line_number << "  " << MO_energy[MO_energy.size()-1] << endl;         
                        }
			else
			{
				// We need some check for when we reach the Beta MOs
				// Let's use the same check as for entering this loop

				if (cmpStr(tokens[0],"DFT") &&
			                cmpStr(tokens[1],"Final") &&
			                cmpStr(tokens[3],"Molecular") &&
			                cmpStr(tokens[4],"Orbital") &&
			                cmpStr(tokens[5],"Analysis"))
			        {	
					j = nwchem_output->size();
				}
			}
                }

                // Getting minimum value of Anion Orbitals
                // Check if this is the first or second loop through
                // And get spread inclusive of alpha and beta values
		if (spin_counter == 0)
		{
                	g->orbital_spread[0].min = MO_energy[region_1_anion_offset];
		}
		else
		{
			// Compare to previous value
			if (g->orbital_spread[0].min > MO_energy[region_1_anion_offset])
			{
				// Save if lower
				g->orbital_spread[0].min = MO_energy[region_1_anion_offset];

			}
		}

                // Getting minimum value of Anion S Orbital 
                // Check if this is the first or second loop through
                // And get spread inclusive of alpha and beta values
                if (spin_counter == 0)
                {
                        g->orbital_spread[0].max = MO_energy[region_1_anion_offset+(g->orbital_spread[0].quantity)-1];
                }
                else
                {
                        // Compare to previous value
                        if (g->orbital_spread[0].max < MO_energy[region_1_anion_offset+(g->orbital_spread[0].quantity)-1])
                        {
                                // Save if higher
                                g->orbital_spread[0].max = MO_energy[region_1_anion_offset+(g->orbital_spread[0].quantity)-1];
                        }
                }

                g->orbital_spread[0].spread = g->orbital_spread[0].max - g->orbital_spread[0].min;

                // Lets work out the average of the Anion Orbital 
                for (int a = region_1_anion_offset; a < (g->orbital_spread[0].quantity)+region_1_anion_offset; a++)
                {
                        g->orbital_spread[0].average += MO_energy[a];
                }

                // Check spin_counter here to set anion_average
                if (spin_counter == 1)
                {
                        // Divide through by number of alpha and beta electrons
                        g->orbital_spread[0].average /= (g->orbital_spread[0].quantity*2);
			//cout << g->orbital_spread[0].average << endl;
                }

		// Getting LUMO value
		// Check if this is the first or second loop through
		// We want the lowest value we can find
		if (spin_counter == 0)
		{
			g->LUMO_value = MO_energy[en];
		}
		else
		{
			// Compare to previous value
			if (g->LUMO_value > MO_energy[en])
			{
				// Save if higher
				g->LUMO_value = MO_energy[en];
			}
		}
		
		// Getting HOMO value
		// Check if this is the first or second loop through
		// We want the highest value we can find
		if (spin_counter == 0)
		{
			g->HOMO_value = MO_energy[en-1];
		}
		else
		{
			// Compare to previous value
			if (g->HOMO_value < MO_energy[en-1])
			{
				// Save if higher
				g->HOMO_value = MO_energy[en-1];
			}
		}
		
		// And we're done! Now to print these an see if they are right....
		// We'll do that outside the loop
		
		// We need this to check if we have been through twice
		// And picked alpha and beta values
		
		// Increment spin_counter
		spin_counter++;
	}
	
	// Haven't found all the results, return spin_counter
	return spin_counter;
}

/*
 Method designed to take aprt the DMA output from NWChem 
 
 @param[in] tokens Current line contents
 @param[in/out] nwchem_output Vector of output content
 @param[in/out] DMA_line_number The position in the file where the DMAs start
 @param[in/out] g Updating gaussin ECP
 @param[in] i Current index in the nwchem_output
 
 @return int Index for continued search
 */
int Nwchem::digest_dma(vector<string> tokens, vector<string> *nwchem_output, int *DMA_line_number,
                                                  gaussian *g, int i)
{
	// This needs testing for a UHF system

	// Now we are looking for distributed multipole analysis.
        // Keyphrase: distributed multipole analysis module
        if (cmpStr(tokens[0],"distributed") &&
            cmpStr(tokens[1],"multipole") &&
            cmpStr(tokens[2],"analysis") &&
            cmpStr(tokens[3],"module"))
        {

		// cout << "FOUND DMA" << endl;
 
                *DMA_line_number = i + 16; // Skip forward to where MOs start
		int DMA_step = 0;
		double value = 0.0;

		tokens.clear();
                Tokenize(nwchem_output->at(*DMA_line_number),tokens," \n\t");

		// cout << nwchem_output->at(*DMA_line_number) << endl;

		while ( cmpStr(tokens[4].substr(tokens[4].length()-1),"1") )
		{
			
/**
Here is the setup I used and the values we are trying to fill:
        g.dma_spread.resize(r1_species.size());
        for (vector<min_max_spread>::size_type i = 0; i < g.dma_spread.size(); i++)
        {
                g.dma_spread[i].max = -1000.0;
                g.dma_spread[i].min =  1000.0;
                g.dma_spread[i].average = 0.0;
                g.dma_spread[i].spread  = 0.0;
                g.dma_spread[i].quantity =  0;
                g.dma_spread[i].label    = r1_species[i];
        }
**/			
			int int_match = -1;

			// Find the matching species
			for (vector<min_max_spread>::size_type j = 0; j < g->dma_spread.size(); j++)
			{
				bool match = true;

				if (tokens[4].length() == g->dma_spread[j].label.length()+1)
                               	{
                               		for (size_t ichar = 0; ichar < g->dma_spread[j].label.length(); ichar++ )
                                       	{
                                       		if ( toupper(tokens[4][ichar]) != g->dma_spread[j].label[ichar])
                                               	{
                                                       	match = false;
                                               	}
                                       	}

                                       	if (match)
                                       	{
                                 	 	int_match = j;
						g->dma_spread[j].quantity++;
						// Quick exit?
						j = g->dma_spread.size();
                                        }
                                }	
			}

			// cout << DMA_step << " " << int_match << " " << tokens[4];
	
			tokens.clear();
			DMA_step++;
			Tokenize(nwchem_output->at(*DMA_line_number+DMA_step),tokens," \n\t");

			// Find modulus of dipolar expansion
                        while ( !cmpStr(tokens[0],"|q2|") )
                        {
                                tokens.clear();
                                DMA_step++;
                                Tokenize(nwchem_output->at(*DMA_line_number+DMA_step),tokens," \n\t");

				// Catch scenario where we have a blank line
                                while ( tokens.size() == 0 )
                                {
                                        tokens.clear();
                                        DMA_step++;
                                        Tokenize(nwchem_output->at(*DMA_line_number+DMA_step),tokens," \n\t");
                                }
                        }

			// cout << " : " << tokens[0] << " " << tokens[1] << " " << " " << tokens[2] << endl;

			// Get DMA value, and store information for later analysis
		        StringToNumber(tokens[2],value);

			g->dma_spread[int_match].average += value;
			
			if ( value > g->dma_spread[int_match].max )	
			{
				g->dma_spread[int_match].max = value;
			}
			
			if ( value < g->dma_spread[int_match].min )
			{
				g->dma_spread[int_match].min = value;
			}

			// Find next start point
			while ( !cmpStr(tokens[0],"site") )
                	{
				tokens.clear();
				DMA_step++;
				Tokenize(nwchem_output->at(*DMA_line_number+DMA_step),tokens," \n\t");

				// Catch scenario where we have a blank line
				while ( tokens.size() == 0 )
				{
					tokens.clear();
					DMA_step++;
					Tokenize(nwchem_output->at(*DMA_line_number+DMA_step),tokens," \n\t");
				}
			}				
		}

		for (vector<min_max_spread>::size_type j = 0; j < g->dma_spread.size(); j++)
                {
			g->dma_spread[j].spread = g->dma_spread[j].max - g->dma_spread[j].min;
			g->dma_spread[j].average /= g->dma_spread[j].quantity;
                }
                
		// And we're done! Now to print these an see if they are right....
                // We'll do that outside the loop

                // This will break us out of the loop now we have the information we want
                return nwchem_output->size();
        }

        // No luck finding results, return i
        return i;
}
