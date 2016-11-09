/*
 *  @file Functions.cpp
 *  fit_my_ecp
 *
 *  Created by Andrew Logsdail on 01/06/2012.
 *  Copyright 2012 University of Birmingham. All rights reserved.
 *
 */

#include "Functions.h"

using namespace std;

/*
 Constructor
 
 No params
 */
Functions::Functions()
{
	// We are going to put in the Functions weights here
	// These will be put to defaults in the appropriate region
	// Initiate to -1 as 0 is acceptable
	weight_gnorm_region1 = -1;
	weight_gnorm_region2 = -1;
	weight_gnorm_region3 = -1;
	weight_anion_spread = -1;
	weight_eigenvalues = -1;
	weight_dma_spread = -1;	
	// Target energies
	// target_homo = 0;	
        // target_lumo = 0;
	// Other targets
	// target_gnorm_region1 = 0;
        // target_gnorm_region2 = 0;
	// target_gnorm_region3 = 0;
	// target_anion_spread = 0;
	// target_dma_spread = 0;
	// Counter keeps track of the length of the arrays
	targets_length = 0;
}

/*
 *  Set parameters for calculating the RMS function
 *   
 *  @param[in] tgr1 Target for gradient in region 1
 *  @param[in] tgr2 Target for gradient in region 2
 *  @param[in] tgr3 Target for gradient in region 3
 *  @param[in] tgas Target for anion spread
 *  @param[in] tghomo Target for the HOMO
 *  @param[in] tglumo Target for the LUMO
 *  @param[in] tgdma Target for the Distributed Multipole Analysis
 *            
 */
void Functions::set_targets(vector<double> tgr1, vector<double> tgr2, vector<double> tgr3, vector<double> tgas,
			    vector<double> tghomo, vector<double> tglumo, vector<double> tgdma)
{
	target_gnorm_region1 = tgr1;
	target_gnorm_region2 = tgr2;
	target_gnorm_region3 = tgr3;
	target_anion_spread = tgas;
	target_homo = tghomo;
	target_lumo = tglumo;
	target_dma_spread = tgdma;

	// Let's make sure they are all the same sizes here. Prevents any segfaults further on
	// All unassigned values will be assigned to zero
	
	// Firstly, what's the longest vector?
	// targets_length = target_gnorm_region1.size();
	// This is now fixed against the number of punch template files

	// Compare this value to other vectors
	// if (target_gnorm_region2.size() > targets_length)
	// {
	//	targets_length = target_gnorm_region2.size();
	//}

        // if (target_gnorm_region3.size() > targets_length)
        // {
	//        targets_length = target_gnorm_region3.size();
        // }	

        //if (target_anion_spread.size() > targets_length)
        //{
        //        targets_length = target_anion_spread.size();
        //}

        //if (target_homo.size() > targets_length)
        //{
        //        targets_length = target_homo.size();
        //}

        //if (target_lumo.size() > targets_length)
        //{
        //        targets_length = target_lumo.size();
        //}

        //if (target_dma_spread.size() > targets_length)
        //{
        //        targets_length = target_dma_spread.size();
        //}

	// Now length is a maximum value. Resize all other vectors
        target_gnorm_region1.resize(targets_length,0.0);
        target_gnorm_region2.resize(targets_length,0.0);
        target_gnorm_region3.resize(targets_length,0.0);
        target_anion_spread.resize(targets_length,0.0);
        target_homo.resize(targets_length,0.0);
        target_lumo.resize(targets_length,0.0);
        target_dma_spread.resize(targets_length,0.0);
}

/*
 Set parameters for calculating the RMS function
 
 @param[in] wgr1 Weight for gradient in region 1
 @param[in] wgr2 Weight for gradient in region 2
 @param[in] wgr3 Weight for gradient in region 3
 @param[in] wOs Weight for anion spread
 @param[in] weigens Weight for the Eigenvalues
 @param[in] wdma Weight on the Distributed Multipole Analysis
 */
void Functions::set_weights(double wgr1, double wgr2, double wgr3, double was, double weigens, double wdma)
{
	// These set weights for linear Functions
	if (wgr1 == -1)
	{
		cout << "Using default weight of Region 1 Gnorm Average in Functions: 70.0" << endl;
		weight_gnorm_region1 = 70.0;
	}
        else
        {
		weight_gnorm_region1 = wgr1;
        }

	if (wgr2 == -1)
	{
		cout << "Using default weight of Region 2 Gnorm Average in Functions: 20.0" << endl;
		weight_gnorm_region2 = 20.0;
	}
	else
        {
                weight_gnorm_region2 = wgr2;
        }

	if (wgr3 == -1)
	{
		cout << "Using default weight of Region 3 Gnorm Average in Functions: 0.1" << endl;
		weight_gnorm_region3 = 0.1;
	}
	else
        {
                weight_gnorm_region3 = wgr3;
        }

	if (was == -1)
	{
		cout << "Using default weight of Anion S-Orbital Spread in Functions: 1.0" << endl;
		weight_anion_spread = 1.0; 
	}
	else
        {
                weight_anion_spread = was;
        }

	// Give error if HOMO is not defined
	// And we'll set this weighting to zero
	//if (target_homo == 0)
	//{
	//	cout << "** Target energy for Highest Occupied Molecular Orbital is not defined **" << endl;
	//	cout << "Using default weight of difference between HOMO values: 0.0" << endl; 
	//	weight_homo = 0;
	//}
	//else 
	//{
		
	if (weigens == -1)
	{
		cout << "Using default weight of difference between Eigenvalues: 0.05" << endl;
		weight_eigenvalues = 0.05;
	}
	else
        {
                weight_eigenvalues = weigens;
	}

        // Give error if LUMO is not defined
        // And we'll set this weighting to zero
        //if (target_lumo == 0)
        //{
        //        cout << "** Target energy for Lowest Unoccupied Molecular Orbital is not defined **" << endl;
        //        cout << "Using default weight of difference between LUMO values: 0.0" << endl;
        //        weight_lumo = 0;
        //}
        //else
        //{
        //        if (wlumo == -1)
        //        {
        //                cout << "Using default weight of difference between LUMO values: 0.05" << endl;
        //                weight_lumo = 0.05;
        //        }
        //        else
        //        {
        //                weight_lumo = whomo;
        //        }
        //}
        //
        if (wdma == -1)
        {
                cout << "Using default weight of Atomic Spread of Distributed Multipole Analysis in Function: 0.1" << endl;
                weight_dma_spread = 0.2; 
        }
        else
        {
                weight_dma_spread = wdma;
        }
	
} 

/*
 Function to see if the weights are the same as taken from the history files

 @param[in] old_weights Vector containing the old weights in order
 */
bool Functions::check_weights(vector<double> old_weights)
{
        //cout << old_weights[0] << " 1 " << weight_gnorm_region1 << endl;
        if (old_weights[0] != weight_gnorm_region1)
        {
		return false;
        }
        
        //cout << old_weights[1] << " 2 " << weight_gnorm_region2 << endl;
        if (old_weights[1] != weight_gnorm_region2)
        {
		return false;
        }

        //cout << old_weights[2] << " 3 " << weight_gnorm_region3 << endl;
        if (old_weights[2] != weight_gnorm_region3)
        {
		return false;
        }

        //cout << old_weights[3] << " 4 " << weight_anion_spread << endl;
        if (old_weights[3] != weight_anion_spread)
	{
		return false;
	}

        //cout << old_weights[4] << " 5 " << weight_eigenvalues << endl;
        if (old_weights[4] != weight_eigenvalues)
	{
		return false;
	}

        //This is repeated for the 5th position (LUMO), and 6th is RANK

	//cout << old_weights[7] << " 8 " << weight_dma_spread << endl;
        if (old_weights[7] != weight_dma_spread)
	{
		return false;
	}

        //cout << "Finished OK" << endl;

	return true;
}

/*
 Function to see if the targets are the same as taken from the history files

 @param[in] old_targets Vector containing the old targets in order
 */
bool Functions::check_targets(vector<double> old_targets, int dataset)
{
	// Loop over both target sets to compare. This needs to be made
	// fire proof if the arrays are different sizes
	// Lets put an if statement, and just return false should they not match up.
	// All target arrays are the same length at this point

	// We have used an int as these arrays are different shapes
	//cout << old_targets[0] << " 1 " << target_gnorm_region1 << endl;
       	if (old_targets[0] != target_gnorm_region1[dataset])
        {
       	        return false;
        }

       	//cout << old_targets[1] << " 2 " << target_gnorm_region2 << endl;
        if (old_targets[1] != target_gnorm_region2[dataset])
       	{
               	return false;
        }

       	//cout << old_targets[2] << " 3 " << target_gnorm_region3 << endl;
        if (old_targets[2] != target_gnorm_region3[dataset])
       	{
               	return false;
        }

        //cout << old_targets[3] << " 4 " << target_anion_spread << endl;
       	if (old_targets[3] != target_anion_spread[dataset])
        {
       	        return false;
        }

       	//cout << old_targets[4] << " 5 " << target_homo << endl;
        if (old_targets[4] != target_homo[dataset])
       	{
               	return false;
        }

        //cout << old_targets[5] << " 6 " << target_lumo << endl;
       	if (old_targets[5] != target_lumo[dataset])
        {
		return false;
        }

        //cout << old_targets[6] << " 7 " << target_dma_spread << endl;
       	if (old_targets[6] != target_dma_spread[dataset])
       	{
               	return false;
        }

        //cout << "Finished OK" << endl;

        return true;
}

/*
 Function to calculate weighted fitness
 This is simply a weighted sum of squares
 
 @param[in] gnorm_1 Gradient in region 1
 @param[in] gnorm_2 Gradient in region 2
 @param[in] gnorm_3 Gradient in region 3
 @param[in] r1_anion The Anion S-Orbital spread
 @param[in] HOMO Energy of the HOMO orbital (Used tentatively to compare to IP)
 @param[in] LUMO Energy of the LUMO orbital (Used for tentative comparison to all electron calculations)
 @param[in] dma_spread Vector containing DMA Spreads
 @param[in] dataset Gives an index for comparing the target values to
 
 @return double Function value
 */
double Functions::calculate_linear_function(double gnorm_1, double gnorm_2, double gnorm_3, double r1_anion, double HOMO, double LUMO, vector<double> dma_spread, int dataset)
{
	// Calculated weighted Functions
	// Weights? We'll make these dynamic. Weights defined at top
	// Johns Functions deals with squares, but here we have square roots (i.e. true norm)
	// These are taken as is from the gradient directory
	// Added ability to catch invalid numbers
	
	gnorm_1 -= target_gnorm_region1[dataset];;
	double weighted_gnorm_region1 = 0;
	if (gnorm_1 != -1) weighted_gnorm_region1 = weight_gnorm_region1*(gnorm_1*gnorm_1);
	//weighted_gnorm_region1 *= weighted_gnorm_region1;
	
	gnorm_2 -= target_gnorm_region2[dataset];
	double weighted_gnorm_region2 = 0;
	if (gnorm_2 != -1) weighted_gnorm_region2 = weight_gnorm_region2*(gnorm_2*gnorm_2);
	//weighted_gnorm_region2 *= weighted_gnorm_region2;
	
	gnorm_3 -= target_gnorm_region3[dataset];
	double weighted_gnorm_region3 = 0;
	if (gnorm_3 != -1) weighted_gnorm_region3 = weight_gnorm_region3*(gnorm_3*gnorm_3);
	//weighted_gnorm_region3 *= weighted_gnorm_region3;
	// r1_anion spread needs converting to eV
	r1_anion *= hartree_to_eV;
	r1_anion -= target_anion_spread[dataset];
	double weighted_anion_spread = weight_anion_spread*r1_anion*r1_anion;
	//weighted_anion_spread *= weighted_anion_spread;
	
	// Here we compare the eigenvalues. Sign IS important
 	double weighted_HOMO = HOMO;
	// Check a Target has been set
        if (target_homo[dataset] != 0)
        {
                // weighted_HOMO = sqrt(HOMO*HOMO)*hartree_to_eV;
                weighted_HOMO *= hartree_to_eV;
		// So we have absolute value of ionisation potential. Next calculate the difference from experiment
		weighted_HOMO -= target_homo[dataset];
		// cout << weighted_HOMO << endl;
		// And now we'll get this as an absolute value
		// Change of plan, we're squaring it to match John's Functions
       		// And now we just need to multiply by the weight
		weighted_HOMO *= weighted_HOMO * weight_eigenvalues;
        }
        else
        {
		weighted_HOMO = 0;
        }

        // Repeat for LUMO
        double weighted_LUMO = LUMO;
        if (target_lumo[dataset] != 0)
        {
                // weighted_LUMO = sqrt(LUMO*LUMO)*hartree_to_eV;
                weighted_LUMO *= hartree_to_eV;
	        // So we have absolute value of ionisation potential. Next calculate the difference from experiment
	        weighted_LUMO -= target_lumo[dataset];
		// cout << weighted_LUMO << endl;
	        // And now we'll get this as an absolute value
        	// Change of plan, we're squaring it to match John's Functions
	        // And now we just need to multiply by the weight
	        weighted_LUMO *= weighted_LUMO * weight_eigenvalues;
        }
	else
	{
		weighted_LUMO = 0;
	}

	//Now to calculate the weighted values of dma_spread
	double weighted_dma_spread = 0;
	//We'll sum these values in to one variable
	for (vector<double>::size_type i = 0; i < dma_spread.size(); i++)
	{
		dma_spread[i] -= target_dma_spread[dataset];
		weighted_dma_spread += weight_dma_spread*(dma_spread[i]*dma_spread[i]);
	}
	// This is actually a recalculation if the gaussian is taken from history
	// but hopefully that is not an expensive calculation
	// and it can always be disabled with an if loop
	return weighted_gnorm_region1 + weighted_gnorm_region2 + weighted_gnorm_region3 + weighted_anion_spread + weighted_HOMO + weighted_LUMO + weighted_dma_spread;
}

/*
 Method just to return the header for the output log file
 
 No params
 */
string Functions::get_header()
{
	const int precision = 4;
	string sentence = "";
        string temp = "";
	// sentence =  "#\tPre.\t\tA\t     zeta\t\t | \t Gnorm 1 (" + NumberToString(weight_gnorm_region1); //UPDATE 1
	// sentence =  "#\t\tA\t     zeta\t\t | \t Gnorm 1 (" + NumberToString(weight_gnorm_region1); // UPDATE 2
	sentence =  "Entry       | Gnorm.1 ";
        temp = NumberToString(weight_gnorm_region1);

        for (int i = temp.length(); i < precision; i++) {
                sentence += " ";
        }
        sentence += "(" + temp;

	sentence += ") Gnorm.2 ";
        temp = NumberToString(weight_gnorm_region2);

        for (int i = temp.length(); i < precision; i++) {
                sentence += " ";
        }
        sentence += "(" + temp;

	sentence += ") Gnorm.3 ";
        temp = NumberToString(weight_gnorm_region3);

        for (int i = temp.length(); i < precision; i++) {
                sentence += " ";
        }
        sentence += "(" + temp;

	sentence += ") Anion.Spread ";
        temp = NumberToString(weight_anion_spread);

        for (int i = temp.length(); i < precision; i++) {
                sentence += " ";
        }
        sentence += "(" + temp;

	sentence += ") HOMO ";
        temp = NumberToString(weight_eigenvalues);

        for (int i = temp.length(); i < precision; i++) {
                sentence += " ";
        }
        sentence += "(" + temp;

        sentence += ") LUMO "; 
        temp = NumberToString(weight_eigenvalues);

        for (int i = temp.length(); i < precision; i++) {
                sentence += " ";
        }
        sentence += "(" + temp;

        sentence += ")   |      Function (Rank) ";
        sentence += "    | DMA.Spread "; 
	temp = NumberToString(weight_dma_spread);

        for (int i = temp.length(); i < precision; i++) {
                sentence += " ";
        }
        sentence += "(" + temp;

        sentence += ")";
	return sentence;
}

/*
 * Method to return the targets in line with the weights above
 *
 * @param[in] dataset Which dataset are we creating a header for?
 */
//string Functions::get_targets_header()
string Functions::get_targets_header(int dataset)
{
	const int precision = 8;
	// Needs expanding
	
	string sentence = "";
        string temp = "";

	sentence += "Targets     |     ";

        temp = NumberToString(target_gnorm_region1[dataset],precision);
        
        for (int i = temp.length(); i <= precision; i++) {
		sentence += " ";
	}
	sentence += temp;

	sentence += "      ";
        temp = NumberToString(target_gnorm_region2[dataset],precision);

        for (int i = temp.length(); i <= precision; i++) {
                sentence += " ";
        }
        sentence += temp;

	sentence += "      ";
        temp = NumberToString(target_gnorm_region3[dataset],precision);

        for (int i = temp.length(); i <= precision; i++) {
                sentence += " ";
        }
        sentence += temp;

	sentence += "         ";
        temp = NumberToString(target_anion_spread[dataset],precision);
   
        for (int i = temp.length(); i <= precision; i++) {
                sentence += " ";
        }
        sentence += temp;

	sentence += "   ";
        temp = NumberToString(target_homo[dataset],precision);

        for (int i = temp.length(); i <= precision; i++) {
                sentence += " ";
        }
        sentence += temp;

	sentence += "   ";
        temp = NumberToString(target_lumo[dataset],precision);

        for (int i = temp.length(); i <= precision; i++) {
                sentence += " ";
        }
        sentence += temp;

	sentence += "    |      **************      |       "; 
        temp =  NumberToString(target_dma_spread[dataset],precision);

        for (int i = temp.length(); i <= precision; i++) {
                sentence += " ";
        }
        sentence += temp;

	return sentence;
}
