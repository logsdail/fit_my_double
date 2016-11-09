/*
 *  @Functions.h
 *  fit_my_ecp
 *
 *  @brief Class to calculate the RMS function value
 *
 *  Created by Andrew Logsdail on 01/06/2012.
 *  Copyright 2012 University College London. All rights reserved.
 *
 */

#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <iostream>
// Personal headers
#include "Utils.h"
#include "Structures.h"

class Functions {
	
public:
	
	Functions();
	
	/*
	 Deconstructor
	 
	 No params
	 */
	~Functions() {;}
	
	void set_weights(double wgr1, double wgr2, double wgr3, double was, double weigens, double wdma);
	// void set_parameters(double wgr1, double wgr2, double wgr3, double was, double weigens, double wdma, double homo, double lumo);
	
	void set_targets(std::vector<double> tgr1, std::vector<double> tgr2, std::vector<double> tgr3, 
                         std::vector<double> tgas, std::vector<double> tghomo, std::vector<double> tglumo, std::vector<double> tgdma);

        bool check_weights(std::vector<double> old_weights);

	bool check_targets(std::vector<double> old_targets, int dataset);
	
	/*
	 Calculate the function value and return it!
	 
	 @param[in] gnorm_1 Gradient in region 1
	 @param[in] gnorm_2 Gradient in region 2
	 @param[in] gnorm_3 Gradient in region 3
	 @param[in] r1_anion The Anion S-Orbital spread
	 @param[in] HOMO Energy of the HOMO orbital (Used tentatively to compare to IP
         @param[in] dma_spread Vector containing the values of the DMA spread
	 @return double Function value
	 */
	double calculate_function(double gnorm_1, double gnorm_2, double gnorm_3, double r1_anion, double HOMO, double LUMO, std::vector<double> dma_spread, int dataset)
	{
		return calculate_linear_function(gnorm_1, gnorm_2, gnorm_3, r1_anion, HOMO, LUMO, dma_spread, dataset);
	}
	
	std::string get_header();

	//std::string get_targets_header();
	std::string get_targets_header(int dataset);

	/*
 	 Return the length of the Target vectors

 	 @return targets_length Length of all targets vectors
	 */
	int get_targets_length()
	{
		return targets_length;
	}

	void set_targets_length(int i)
	{
		targets_length = i;
	}
	
private:
	
	// We are going to put in the Functionstion weights here
	double weight_gnorm_region1;
	double weight_gnorm_region2;
	double weight_gnorm_region3;
	double weight_anion_spread;
	double weight_eigenvalues;
	double weight_dma_spread;
	// Experimentally defined ionisation energy
	//double target_homo;
        //double target_lumo;
        // Other targets
        //double target_gnorm_region1;
        //double target_gnorm_region2;
        //double target_gnorm_region3;
        //double target_anion_spread;
        //double target_dma_spread;
        // Vector replacements...
        std::vector<double> target_homo;
	std::vector<double> target_lumo;
	std::vector<double> target_gnorm_region1;
	std::vector<double> target_gnorm_region2;
	std::vector<double> target_gnorm_region3;
	std::vector<double> target_anion_spread;
        std::vector<double> target_dma_spread;

	// Keep track of targets length
	unsigned int targets_length;
	
	double calculate_linear_function(double gnorm_1, double gnorm_2, double gnorm_3, double r1_anion, double HOMO, double LUMO, std::vector<double> dma_spread, int dataset);
};

#endif
