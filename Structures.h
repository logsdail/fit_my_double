/*
 *  @Structures.h
 *  fit_my_ecp
 *
 *  @brief Some customised structures for helping our cause.
 *  These are almost exclusively used in vectors to make things easier.
 *
 *  Created by Andrew Logsdail on 01/06/2012.
 *  Copyright 2012 University of Birmingham. All rights reserved.
 *
 */

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>
#include <string>

// Some important constants we need to define
#define hartree_to_eV 27.21138386
#define bohr_to_angstrom 0.52917721092
#define gradient_converter hartree_to_eV/bohr_to_angstrom
#define spacer "================================================================================================"

struct regions_data
{
	double gnorm;
	double gnorm_max;
	double gradx_max;
	double grady_max;
	double gradz_max;
};

struct gaussian_info
{
	int line_number;
	double value;
	int type;
};

struct min_max_spread
{
        double max;
        double min;
        double average;
        double spread;
        int quantity;
        std::string label;
};

struct gaussian
{	
	// unsigned int line_number;
	// std::string preamble;
	// std::vector<double> values; 
	std::vector<gaussian_info> values;
	std::vector<regions_data> regions;
        std::vector<min_max_spread> orbital_spread;
        std::vector<min_max_spread> dma_spread;
	double HOMO_value;
	double LUMO_value;
	//double anion_max;
	//double anion_min;
	//double anion_average;
	//double anion_spread;
	double function;
	int rank;
	int index;
	bool failed;
};

#endif
