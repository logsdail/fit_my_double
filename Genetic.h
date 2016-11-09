/*
 *  @file Genetic.h
 *  fit_my_ecp
 *
 *  @brief This class controls the genetic algorithm (GA) for fitting ECPs.
 *  Currently very inefficient as the problem size is too small
 *
 *  Created by Andrew Logsdail on 02/06/2012.
 *  Copyright 2012 University College London. All rights reserved.
 *
 *
 */


#ifndef GENETIC_H
#define GENETIC_H

#include "Powells.h"

class Genetic : public Powells {
	
public:
	
	/*
	 Constructor
	 
	 No params
	 */
	Genetic(){}
	
	Genetic(int *seed);
	
	/*
	 Deconstructor
	 
	 No params
	 */
	~Genetic(){}
	
	void set_ga_parameters(int ps,
						   int ms,
						   int os,
						   int cc,
						   bool md);
	
	/*
	 Method to print the search type being conducted, common with all other search methods
	 
	 No Param
	 */
	void print_type()
	{
		std::cout << "Performing a Genetic Algorithm search" << std::endl;
	}
	
protected:
	
	std::vector<gaussian> population;
	// std::vector<gaussian> mutants_offspring;
	
	int calculate_ecps_to_test_counter;
	unsigned int population_size;
	unsigned int mutations_size;
	unsigned int offspring_size;
	int convergence_criteria;
	bool mutation_dynamic;

	void calculate_ecps_to_test();

	void check_converged();	

	gaussian get_random();

	gaussian get_offspring();
	
	gaussian get_mutant();
	
	int get_population_worst_option();
	
	double fitness(const int pos);
};

#endif
