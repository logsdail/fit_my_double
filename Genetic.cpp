/*
 *  @file Genetic.cpp
 *  fit_my_ecp
 *
 *  Created by Andrew Logsdail on 02/06/2012.
 *  Copyright 2012 University College London. All rights reserved.
 *
 */

#include "Genetic.h"

using namespace std;

/*
 Constructor
 
 @param[in/out]	*seed Pointer to random number seed
 */
Genetic::Genetic(int *seed) : Powells(seed)
{
	// idum = seed;
	calculate_ecps_to_test_counter = 0;
	population_size = 0;
	mutations_size = 0;
	offspring_size = 0;
	convergence_criteria = 0;
	mutation_dynamic = false;
}

/*
 Method to set the initial GA parameters.
 Requirement if we are running a quick GA search.
 
 @param[in] ps The size of the GA population
 @param[in] ms The size of the mutant population
 @param[in] os The size of the offspring population
 @param[in] cc Convergence criteria to terminate a GA search
 @param[in] md Check for dynamic mutation
 */
void Genetic::set_ga_parameters(int ps,
				int ms,
				int os,
				int cc,
				bool md)
{
	population_size = ps;
	mutations_size = ms;
	offspring_size = os;
	convergence_criteria = cc;
	mutation_dynamic = md;
	
	// - population size
	if (population_size == 0)
	{
		cout << "Using default population size for GA : 4" << endl;
		population_size = 4;
	}
	// - mutations size
	if (mutations_size == 0)
	{
		cout << "Using default mutations size for GA : 2" << endl;
		mutations_size = 2;
	}
	// - offspring size
	if (offspring_size == 0)
	{
		cout << "Using default offspring size for GA : 2" << endl;
		offspring_size = 2;
	}
	// - convergence criteria
	if (convergence_criteria == 0)
	{
		cout << "Using default convergence criteria for GA : 3" << endl;
		convergence_criteria = 3;
	}
	// - dynamic mutation
	if (!mutation_dynamic)
	{
		cout << "Using default mutation type for GA : Static" << endl;
	}
}

/*
 Method to calculate next ECPs.
 Internal workings call other internal functions
 */
void Genetic::calculate_ecps_to_test()
{
	cout << "Calculating ECPs to Test" << endl;
	ecps_to_test.clear();
	// First off insert random numbers 
	if (calculate_ecps_to_test_counter == 0)
	{
		cout << "Creating a random population of size " << population_size << endl;
		// Fill the population with random values
		while (ecps_to_test.size() < population_size)
		{
			// Insert random value
			while (ecps_to_test.size() < population_size)
			{
				ecps_to_test.push_back(get_random());
			}
			
			// This double loop means we can remove anything we don't want
			check_boundaries_duplicates();
		}
	}
	else 
	{
		if (calculate_ecps_to_test_counter%2 == 1)
		{
			// We've calculated the population.
			population = ecps_tested;
			// Calculate the mutants and offspring
			cout << "Calculating a new set of offspring and mutants" << endl;
  			// cout << mutants_offspring.size() << " " << mutations_size + offspring_size << endl;
			// We need an exit strategy if we get stuck in one of these loops
			int loop_counter = 0;
			// Enter while loop
			while (ecps_to_test.size() < offspring_size)
			{
				// cout << ecps_to_test.size() << " " << offspring_size << endl;
				// This double loop means we can check for duplicates at the end, which is computationally more efficient
				while (ecps_to_test.size() < offspring_size)
				{
					// Firstly fill up with the offspring we need
					// cout << "Creating Offspring " << ecps_to_test.size()+1 << " of " << offspring_size << endl;
					ecps_to_test.push_back(get_offspring());
				}
				// And check boundary conditions
				// Give it ten attempts otherwise just accept that some might be the same as we've approached the minima.
				if (loop_counter < 10 )
				{
					check_boundaries_duplicates();
				}
				loop_counter++;
			}
			// mutants_offspring.clear();
			// Use out exit strategy again
			loop_counter = 0;
                        while (ecps_to_test.size() < (mutations_size + offspring_size))
                        {
                                // cout << ecps_to_test.size() << " " << mutations_size + offspring_size << endl;
                                // This double loop means we can check for duplicates at the end, which is computationally more efficient
                                while (ecps_to_test.size() < (mutations_size + offspring_size))
                                {
                                        // cout << "Creating Mutant " << ecps_to_test.size()+1-offspring_size << " of " << mutations_size << endl;
                                        // Fill up with mutants of the population
                                        ecps_to_test.push_back(get_mutant());
                                }
                                // And check boundary conditions
				// Give it ten attempts otherwise just accept that some might be the same as we've approached the minima.
                                if (loop_counter < 10 )
                                {
                                        check_boundaries_duplicates();
                                }
                                loop_counter++;
                        }
		}
		else
		{
			cout << "Creating new population of size " << population_size << endl;
			// We've just calculated the mutants and offspring
			// Form a new population from all the values we have
			// Let's put them all together
			// population.insert(population.end(),ecps_tested.begin(),ecps_tested.end());
			ecps_to_test = population;
			ecps_to_test.insert(ecps_to_test.end(),ecps_tested.begin(),ecps_tested.end());
			population.clear();
			// All together. Let's remove duplicates
			check_boundaries_duplicates();			
			// Shrink the population down
			// cout << ecps_to_test.size() << endl;
			while (ecps_to_test.size() > population_size)
			{
				// Remove the worst option
				int i = get_population_worst_option();
				// cout << "Removing member: " << i << endl;
				ecps_to_test.erase(ecps_to_test.begin()+i);
			} 
			// And we are done
		}
	}

	calculate_ecps_to_test_counter++;
}


/*
 Method to check convergence
 */
void Genetic::check_converged()
{
	// Check if we are calculating mutants and offspring, or population
	// Compare results to see if we've converged

	// We must only check these results when we have a population in ecps_tested
	if (calculate_ecps_to_test_counter%2 == 1)
	{
		cout << endl;
		cout << "Checking convergence for population : " << calculate_ecps_to_test_counter/2 << endl;
		cout << "Current loops without change : " << minimisation_count << ", Convergence criteria : " << convergence_criteria << endl;
		cout << endl;

		if (compare_ecps(ecps_tested[number_one_ranked],previous_minimum))
		{
			// Number of times we have been through the minima without change
			minimisation_count++;
		
			// We multiply convergence by two
			// As we do a loop for the population and the mutants/offspring
			if (minimisation_count == convergence_criteria)
			{
				converged = true;
			}
		}
		else
		{
			converged = false;
			previous_minimum = ecps_tested[number_one_ranked];
			minimisation_count = 0;
		}

	}
	
	// Calculate new ecps if we are not converged
	if (!converged)
	{
		calculate_ecps_to_test();
		
	}
}	

/*
 Returns random ECP within search limits
 
 @return Random ECP
 */
gaussian Genetic::get_random()
{
	// Get a template
	gaussian g = previous_minimum;
	
	// Dynamic initialisation
	for (vector<gaussian_info>::size_type i = 0; i < g.values.size(); i++)
	{
		g.values[i].value = minimums[g.values[i].type] + randomNumberF((maximums[g.values[i].type]-minimums[g.values[i].type]), idum);
	}
	
	return g;
}

/*
 Figures out the worst member of the current population
 We should couple this with finding the best in the long term
 As finding the best needs moving out of the Main.cpp
 
 @param int Index of worst ranked option in population
 */
int Genetic::get_population_worst_option()
{
	int worst_option = 0;
	/** LOOP OVER THE VECTOR AND RETURN WORST OPTIONS **/
	for (vector<gaussian>::size_type i = 1; i < ecps_to_test.size(); i++)
	{
		if (ecps_to_test[i].function > ecps_to_test[worst_option].function)
		{
			worst_option = i;
		}
	}
	// Return the position of the highest function
	return worst_option;
}

/*
 Generate a new offpsring member
 Bases decision making on members of current population
 
 @return gaussian Offspring ECP
 */
gaussian Genetic::get_offspring()
{
	int i = 0;
	// This is is essentially roulette selection
	while (randomNumber(idum) > fitness(i))
	{
		i = randomNumber(population_size,idum);
	}
	
	int j = 0;
	// Select a second parent. Roulette selection
	// Make sure it isn't the same as the other parent!
	while ((randomNumber(idum) > fitness(j)) && (j != i))
	{
		j = randomNumber(population_size,idum);
	}
	
	// We'll do uniform crossover
	gaussian g = population[i];
	gaussian second = population[j];

	if (g.values.size() < 3)
	{
		// The number of values is small. We'll force a mutation
                i = randomNumber(g.values.size(), idum);
                // Picked a value, now copy in from second parent
                g.values[i] = second.values[i];
		// Set function to big value
        	g.function = 999999;
	}
	else
	{
		// Formal method but doesn't work well with small number of values
		vector<double>::size_type crossover_counter = 0;
	 
		for (vector<double>::size_type k = 0; k < g.values.size(); k++)
		{
			// If random number is greater than 0.5, copy in value from second
			if (randomNumber(idum) > 0.5)
			{
				g.values[k] = second.values[k];
				// Set function to big value
				g.function = 999999;
				// Crossover counter
				crossover_counter++;
			}
		}
	
		// If we've copied everything across, just return the second gaussian as that saves a duplicate calculation
		if (crossover_counter == g.values.size())
		{
			return second;
		}
		// Else return the updated gaussian
		else 
		{
			return g;
		} 
        }
	
	// Fireproof escape
	return g;
}

/*
 Generates a mutant population member
 Bases decision on current population
 
 @return gaussian Mutant ECP
 */
gaussian Genetic::get_mutant()
{
	int i = 0;
	// This is is essentially roulette selection
	while (randomNumber(idum) > fitness(i))
	{
		i = randomNumber(population_size,idum);
	}
	
	// Get some gaussian values
	gaussian g = population[i];
	gaussian random = get_random();
	
	// We are going to mutate in just one direction
	i = randomNumber(g.values.size(), idum);

	// Dynamic mutation if we are not seeing much variation
	if (mutation_dynamic)
	{
		// Get a weighted new value close to previous value
		g.values[i].value *= minimisation_count;
		g.values[i].value += random.values[i].value;
		g.values[i].value /= (minimisation_count + 1);
	}
	else
	{
		// Picked a value, now copy in random variable
		g.values[i] = random.values[i];
	}
	// Set function and values to zero
	g.function = 999999;
	
	return g;
}


/*
 Calculates the fitness of the chosen population member, for mating/mutation purposes
 Lifted from zContrast
 
 @param[in] pos Position in population
 @return double Fitness relative to other members of the population
 */
double Genetic::fitness(const int pos)
// Work out fitness of a value compared to overall population
// Input: Current Selection of Points
// Output : double of fitness value
{
	int worst = get_population_worst_option(); // MAX
	int best = number_one_ranked; // MIN
	double p = ((population[pos].function - population[best].function)/(population[worst].function - population[best].function)); // Normalise
	// EXPONENTIAL FACTOR //
	const double alpha = 3;
	// Should we soft code this? Probably //
	////////////////////////
	
	// if (bFLinear) 
	// {
	// 	return (1-(0.7*p)); // Linear function
	// }
	//else if (bFExponential) 
	//{
		return exp(-alpha*p); // Exponential function
	//}
	//else if (bFTanh) 
	//{
	//	return 0.5*(1-tanh((2*p)-1)); // Hyperbolic Tangent function
	//}
	//else 
	//{
	//	return p;
	//}
}
