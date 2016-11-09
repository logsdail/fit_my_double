/*
 *  @History.h
 *  fit_my_ecp
 *
 *  @brief Class to store all the historical information
 *
 *  Created by Andrew Logsdail on 02/06/2012.
 *  Copyright 2012 University College London. All rights reserved.
 *
 */

#ifndef HISTORY_H
#define HISTORY_H

#include <vector>
// Personal headers
#include "Utils.h"
#include "Structures.h"
#include "Functions.h"

class History {
	
public:
	
	/*
	 Constructor
	 
	 No params
	 */
	History() {;}
	
	/*
	 Deconstructor
	 
	 No params
	 */
	~History() {;}
	
	void insert_old_data(std::vector<std::string> input, std::vector<std::string> regions_input);
	
	int check_history(gaussian g);
	
	/*
	 Return the size of ECP history
	 
	 @return int Size of array. We should get rid of this method
	 */
	int size()
	{
		return ecp_history.size();
	}
	
	/*
	 Returns the entire history vector
	 
	 @return Vector of type gaussian with all ECPs from history
	 */
	std::vector<gaussian> get_history()
	{
		return ecp_history;
	}

	 /*
         Returns the numer of entries last added to the history

         @return Integer counter
         */
        int get_number_of_entries_last_added()
        {
                return number_of_entries_last_added;
        }

        /*
         Sets the entire history vector

         @input Vector of type gaussian with ECPs for history
         */
        void set_history(std::vector<gaussian> v)
        {
		for (std::vector<gaussian>::size_type i = 0; i < v.size(); i++)
		{
			add_to_history(v[i]);
		}
        }
	
	/*
	 Method to get data from the history vector
	 
	 @param[in] i Index to be retrieved
	 @return gaussian ECP value at i
	 */
	gaussian get(int i)
	{
		return ecp_history[i];
	}

	/*
	 Add new data into history
	 
	 @param[in] v Vector containing all gaussians (Should make into a pointer)
	 @param[in] v2 Vector with 1's if the gaussian values are new
	 */
        void append(std::vector<gaussian> v, std::vector<int> v2)
	{
		number_of_entries_last_added = 0;

		//std::cout << "v size " << v.size() << std::endl;
		//std::cout << "v2 size " << v2.size() << std::endl;
		std::vector<gaussian> v_new;

		for (std::vector<int>::size_type a = 0; a < v2.size(); a++)
		{
			// std::cout << "v2 " << v2[a] << std::endl;
			// Check if this is a new addition to the history
			if (v2[a] == 0)
			{
				v_new.push_back(v[a]);

				number_of_entries_last_added++;
			}
		}

		// std::cout << "v_new size " << v_new.size() << std::endl;
                // std:: cout << "ecp_history size " << ecp_history.size() << std::endl;
		append(v_new);		
		// std:: cout << "ecp_history size " << ecp_history.size() << std::endl;
	}

	void recalc_function(Functions *func_calc, int dataset, bool force_recalc);

        /*
         Method to delete duds from the history array 
        */
	void remove_duds()
	{
		std::vector<gaussian>::size_type i = 0;
		while (i < ecp_history.size())
                {
        	       	if (ecp_history[i].function == 888888)
			{
				ecp_history.erase(ecp_history.begin()+i);
			}
			else
			{
				i++;
			}
	        }
	}
	
private:

	/*
	 Method to insert subarray in the history array
	 
	 @param[in] v Vector of new ECPs
	 */
	void append(std::vector<gaussian> v)
        {
              if (v.size() > 0)
              {
                      ecp_history.insert(ecp_history.end(), v.begin(), v.end());
              }
        }

	
	void add_to_history(gaussian g);
	
	std::vector< std::vector<gaussian_info> > decompose_regions(std::vector<std::string> regions_input);
	
	std::vector<gaussian> ecp_history;

	int number_of_entries_last_added;
};

#endif

