/*
 *  @file io.cpp
 *  fit_my_ecp
 *
 *  Created by Andrew Logsdail on 01/06/2012.
 *  Copyright 2012 University of Birmingham. All rights reserved.
 *
 */

#include "IO.h"

using namespace std;

/*
 This is a function to update the log file
 
 @param[in] output Filename
 @param[in] v Vector of gaussians, containing all the ECPs to be written
 @param[in] critical Catch for critical errors to terminate program
 */
void update_log_file(string output, vector<gaussian> v, bool critical)
{
	ofstream o;
        const int precision = 9;
	string temp = "";
	double temp_double = 0;
	o.open( output.c_str(), ios::app );
	
	if (o)
	{
		for (vector<gaussian>::size_type i = 0; i != v.size(); i++)
		{
			o.precision(precision);
			temp = NumberToString(v[i].index,precision);
			o << temp;

                	for (int j = temp.length(); j <= precision; j++) {
                        	o << " ";
        		}
                        
			o << "  |    ";
                        //o << setw(12) << v[i].regions[0].gnorm;
                        temp = NumberToString(v[i].regions[0].gnorm);

                        for (int j = temp.length(); j <= precision; j++) {
                                o << " ";
                        }

                        o << temp << "     ";
			//o << setw(12) << v[i].regions[1].gnorm;
                        temp = NumberToString(v[i].regions[1].gnorm);

                        for (int j = temp.length(); j <= precision; j++) {
                                o << " ";
                        }

                        o << temp << "     ";
                        //o << setw(13) << v[i].regions[2].gnorm;
                        temp = NumberToString(v[i].regions[2].gnorm);

                        for (int j = temp.length(); j <= precision; j++) {
                                o << " ";
                        }

                        o << temp << "          ";
			//o << setw(14) << v[i].orbital_spread[0].spread*hartree_to_eV;
			temp_double = v[i].orbital_spread[0].spread*hartree_to_eV;
			temp = NumberToString(temp_double);

                        for (int j = temp.length(); j <= precision; j++) {
                                o << " ";
                        }

                        o << temp << "  ";
                        //o << setw(19) << v[i].HOMO_value*hartree_to_eV;
                        temp_double = v[i].HOMO_value*hartree_to_eV;
                        temp = NumberToString(temp_double);

                        for (int j = temp.length(); j <= precision; j++) {
                                o << " ";
                        }

                        o << temp << "  ";
                        //o << setw(11) << v[i].LUMO_value*hartree_to_eV;
                        temp_double = v[i].LUMO_value*hartree_to_eV;
                        temp = NumberToString(temp_double);

                        for (int j = temp.length(); j <= precision; j++) {
                                o << " ";
                        }

                        o << temp << "    |  ";
			// o << setw(8) << v[i].function;
			temp = NumberToString(v[i].function);

                        for (int j = temp.length(); j <= precision; j++) {
                                o << " ";
                        }

                        o << temp << " ";	

                        temp = NumberToString(v[i].rank);

                        for (int j = temp.length(); j < 4; j++) {
                                o << " ";
                        }

                        o << "(" << temp << ")";
                        o << "       |      ";
                        // Add in DMA Spread at the end for now
                        if (v[i].dma_spread.size() > 0)
			{
				double sum_total = 0;
				for (vector<min_max_spread>::size_type j = 0; j < v[i].dma_spread.size(); j++)
				{
					sum_total += v[i].dma_spread[j].spread;
				}

				temp = NumberToString(sum_total);

	                        for (int j = temp.length(); j <= precision; j++) {
        	                        o << " ";
                	        }

                        	o << temp;
				 
			}
			o << endl;
		}
		o << spacer << spacer << endl;
		
	}
	else
	{
		string error = "Could not append to log file: " + output;
		cout << error << endl;
		if (critical)
		{
			cout << "Critical Error" << endl;
			exit(EXIT_FAILURE);
		}
	}
	// Close output
	o.close();
}

/*
 This is a function to update the regions file
 
 @param[in] output Filename
 @param[in] v Vector of gaussians, containing all the ECPs to be written
 @param[in] critical Catch for critical errors to terminate program
 */
void update_regions_file(string output, vector<gaussian> v, vector<int> v_history, bool critical)
{
	ofstream o;
	o.open( output.c_str(), ios::app );
	
	if (o)
	{
		// THIS NEEDS UPDATING FOR NEW STRUCTURE
		for (vector<gaussian>::size_type i = 0; i != v.size(); i++)
		{
			if (v_history[i] != 1)
			{
				o.precision(8);
				// Loop over all values and print
				for (vector<gaussian_info>::size_type j = 0; j < v[i].values.size(); j++)
				{
					if (j == 0)
					{
						o << v[i].index;
					}
					else
					{
						o << " ";
					}

					o << "\t|\t" << setw(10) << v[i].values[j].line_number << "\t\t" << v[i].values[j].type << "\t\t" << v[i].values[j].value << endl;
				}
				// Separate from other values
				o << spacer << endl;
			}
		}
	}
	else
	{
		string error = "Could not append to log file: " + output;
		cout << error << endl;
		if (critical)
		{
			cout << "Critical Error" << endl;
			exit(EXIT_FAILURE);
		}
	}
	// Close output
	o.close();
}

/*
 Generic function to read in a TEXT file and return it in a vector<string>
 
 @param[in] input Filename
 @param[in] critical Flag to terminate if there is a critical problem
 @return Vector with string contents of file in question.
  */
vector<string> read_in_lines(string input, bool critical)
{
	ifstream inData;
	vector<string> outData;
	// Open file
	inData.open( input.c_str() );
	if (inData) 
	{
		string buffer = "";
		// Get all lines
		while (getline(inData,buffer))
		{
			outData.push_back(buffer);
		}
	}
	else
	{
		// Else throw error
		string error = "Could not find file: " + input;
		cout << error << endl;
		if (critical)
		{
			cout << "Critical Error" << endl;
			exit(EXIT_FAILURE);
		}
	}
	inData.close();	
	
	return outData;
}

/*
 Generic function to read in a BINARY file and return it in a vector<string>

 @param[in] input Filename
 @param[in] critical Flag to terminate if there is a critical problem
 @return Vector of Gaussians with contents of file in question.
  */
vector<gaussian> read_in_binary(string input, bool critical)
{
        ifstream inData;

        // cout << "OPENING: " << input.c_str() << endl;

	inData.open(input.c_str(), ios::in | ios::binary);

        vector<gaussian> outData;

        if (inData)
        {
	        // Temp structure;
	        gaussian g;
                size_t sz;

                while (inData.good())
		{
                        inData.read(reinterpret_cast<char*>(&g.index), sizeof(int));
	                inData.read(reinterpret_cast<char*>(&g.rank), sizeof(int));
                        inData.read(reinterpret_cast<char*>(&g.HOMO_value), sizeof(double));
                        inData.read(reinterpret_cast<char*>(&g.LUMO_value), sizeof(double));
                        inData.read(reinterpret_cast<char*>(&g.function), sizeof(double));
                        inData.read(reinterpret_cast<char*>(&g.failed), sizeof(bool));

	                inData.read(reinterpret_cast<char*>(&sz), sizeof(sz));
	                g.values.resize(sz);
        	        for (vector<gaussian_info>::size_type sz_t = 0; sz_t < sz; sz_t++)
                	{
	                        inData.read(reinterpret_cast<char*>(&g.values[sz_t]), sizeof(gaussian_info));
                	}

	                inData.read(reinterpret_cast<char*>(&sz), sizeof(sz));
        	        g.regions.resize(sz);
                	for (vector<regions_data>::size_type sz_t = 0; sz_t < sz; sz_t++)
	                {
                	        inData.read(reinterpret_cast<char*>(&g.regions[sz_t]), sizeof(regions_data));
	                }


        	        inData.read(reinterpret_cast<char*>(&sz), sizeof(sz));
                	g.orbital_spread.resize(sz);
                        // Generate character buffer for reading string 
                        char buffer [11];
                        // This is necessary as we cannot write a string directly to binary
	                for (vector<min_max_spread>::size_type sz_t = 0; sz_t < sz; sz_t++)
        	        {
//                                inData.read(reinterpret_cast<char*>(&g.orbital_spread[sz_t]), sizeof(min_max_spread));
                                inData.read(reinterpret_cast<char*>(&g.orbital_spread[sz_t].max), sizeof(double));
                                inData.read(reinterpret_cast<char*>(&g.orbital_spread[sz_t].min), sizeof(double));
                                inData.read(reinterpret_cast<char*>(&g.orbital_spread[sz_t].average), sizeof(double));
                                inData.read(reinterpret_cast<char*>(&g.orbital_spread[sz_t].spread), sizeof(double));
                                inData.read(reinterpret_cast<char*>(&g.orbital_spread[sz_t].quantity), sizeof(int));
                                // Read characters
                                inData.read(buffer, 10);
                                // Set null point to signify end of string
                                buffer[10] = 0;
                                // Assign
                                g.orbital_spread[sz_t].label = buffer;
        	        }
	
        	        inData.read(reinterpret_cast<char*>(&sz), sizeof(sz));
                	g.dma_spread.resize(sz);
                        // Generate character buffer for reading string
                        char buffer2 [11];
                        // This is necessary as we cannot write a string directly to binary
	                for (vector<min_max_spread>::size_type sz_t = 0; sz_t < sz; sz_t++)
        	        {
//                                inData.read(reinterpret_cast<char*>(&g.dma_spread[sz_t]), sizeof(min_max_spread));
                                inData.read(reinterpret_cast<char*>(&g.dma_spread[sz_t].max), sizeof(double));
                                inData.read(reinterpret_cast<char*>(&g.dma_spread[sz_t].min), sizeof(double));
                                inData.read(reinterpret_cast<char*>(&g.dma_spread[sz_t].average), sizeof(double));
                                inData.read(reinterpret_cast<char*>(&g.dma_spread[sz_t].spread), sizeof(double));
                                inData.read(reinterpret_cast<char*>(&g.dma_spread[sz_t].quantity), sizeof(int));
                                // Read characters
                                inData.read(buffer2, 10);
                                // Set null point to signify end of string
                                buffer2[10] = 0;
                                // Assign
                                g.dma_spread[sz_t].label = buffer2;
        	        }

			// Horrible segfaults from trying to put the Gaussian straight into the vector.
			// i.e. outData.push_pack(g);
			// Trying to do it manually now... 

                        // cout << "READ ONE ENTRY" << endl;
		
			if (!inData.eof())
			{
                                // cout << "RESIZING OUTDATA. OLD SIZE:" << outData.size();
	                        outData.resize(outData.size()+1);
                                // cout << "; NEW SIZE: " << outData.size() << endl;
                                // cout << "ADDING INDEX: " << g.index << endl;
				outData[outData.size()-1].index = g.index;
                                // cout << "ADDING RANK: " << g.rank << endl;
                	        outData[outData.size()-1].rank = g.rank;
                                // cout << "ADDING HOMO: " << g.HOMO_value << endl;
                        	outData[outData.size()-1].HOMO_value = g.HOMO_value;
                                // cout << "ADDING LUMO: " << g.LUMO_value << endl;
	                        outData[outData.size()-1].LUMO_value = g.LUMO_value;
                                // cout << "ADDING FUNCTION: " << g.function << endl;
        	                outData[outData.size()-1].function = g.function;
                                // cout << "ADDING FAILED: " << g.failed << endl;
                	        outData[outData.size()-1].failed = g.failed;
                                // cout << "ADDING VALUES: " << g.values.size() << endl;
                        	outData[outData.size()-1].values = g.values;
                                // cout << "ADDING REGIONS: " << g.regions.size() << endl;
	                        outData[outData.size()-1].regions = g.regions;
                                // cout << "ADDING ORBITAL SPREAD: " << g.orbital_spread.size() << endl;
                                // cout << "ORBITAL SPREAD VALUES: " << endl;
                                // cout << "MAX: " << g.orbital_spread[0].max << endl;
                                // cout << "MIN: " << g.orbital_spread[0].min << endl;
                                // cout << "AVERAGE: " << g.orbital_spread[0].average << endl;
                                // cout << "SPREAD: " << g.orbital_spread[0].spread << endl;
                                // cout << "QUANTITY: " << g.orbital_spread[0].quantity << endl;
                                // cout << "LABEL: " << g.orbital_spread[0].label << endl; 
        	                outData[outData.size()-1].orbital_spread = g.orbital_spread;
                                // cout << "ADDING DMA SPREAD: " << g.dma_spread.size() << endl;
                	        outData[outData.size()-1].dma_spread = g.dma_spread;
			}
                        
                        // cout << "ASSIGNED ONE ENTRY" << endl;
		}
        }
        else
        {
                // Else throw error
                string error = "Could not find file: " + input;
                cout << error << endl;
                if (critical)
                {
                        cout << "Critical Error" << endl;
                        exit(EXIT_FAILURE);
                }
        }

        inData.close();

        return outData;
}


/*
 Generic function take a vector and write it to file AS TEXT.
 
 @param[in] output Filename
 @param[in] content Pointer to Vector containing output data (Can't remember why I used a pointer?)
 @param[in] critical Error flag if there is a problem
 */
void write_out_lines(string output, vector<string> *content, bool critical)
{
	ofstream outData;
	// Open file
	outData.open( output.c_str() );
	if (outData) 
	{
		// Output all the data
		for (vector<string>::size_type a = 0; a != content->size(); a++)
		{
			outData << content->at(a) << endl;
		}
	}
	else
	{
		// Else throw error
		string error = "Could not open output file: " + output;
		cout << error << endl;
		if (critical)
		{
			cout << "Critical Error" << endl;
			exit(EXIT_FAILURE);
		}
	}

        // Empty the vector before returning
        content->clear();

	outData.close();	
}

/*
 Generic function take a vector of Gaussians and write it to file AS BINARY.

 @param[in] output Filename
 @param[in] content Pointer to Vector containing output data (Can't remember why I used a pointer?)
 @param[in] critical Error flag if there is a problem
 */
void write_out_binary(string output, vector<gaussian> content, bool critical)
{
        ofstream outData;
	outData.open(output.c_str(), ios::out | ios::binary);

        if (outData)
        {
                // Output all the data
                for (vector<gaussian>::size_type a = 0; a < content.size(); a++)
                {
                        outData.write(reinterpret_cast<char*>(&content[a].index), sizeof(int));
                        outData.write(reinterpret_cast<char*>(&content[a].rank), sizeof(int));
                        outData.write(reinterpret_cast<char*>(&content[a].HOMO_value), sizeof(double));
                        outData.write(reinterpret_cast<char*>(&content[a].LUMO_value), sizeof(double));
                        outData.write(reinterpret_cast<char*>(&content[a].function), sizeof(double));
                        outData.write(reinterpret_cast<char*>(&content[a].failed), sizeof(bool));

                        size_t sz = content.at(a).values.size();
			outData.write(reinterpret_cast<char*>(&sz), sizeof(sz));
			// Write contents of values
			for (vector<gaussian_info>::size_type sz_t = 0; sz_t < sz; sz_t++)
			{
				outData.write(reinterpret_cast<char*>(&content[a].values[sz_t]), sizeof(gaussian_info));
			}

                        sz =  content.at(a).regions.size(); 
			outData.write(reinterpret_cast<char*>(&sz), sizeof(sz));
                        // Write contents of regions
                        for (vector<regions_data>::size_type sz_t = 0; sz_t < sz; sz_t++)
                        {
                                outData.write(reinterpret_cast<char*>(&content[a].regions[sz_t]), sizeof(regions_data));
                        }

			sz =  content.at(a).orbital_spread.size();
			outData.write(reinterpret_cast<char*>(&sz), sizeof(sz));
                        // Can we write this all at once? The answer is no due to the string 
                        // Write contents of orbital_spread
                        for (vector<min_max_spread>::size_type sz_t = 0; sz_t < sz; sz_t++)
                        {
//                                outData.write(reinterpret_cast<char*>(&content[a].orbital_spread[sz_t]), sizeof(min_max_spread));
                                  outData.write(reinterpret_cast<char*>(&content[a].orbital_spread[sz_t].max), sizeof(double));
                                  outData.write(reinterpret_cast<char*>(&content[a].orbital_spread[sz_t].min), sizeof(double));
                                  outData.write(reinterpret_cast<char*>(&content[a].orbital_spread[sz_t].average), sizeof(double));
                                  outData.write(reinterpret_cast<char*>(&content[a].orbital_spread[sz_t].spread), sizeof(double));
                                  outData.write(reinterpret_cast<char*>(&content[a].orbital_spread[sz_t].quantity), sizeof(int));
//                                  outData << &content[a].orbital_spread[sz_t].label << std::endl;
//                                  cout << content[a].orbital_spread[sz_t].label.c_str() << endl;
                                  outData.write(content[a].orbital_spread[sz_t].label.c_str(), 10);
                        }

			sz =  content.at(a).dma_spread.size();
			outData.write(reinterpret_cast<char*>(&sz), sizeof(sz));
                        // Write contents of dma_spread
                        for (vector<min_max_spread>::size_type sz_t = 0; sz_t < sz; sz_t++)
                        {
//                                outData.write(reinterpret_cast<char*>(&content[a].dma_spread[sz_t]), sizeof(min_max_spread));
                                  outData.write(reinterpret_cast<char*>(&content[a].dma_spread[sz_t].max), sizeof(double));
                                  outData.write(reinterpret_cast<char*>(&content[a].dma_spread[sz_t].min), sizeof(double));
                                  outData.write(reinterpret_cast<char*>(&content[a].dma_spread[sz_t].average), sizeof(double));
                                  outData.write(reinterpret_cast<char*>(&content[a].dma_spread[sz_t].spread), sizeof(double));
                                  outData.write(reinterpret_cast<char*>(&content[a].dma_spread[sz_t].quantity), sizeof(int));
//                                  outData << &content[a].dma_spread[sz_t].label << std::endl;
//                                  cout << content[a].orbital_spread[sz_t].label.c_str() << endl;
                                  outData.write(content[a].dma_spread[sz_t].label.c_str(), 10);
                        } 

			// Write all data
                        //outData.write(reinterpret_cast<char*>(&content.at(a)), sizeof(gaussian));
                }
        }
        else
        {
                // Else throw error
                string error = "Could not open output file: " + output;
                cout << error << endl;
                if (critical)
                {
                        cout << "Critical Error" << endl;
                        exit(EXIT_FAILURE);
                }
        }

        outData.close();
}

