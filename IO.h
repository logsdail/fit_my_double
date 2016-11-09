/*
 *  @io.h
 *  fit_my_ecp
 *
 *  @brief Handles IO interactions, reading in files and writing text output
 *
 *  Created by Andrew Logsdail on 01/06/2012.
 *  Copyright 2012 University College London. All rights reserved.
 *
 */

#ifndef IO_H
#define IO_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
// Personal headers
#include "Utils.h"
#include "Structures.h"

// Function to update log file
void update_log_file(std::string output, std::vector<gaussian> v, bool critical = true);
// Function to update regions file
// void update_regions_file(std::string output, std::vector<gaussian> v, bool critical = true);
void update_regions_file(std::string output, std::vector<gaussian> v, std::vector<int> v_history, bool critical = true);
// Generic function read in a file and return it in a vector
std::vector<std::string> read_in_lines(std::string input, bool critical = true);
std::vector<gaussian> read_in_binary(std::string input, bool critical = true);
// Generic function take a vector and write it to file
void write_out_lines(std::string output, std::vector<std::string> *content, bool critical = true);
void write_out_binary(std::string output, std::vector<gaussian> content, bool critical = true);
#endif
