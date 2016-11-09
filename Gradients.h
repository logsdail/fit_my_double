/*
 *  @file Gradients.h
 *  fit_my_ecp
 *
 *  @brief Class specifically designed for understanding the ChemShell gradient outputs
 *  which are normally written in punch card format
 *
 *  Created by Andrew Logsdail on 03/06/2012.
 *  Copyright 2012 University College London. All rights reserved.
 *
 */

#ifndef GRADIENTS_H
#define GRADIENTS_H
#include <iostream>
// Personal headers
#include "Utils.h"
#include "Structures.h"

std::vector<regions_data> digest_gradients(std::vector<std::string> gradient_input, int total_centres, 
					   std::vector<int> centre_regions_total, std::vector<int> centre_regions,
					   bool absolute_gradients);

#endif
