/*
 *  @file Main.cpp
 *  fit_my_ecp
 *
 *  @brief This class brings everything together. Currently it decomposes the inputs,
 *  but long term this needs to be outsourced to a variable class
 *
 *  Created by Andrew Logsdail on 14/04/2012.
 *  Copyright 2012 University College London. All rights reserved.
 *
 *
 */

// TO DO:
// - Forces on Shells. Editing of ChemShell
// - Spin Polarised ECPs

// System headers
#include <iostream>
#include <iomanip>
// Personal headers
#include "Utils.h"
#include "Structures.h"
#include "IO.h"
#include "Functions.h"
#include "Genetic.h"
#include "Linear.h"
#include "Powells.h"
#include "Outputs.h"
#include "Newton_Raphson.h"
#include "History.h"
#include "Gradients.h"
#include "Punch.h"
#include "DFT_Program.h"
#include "Gamess_UK.h"
#include "Nwchem.h"

using namespace std;

/*
 Shows the "help" contents for running the program
 
 No params
 */
void inputs()
{
	// Current options and descriptor
	// This is going to need updating at times
	cout << "*** System Options ***" << endl;
	cout << endl;
	cout << "-d,--dryrun : Check everything runs without submitting ChemShell calculations" << endl;
	cout << "-h,--help   : Display this message" << endl;
	cout << endl;
	cout << "*** Character Values ***" << endl;
	cout << endl;
	cout << "-f,--function=outputs : Just perform analysis of output files.               Requires -pt" << endl;
	cout << "              linear  : Performs Linear scan between mins and maxes.         Requires -pt, -ef, -et, -pf and -cf" << endl;
	cout << "              powells : Performs downhill minimisation to best ECP.          Requirements as above" << endl;
	cout << "              newton  : Performs quasi-newtonian minimisation.               Requirements as above" << endl;
        cout << "              lbfgs   : Performs minimimisation.                             Requirements as above" << endl; 
	cout << "              ga      : Performs global optimisation with genetic algorithm. Requirements as above" << endl;
	cout << endl;
	cout << "-ef,--ecpfile=ECP_FILENAME              : Output location of ECP file, as defined in CHM_FILE" << endl;
	cout << "-et,--ecptemplate=ECP_FILENAME          : Locaton of input ECP template file" << endl;
	cout << "-cf,--chmfile=CHM_FILENAME              : Location of CHM file" << endl;
	cout << "-pf,--punchfile=PUNCH_FILENAME          : Output location of punch file, as defined in CHM_FILE" << endl;
	cout << "-pt,--punchtemplate=PUNCH_FILENAME      : Input punch template file" << endl;
	cout << "-qmo,--qmoutput=QUANTUM_OUTPUT_FILENAME : QM Calculation output. Default: gamess1.out.1" << endl;
	cout << "-go,--gradientoutput=GRADIENT_FILENAME  : Gradient output. Default: gradient" << endl;
	cout << "-of,--outputfolder=OUTPUT_FOLDER        : Folder to move results in to. Default: results" << endl;
        cout << "-e,--executable=EXECUTABLE              : ChemShell executable" << endl;
        cout << "-a,--anion=CHEMICAL_SYMBOL              : Species for which we are comparing deep lying orbitals (e.g. O 1s)" << endl;
	cout << endl;
	cout << "*** Numeric Values ***" << endl;
	cout << endl;
	cout << "--seed=NUMBER                : Intial seed for random number generator. Default: 0" << endl;
	cout << "--stepsize(A/Z)=NUMBER       : Intial step size for search of A/Z. Default: 1" << endl;
	cout << "--stepreduction(A/Z)=NUMBER  : Intial step size reduction when at minimum of A/Z. Default: 0.5" << endl;
	cout << "--stepmin(A/Z)=NUMBER        : Convergence requirement which must be bettered to exit search of A/Z. Default: 1" << endl;	
	cout << "--min(A/Z)=NUMBER            : Minimum value of A/Z. Default: 0.0" << endl;
	cout << "--max(A/Z)=NUMBER            : Maximum value of A/Z. Default: 1000" << endl;
	cout << "--maxsteps=NUMBER            : Maximum number of steps in one direction. Default: OFF" << endl;
	cout << "--chemshellmax=NUMBER        : Alternative criteria to exit search. Default: 5" << endl;		
	cout << "--anion_offset=NUMBER        : Offset from top of eigenvalue list for spread of interest. Default: 0" << endl;
// To be considered implementing. Will require something nifty for different punch files.
//        cout << "--anion_spread=NUMBER      : Number of eigenvalues to include in spread of interest. Default: Number of anion species in Region 1)" << endl;
        cout << "--processors=NUMBER          : Total number of processors (HECToR)" << endl;
        cout << "--processors_per_node=NUMBER : Number of processors per node (HECToR)" << endl;
        cout << endl;
	cout << "*** Function Weights ***" << endl;
	cout << endl;
	cout << "--wt_gnorm1=NUMBER       : Weighting for region 1 gnorm average in quadratic. Default: 70.0" << endl;
	cout << "--wt_gnorm2=NUMBER       : Weighting for region 2 gnorm average in quadratic. Default: 20.0" << endl;
	cout << "--wt_gnorm3=NUMBER       : Weighting for region 3 gnorm average in quadratic. Default: 0.1" << endl;
	cout << "--wt_anion_spread=NUMBER : Weighting for Anionic Spread of S-Orbitals in quadratic. Default: 1.0" << endl;
	cout << "--wt_eigenvalues=NUMBER  : Weighting for difference between eigenvalues in quadratic. Default: 0.05" << endl;
        cout << "--wt_dma_spread=NUMBER   : Weighting for Atomic Spread of Distributed Multipole Analysis in quadratic. Default: 0.2" << endl;
        cout << endl;
        cout << "*** Function Targets ***" << endl;
	cout << endl;
        cout << "--tg_gnorm1=NUMBER       : Target for region 1 gnorm average in quadratic. Default: 0.0" << endl;
        cout << "--tg_gnorm2=NUMBER       : Target for region 2 gnorm average in quadratic. Default: 0.0" << endl;
        cout << "--tg_gnorm3=NUMBER       : Target for region 3 gnorm average in quadratic. Default: 0.0" << endl;
        cout << "--tg_anion_spread=NUMBER : Target for Anionic Spread of S-Orbitals in quadratic. Default: 0.0" << endl;
        cout << "--tg_homo=NUMBER         : Target for highest occupied molecular orbital in quadratic. Default: 0.0" << endl;
        cout << "--tg_lumo=NUMBER         : Target for lowest unoccupied molecular orbital in quadratic. Default: 0.0" << endl;
        cout << "--tg_dma_spread=NUMBER   : Target for Atomic Spread of Distributed Multipole Analysis in quadratic. Default: 0.0" << endl;
        cout << endl;
        cout << "*** GA Settings ***" << endl;
	cout << endl;
	cout << "--ga_population=NUMBER   : Population Size for GA run. Default: 4" << endl;
	cout << "--ga_offspring=NUMBER    : Offspring Size for GA run. Default: 2" << endl;
	cout << "--ga_mutations=NUMBER    : Mutations Size for GA run. Default: 2" << endl;
	cout << "--ga_convergence=NUMBER  : Convergence Criteria for GA run. Default: 3" << endl;
	cout << endl;
	cout << "*** Boolean Options ***" << endl;
	cout << endl;
	cout << "--ga_mutation_dynamic    : Use dynamic mutation in GA Search" << endl;
        cout << "--force_history_recalc   : All functions read in from History will be recalculated, to account for lost accuracy in outputs" << endl;
//        cout << "--remove_history_duds    : Remove all Gaussians from History that are duds (888888), thus forcing their re-run" << endl;
	cout << "--not_absolute_gradients : Do not use absolute gradients, but just as-read values (for 1D systems)" << endl;
}

/*
 Error function which exits the program if we have a critical error
 
 @param[in] crit Defines if this is a critical error
 */
void critical_error(bool crit)
{
	if (crit)
	{
		cout << "Critical Error" << endl;
		// Exit as we won't manage to run
		exit(EXIT_FAILURE);	
	}
}

/*
 Function to display standard "ChemShell" message as calculation proceeds
 
 @param[in] g Current ECP that is going to be used
 @param[in] cc Calculation number
 @param[in] ci Current Index for thsi Gaussian
 */
void print_chemshell_message(gaussian g, int cc, int ci)
{
	cout << endl;
	cout << "Calculation   : " << cc << endl;
	cout << "Current Index : " << ci << endl;
	// Print to screen so we know what is being fitted
	for (vector<gaussian_info>::size_type a = 0; a < g.values.size(); a++)
	{
		// Default is a coefficient
		string temp = "Coefficient";
		// Coefficient Type == 1
		if (g.values[a].type == 1)
		{
			temp = "Exponent   ";
		}
		cout << temp << "   : " << g.values[a].value << "\t ( Line : " << g.values[a].line_number << ")" << endl;
	}
	cout << endl;	
}

/*
 Main method. Here we read in, organise and perform the ECP minimisation
 Most of the IO is outsourced, as is managing which ECPs to calculate with
 
 @param[in] argc Number of run time arguments
 @param[in] argv[] Runtime arguments
 @return Exit success integer
 */
int main (int argc, char * const argv[]) {
	
	// Initial mark to ensure we get ALL the detail we want
	cout.precision(8);
	// Open up with some spacing
	cout << endl;
	
	// Integers
	int random_seed = 0;
	int current_index = 0;
	int chemshell_counter = 0;
	int chemshell_counter_max = 0;
	int region_1_anion_offset = 0;
	// Numbers for the processors
	// To be used later with parallelisation
	int processors_per_node = 0;
	int processors = 0;
	// Strings
	string function = "";
	string ecp_file = "";
	string punch_file = "";
	string chm_file = "";
        string executable = "";
	string ecp_template_file = "";
	//string punch_template_file = "";
	string qm_output_file = "";
	string log_output_file = "output.log";
	string regions_output_file = "regions.log";
        string dma_output_file = "dma.log";
	string gradient_output_file = "";
	string output_folder = "";
        string anion_species = "";
	// string command_line = "";
	// Vectors
	vector<string> outData;
	// Added for multiple fits
	vector<string> punch_template_files;
        vector<string> log_output_files;
	// Boolean
	bool dry_run = false;
	bool outputs_only = false;
	bool punch_output_check = false;
	bool force_recalc = false;
//        bool remove_duds = false;
        bool absolute_gradients = true;
	// Some classes to do the important stuff
	Outputs *ecp_searcher = NULL;
	Functions *func_calc = new Functions();
	// History *ecps_history = new History();
        vector<History *> ecps_history;

	//Gamess_UK *DFT_program = new Gamess_UK();
        DFT_Program *qm_program = NULL;
	// Punch punch;
	vector<Punch> punch;

	cout << "Parsing inputs:" << endl;
	
	// Read input parameters.
	// These can be taken from argv.
	// Remember that argc by default has length 1, with nothing in the first part of argv i.e. argv[0] == ''
	if (argc > 1) // Check if we have input parameters
	{
		// Variables we might need just in here
		int max_steps_in_one_direction = 0;
		// This will save a little on memory
		vector< vector<double> > step_size(2);
		vector< vector<double> > step_reduction(2);
		vector< vector<double> > step_size_min(2);
		vector<double> minimums(2,0.0);
		vector<double> maximums(2,0.0);
		// Initialise GA parameters
		int population_size = 0;
		int mutations_size = 0;
		int offspring_size = 0;
		int convergence_criteria = 0;
		// We are going to put in the function weights here
		// These will be put to defaults in the appropriate region
		// Initiate to -1 as 0 is acceptable
		double weight_gnorm_region1 = -1;
		double weight_gnorm_region2 = -1;
		double weight_gnorm_region3 = -1;
		double weight_anion_spread = -1;
		double weight_eigenvalues = -1;
		double weight_dma_spread = -1;	
		// Defined HOMO/LUMO targets
		// double target_homo = 0;
                // double target_lumo = 0;
		// Define other targets
		// double target_gnorm_region1 = 0;
		// double target_gnorm_region2 = 0;
		// double target_gnorm_region3 = 0;
		// double target_anion_spread = 0;
		// double target_dma_spread = 0;
                // Redefine for simultaneous fitting
                vector<double> target_homo;
                vector<double> target_lumo;
                vector<double> target_gnorm_region1;
                vector<double> target_gnorm_region2;
                vector<double> target_gnorm_region3;
                vector<double> target_anion_spread;
                vector<double> target_dma_spread;
		// Strings
		string argv_string = "";
		size_t argv_splitter;
		string argv_value = "";
		string argv_variable = "";
		// Boolean value
		bool mutation_dynamic = false;
		
		for (int argc_counter = 1; argc_counter < argc; argc_counter++)
		{
			argv_string = argv[argc_counter];
                        // Remove leading hyphen. We don't need it
			// This is an update from the initial design
			// Hopefully that'll allow me to pass files to the program as input
			// if (argv_string[0] == '-')
			while (argv_string[0] == '-')
			{
				string::iterator it = argv_string.begin();
				argv_string.erase(it);
			}
			// Separate value from parameter
			argv_splitter = argv_string.find_first_of('=');
			if (argv_splitter != string::npos)
			{
				argv_variable = argv_string.substr(0,argv_splitter);
				argv_value = argv_string.substr(argv_splitter+1);
				
				//cout << argv_variable << " " << argv_value;
				// This seems to split the input fine.
				// I'm going to include my favourite file Utils.cpp here, and use the word comparison
				
				if (cmpStr("function",argv_variable) || cmpStr("f",argv_variable))
				{
					function = argv_value;	
				} 
				else if (cmpStr("ecpfile",argv_variable) || cmpStr("ef",argv_variable))
				{
					ecp_file = argv_value;	
				} 
				else if (cmpStr("ecptemplate",argv_variable) || cmpStr("et",argv_variable))
				{
					ecp_template_file = argv_value;
				}
				else if (cmpStr("chmfile",argv_variable) || cmpStr("cf",argv_variable))
				{
					chm_file = argv_value;
				}
				else if (cmpStr("punchfile",argv_variable) || cmpStr("pf",argv_variable))
				{
					punch_file = argv_value;
				}
				else if (cmpStr("punchtemplate",argv_variable) || cmpStr("pt",argv_variable))
				{
					// These need separating
					vector<string> tokens;
					//cout << argv_value << " " << argv_value.length() << endl;
                                        Tokenize(argv_value,tokens,", ");
					// Push these values back
					for (vector<string>::size_type a = 0; a < tokens.size(); a++)
                                        {
						//cout << tokens[a] << " " << tokens[a].length() << endl;
                                                punch_template_files.push_back(tokens[a]);
                                       	}
					//punch_template_file = argv_value;
				}
				else if (cmpStr("qmoutput",argv_variable) || cmpStr("qmo",argv_variable))
				{
					qm_output_file = argv_value;
				}
				else if (cmpStr("gradientoutput",argv_variable) || cmpStr("go",argv_variable))
				{
					gradient_output_file = argv_value;
				}
				else if (cmpStr("outputfolder",argv_variable) || cmpStr("of",argv_variable))
				{
					output_folder = argv_value;
				}
				else if (cmpStr("executable",argv_variable) || cmpStr("e",argv_variable))
				{
					executable = argv_value;
				}
				else if (cmpStr("anion",argv_variable) || cmpStr("a",argv_variable))
                                {
                                        anion_species = argv_value;
                                }
				// And now we'll collect numbers
				else if (cmpStr("seed",argv_variable))
				{
					StringToNumber(argv_value,random_seed);
					// Set this to negative as we need to make that happen
					random_seed -= random_seed;
				}
                                else if (cmpStr("step",argv_variable.substr(0,4)))
				{
                                        // Split this up as it may contain multiple values
                                        // For now we assume they are separated by commas (,)
                                        vector<string> tokens;
                                        Tokenize(argv_value,tokens,",");
                                        double temp_double = 0;
                                        // Add values to arrays
					if  (cmpStr("stepsizeA",argv_variable))
					{
	                                        for (vector<string>::size_type a = 0; a < tokens.size(); a++)
        	                                {
                	                                StringToNumber(tokens[a],temp_double);
                        	                        step_size[0].push_back(temp_double);
                                	        }
					}
					else if (cmpStr("stepsizeZ",argv_variable))
					{
                                                for (vector<string>::size_type a = 0; a < tokens.size(); a++)
                                                {
                                                        StringToNumber(tokens[a],temp_double);
                                                        step_size[1].push_back(temp_double);
                                                }
                                        }
					else if (cmpStr("stepreductionA",argv_variable))
	                                {
                                                for (vector<string>::size_type a = 0; a < tokens.size(); a++)
                                                {
                                                        StringToNumber(tokens[a],temp_double);
                                                        step_reduction[0].push_back(temp_double);
                                                }
                                        }
                                        else if (cmpStr("stepreductionZ",argv_variable))
                                        {
                                                for (vector<string>::size_type a = 0; a < tokens.size(); a++)
                                                {
                                                        StringToNumber(tokens[a],temp_double);
                                                        step_reduction[1].push_back(temp_double);
                                                }
                                        }
                                        else if (cmpStr("stepminA",argv_variable))
                                        {
                                                for (vector<string>::size_type a = 0; a < tokens.size(); a++)
                                                {
                                                        StringToNumber(tokens[a],temp_double);
                                                        step_size_min[0].push_back(temp_double);
                                                }
                                        }
                                        else if (cmpStr("stepminZ",argv_variable))
                                        {
                                                for (vector<string>::size_type a = 0; a < tokens.size(); a++)
                                                {
                                                        StringToNumber(tokens[a],temp_double);
                                                        step_size_min[1].push_back(temp_double);
                                                }
                                        }
				}
//				else if (cmpStr("stepreductionA",argv_variable))
//				{
//					StringToNumber(argv_value,step_reduction[0]);
//				}
//				else if (cmpStr("stepminA",argv_variable))
//				{
//					StringToNumber(argv_value,step_size_min[0]);
//				}
//				else if (cmpStr("stepsizeZ",argv_variable))
//				{
//					StringToNumber(argv_value,step_size[1]);
//				}
//				else if (cmpStr("stepreductionZ",argv_variable))
//				{
//					StringToNumber(argv_value,step_reduction[1]);
//				}
//				else if (cmpStr("stepminZ",argv_variable))
//				{
//					StringToNumber(argv_value,step_size_min[1]);
//				}
				else if (cmpStr("minA",argv_variable))
				{
					StringToNumber(argv_value,minimums[0]);
				}
				else if (cmpStr("maxA",argv_variable))
				{
					StringToNumber(argv_value,maximums[0]);
				}
				else if (cmpStr("minZ",argv_variable))
				{
					StringToNumber(argv_value,minimums[1]);
				}
				else if (cmpStr("maxZ",argv_variable))
				{
					StringToNumber(argv_value,maximums[1]);
				}
				else if (cmpStr("maxsteps",argv_variable))
				{
					StringToNumber(argv_value,max_steps_in_one_direction);
				}
				else if (cmpStr("chemshellmax",argv_variable))
				{
					StringToNumber(argv_value,chemshell_counter_max);
				}
				else if (cmpStr("anion_offset",argv_variable))
				{
					StringToNumber(argv_value,region_1_anion_offset);
				}
                                else if (cmpStr("tg_",argv_variable.substr(0,3)))
				{
					// Split this up as it may contain multiple values
					// For now we assume they are separated by commas (,)
					vector<string> tokens;
			                Tokenize(argv_value,tokens,",");
					double temp_double = 0;
					// Add values to arrays
					if (cmpStr("tg_homo",argv_variable))
                                	{
						for (vector<string>::size_type a = 0; a < tokens.size(); a++)
                                                {
							StringToNumber(tokens[a],temp_double);
                                                        target_homo.push_back(temp_double);
                                                }
                                        	// StringToNumber(argv_value,target_homo);
                                        	// Let's make sure this number is positive
                                        	// target_homo = sqrt(target_homo*target_homo);
                                        }
	                                else if (cmpStr("tg_lumo",argv_variable))
        	                        {
						for (vector<string>::size_type a = 0; a < tokens.size(); a++)
						{
							StringToNumber(tokens[a],temp_double);
                                        		target_lumo.push_back(temp_double);
						}
						// StringToNumber(argv_value,target_lumo);
                                        	// Let's make sure this number is positive
                                        	// target_lumo = sqrt(target_lumo*target_lumo);
	                                }
					else if (cmpStr("tg_gnorm1",argv_variable))
                	                {
						for (vector<string>::size_type a = 0; a < tokens.size(); a++)
                                                {
							StringToNumber(tokens[a],temp_double);
                                                        target_gnorm_region1.push_back(temp_double);
                                                }
                        	                //StringToNumber(argv_value,target_gnorm_region1);
	                                }
        	                        else if (cmpStr("tg_gnorm2",argv_variable))
                	                {
						for (vector<string>::size_type a = 0; a < tokens.size(); a++)
                                                {
							StringToNumber(tokens[a],temp_double);
                                                        target_gnorm_region2.push_back(temp_double);
                                                }
                        	                //StringToNumber(argv_value,target_gnorm_region2);
                                	}
					else if (cmpStr("tg_gnorm3",argv_variable))
        	                        {
						for (vector<string>::size_type a = 0; a < tokens.size(); a++)
                                                {
							StringToNumber(tokens[a],temp_double);
                                                        target_gnorm_region3.push_back(temp_double);
                                                }
                	                        //StringToNumber(argv_value,target_gnorm_region3);
                        	        }
					else if (cmpStr("tg_anion_spread",argv_variable))
        	                        {
						for (vector<string>::size_type a = 0; a < tokens.size(); a++)
                                                {
							StringToNumber(tokens[a],temp_double);
                                                        target_anion_spread.push_back(temp_double);
                                                }
                	                        //StringToNumber(argv_value,target_anion_spread);
                        	        }
	                                else if (cmpStr("tg_dma_spread",argv_variable))
        	                        {
						for (vector<string>::size_type a = 0; a < tokens.size(); a++)
                                                {
							StringToNumber(tokens[a],temp_double);
                                                        target_dma_spread.push_back(temp_double);
                                                }
                	                        //StringToNumber(argv_value,target_dma_spread);
                        	        }
                                }
				else if (cmpStr("wt_gnorm1",argv_variable))
				{
					StringToNumber(argv_value,weight_gnorm_region1);
				}
				else if (cmpStr("wt_gnorm2",argv_variable))
				{
					StringToNumber(argv_value,weight_gnorm_region2);
				}
				else if (cmpStr("wt_gnorm3",argv_variable))
				{
					StringToNumber(argv_value,weight_gnorm_region3);
				}
				else if (cmpStr("wt_anion_spread",argv_variable))
				{
					StringToNumber(argv_value,weight_anion_spread);
				}
				else if (cmpStr("wt_eigenvalues",argv_variable))
				{
					StringToNumber(argv_value,weight_eigenvalues);
				}
                                else if (cmpStr("wt_dma_spread",argv_variable))
                                {
                                        StringToNumber(argv_value,weight_dma_spread);
                                }
				else if (cmpStr("processors",argv_variable) || 
					 cmpStr("mppwidth",argv_variable)) // This dates to the previous key word
				{
					StringToNumber(argv_value,processors);
				}
				else if (cmpStr("processors_per_node",argv_variable) ||
					 cmpStr("mppnppn",argv_variable)) // This dates to the previous key word
				{
					StringToNumber(argv_value,processors_per_node);
				}
				else if (cmpStr("ga_population",argv_variable))
				{
					StringToNumber(argv_value,population_size);
				}
				else if (cmpStr("ga_offspring",argv_variable))
				{
					StringToNumber(argv_value,offspring_size);
				}
				else if (cmpStr("ga_mutations",argv_variable))
				{
					StringToNumber(argv_value,mutations_size);
				}
				else if (cmpStr("ga_convergence",argv_variable))
				{
					StringToNumber(argv_value,convergence_criteria);
				}
				else
				{
					// Insert error about correct usage
					cout << "Variable " << argv_variable << " is not recognised. Ignoring." << endl;;
					// But non-fatal
				}
			}
			else if (cmpStr("force_history_recalc",argv_string))
                        {
                                force_recalc = true;
                        }
//			else if (cmpStr("remove_history_duds",argv_string))
//                        {
//                                remove_duds = true;
//                        }
			else if (cmpStr("not_absolute_gradients",argv_string))
                        {
                                absolute_gradients = false;
                        }
			else if (cmpStr("ga_mutation_dynamic",argv_string))
			{
				mutation_dynamic = true;
			}
			else if (cmpStr("dryrun",argv_string) || cmpStr("d",argv_variable))
			{
				dry_run = true;
			}
			else if (cmpStr("help",argv_string) || cmpStr("h",argv_variable))
			{
				inputs();
				return EXIT_FAILURE;
			}
			else
			{
				// Insert error about correct usage
				cout << "Variable " << argv_string << " is not recognised. Ignoring." << endl;
				// But non-fatal
			}
		}
		
		if (cmpStr(function,"outputs"))
		{
			outputs_only = true;
			ecp_searcher = new Outputs(&random_seed);
		}
		else if (cmpStr(function,"linear"))
		{
			ecp_searcher = new Linear(&random_seed);
		} 
		else if (cmpStr(function,"powells"))
		{
			ecp_searcher = new Powells(&random_seed);
		}
		else if (cmpStr(function,"ga"))
		{
			ecp_searcher = new Genetic(&random_seed);
			
		}
                else if (cmpStr(function,"newton"))
    		{
			ecp_searcher = new Newton_Raphson(&random_seed);
		}
                else if (cmpStr(function,"lbfgs"))
                {
                        ecp_searcher = new Newton_Raphson(&random_seed,true);
                }
		else
		{
			cout << "Function is not defined. Please address this." << endl;
			critical_error(true);
		}

		// Check what we are doing
		ecp_searcher->print_type();
		
		// Set search direction parameters
		ecp_searcher->set_parameters(step_size,step_reduction,step_size_min,minimums,maximums,max_steps_in_one_direction);
		
		// Set GA parameters
		if (cmpStr(function,"ga"))
		{
			ecp_searcher->set_ga_parameters(population_size,mutations_size,offspring_size,
											convergence_criteria,mutation_dynamic);
		}
		
		// Set function calculation parameters
		// This is split so I can set the targets beforehand. They are normalised within this subroutine
		func_calc->set_targets_length(punch_template_files.size());

		func_calc->set_targets(target_gnorm_region1,target_gnorm_region2,target_gnorm_region3,
				       target_anion_spread,target_homo,target_lumo,target_dma_spread);
		//
		func_calc->set_weights(weight_gnorm_region1,weight_gnorm_region2,weight_gnorm_region3,
				       weight_anion_spread,weight_eigenvalues,weight_dma_spread);

		
	}
	else
	{
		// We need to put in something here for correct usage as we have no inputs
		inputs();
		return EXIT_FAILURE;
	}
	
	// Lets do the validity checks here;
	// This needs some more work for dryrun and outputs_only
	if (!outputs_only && (ecp_file.length() == 0 ||
		ecp_template_file.length() == 0 ||
		chm_file.length() == 0 ||
		punch_file.length() == 0 ||
		punch_template_files.size() == 0))
	{
		// Insert error about correct usage
		// But non-fatal unless doing doing a full run
		cout << "We are missing a required filename. Please check the requirements:" << endl;
		if (ecp_file.length() == 0)
		{
			cout << "ECP File" << endl;
		}
                if (ecp_template_file.length() == 0)
                {
                        cout << "ECP Template" << endl;
                }
                if (chm_file.length() == 0)
                {
                        cout << "CHM File" << endl;
                }
                if (punch_file.length() == 0)
                {
                        cout << "PUNCH File" << endl;
                }
                if (punch_template_files.size() == 0)
                {
                        cout << "PUNCH Template" << endl;
                }
		
		cout << endl;
		// Print out error message about inputs
		inputs();
		cout << endl;
		critical_error(!dry_run);
	}
	
	// We need something here to set defaults for values like
	// - qm output
	if (qm_output_file.length() == 0.0)
	{
		cout << "Using default GAMESS-UK output file: gamess1.out.1" << endl;
		qm_output_file = "gamess1.out.1";
	}
	// - gradient output
	if (gradient_output_file.length() == 0.0)
	{
		cout << "Using default gradient output file: gradient" << endl;
		gradient_output_file = "gradient";
	}
	// - output_folder
	if (output_folder.length() == 0.0)
	{
		if (!outputs_only)
		{
			cout << "Using default output folder: results" << endl;
		}
		output_folder = "results";
	}
        if (executable.length() == 0.0)
	{
		if (!outputs_only)
		{
			cout << "Using default executable: chemsh.x" << endl;
		}
		executable = "chemsh.x";
	}
        if (anion_species.length() == 0.0)
        {
                if (!outputs_only)
                {
                        cout << "Using default anionic species: O" << endl;
                }
                anion_species = "O";
        }
	// - chemshell counter
	if (chemshell_counter_max == 0)
	{
		if (!outputs_only)
		{
			cout << "Using default max Chemshell calculations: 5" << endl;
		}
		chemshell_counter_max = 5;
	}
	// - random seed
	if (random_seed == 0)
	{
		if (!outputs_only)
		{
			cout << "Using default Random seed: 1" << endl;
		}
		random_seed = 1;
	}
	// Check the number of processors in use
	if (processors == 0)
	{
		if (!outputs_only)
		{
			cout << "Total Number of Processors is not defined. Executing ChemShell in serial" << endl;
		}
		//critical_error(!outputs_only || dry_run);
	}
	if (processors_per_node == 0)
	{
		if (!outputs_only)
		{
			cout << "Number of Processors per node is not defined. Executing ChemShell in serial" << endl;
		}
		//critical_error(!outputs_only || dry_run);
	}
	// Defaults are set
	cout << endl;

        // Decide what type of QM calculation we are looking at
        string qm_type = qm_output_file.substr(0,6);

        if (cmpStr(qm_type,"nwchem"))
        {
		qm_program = new Nwchem();
        }
        else
        { 
        	qm_program = new Gamess_UK();
        }
	
	// Set anion offset
	qm_program->set_anion_offset(region_1_anion_offset);
        qm_program->set_anion_species(anion_species);

	if (!outputs_only)
	{
		// Firstly we will read in the ECP input file
		qm_program->set_ecp_template(read_in_lines(ecp_template_file,!dry_run));
		
		// Identify ECPs to be fitted.
		gaussian sg = qm_program->get_starting_gaussian_ecps();
		
		if (sg.values.size()> 0)
		{
			ecp_searcher->set_starting_gaussians(sg);
		}
		else
		{
			chemshell_counter_max = 1;
		}
	}
	else
	{
		chemshell_counter_max = 1;
	}
	
	// Next we will asign data from the punch template file
	punch.resize(punch_template_files.size());
	// Put in all the important information
	for (vector<string>::size_type a = 0; a < punch_template_files.size(); a++)
        {
		cout << "Punch: " << punch_template_files[a] << " (#" << a << ")" << endl;;
        	punch[a].set_punch_template(read_in_lines(punch_template_files[a]),anion_species);
		//cout << "Punch: " << punch_template_files[a] << " (#" << a << ")";
		punch[a].print_regions();

	        if (outputs_only)
        	{
                	// Print bond lengths on analysis
                	punch[a].print_bond_data();
                	// punch[a].print_region_1_species();
                }
        } 

	//punch.set_punch_template(read_in_lines(punch_template_file),anion_species);
	//punch.print_regions();
	
	// We need to resize ecps_history to account for multiple inputs
	ecps_history.resize(punch_template_files.size());
	// Assign new functions to this vector
	for (vector<History *>::size_type a = 0; a < ecps_history.size(); a++)
	{
		ecps_history[a] = new History();
		// At the same time we can set up the output log files
		string temp_output_log = log_output_file;
		temp_output_log += ".";
		temp_output_log += NumberToString(a);
		log_output_files.push_back(temp_output_log);
	}
	

	if (!outputs_only)
	{
		// Gather in all the different old outputs for restart
		for (vector<History *>::size_type i_history = 0; i_history < ecps_history.size(); i_history++)
        	{
			// Let's check if there are any old inputs we add to our history
                        vector<gaussian> temp = read_in_binary(log_output_files[i_history]+".restart",false);

			// If the read of the binary file fails then we copy from the text outputs
			if (temp.size() > 0)
                        {
                                // cout << "SUCCESS, temp.size = " << temp.size() << endl;
				ecps_history[i_history]->set_history(temp);
			}
			else
			{
                                // cout << "FAILURE, temp.size = " << temp.size() << endl;
				ecps_history[i_history]->insert_old_data(read_in_lines(log_output_files[i_history], false),read_in_lines(regions_output_file, false));
                                // cout << "SUCCESS?, temp.size = " << temp.size() << endl;
                        }

                        // We need to add something here to do the function recalculations.
                        // Easiest method is going to be.... separate routine in history.
			ecps_history[i_history]->recalc_function(func_calc,i_history,force_recalc);

			// Also need to remove the duds
			//if (remove_duds)
			//{
			//	ecps_history[i_history]->remove_duds();
			//}
			
			// Print Header
			string sentence = func_calc->get_header();
			outData.push_back(sentence);
			sentence = spacer;
			sentence += spacer;
			outData.push_back(sentence);

			// Print Targets
			sentence = func_calc->get_targets_header(i_history);
			outData.push_back(sentence);
			sentence = spacer;
			sentence += spacer;
			outData.push_back(sentence);

			write_out_lines(log_output_files[i_history],&outData,!dry_run);
                        write_out_binary(log_output_files[i_history]+".restart",ecps_history[i_history]->get_history(),!dry_run);
		}
		
		// Set up regions file
		string sentence = "Entry\t|\t\tLine\t\tType\t\tValue";
		outData.clear();
		outData.push_back(sentence); // Change the header
		sentence = spacer;
		outData.push_back(sentence);		
		write_out_lines(regions_output_file,&outData,!dry_run);

		// Rewrite all gathered information back to the log files
		for (vector<History *>::size_type i_history = 0; i_history < ecps_history.size(); i_history++)
		{

			// If there was anything in the old log file add it to the new one
			if (ecps_history[i_history]->size() > 0)
			{
				// Update results 
				update_log_file(log_output_files[i_history],ecps_history[i_history]->get_history(),!dry_run);
			
				// Update regions
				vector<int> v_history(ecps_history[i_history]->size(),0);
				update_regions_file(regions_output_file, ecps_history[i_history]->get_history(),v_history,!dry_run);

				// Check the current size of ecp_history; add to the total as each is saved locally
				current_index += ecps_history[i_history]->size();
			}
		}

		// Add one to current index so we are at a new value 
		current_index += 1;
	}
	
	// So this will loop until we get the step size small enough or we just run too many calculations
	while (!ecp_searcher->get_converged() &&
		   (chemshell_counter < chemshell_counter_max))
	{
		vector<gaussian> ecps_to_test = ecp_searcher->get_ecps_to_test();
		vector<gaussian> ecps_tested(ecps_to_test.size());

		// Let's duplicate these structures for the varying number of punch templates we are testing
		vector< vector<gaussian> > ecps_to_test_vector;
		vector< vector<gaussian> > ecps_tested_vector;
		ecps_to_test_vector.resize(punch.size());
		ecps_tested_vector.resize(punch.size());

		// Assign values
		for (vector<gaussian>::size_type i_ecps = 0; i_ecps < ecps_to_test_vector.size(); i_ecps++)
		{
			ecps_to_test_vector[i_ecps] = ecps_to_test;
			ecps_tested_vector[i_ecps].resize(ecps_to_test.size());
		}
		
		// Check history. If they've already been tested, put them in the results section
		// For some reason I couldn't get find() to work here. Probably needs some attention in the long term
		
		for (vector<gaussian>::size_type i_punch = 0; i_punch < punch.size(); i_punch++)
	        {

	                // We are going to explicitly define a vector to show pulled from history
        	        // With values set to 0 for false, 1 for true
 	                vector<int> from_history(ecps_to_test_vector[i_punch].size(),0);

        	        // Define a counter for failures
	                unsigned int failures = 0;
	
			// Loop through our ecps to test
			for (vector<gaussian>::size_type a = 0; a < ecps_to_test_vector[i_punch].size(); a++)
			{

				// Loop through the history file. We'll go in reverse order as most recent ecps will be at the back
				int exists = ecps_history[i_punch]->check_history(ecps_to_test_vector[i_punch][a]);
			
				if (exists != -1)
				{

					// ecps_history[b] contains ecps_to_test[a]
					ecps_tested_vector[i_punch][a] = ecps_history[i_punch]->get(exists);

					// This is our note that we've pulled this from history
					from_history[a] = 1;

					// We need to make an escape route if we are doing a dry_run
					if (dry_run)
					{
						chemshell_counter++;
						print_chemshell_message(ecps_tested_vector[i_punch][a],chemshell_counter,current_index);
						cout << "We are in dry run, but have picked up a value from history. Increasing chemshell counter." << endl;
					}
				}
			}

			// This would be the start of our loop function to test current ecps
			for (vector<gaussian>::size_type a = 0; a != ecps_to_test_vector[i_punch].size(); a++)
			{
				// Check that we have not copied the result from history
				if (from_history[a] != 1)
				{
					// This is going to be our carrier, and written directly to ecps_tested at the end
					gaussian g = ecps_to_test_vector[i_punch][a]; 

					// Set marker to see if calculation runs ok
					g.failed = false;

					chemshell_counter++;
					cout << spacer << endl;
					print_chemshell_message(g,chemshell_counter,current_index);

					// Create folder name for moving around results
					string current_counter = "";
					NumberToString(current_index,current_counter);
	                		string current_folder = output_folder;
        	        		current_folder += "_";
					current_folder += current_counter;

					// Write ECP and Punch file for this run 
					if (!outputs_only)
					{
						outData = qm_program->get_ecp_template(g);
						write_out_lines(ecp_file,&outData,!dry_run);
				
						outData = punch[i_punch].get_punch_template();
						write_out_lines(punch_file,&outData,!dry_run);
					}

					// Run Chemshell QM/MM calculator, using predefined setup.
					// command_line="aprun -n $NPROC -N $NTASK chemsh.x " + chm_file; // We should softcode this; Done a bit for now
					// Having processors and processors_per_node in the code means we can parallelise easily
					if (!dry_run && !outputs_only)
					{
						string command_line = "";
						if (processors > 0 && processors_per_node > 0)
						{
							command_line = "aprun -n " + NumberToString(processors) + " -N " + NumberToString(processors_per_node) + " " + executable + " " + chm_file;
						}
						else
						{
							command_line = executable + " " + chm_file;
						}
						cout << "Running Chemshell" << endl;
						cout << spacer << endl;
						system(command_line.c_str());
						cout << spacer << endl;
						cout << endl;
					}

					// Read in data from output once finished. We need a check here in case something hasn't converged.
					// Check if punch_output is different size to punch_template. This will only be performed after the first calculation
					if (outputs_only)
					{

						// In this case we should need to reread, as the calculation-used punch file will be the template
						punch_output_check = true;
					}
					else if (!punch_output_check)
					{
						// We are going to reread the punch output file and check nothing has changed
						// In a defected system this will have changed.
						punch[i_punch].compare(read_in_lines(punch_file, false));
						
						// Let's check if the punch output is different from the original
						// Note! If the system in question has a defect this will introduce errors.
						punch_output_check = true;
					}
			
					// Read in gradients to see if they are defined
					g.regions = digest_gradients(read_in_lines(gradient_output_file, outputs_only),punch[i_punch].get_total_centres(),
  								     punch[i_punch].get_centre_regions_total(),punch[i_punch].get_centre_regions(),absolute_gradients);
				
					if ((g.regions[0].gnorm_max == 0) &&
						(g.regions[1].gnorm_max == 0) &&
						(g.regions[2].gnorm_max == 0) &&
						(g.regions[3].gnorm_max == 0) &&
						(g.regions[4].gnorm_max == 0))
					{
						g.failed = true;
					}
			
					// Now to calculate electronic information from DFT output
					g = qm_program->digest_electronic(read_in_lines(qm_output_file, outputs_only), g, punch[i_punch].get_region_1_anions(), punch[i_punch].get_region_1_species());
				
					// Copy output to temporary location in case we want to check it.
					// This should be optional otherwise we'll end up with lots of datafiles.
					if (!dry_run && !outputs_only)
					{
						// Now lets create the command line to move files around
						string command_line = " mkdir ";
						command_line += current_folder;
						command_line += " ; cp ";
						command_line += qm_type;
						command_line += "* ";
						command_line += current_folder;
						command_line += " ; mv gulp* hybrid* ";
						command_line += ecp_file;
						command_line += " ";
						command_line += gradient_output_file;
						command_line += " *nergy *xyz ";
						command_line += current_folder;
						// cout << command_line << endl;
						system(command_line.c_str());
					}	
			
					// Calculate function value
					if (!g.failed && g.function != 888888)
					{

						// Create a temporary vector to hold the DMA spreads
						vector<double> temp(g.dma_spread.size());
						for (vector<double>::size_type i = 0; i < temp.size(); i++)
						{
							temp[i] = g.dma_spread[i].spread;
						}

						// Calculate function
						g.function = func_calc->calculate_function(g.regions[0].gnorm, g.regions[1].gnorm, g.regions[2].gnorm, g.orbital_spread[0].spread, g.HOMO_value, g.LUMO_value, temp, i_punch);

						// Print the function value to screen
						cout << "Function value : " << g.function << endl;
						cout << endl;
					}
					// If not set the function value arbitrarily high
					else 
					{
						g.function = 888888;
						failures++;
					}
				
					// Add index to reflect this calculation
					g.index = current_index;
					current_index++;
				
					// Copy this complete ECP to ECP_tested
					ecps_tested_vector[i_punch][a] = g;
				}
				else
				{
					if (ecps_tested_vector[i_punch][a].failed || ecps_tested_vector[i_punch][a].function == 888888)
					{
						failures++;
					}
				}
			}
		
			// At this point we have done all the ECPs_to_test, and must now compare
			// We need to see whether a compund function fits the required goal.
			// Perhaps Spearman's Rank Correlation Coefficient: http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient
			// And except uphill movements on a probablistic nature (i.e. MC)
		
			#define DELTA (0.000000001)
			// We've calculated all the function values, so now we just need to rank all of them
			// We'll do this in a simple double loop, as it should be a quick calculation
			for (vector<gaussian>::size_type i = 0; i < ecps_tested_vector[i_punch].size(); i++)
			{
				// Check if calculation failed
				if (ecps_tested_vector[i_punch][i].failed)
				{
					// If so it's rank is maximum possible
					ecps_tested_vector[i_punch][i].rank = ecps_tested_vector[i_punch].size();
				}
				else
				{
					ecps_tested_vector[i_punch][i].rank = 1;
					for (vector<gaussian>::size_type j = 0; j < ecps_tested_vector[i_punch].size(); j++)
					{
						// Originally this had an a != i clause, but that shouldn't make a difference as it >, not >= clause
						if ((ecps_tested_vector[i_punch][i].function - ecps_tested_vector[i_punch][j].function) > DELTA)
						{
							// Increment rank number if the function of another gaussian
							// is lower, and we are not looking at the same gaussian
							ecps_tested_vector[i_punch][i].rank++;
						}
					}
				}
			}
			#undef DELTA
		
			if (!outputs_only)
			{
				// Now we have everything ranked lets append to the output file
				update_log_file(log_output_files[i_punch],ecps_tested_vector[i_punch],!dry_run);
			
				// This needs to be done on every loop, as otherwise restart won't work.
				update_regions_file(regions_output_file,ecps_tested_vector[i_punch],from_history,!dry_run);

				// Print to screen to check failures counter
				cout << failures << " of the " << ecps_tested.size() << " " << qm_program->type() << " calculations have failed." << endl;
			
				if (!dry_run && (failures == ecps_tested_vector[i_punch].size()))
	        	        // All calculations have failed so we need to exit otherwise we're wasting CPU time
				{
					string error = "All of the "; 
	                                error += qm_program->type();
        	                        error += " calculations have failed! Quiting.";
					cout << error << endl;
					critical_error(true);
				}
			
				// Save any results from the post calculation strip down. This should be placed in "history" file so we can look at previous results
				// Due to the nature of this search, the history will be unordered. But that shouldn't be to big a problem as we won't be running lots of calculations
				ecps_history[i_punch]->append(ecps_tested_vector[i_punch],from_history);

                                //Update the binary file as well with new history (Append does not work nicely at present...)
				if (ecps_history[i_punch]->get_number_of_entries_last_added() > 0)
				{
                                	write_out_binary(log_output_files[i_punch]+".restart",ecps_history[i_punch]->get_history(),!dry_run);
				}
			}
		}

		// Need to put the ecp_searcher command in here. Then we have finished and need to test
		// We need to firstly work out which is the number_one_ranked overall!
		if (!outputs_only)
		{
			vector<int> summed_ranks(ecps_tested.size(),0);
			vector<double> summed_functions(ecps_tested.size(),0.0);

			for (vector<gaussian>::size_type i_punch = 0; i_punch < punch.size(); i_punch++)
                	{
				for (vector<gaussian>::size_type i = 0; i < ecps_tested_vector[i_punch].size(); i++)
                        	{
					summed_ranks[i] += ecps_tested_vector[i_punch][i].rank;
					summed_functions[i] += ecps_tested_vector[i_punch][i].function;
				}
			}
			
			// Counter for the best result
			int number_one_ranked = 0;
			
			// We've got all the summed ranks and functions, now to find the best.
			for (vector<gaussian>::size_type i = 0; i < ecps_tested.size(); i++)
			{
				//cout << i << " " << number_one_ranked << endl;
				cout << endl;
				cout << "ECP Tested       : " << i << endl;
				cout << "Combined Rank    : " << summed_ranks[i] << endl;
				cout << "Combined Function: " << summed_functions[i] << endl;
				
				// I need to make a decision which of these is better for the ranking procedure
				// Or if I should combine them. For now we are using the summed functions
				if (summed_functions[i] < summed_functions[number_one_ranked])
				{
					number_one_ranked = i;
				}

				//if (summed_ranks[i] < summed_ranks[number_one_ranked])
				//{
				//	number_one_ranked = i;
				//}
				
				//cout << i << " " << number_one_ranked << endl;
			}
			
			cout << endl;
			cout << "Number one ranked ECP: " << number_one_ranked << endl;
			cout << endl;
			
			ecp_searcher->set_ecps_tested(ecps_tested_vector[0],number_one_ranked);
		}
	}
	
	// Confirm we've converged
	if (!outputs_only)
	{
		if (chemshell_counter >= chemshell_counter_max)
		{
			cout << endl;
			cout << "Max Chemshell calculations reached" << endl;
			cout << endl;
		}
		else 
		{
			cout << endl;
			cout << "Convergence achieved" << endl;
			cout << endl;	
		}
		
	}
	
	return EXIT_SUCCESS;
}
