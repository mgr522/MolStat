/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file main-simulator-inputparse.cc
 * \brief Definition of the SimulatorInputParse class for parsing the simulator
 *    input deck.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#include "main-simulator.h"
#include <iomanip>
#include <general/simulator_tools/simulator_exceptions.h>
#include <electron_transport/simulator_models/transport_simulate_module.h>

using namespace std;

void SimulatorInputParse::printError(std::ostream &output, std::size_t lineno,
	std::string message)
{
	output << "Error on line " << setw(2) << lineno << ": " << message << endl;
}

void SimulatorInputParse::readInput(std::istream &input)
{
	std::size_t lineno{ 1 }; // line number

	// process the input deck
	// input evaluates to "false" when we hit EOF
	while(input)
	{
		string line;
		getline(input, line);

		// tokenize the string
		molstat::TokenContainer tokens = molstat::tokenize(line);
		if(tokens.size() == 0) // empty line
			continue;

		// the first token is the command name, pop it off and then process the
		// rest of the tokens
		string command { molstat::to_lower(tokens.front()) };
		tokens.pop();

		// go through the supported commands
		if(command == "bin" || command == "bin_x" || command == "bin_y")
		{
			printError(cout, lineno, "bin handler not implemented.");
		}
		else if(command == "model")
		{
			// the total lineno needs to be passed to the model reader
			const size_t mylineno{ lineno };

			// enter the model reader to process input lines until the
			// appropriate "endmodel" command is found
			ModelInformation model = readModel(input, ++lineno);

			// make sure there is a name (model type) specified
			if(tokens.size() == 0)
			{
				printError(cout, mylineno, "No model type specified.");
			}
			else
			{
				// store the name of the model type
				model.name = molstat::to_lower(tokens.front());

				// set this to the parser's top model
				top_model = move(model);
			}
		}
		else if(command == "observable" || command == "observable_x" ||
			command == "observable_y")
		{
			// make sure we have a name for the observable
			if(tokens.size() == 0)
			{
				printError(cout, lineno, "No observable specified.");
			}
			else {
				// store the name of the observable for later
				if(command == "observable_y")
					used_observables.emplace(1, tokens.front());
				else
					used_observables.emplace(0, tokens.front());
			}
		}
		else if(command == "output")
		{
			if(tokens.size() == 0)
			{
				printError(cout, lineno, "No output file name specified.");
			}
			else
			{
				histfilename = tokens.front();
			}
		}
		else if(command == "trials")
		{
			printError(cout, lineno, "trials handler not implemented.");
		}
		else
		{
			printError(cout, lineno, "Unknown command: \"" + command + "\".");
		}

		// move to the next line
		++lineno;
	}

#if 0
	// construct the Simulator
	unique_ptr<molstat::Simulator> sim{ new molstat::Simulator(model) };

	// load the observables
	for(const auto obs : used_observables)
	{
		try
		{
			sim->setObservable(obs.first, obs.second);
		}
		catch(const out_of_range &e)
		{
			throw runtime_error(
				"Non-contiguous observable numbering is not allowed.");
		}
		catch(const molstat::IncompatibleObservable &e)
		{
			throw runtime_error("The model is incompatible with observable "
				+ to_string(obs.first) + '.');
		}
	}
#endif
}

SimulatorInputParse::ModelInformation SimulatorInputParse::readModel(
	std::istream &input, std::size_t &lineno)
{
	ModelInformation ret;

	while(input)
	{
		string line;
		getline(input, line);

		// tokenize the string
		molstat::TokenContainer tokens = molstat::tokenize(line);
		if(tokens.size() == 0) // empty line
			continue;

		// the first token is the command name, pop it off and then process the
		// rest of the tokens
		string command { molstat::to_lower(tokens.front()) };
		tokens.pop();

		// go through the supported commands
		if(command == "endmodel")
		{
			// we're done here
			return ret;
		}
		else if(command == "model")
		{
			// the total lineno needs to be passed to the (sub)model reader
			const size_t mylineno{ lineno };

			// enter the (sub)model reader to process input lines until the
			// appropriate "endmodel" command is found
			ModelInformation model = readModel(input, ++lineno);

			// make sure there is a name (model type) specified
			if(tokens.size() == 0)
			{
				printError(cout, mylineno, "No submodel type specified.");
			}
			else
			{
				// store the name of the model type
				model.name = molstat::to_lower(tokens.front());

				// set this to the parser's top model
				ret.submodels.emplace_back( move(model) );
			}
		}
		else if(command == "distribution")
		{
			// make sure there are tokens, if so, push the tokens
			if(tokens.size() == 0)
			{
				printError(cout, lineno, "No distribution type specified.");
			}
			else
			{
				ret.dists.emplace_back(tokens);
			}
		}
		else
		{
			printError(cout, lineno,
				"Unknown model command: \"" + command + "\".");
		}

		// move to the next line
		++lineno;
	}

	// if we're here, we hit EOF before finding the appropriate endmodel command
	throw runtime_error("Missing \"endmodel\" command.");
}

std::unique_ptr<molstat::Simulator> SimulatorInputParse::createSimulator()
{
	// the map of model instantiators
	map<string, molstat::SimulateModelFactoryFunction> models;

	// the map of observable indexes
	map<string, molstat::ObservableIndex> observables;

	// load model and observable names
	// the general syntax for loading models is
	// models.emplace( to_lower(name),
	//                 molstat::GetSimulateModelFactory<model_type>() );
	//
	// and likewise for observables,
	// observables.emplace( to_lower(name),
	//                      GetObservableIndex<observable_type>() );
	molstat::transport::load_models(models);
	molstat::transport::load_observables(observables);

	// make the model

	// set the observables

	return nullptr;
}

#if 0
std::shared_ptr<molstat::SimulateModel> processModel(
	std::queue<std::string> &&tokens,
	std::istream &input,
	const std::map<std::string,
	               molstat::SimulateModelFactoryFunction> &models)
{
	if(tokens.size() == 0)
		throw runtime_error("Model name not specified.");

	// get the factory
	molstat::SimulateModelFactoryFunction ffunc;
	try
	{
		// the next token is the name of the model
		// use "at" to check existence of the key in the models map
		ffunc = models.at(molstat::to_lower(tokens.front()));
	}
	catch(const out_of_range &e)
	{
		// model not found
		throw runtime_error("Model \"" + tokens.front() + "\" not found.");
	}
	molstat::SimulateModelFactory factory{ ffunc() };

	// start processing information from input...
	string command{ "not endmodel" };
	while(input && command != "endmodel")
	{
		string line;
		getline(input, line);

		// tokenize the string
		tokens = molstat::tokenize(line);
		if(tokens.size() == 0) // empty line
			continue;

		// the first token is the command name, pop it off and then process the
		// rest of the tokens
		command = molstat::to_lower(tokens.front());
		tokens.pop();

		// process the two available commands (get a submodel or set a
		// distribution)
		if(command == "distribution")
		{
			// get the name
			if(tokens.size() == 0) // no name specified
			{
				cout << "   No name specified on \"distribution\"." << endl;
			}
			else
			{
				const string name{ molstat::to_lower(tokens.front()) };
				tokens.pop();
				shared_ptr<const molstat::RandomDistribution> dist{ nullptr };

				// process the remaining tokens as a distribution
				try
				{
					dist = molstat::RandomDistributionFactory(move(tokens));

					// set the distribution
					factory.setDistribution(name, dist);
				}
				catch(const invalid_argument &e)
				{
					// print out a message, but don't error out
					cout << "   Error constructing random distribution: " <<
						e.what() << endl;
				}
			}
		}
		else if(command == "endmodel")
		{
			// don't need to do anything; we'll exit the input reading loop after
			// this cycle. we just don't want to print the error message.
		}
		else if(command == "model")
		{
			// construct a submodel... recursion!
			shared_ptr<molstat::SimulateModel> model{ nullptr };
			try
			{
				model = processModel(move(tokens), input, models);

				// add the submodel
				factory.addSubmodel(model);
			}
			catch(const molstat::NotCompositeSimulateModel &e)
			{
				cout << "Error setting submodel: invalid model type." << endl;
			}
			catch(const molstat::IncompatibleObservable &e)
			{
				cout << "Error setting submodel: invalid model type." << endl;
			}
			catch(const runtime_error &e)
			{
				cout << "Error constructing submodel: " << e.what() << endl;
			}
		}
		else
		{
			cout << "   Unrecognized model command: \"" + command + "\"." << endl;
		}
	}

	// make sure the last command was endmodel; that is, we didn't hit EOF
	if(command != "endmodel")
		throw runtime_error("EOF encountered in \"model\" block.");

	return factory.getModel(); // let any exceptions pass up
}

molstat::ObservableIndex processObservable(
	std::queue<std::string> &&tokens,
	const std::map<std::string,
	          molstat::ObservableIndex> &observables)
{
	// make sure there is a token left, otherwise throw a runtime_error
	if(tokens.size() == 0)
		throw runtime_error("Observable name not specified.");

	molstat::ObservableIndex obsindex{ typeid(void*) };

	// look up the specified observable
	try
	{
		// the next token is the name of the observable
		// use "at" to check existence of the key in the observables map
		obsindex = observables.at(molstat::to_lower(tokens.front()));
	}
	catch(const out_of_range &e)
	{
		// observable not found
		throw runtime_error("Observable \"" + tokens.front() + "\" not found.");
	}

	return obsindex;
}
#endif

std::string SimulatorInputParse::ModelInformation::to_string() const
{
	// first put in the name
	string ret{ name };
	ret += "\n   Distributions:";

	// load in the distributions
	for(auto dist : dists)
	{
		ret += "\n      " + dist.front();
	}

	// submodel information
	for(auto submodel : submodels)
	{
		ret += "\n   Submodel type: ";

		// get the submodel info; indent it
		string submodel_string{ submodel.to_string() };
		size_t start_pos{ 0 }, newl_pos{ 0 };
		while((newl_pos = submodel_string.find("\n", start_pos)) != string::npos)
		{
			// move up to (and including) the newline to the return string
			ret += submodel_string.substr(start_pos, newl_pos - start_pos + 1);
			ret += "   "; // the indent
			start_pos = newl_pos + 1;
		}

		// copy the last segment
		ret += submodel_string.substr(start_pos, string::npos);
	}

	return ret;
}

void SimulatorInputParse::printState(std::ostream &output) const
{
	output << "--------------------------------------------------\n" <<
		"State of input parser.\n\n";

	output << "Model type: " << top_model.to_string() << "\n\n";

	output << "Observables:\n";
	for(auto obspair : used_observables)
	{
		output << obspair.first << " -> " << obspair.second << '\n';
	}
	output << '\n';

	output << "Histogram Output File: " << histfilename << "\n";

	output << "--------------------------------------------------" << endl;
}