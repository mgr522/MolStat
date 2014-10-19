/* This file is a part of MolStat, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.
   MolStat (c) 2014, Northwestern University. */
/**
 * \file simulate_model_factory.h
 * \brief Declares the molstat::SimulateModelFactory class for constructing
 *    molstat::SimulateModel objects.
 *
 * \author Matthew G.\ Reuter
 * \date October 2014
 */

#ifndef __simulate_model_factory_h__
#define __simulate_model_factory_h__

#include <memory>
#include <string>
#include <vector>
#include <set>
#include <general/string_tools.h>
#include "simulate_model.h"

namespace molstat {

// forward declarations ... we only need pointers in this header file
class RandomDistribution;

class SimulateModelFactory {
private:
	SimulateModelFactory() = default;

	/**
	 * \brief Pointer to the model being constructed.
	 */
	std::shared_ptr<SimulateModel> model;

	/**
	 * \brief The set of distribution names for the model that still need to
	 *    be specified.
	 */
	std::set<std::string> remaining_names;

	/**
	 * \brief Cache of the names of distributions required by the model.
	 */
	std::vector<std::string> model_names;

public:
	SimulateModelFactory(SimulateModelFactory &&) = default;
	SimulateModelFactory &operator=(SimulateModelFactory &&) = default;

	/**
	 * \brief Creates a molstat::SimulateModelFactory with the underlying model
	 *    of type T.
	 *
	 * \tparam T The type of molstat::SimulateModel to build with this factory.
	 * \return The factory.
	 */
	template<typename T>
	static SimulateModelFactory makeFactory();

	/**
	 * \brief Adds a random distribution to the model.
	 *
	 * \param[in] name The name of the distribution being added.
	 * \param[in] dist The distribution being added.
	 * \return The factory.
	 */
	 SimulateModelFactory &setDistribution(std::string name,
	 	std::shared_ptr<const RandomDistribution> dist);

	 /**
	  * \brief Returns the constructed model.
	  *
	  * Some runtime error checking, such as making sure all distributions are
	  * specified is performed here.
	  *
	  * \throw molstat::MissingDistribution if one of the required distributions
	  *    has not been specified.
	  *
	  * \return Pointer to the constructed model.
	  */
	 std::shared_ptr<SimulateModel> getModel();
};

// templated definitions
template<typename T>
SimulateModelFactory SimulateModelFactory::makeFactory() {
	using namespace std;

	SimulateModelFactory factory;

	factory.model = make_shared<T>();

	// get the names of required distributions
	factory.model_names = factory.model->get_names();

	// convert the ordered vector to just a set
	for(const std::string &iter : factory.model_names)
		factory.remaining_names.emplace(to_lower(iter));

	// set the size of the model's vector of distributions
	factory.model->dists.resize(factory.model->get_num_parameters());

	return factory;
}

} // namespace MolStat

#endif