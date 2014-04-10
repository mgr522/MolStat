/**
 * \file model_interface.cc
 * \brief Implementation for the interface for various tight-binding models
 *        and their transmission functions/conductances.
 *
 * \author Matthew G.\ Reuter
 * \date April 2014
 */

#include "model_interface.h"

ConductanceModel::ConductanceModel(shared_ptr<const RandomDistribution> eta)
	: dist_eta(eta) {}
