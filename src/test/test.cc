#include <cstdio>
#include "../general/simulator_tools/simulate_model_interface.h"

#include <general/random_distributions/constant.h>

using namespace std;
using namespace molstat;

template<size_t MPs>
class ZBConductance : Observable<MPs> {
public:
	virtual ~ZBConductance() = default;
	virtual double ZBG(const array<double, MPs> &vals) const = 0;

	virtual double operator()(const array<double, MPs> &vals) const override final {
		return ZBG(vals);
	}
};

template<size_t MPs>
class DiffConductance : Observable<MPs> {
public:
	virtual ~DiffConductance() = default;
	virtual double DiffG(const array<double, MPs> &vals) const = 0;

	virtual double operator()(const array<double, MPs> &vals) const override final {
		return DiffG(vals);
	}
};

class Tester : public SimulateModel<1>, public ZBConductance<1>, public DiffConductance<1> {
public:
	Tester(const map<string, shared_ptr<RandomDistribution>> &avail)
		: SimulateModel<1>(avail, {{"a"}}) {

	}

	double ZBG(const array<double, 1> &vals) const override {
		return vals[0];
	}

	double DiffG(const array<double, 1> &vals) const override {
		return 2.;
	}
};

int main(int argc, char **argv) {

	map<string, shared_ptr<RandomDistribution>> dists;
	dists["a"] = make_shared<ConstantDistribution>(5.);
	gsl_rng_ptr r(nullptr, &gsl_rng_free);

	shared_ptr<SimulateObservables<3, 1>> obs = SimulateObservables<3, 1>::Factory<Tester>(dists);

	obs->setObservable<ZBConductance>(0);
	obs->setObservable<DiffConductance>(2);

	array<double, 3> vals = obs->simulate(r);

	printf("%f %f %f\n", vals[0], vals[1], vals[2]);

	return 0;
}