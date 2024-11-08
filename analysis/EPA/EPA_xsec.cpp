#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <map>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_interp.h>

//*****************************************************************
// CONSTANTS
//*****************************************************************
const int NDIM = 3; // Dimensionality of integral

const double pi = 3.1415926535897;

const double GF = 0.0000116637; // Fermi constant [1/GeV^2]
const double aem = 1/137.036; // em coupling [1]

const double mtau = 1.778; // Tau lepton mass [GeV]
const double mtau2 = mtau*mtau; // Tau lepton mass squared [GeV^2]

const double MArgon = 39.95*0.9315; // Mass of argon nucleus
const double Mproton = 0.938272; // Mass of proton
const double Mneutron = 0.939565; // Mass of neutron

const int Z_Ar = 18;

//*****************************************************************
// RELEVANT KINEMATIC DOT PRODUCTS
//*****************************************************************

double p1q(double s){
    double result = s / 2;
    return result;
}

double p1p3(double cos34, double cosCM, double l, double s) {
    double result = (mtau2*((s-l)*cosCM+l+s) - l*((1-mtau2/l)*(cos34*(cosCM*((sqrt(l)-sqrt(s))*(sqrt(l)-sqrt(s))*cosCM-l+s)+2*sqrt(l*s)) + sqrt(1-cos34*cos34)*sqrt(1-cosCM*cosCM)*((sqrt(l)-sqrt(s))*(sqrt(l)-sqrt(s))*cosCM-l+s)) + (s-l)*cosCM + l + s))/(8*l);
    return result;
}

double p1p4(double cos34, double cosCM, double l, double s) {
    double result = (l*((1-mtau2/l)*(cos34*(cosCM*((sqrt(l)-sqrt(s))*(sqrt(l)-sqrt(s))*cosCM-l+s) + 2*sqrt(l*s)) + sqrt(1-cos34*cos34)*sqrt(1-cosCM*cosCM)*((sqrt(l)-sqrt(s))*(sqrt(l)-sqrt(s))*cosCM - l + s)) + (s-l)*cosCM + l + s) - mtau2*((s-l)*cosCM + l + s))/(8*l);
    return result;
}

double p2q(double cosCM, double l, double s) {
    double result = (s-l)*(cosCM+1)/4;
    return result;
}

double p2p3(double cos34, double cosCM, double l, double s) {
    double result = ((l-s)*((l-mtau2)*(sqrt(1-cos34*cos34)*sqrt(1-cosCM*cosCM)+cos34*cosCM)-l-mtau2))/(4*l);
    return result;
}

double p2p4(double cos34, double cosCM, double l, double s) {
    double result = -((l-s)*(l-mtau2)*(sqrt(1-cos34*cos34)*sqrt(1-cosCM*cosCM) + cos34*cosCM + 1))/(4*l);
    return result;
}

double p3q(double cos34, double cosCM, double l, double s) {
    double result = (mtau2*((l-s)*cosCM+l+s) + (l-mtau2)*(cos34*(cosCM*((sqrt(l)-sqrt(s))*(sqrt(l)-sqrt(s))*cosCM+l-s)+2*sqrt(l*s)) + sqrt(1-cos34*cos34)*sqrt(1-cosCM*cosCM)*((sqrt(l)-sqrt(s))*(sqrt(l)-sqrt(s))*cosCM+l-s)) + l*(l-s)*cosCM + l*l + l*s)/(8*l);
    return result;
}

double p3p4(double l) {
    double result = (l-mtau2)/2;
    return result;
}

double p4q(double cos34, double cosCM, double l, double s) {
    double result = -(mtau2*((l-s)*cosCM+l+s) + (l-mtau2)*(cos34*(cosCM*((sqrt(l)-sqrt(s))*(sqrt(l)-sqrt(s))*cosCM + l - s) + 2*sqrt(l*s)) + sqrt(1-cos34*cos34)*sqrt(1-cosCM*cosCM)*((sqrt(l)-sqrt(s))*(sqrt(l)-sqrt(s))*cosCM + l - s)) - l*((l-s)*cosCM + l + s))/(8*l);
    return result;
}

//*****************************************************************
// AMPLITUDE TERMS
//*****************************************************************

double amp1(double p1p3_s, double p2p4_s, double p3q_s) {
    double result = 128*4*pi*aem*GF*GF*mtau2*p1p3_s*p2p4_s/(p3q_s*p3q_s);
    return result;
}

double amp2(double p1q_s, double p2p4_s, double p3q_s) {
    double result = -128*4*pi*aem*GF*GF*mtau2*p1q_s*p2p4_s/(p3q_s*p3q_s);
    return result;
}

double amp3(double p1p3_s, double p2p4_s, double p3q_s) {
    double result = 128*4*pi*aem*GF*GF*p1p3_s*p2p4_s/p3q_s;
    return result;
}

double amp4(double p1q_s, double p2p4_s, double p3q_s) {
    double result = -128*4*pi*aem*GF*GF*p1q_s*p2p4_s/p3q_s;
    return result;
}

double amp5(double p1p3_s, double p2p4_s, double p4q_s) {
    double result = 128*4*pi*aem*GF*GF*p1p3_s*p2p4_s/p4q_s;
    return result;
}

double amp6(double p1p3_s, double p2q_s, double p4q_s) {
    double result = -128*4*pi*aem*GF*GF*p1p3_s*p2q_s/p4q_s;
    return result;
}

double amp7(double p1p3_s, double p2p4_s, double p3p4_s, double p3q_s, double p4q_s) {
    double result = -256*4*pi*aem*GF*GF*p1p3_s*p2p4_s*p3p4_s/(p3q_s*p4q_s);
    return result;
}

double amp8(double p1q_s, double p2p4_s, double p3p4_s, double p3q_s, double p4q_s) {
    double result = 128*4*pi*aem*GF*GF*p1q_s*p2p4_s*p3p4_s/(p3q_s*p4q_s);
    return result;
}

double amp9(double p1p3_s, double p2q_s, double p3p4_s, double p3q_s, double p4q_s) {
    double result = 128*4*pi*aem*GF*GF*p1p3_s*p2q_s*p3p4_s/(p3q_s*p4q_s);
    return result;
}

double amp10(double p1p3_s, double p2p3_s, double p3q_s) {
    double result = -128*4*pi*aem*GF*GF*p1p3_s*p2p3_s/p3q_s;
    return result;
}

double amp11(double p1p4_s, double p2p4_s, double p4q_s) {
    double result = -128*4*pi*aem*GF*GF*p1p4_s*p2p4_s/p4q_s;
    return result;
}

//*****************************************************************
// FUNCTIONS
//*****************************************************************

// Final amplitude
double final_amplitude(double *k, size_t dim, void *params) {
    (void)(dim); /* avoid unused parameter warning */

    double s = *(double *)params;

    double p1q_s = p1q(s);
    double p1p3_s = p1p3(k[0], k[1], k[2], s);
    double p1p4_s = p1p4(k[0], k[1], k[2], s);
    double p2q_s = p2q(k[1], k[2], s);
    double p2p3_s = p2p3(k[0], k[1], k[2], s);
    double p2p4_s = p2p4(k[0], k[1], k[2], s);
    double p3q_s = p3q(k[0], k[1], k[2], s);
    double p3p4_s = p3p4(k[2]);
    double p4q_s = p4q(k[0], k[1], k[2], s);

    double amp1_s = amp1(p1p3_s, p2p4_s, p3q_s);
    double amp2_s = amp2(p1q_s, p2p4_s, p3q_s);
    double amp3_s = amp3(p1p3_s, p2p4_s, p3q_s);
    double amp4_s = amp4(p1q_s, p2p4_s, p3q_s);
    double amp5_s = amp5(p1p3_s, p2p4_s, p4q_s);
    double amp6_s = amp6(p1p3_s, p2q_s, p4q_s);
    double amp7_s = amp7(p1p3_s, p2p4_s, p3p4_s, p3q_s, p4q_s);
    double amp8_s = amp8(p1q_s, p2p4_s, p3p4_s, p3q_s, p4q_s);
    double amp9_s = amp9(p1p3_s, p2q_s, p3p4_s, p3q_s, p4q_s);
    double amp10_s = amp10(p1p3_s, p2p3_s, p3q_s);
    double amp11_s = amp11(p1p4_s, p2p4_s, p4q_s);

    double final_amp = amp1_s + amp2_s + amp3_s + amp4_s + amp5_s + amp6_s + amp7_s + amp8_s + amp9_s + amp10_s + amp11_s;
    return final_amp;
}

// Form factor
// Structure to hold interpolation parameters
struct InterpolationParams {
    gsl_interp *interp;
    gsl_interp_accel *acc;
    std::vector<double> x;
    std::vector<double> y;
};

// Function to evaluate the interpolated function
double interpolated_function(double x, void *params) {
    InterpolationParams *p = (InterpolationParams *)params;
    if (x < p->x.front() || x > p->x.back()) {
        // std::cerr << "Interpolation error: x = " << x << " is out of bounds." << std::endl;
        return 0.0; // Or handle error appropriately
    }
    return gsl_interp_eval(p->interp, p->x.data(), p->y.data(), x, p->acc);
}

// Function to read CSV file
void read_csv_and_process(const std::string& filename, std::vector<double>& Q2, std::vector<double>& integrand) {
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        double q, FF;
        char comma;

        ss >> q >> comma >> FF;
	double q2 = q*q;
	double intgrd = FF*FF/(q2);

	if (!Q2.empty() and (Q2.back() == q2)) {
	    integrand.back() = intgrd;
	    continue;
	}
	Q2.push_back(q2);
	integrand.push_back(intgrd);
    }
}


// Function to perform form factor integration using GSL
double integrate_form_factor(double lower_limit, double upper_limit, InterpolationParams& params) {

    // Check if integration has been done before
//    auto key = std::make_pair(lower_limit, upper_limit);
//    auto it = memo.find(key);
//    if (it != memo.end()) {
//        return it->second;
//    }


    if (upper_limit == 0.0 || upper_limit < lower_limit) {
	std::cout << "Form factor integral for the following integration limits will be set to 0.\n";
	std::cout << "Lower limit = " << lower_limit << " GeV.\n";
	std::cout << "Upper limit = " << upper_limit << " GeV.\n";
        return 0.0;
    }

    else {
	double result = gsl_interp_eval_integ(params.interp, params.x.data(), params.y.data(), lower_limit, upper_limit, params.acc);
    	std::cout << "Form factor integral calculated for lower limit " << lower_limit << " GeV^2 and upper limit " << upper_limit << " GeV^2.\n";
	return result;
    }

    // Set up GSL integration
//    gsl_function F;
//    F.function = &interpolated_function;
//    F.params = &params;

//    double result, error;

    // Relaxed tolerances
//    double epsabs = 1e-3;
//    double epsrel = 1e-3;

    // Perform the integration over the range of x values
//    int status = gsl_integration_qags(&F, lower_limit, upper_limit, epsabs, epsrel, 1000, workspace, &result, &error);
//    if (status != 0) {
//        std::cerr << "GSL integration error: " << gsl_strerror(status) << std::endl;
//    }

    // Store result in memo
//    memo[key] = result;
}

//*****************************************************************
// INTEGRATORS
//*****************************************************************

// Function to perform the integration neutrino-photon cross section
double neutrino_photon_xsec(double s, double xl[], double xu[], size_t calls) {
    gsl_monte_function G = { &final_amplitude, NDIM, &s };

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    double res, err;

    gsl_monte_vegas_state *state = gsl_monte_vegas_alloc(NDIM);
    gsl_monte_vegas_integrate(&G, xl, xu, NDIM, 10000, r, state, &res, &err);

    do {
        gsl_monte_vegas_integrate(&G, xl, xu, 3, calls/5, r, state, &res, &err);
    } while (fabs(gsl_monte_vegas_chisq(state)-1.0) > 0.5);

    gsl_monte_vegas_free(state);
    gsl_rng_free(r);

    return res;
}

// Wrapper function for the outer integral
double neutrino_nucleus_xsec(double *w, size_t dim, void *params) {
    (void)(dim); // avoid unused parameter warning

    auto *params_tuple = (std::tuple<size_t, double, InterpolationParams*> *)params;
    size_t calls = std::get<0>(*params_tuple);
    double Ev = std::get<1>(*params_tuple);
    InterpolationParams* interp_params = std::get<2>(*params_tuple);

    double max_Q2 = interp_params->x.back();
    double min_Q2 = interp_params->x.front();

    double xl3[3] = { -1.0, -1.0, mtau2};
    double xu3[3] = { 1.0, 1.0, w[0]};

    // Set up lower and upper limit
    double lower_limit = (w[0]*w[0])/(4*Ev*Ev);
    double upper_limit = w[0] - mtau2;

    std::cout << "Neutrino energy " << Ev << "\n";
    std::cout << "Neutrino-photon center-of-mass energy " <<  w[0] << "\n";
    std::cout << "Initial lower limit " << lower_limit << " and upper limit " << upper_limit << "\n";
   
    // Ensure the upper limit does not exceed the maximum x value
    if (upper_limit > max_Q2) {
        upper_limit = max_Q2;
	std::cout << "Upper limit exceeded the maximum Q2 " << max_Q2 << ", upper limit is now set to " << upper_limit << "\n";
    }
    if (lower_limit < min_Q2) {
        lower_limit = min_Q2;
	std::cout << "Lower limit exceeded the minimum Q2 " << min_Q2 << ", lower limit is now set to " << lower_limit << "\n";
    }

    if (lower_limit == upper_limit) {
        return 0.0;
    }
//    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(1000);
//    if (!workspace) {
//        std::cerr << "Failed to allocate GSL integration workspace." << std::endl;
//        return 0.0;
//    }

    double form_factor_integral = integrate_form_factor(lower_limit, upper_limit, *interp_params);
//    gsl_integration_workspace_free(workspace);

    return (Z_Ar*Z_Ar*aem/pi) * (-1/(4*w[0]*w[0])) * neutrino_photon_xsec(w[0], xl3, xu3, calls) * form_factor_integral;
}

// Function for upper limit smax
std::pair<double, double> calculate_slimits(double Ev, double MTarget) {
    double smax = (Ev*(mtau2 + 2*Ev*MTarget + sqrt((2*Ev*MTarget - mtau2)*(2*Ev*MTarget - mtau2) - 4*MTarget*MTarget*mtau2)))/(2*Ev + MTarget);
    double smin = (Ev*(mtau2 + 2*Ev*MTarget - sqrt((2*Ev*MTarget - mtau2)*(2*Ev*MTarget - mtau2) - 4*MTarget*MTarget*mtau2)))/(2*Ev + MTarget);
    return std::make_pair(smin, smax);
}

double smax_threshold(double Ev, double max_Q2) {
    return 2*Ev*sqrt(max_Q2);
}

double Ethresh(double m1, double m2, double MTarget) {
    double result = ((m1 + m2 + MTarget)*(m1 + m2 + MTarget) - MTarget*MTarget)/(2*MTarget);
    return result;
}

int main() {
    // Number of Monte Carlo samples
    size_t calls = 500000;

    // Input value for the neutrino energy
    double Ev;
    std::cout << "Enter the neutrino energy (GeV)\n\n";
    std::cin >> Ev;
    double Emin = Ethresh(mtau, 0.0, MArgon);
    if(Ev < Emin){
        std::cout << "Invalid value; Ev must be greater than the threshold energy E = " << Emin << " GeV.";
	return 0.0;
    }

    std::string filename = "../../csv/form_factors/FF_Ar_3Fp-red_Alt.csv";

    // Initilize argon form factor data
    std::vector<double> Q2, integrand;
    read_csv_and_process(filename, Q2, integrand);
    double max_Q2 = Q2.back();
    double min_Q2 = Q2.front();
    
    // Define limits of integration for neutrino-photon center-of-mass integral.
    std::pair<double, double> slimit = calculate_slimits(Ev, MArgon);
    double smin = slimit.first;
    double smax = slimit.second;
    std::cout << "s limits calculated: smin = " << smin << ", smax = " << smax << "\n";

//    double smax_thresh = smax_threshold(Ev, max_Q2);
//    std::cout << "s max threshold is " << smax_thresh << "\n";
//    if (smax > smax_thresh) {
//        smax = smax_thresh;
//	std::cout << "new s upper limit is " << smax << "\n";
//    }

    double xl4[1] = { smin };
    double xu4[1] = { smax };

    // Initialize GSL interpolation
    gsl_interp *interp = gsl_interp_alloc(gsl_interp_steffen, Q2.size());
    if (!interp) {
        std::cerr << "Failed to allocate GSL interpolation object." << std::endl;
        return 0.0;
    }
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    if (!acc) {
        gsl_interp_free(interp);
        std::cerr << "Failed to allocate GSL interpolation accelerator." << std::endl;
        return 0.0;
    }
    if (gsl_interp_init(interp, Q2.data(), integrand.data(), Q2.size()) != 0) {
        gsl_interp_free(interp);
        gsl_interp_accel_free(acc);
        std::cerr << "Failed to initialize GSL interpolation." << std::endl;
        return 0.0;
    }

    InterpolationParams interp_params = {interp, acc, Q2, integrand};

    std::tuple<size_t, double, InterpolationParams*> data = std::make_tuple(calls, Ev, &interp_params);

    gsl_monte_function G = { &neutrino_nucleus_xsec, 1, &data };

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    double res, err;

    std::cout << "Vegas integration running...\n";
    std::cout << "\ts lower limit = " << xl4 << " GeV^2\n";
    std::cout << "\ts upper limit = " << xu4 << " GeV^2\n";

    gsl_monte_vegas_state *state = gsl_monte_vegas_alloc(1);
    gsl_monte_vegas_integrate(&G, xl4, xu4, 1, 10000, r, state, &res, &err);

    do {
        gsl_monte_vegas_integrate(&G, xl4, xu4, 1, calls/5, r, state, &res, &err);
    } while (fabs(gsl_monte_vegas_chisq(state)-1.0) > 0.5);


    // Free GSL resources
    gsl_monte_vegas_free(state);
    gsl_rng_free(r);

    gsl_interp_free(interp);
    gsl_interp_accel_free(acc);

    std::cout << "The integral result is: " << res << " with error " << err << std::endl;
   
    return 0;
}
