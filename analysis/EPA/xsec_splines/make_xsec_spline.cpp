#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <map>
#include <algorithm>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

//*****************************************************************
// CONSTANTS
//*****************************************************************
const double pi = 3.1415926535897;

const double GF = 0.0000116637; // Fermi constant [1/GeV^2]
const double aem = 1/137.036; // em coupling [1]
const double sinW2 = 0.23129; // MSbar (sinThetaW)^2 at electroweak scale

const double me = 5.109e-4; // Electron lepton mass [GeV]
const double me2 = me*me;

const double mmu = 0.105; // Muon lepton mass [GeV]
const double mmu2 = mmu*mmu;

const double mtau = 1.778; // Tau lepton mass [GeV]
const double mtau2 = mtau*mtau; // Tau lepton mass squared [GeV^2]

//*****************************************************************
// GLOBAL VARIABLES
//*****************************************************************
int process;
std::string command;
std::string type;

double m3sqr;
double m4sqr;
double ERange_start;
double ERange_end;
double nXsec;
std::vector<std::tuple<double,double,double>> xsections;

std::string xsec_file_out;

// Number of Monte Carlo samples
size_t calls = 500000;

//*****************************************************************
// KINEMATIC DOT PRODUCTS
//*****************************************************************

double p1q(double s){
    double result = s / 2;
    return result;
}

double p1p2(double cosCM, double l, double s) {
    double result = (s-l)/4*(1-cosCM);
    return result;
}

double p1p3(double cos34, double cosCM, double l, double s, double m3sqr, double m4sqr, double beta34) {
    double result = (m3sqr*((s-l)*cosCM+l+s) + m4sqr*((l-s)*cosCM-l-s) - beta34*l*(cos34*(cosCM*((sqrt(l)-sqrt(s))*(sqrt(l)-sqrt(s))*cosCM-l+s)+2*sqrt(l*s)) + sqrt(1-cos34*cos34)*sqrt(1-cosCM*cosCM)*((sqrt(l)-sqrt(s))*(sqrt(l)-sqrt(s))*cosCM-l+s)) + l*((s-l)*cosCM + l + s))/(8*l);
    return result;
}

double p1p4(double cos34, double cosCM, double l, double s, double m3sqr, double m4sqr, double beta34) {
    double result = 1/(8*l)*((m4sqr-m3sqr)*((s-l)*cosCM+l+s) + l*(beta34*(cos34*(cosCM*((sqrt(l)-sqrt(s))*(sqrt(l)-sqrt(s))*cosCM-l+s) + 2*sqrt(l*s)) + sqrt(1-cos34*cos34)*sqrt(1-cosCM*cosCM)*((sqrt(l)-sqrt(s))*(sqrt(l)-sqrt(s))*cosCM-l+s)) + (s-l)*cosCM + l + s));
    return result;
}

double p2q(double cosCM, double l, double s) {
    double result = (s-l)*(cosCM+1)/4;
    return result;
}

double p2p3(double cos34, double cosCM, double l, double s, double m3sqr, double m4sqr, double beta34) {
    double result = ((l-s)*(l*(beta34*(sqrt(1-cos34*cos34)*sqrt(1-cosCM*cosCM)+cos34*cosCM)-1)-m3sqr+m4sqr))/(4*l);
    return result;
}

double p2p4(double cos34, double cosCM, double l, double s, double m3sqr, double m4sqr, double beta34) {
    double result = -((l-s)*(beta34*l*(sqrt(1-cos34*cos34)*sqrt(1-cosCM*cosCM) + cos34*cosCM) + l - m3sqr + m4sqr))/(4*l);
    return result;
}

double p3q(double cos34, double cosCM, double l, double s, double m3sqr, double m4sqr, double beta34) {
    double result = ((m3sqr-m4sqr)*((l-s)*cosCM+l+s) + l*(beta34*(cos34*(cosCM*((sqrt(l)-sqrt(s))*(sqrt(l)-sqrt(s))*cosCM+l-s)+2*sqrt(l*s)) + sqrt(1-cos34*cos34)*sqrt(1-cosCM*cosCM)*((sqrt(l)-sqrt(s))*(sqrt(l)-sqrt(s))*cosCM+l-s)) + (l-s)*cosCM + l + s))/(8*l);
    return result;
}

double p3p4(double l, double m3sqr, double m4sqr) {
    double result = (l-m3sqr-m4sqr)/2;
    return result;
}

double p4q(double cos34, double cosCM, double l, double s, double m3sqr, double m4sqr, double beta34) {
    double result = ((m4sqr-m3sqr)*((l-s)*cosCM+l+s) - beta34*l*(cos34*(cosCM*((sqrt(l)-sqrt(s))*(sqrt(l)-sqrt(s))*cosCM + l - s) + 2*sqrt(l*s)) + sqrt(1-cos34*cos34)*sqrt(1-cosCM*cosCM)*((sqrt(l)-sqrt(s))*(sqrt(l)-sqrt(s))*cosCM + l - s)) + l*((l-s)*cosCM + l + s))/(8*l);
    return result;
}

//*****************************************************************
// AMPLITUDES
//*****************************************************************

double amp_W(double m3sqr, double m4sqr, double p1p2_s, double p1p3_s, double p1p4_s, double p1q_s, double p2p3_s, double p2p4_s, double p2q_s, double p3p4_s, double p3q_s, double p4q_s) {
    double amp1 = 512*pi*aem*GF*GF*m4sqr*p1p3_s*p2p4_s/p4q_s/p4q_s;
    double amp2 = -512*pi*aem*GF*GF*m4sqr*p1p3_s*p2q_s/p4q_s/p4q_s;
    double amp3 = 512*pi*aem*GF*GF*m3sqr*p1p3_s*p2p4_s/p3q_s/p3q_s;
    double amp4 = -512*pi*aem*GF*GF*m3sqr*p1q_s*p2p4_s/p3q_s/p3q_s;
    double amp5 = 512*pi*aem*GF*GF*p1p3_s*p2p4_s/p3q_s;
    double amp6 = -512*pi*aem*GF*GF*p1q_s*p2p4_s/p3q_s;
    double amp7 = 512*pi*aem*GF*GF*p1p3_s*p2p4_s/p4q_s;
    double amp8 = -512*pi*aem*GF*GF*p1p3_s*p2q_s/p4q_s;
    double amp9 = -1024*pi*aem*GF*GF*p1p3_s*p2p4_s*p3p4_s/p3q_s/p4q_s;
    double amp10 = 512*pi*aem*GF*GF*p1q_s*p2p4_s*p3p4_s/p3q_s/p4q_s;
    double amp11 = 512*pi*aem*GF*GF*p1p3_s*p2q_s*p3p4_s/p3q_s/p4q_s;
    double amp12 = -512*pi*aem*GF*GF*p1p3_s*p2p3_s/p3q_s;
    double amp13 = -512*pi*aem*GF*GF*p1p4_s*p2p4_s/p4q_s;
    return amp1+amp2+amp3+amp4+amp5+amp6+amp7+amp8+amp9+amp10+
	   amp11+amp12+amp13;
}

double amp_WZ(double ml2, double p1p2_s, double p1p3_s, double p1p4_s, double p1q_s, double p2p3_s, double p2p4_s, double p2q_s, double p3p4_s, double p3q_s, double p4q_s) {
    double amp1 = 512*pi*aem*GF*GF*sinW2*sinW2*ml2*ml2*p1p2_s/p3q_s/p3q_s;
    double amp2 = 512*pi*aem*GF*GF*sinW2*sinW2*ml2*ml2*p1p2_s/p4q_s/p4q_s;
    double amp3 = 256*pi*aem*GF*GF*sinW2*ml2*ml2*p1p2_s/p3q_s/p3q_s;
    double amp4 = 256*pi*aem*GF*GF*sinW2*ml2*ml2*p1p2_s/p4q_s/p4q_s;
    double amp5 = 128*pi*aem*GF*GF*ml2*p1p3_s*p2p4_s/p3q_s/p3q_s;
    double amp6 = -128*pi*aem*GF*GF*ml2*p1q_s*p2p4_s/p3q_s/p3q_s;
    double amp7 = 128*pi*aem*GF*GF*ml2*p1p3_s*p2p4_s/p4q_s/p4q_s;
    double amp8 = -128*pi*aem*GF*GF*ml2*p1p3_s*p2q_s/p4q_s/p4q_s;
    double amp9 = 1024*pi*aem*GF*GF*ml2*sinW2*sinW2*p1q_s*p2q_s/p3q_s/p4q_s;
    double amp10 = -1024*pi*aem*GF*GF*ml2*sinW2*sinW2*p1p2_s*p3p4_s/p3q_s/p4q_s;
    double amp11 = 512*pi*aem*GF*GF*ml2*sinW2*sinW2*p1p4_s*p2p3_s/p3q_s/p3q_s;
    double amp12 = 512*pi*aem*GF*GF*ml2*sinW2*sinW2*p1p3_s*p2p4_s/p3q_s/p3q_s;
    double amp13 = -512*pi*aem*GF*GF*ml2*sinW2*sinW2*p1q_s*p2p4_s/p3q_s/p3q_s;
    double amp14 = -512*pi*aem*GF*GF*ml2*sinW2*sinW2*p1p4_s*p2q_s/p3q_s/p3q_s;
    double amp15 = 512*pi*aem*GF*GF*ml2*sinW2*sinW2*p1p4_s*p2p3_s/p4q_s/p4q_s;
    double amp16 = -512*pi*aem*GF*GF*ml2*sinW2*sinW2*p1q_s*p2p3_s/p4q_s/p4q_s;
    double amp17 = 512*pi*aem*GF*GF*ml2*sinW2*sinW2*p1p3_s*p2p4_s/p4q_s/p4q_s;
    double amp18 = -512*pi*aem*GF*GF*ml2*sinW2*sinW2*p1p3_s*p2q_s/p4q_s/p4q_s;
    double amp19 = 512*pi*aem*GF*GF*ml2*sinW2*p1q_s*p2q_s/p3q_s/p4q_s;
    double amp20 = -512*pi*aem*GF*GF*ml2*sinW2*p1p2_s*p3p4_s/p3q_s/p4q_s;
    double amp21 = 512*pi*aem*GF*GF*ml2*sinW2*p1p3_s*p2p4_s/p3q_s/p3q_s;
    double amp22 = -512*pi*aem*GF*GF*ml2*sinW2*p1q_s*p2p4_s/p3q_s/p3q_s;
    double amp23 = 512*pi*aem*GF*GF*ml2*sinW2*p1p3_s*p2p4_s/p4q_s/p4q_s;
    double amp24 = -512*pi*aem*GF*GF*ml2*sinW2*p1p3_s*p2q_s/p4q_s/p4q_s;
    double amp25 = -128*pi*aem*GF*GF*p1p3_s*p2p3_s/p3q_s;
    double amp26 = 128*pi*aem*GF*GF*p1p3_s*p2p4_s/p3q_s;
    double amp27 = -128*pi*aem*GF*GF*p1q_s*p2p4_s/p3q_s;
    double amp28 = 128*pi*aem*GF*GF*p1p3_s*p2p4_s/p4q_s;
    double amp29 = -128*pi*aem*GF*GF*p1p4_s*p2p4_s/p4q_s;
    double amp30 = -128*pi*aem*GF*GF*p1p3_s*p2q_s/p4q_s;
    double amp31 = -256*pi*aem*GF*GF*p1p3_s*p2p4_s*p3p4_s/p3q_s/p4q_s;
    double amp32 = 128*pi*aem*GF*GF*p1q_s*p2p4_s*p3p4_s/p3q_s/p4q_s;
    double amp33 = 128*pi*aem*GF*GF*p1p3_s*p2q_s*p3p4_s/p3q_s/p4q_s;
    double amp34 = -1024*pi*aem*GF*GF*sinW2*sinW2*p1p3_s*p2p3_s/p3q_s;
    double amp35 = 512*pi*aem*GF*GF*sinW2*sinW2*p1p4_s*p2p3_s/p3q_s;
    double amp36 = 512*pi*aem*GF*GF*sinW2*sinW2*p1p3_s*p2p4_s/p3q_s;
    double amp37 = -512*pi*aem*GF*GF*sinW2*sinW2*p1q_s*p2p4_s/p3q_s;
    double amp38 = -512*pi*aem*GF*GF*sinW2*sinW2*p1p4_s*p2q_s/p3q_s;
    double amp39 = 512*pi*aem*GF*GF*sinW2*sinW2*p1p4_s*p2p3_s/p4q_s;
    double amp40 = -512*pi*aem*GF*GF*sinW2*sinW2*p1q_s*p2p3_s/p4q_s;
    double amp41 = 512*pi*aem*GF*GF*sinW2*sinW2*p1p3_s*p2p4_s/p4q_s;
    double amp42 = -1024*pi*aem*GF*GF*sinW2*sinW2*p1p4_s*p2p4_s/p4q_s;
    double amp43 = -512*pi*aem*GF*GF*sinW2*sinW2*p1p3_s*p2q_s/p4q_s;
    double amp44 = -1024*pi*aem*GF*GF*sinW2*sinW2*p1p4_s*p2p3_s*p3p4_s/p3q_s/p4q_s;
    double amp45 = 512*pi*aem*GF*GF*sinW2*sinW2*p1q_s*p2p3_s*p3p4_s/p3q_s/p4q_s;
    double amp46 = -1024*pi*aem*GF*GF*sinW2*sinW2*p1p3_s*p2p4_s*p3p4_s/p3q_s/p4q_s;
    double amp47 = 512*pi*aem*GF*GF*sinW2*sinW2*p1q_s*p2p4_s*p3p4_s/p3q_s/p4q_s;
    double amp48 = 512*pi*aem*GF*GF*sinW2*sinW2*p1p3_s*p2q_s*p3p4_s/p3q_s/p4q_s;
    double amp49 = 512*pi*aem*GF*GF*sinW2*sinW2*p1p4_s*p2q_s*p3p4_s/p3q_s/p4q_s;
    double amp50 = -512*pi*aem*GF*GF*sinW2*p1p3_s*p2p3_s/p3q_s;
    double amp51 = 512*pi*aem*GF*GF*sinW2*p1p3_s*p2p4_s/p3q_s;
    double amp52 = -512*pi*aem*GF*GF*sinW2*p1q_s*p2p4_s/p3q_s;
    double amp53 = 512*pi*aem*GF*GF*sinW2*p1p3_s*p2p4_s/p4q_s;
    double amp54 = -512*pi*aem*GF*GF*sinW2*p1p4_s*p2p4_s/p4q_s;
    double amp55 = -512*pi*aem*GF*GF*sinW2*p1p3_s*p2q_s/p4q_s;
    double amp56 = -1024*pi*aem*GF*GF*sinW2*p1p3_s*p2p4_s*p3p4_s/p3q_s/p4q_s;
    double amp57 = 512*pi*aem*GF*GF*sinW2*p1q_s*p2p4_s*p3p4_s/p3q_s/p4q_s;
    double amp58 = 512*pi*aem*GF*GF*sinW2*p1p3_s*p2q_s*p3p4_s/p3q_s/p4q_s;

    return amp1+amp2+amp3+amp4+amp5+amp6+amp7+amp8+amp9+amp10+
	   amp11+amp12+amp13+amp14+amp15+amp16+amp17+amp18+amp19+amp20+
	   amp21+amp22+amp23+amp24+amp25+amp26+amp27+amp28+amp29+amp30+
	   amp31+amp32+amp33+amp34+amp35+amp36+amp37+amp38+amp39+amp40+
	   amp41+amp42+amp43+amp44+amp45+amp46+amp47+amp48+amp49+amp50+
	   amp51+amp52+amp53+amp54+amp55+amp56+amp57+amp58;
}

//*****************************************************************
// FUNCTIONS
//*****************************************************************

// Final amplitude
double final_amplitude(double *k, size_t dim, void *params) {
    (void)(dim); /* avoid unused parameter warning */

    auto *params_tuple = (std::tuple<double, double, double, std::string> *)params;
    double s = std::get<0>(*params_tuple);
    double m3sqr = std::get<1>(*params_tuple);
    double m4sqr = std::get<2>(*params_tuple);
    double beta34 = sqrt(1 - 2*(m3sqr+m4sqr)/k[2] + (m4sqr-m3sqr)*(m4sqr-m3sqr)/k[2]/k[2]);
    std::string type = std::get<3>(*params_tuple);

    double p1q_s = p1q(s);
    double p1p2_s = p1p2(k[1], k[2], s);
    double p1p3_s = p1p3(k[0], k[1], k[2], s, m3sqr, m4sqr, beta34);
    double p1p4_s = p1p4(k[0], k[1], k[2], s, m3sqr, m4sqr, beta34);
    double p2q_s = p2q(k[1], k[2], s);
    double p2p3_s = p2p3(k[0], k[1], k[2], s, m3sqr, m4sqr, beta34);
    double p2p4_s = p2p4(k[0], k[1], k[2], s, m3sqr, m4sqr, beta34);
    double p3q_s = p3q(k[0], k[1], k[2], s, m3sqr, m4sqr, beta34);
    double p3p4_s = p3p4(k[2], m3sqr, m4sqr);
    double p4q_s = p4q(k[0], k[1], k[2], s, m3sqr, m4sqr, beta34);

    if (type == "WZ") { // This trident has the same flavor for all neutrinos and leptons. In this instance, ml = mass of the lepton.
        return amp_WZ(m3sqr, p1p2_s, p1p3_s, p1p4_s, p1q_s, p2p3_s, p2p4_s, p2q_s, p3p4_s, p3q_s, p4q_s);
    }
    
    if (type == "W") {
        return amp_W(m3sqr, m4sqr, p1p2_s, p1p3_s, p1p4_s, p1q_s, p2p3_s, p2p4_s, p2q_s, p3p4_s, p3q_s, p4q_s);
    }

    if (type == "Z") {
        return 0.0;
    }

    std::cout << "Error occurred. Only available trident types are WZ, W, and Z. Your type is " << type;
    return 0.0;
}

//*****************************************************************
// INTEGRATORS
//*****************************************************************

// Function to perform the integration neutrino-photon cross section
std::pair<double,double> neutrino_photon_xsec(std::tuple<double, double, double, std::string> data, double xl[], double xu[], size_t calls) {
    gsl_monte_function G = { &final_amplitude, 3, &data };
    double s = std::get<0>(data);

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    double res, err;

    gsl_monte_vegas_state *state = gsl_monte_vegas_alloc(3);
    gsl_monte_vegas_integrate(&G, xl, xu, 3, 10000, r, state, &res, &err);

    do {
        gsl_monte_vegas_integrate(&G, xl, xu, 3, calls/5, r, state, &res, &err);
    } while (fabs(gsl_monte_vegas_chisq(state)-1.0) > 0.5);

    gsl_monte_vegas_free(state);
    gsl_rng_free(r);

    res = -res/(4*s); // For some reason, Ballet et al. and Magill & Plestid have res = - res/(4s). I am omitting minus sign for now.
    err = -err/(4*s);

    return {res,err};
}

std::vector<std::tuple<double,double,double>> ComputeCrossSections(double start, double finish, int total, size_t calls) {
    double h = (double)(finish - start) / (double)(total - 1);
    std::vector<double> v(total);
    std::generate(v.begin(), v.end(), [n = 0, &h, &start]() mutable { return start + n++ * h; });
    std::cout << " \nCross section computation over energy range initialized.\n";
    int nenergy = 1;

    std::vector<std::tuple<double,double,double>> xsections;
    for(const double& s : v){
	double ml = sqrt(m3sqr) + sqrt(m4sqr);
	double ml2 = ml*ml;
	double xl3[3] = { -1.0, -1.0, ml2};
	double xu3[3] = { 1.0, 1.0, s};

	if (s <= ml2) {
	    xsections.push_back(std::make_tuple(s, 0.0, 0.0));
	}
	else {
	    std::pair integral_results = neutrino_photon_xsec(std::make_tuple(s, m3sqr, m4sqr, type), xl3, xu3, calls);
	    double xsec = integral_results.first;
	    double delta = integral_results.second;
	    xsections.push_back(std::make_tuple(s, xsec, delta));
	}
	std::cout << "Cross sections computed: " << nenergy << " / " << total << "\r";
	std::cout.flush();
	++nenergy;
    }
    return xsections;
}

//*************************************************
// Write file with calculated cross sections
//*************************************************

void WriteXSecFile(std::string filename, std::vector<std::tuple<double,double,double>>& xsections){
    
    std::ofstream outfile;
    outfile.open(filename.c_str(), std::ios_base::trunc | std::ios_base::out | std::ios_base::in);

    // format: energy, xsec, deltaxsec
    for(auto const& [energy_i, xsec, delta] : xsections){
        outfile << energy_i << "," << xsec << "," << delta << "\n";
    }

    outfile.close();
}

void SetTridentProcess(){

    if (process == 1) { // nu_e -> nu_e e+ e-
        m3sqr = me2; 
	m4sqr = me2;
	type = "WZ";
    }
    else if (process == 2) { // nu_e -> nu_mu mu+ e-
        m3sqr = mmu2;
	m4sqr = me2;
	type = "W";
    }    
    else if (process == 3) { // nu_e -> nu_tau tau+ e-
        m3sqr = mtau2;
	m4sqr = me2;
	type = "W";
    }         
    else if (process == 4) { // nu_e -> nu_e mu+ mu-
        m3sqr = mmu2;
	m4sqr = mmu2;
	type = "Z";
    }         
    else if (process == 5) { // nu_e -> nu_e tau+ tau-
        m3sqr = mtau2;
	m4sqr = mtau2;
	type = "Z";
    }    
    else if (process == 6) { // nu_mu -> nu_mu mu+ mu-
        m3sqr = mmu2;
	m4sqr = mmu2;
	type = "WZ";
    }    
    else if (process == 7) { // nu_mu -> nu_e e+ mu-
        m3sqr = me2;
	m4sqr = mmu2;
	type = "W";
    }         
    else if (process == 8) { // nu_mu -> nu_tau tau+ mu-
        m3sqr = mtau2;
	m4sqr = mmu2;
	type = "W";
    }         
    else if (process == 9) { // nu_mu -> nu_mu e+ e-
        m3sqr = me2;
	m4sqr = me2;
	type = "Z";
    }         
    else if (process == 10) { // nu_mu -> nu_mu tau+ tau-
        m3sqr = mtau2;
	m4sqr = mtau2;
	type = "Z";
    }
    else if (process == 11) { // nu_tau -> nu_tau tau+ tau-
        m3sqr = mtau2;
	m4sqr = mtau2;
	type = "WZ";
    }
    else if (process == 12) { // nu_tau -> nu_e e+ tau-
        m3sqr = me2;
	m4sqr = mtau2;
	type = "W";
    }
    else if (process == 13) { // nu_tau -> nu_mu mu+ tau-
        m3sqr = mmu2;
	m4sqr = mtau2;
	type = "W";
    } 
    else if (process == 14) { // nu_tau -> nu_tau e+ e-
        m3sqr = me2;
	m4sqr = me2;
	type = "Z";
    }  
    else if (process == 15) { // nu_tau -> nu_tau mu+ mu-
        m3sqr = mmu2;
	m4sqr = mmu2;
	type = "Z";
    }
    return;
  
}
int main() {
    // User input
    std::cout << "\n";
    std::cout << "Select one of the following tridents:\n\n";
    std::cout << " [1] nu_e -> nu_e e+ e-     [6] nu_mu -> nu_mu mu+ mu-     [11] nu_tau -> nu_tau tau+ tau-\n";
    std::cout << " [2] nu_e -> nu_mu mu+ e-   [7] nu_mu -> nu_e e+ mu-       [12] nu_tau -> nu_e e+ tau-\n";
    std::cout << " [3] nu_e -> nu_tau tau+ e- [8] nu_mu -> nu_tau tau+ mu-   [13] nu_tau -> nu_mu mu+ tau-\n";
    std::cout << " [4] nu_e -> nu_e mu+ mu-   [9] nu_mu -> nu_mu e+ e-       [14] nu_tau -> nu_tau e+ e-\n";
    std::cout << " [5] nu_e -> nu_e tau+ tau- [10] nu_mu -> nu_mu tau+ tau-  [15] nu_tau -> nu_tau mu+ mu-\n";
    std::cin >> process;
    if (process != 1 && process != 2 && process != 3 && process != 4 && process != 5 &&
        process != 6 && process != 7 && process != 8 && process != 9 && process != 10 &&
        process != 11 && process != 12 && process != 13 && process != 14 && process != 15) {
        std::cout << "\n Invalid selection \n";
        return 0;
    }
    SetTridentProcess();

    std::cout << "\n";
    std::cout << "You can compute the trident [CrossSection] or [CrossSectionOverRange]\n\n";
    std::cin >> command;

    if (command.compare("CrossSection") != 0 && command.compare("CrossSectionOverRange") != 0){
        std::cout << "\n Invalid selection \n";
	return 0;
    }

    if (command.compare("CrossSectionOverRange") == 0){
        std::cout << "\n";
        std::cout << "Enter range start of neutrino-photon center-of-mass energy squared (GeV^2)\n\n";
	std::cin >> ERange_start;
	std::cout << "\n";
	std::cout << "Enter range end of neutrino-photon center-of-mass energy squared (GeV^2)\n\n";
	std::cin >> ERange_end;
	std::cout << "\n";
	std::cout << "Enter the number of xsec to be calculated \n\n";
	std::cin >> nXsec;
	std::cout << "\n";
	std::cout << "Enter the name of the output file \n\n";
	std::cin >> xsec_file_out;
	std::cout << "\n";
    }


    // Compute cross section
    if (command.compare("CrossSection") == 0){
	double s;
	std::cout << "Enter the neutrino-photon center-of-mass energy squared (GeV^2)\n\n";
	std::cin >> s;

	std::tuple<double, double, double, std::string> data = std::make_tuple(s, m3sqr, m4sqr, type);
	double ml = sqrt(m3sqr) + sqrt(m4sqr);

	double xl3[3] = { -1.0, -1.0, ml*ml};
	double xu3[3] = { 1.0, 1.0, s};

	std::pair xsec_results = neutrino_photon_xsec(data, xl3, xu3, calls);

	double res = xsec_results.first;
	double err = xsec_results.second;
        std::cout << "\n\n";
	std::cout << "The neutrino-photon transverse cross section is  ( " << res << " +- " << err << " ) GeV^-2  \n";
    }

    if(command.compare("CrossSectionOverRange") == 0){
	xsections = ComputeCrossSections(ERange_start, ERange_end, nXsec, calls);
	WriteXSecFile(xsec_file_out.c_str(), xsections);
	std::cout << "\n\n";
	std::cout << "Generated " << nXsec << " trident cross sections for the neutrino-photon center-of-mass energy squared range s=(" << ERange_start << ", " << ERange_end << ") GeV^2\n";
	std::cout << "File was generated and saved as " << xsec_file_out;
    }
  
    return 0;
}
