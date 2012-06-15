//============================================================================
// Name        : semiconductor.cpp
// Author      : Karl Lopker
// Version     :
// Copyright   : none
// Description : Semiconductor simulation for CS240A
//============================================================================

#include <iostream>
#include <cmath>
#include <time.h>
#include <unistd.h>
#include <ios>
#include <fstream>
#include <string>

#include "armadillo"
#include "INIReader.h"
#include "constants.h"
#include "Poisson.h"
#include "FermiIteration.h"
#include "sys/types.h"
#include "sys/sysinfo.h"
using namespace arma;
using namespace std;

// Type of Laplacian Employed
int lap_kind;
// Number of iteration in Poisson Solver
int itertime;
// Maximum Number of iteration for Self-consistent Potential/QuasifermiLevel Calculations
int big_itertime;
// ion=0 implies Complete Ionization; ion=1 implies partial ionization
int ion;
// Tolerance for the Convergence
double ETOL;

// Layer lengths
vec lengths;
//Bandgap in each layer
vec Eg;
//Doping in each layer, atoms/m^2
vec Na_1;
vec Nd_1;
//Activation Energy for Dopants
vec Eac;
//Effective mass of Electrons
vec me;
//Effective mass of Holes
vec mh;
vec mu_p;
vec mu_n;
vec ni;
vec epsilon;

vec Efn; // Electron fermi level
vec Efp; // Hole fermi level
vec Jn; // Current density
vec Jp;
vec V, V_old;

vec Nb;

vec n;
vec p;

int N = 83;
double dx;
double v_bias;
double Efn_left_bound, Efp_right_bound;
double V_left_B, V_right_B;
double main_iteration;

bool checkUsage(int number_of_arguments);
bool setUp(char* config_path);
double main_loop();
void calcJ(const vec* mu_, const vec* v_, const vec* Ef_, double dx, vec* J_);
int getMemoryUsed();
void process_mem_usage(double& vm_usage, double& resident_set);

int main(int argc, char** argv) {

	if (!checkUsage(argc)) {
		return 1;
	}

	if (!setUp(argv[1])) {
		return 1;
	}

	int t1 = clock();
	for (int i = 0; i < main_iteration; ++i) {
		double error = main_loop();
		printf("error: %f at iteration: %d\n", error, i);

		double vm, rss;
		process_mem_usage(vm, rss);
		printf("Memory: %f\n", rss);

		if (error < 1e-0018) {
			break;
		}
	}

	int t2 = clock();
	printf("Time: %f seconds\n", ((float) (t2 - t1)) / 1000000);
	printf("Problem Size: %d\n", N);
	calcJ(&mu_n, &n, &Efn, dx, &Jn);
	calcJ(&mu_p, &p, &Efp, dx, &Jp);

	Efn.save("Efn.txt", raw_ascii);
	Efp.save("Efp.txt", raw_ascii);
	V.save("V.txt", raw_ascii);

	return 0;
}

bool checkUsage(int number_of_arguments) {
	if (number_of_arguments != 2) {
		cout << "Usage: semiconductor <Config_File>" << endl;
		return false;
	}
	return true;
}

bool setUp(char* config_path) {

	INIReader reader(config_path);
	if (reader.ParseError() < 0) {
		cout << "Can't load config file: " << config_path << endl
				<< "Using defaults.\n";
	} else {
		cout << "Config loaded from: " << config_path << endl;
	}

	lengths = vec(reader.Get("Problem", "length", "500e-9 500e-9"));
	vec t_Eg = vec(reader.Get("Problem", "Eg", "3.2 3.2"));
	t_Eg *= q;
	vec t_Na_1 = vec(reader.Get("Problem", "Na_1", "1e23 0"));
	vec t_Nd_1 = vec(reader.Get("Problem", "Nd_1", "0 1e23"));
	vec t_Eac = vec(reader.Get("Problem", "Eac", "0.020 0.020"));
	t_Eac *= q;
	vec t_me = vec(reader.Get("Problem", "me", "0.2 0.2"));
	t_me *= m0;
	vec t_mh = vec(reader.Get("Problem", "mh", "0.8 0.8"));
	t_mh *= m0;
	vec t_mu_n = vec(reader.Get("Problem", "mu_n", "0.044 0.044"));
	vec t_mu_p = vec(reader.Get("Problem", "mu_p", "0.035 0.035"));
	vec t_epsilon = vec(reader.Get("Problem", "epsilon", "10 10"));
	t_epsilon *= e0;

	v_bias = vec(reader.Get("Problem", "V_bias", "1")).at(0);
	main_iteration = vec(reader.Get("Numerical", "big_itertime", "100")).at(0);
	lap_kind = vec(reader.Get("Numerical", "lap_kind", "0")).at(0);
	itertime = vec(reader.Get("Numerical", "itertime", "10000")).at(0);
	ion = vec(reader.Get("Numerical", "ion", "0")).at(0);

	double dx_theo = sqrt(
			10 * e0 * kb * temperature / pow(q, 2) / (t_Na_1.at(0)));
	N = floor(accu(lengths) / dx_theo);
	N += N % 2;
	dx = accu(lengths) / (N - 1);

	vec t_ni = zeros<vec>(lengths.n_elem, 1);
	for (unsigned int i = 0; i < t_ni.n_elem; ++i) {
		t_ni.at(i) = pow(
				2
						* pow(
								t_me.at(i) * kb * temperature
										/ (2 * M_PI * pow(hbar, 2)), 3.0 / 2)
						* (2
								* pow(
										t_mh.at(i) * kb * temperature
												/ (2 * M_PI * pow(hbar, 2)),
										3.0 / 2)
								* exp(-t_Eg.at(i) / (kb * temperature))), 0.5);
	}

//	Initialize arrays
	Eg = zeros<vec>(N, 1);
	Na_1 = zeros<vec>(N, 1);
	Nd_1 = zeros<vec>(N, 1);
	me = zeros<vec>(N, 1);
	mh = zeros<vec>(N, 1);
	ni = zeros<vec>(N, 1);
	epsilon = zeros<vec>(N, 1);
	Eac = zeros<vec>(N, 1);
	mu_p = zeros<vec>(N, 1);
	mu_n = zeros<vec>(N, 1);

	n = zeros<vec>(N, 1);
	p = zeros<vec>(N, 1);

	Efn = zeros<vec>(N, 1); //Electron fermi level
	Efp = zeros<vec>(N, 1); //Hole fermi level
	Jn = zeros<vec>(N, 1); //Electron Current density
	Jp = zeros<vec>(N, 1); //Hole Current density
	V = zeros<vec>(N, 1);

	Efn_left_bound = -(v_bias / 2) * q;
	Efp_right_bound = (v_bias / 2) * q;
	V_left_B = (Efn_left_bound / q)
			+ kb * temperature / q * log(t_Na_1.at(0) / t_ni.at(0));
	V_right_B = (Efp_right_bound / q)
			- kb * temperature / q * log(t_Nd_1.at(1) / t_ni.at(1));

	for (int i = 0; i < N; ++i) {
		V.at(i) = V_left_B + 2 * V_right_B / (N - 1) * i;
	}

	int start = 0;
	for (unsigned int i = 0; i < lengths.n_elem; ++i) {
		double length_sum = sum(lengths);
		int n = floor(lengths.at(i) / length_sum * N);
		for (unsigned int j = start; j < (n + start); ++j) {
			Eg.at(j) = t_Eg.at(i);
			Na_1.at(j) = t_Na_1.at(i);
			Nd_1.at(j) = t_Nd_1.at(i);
			me.at(j) = t_mh.at(i);
			mh.at(j) = t_mh.at(i);
			ni.at(j) = t_ni.at(i);
			epsilon.at(j) = t_epsilon.at(i);
			Eac.at(j) = t_Eac.at(i);
			mu_p.at(j) = t_mu_p.at(i);
			mu_n.at(j) = t_mu_n.at(i);
		}
		start += n;
	}

	Nb = Na_1 - Nd_1;

	return true;
}

double main_loop() {
	V_old = V;

	poisson(N, dx, lap_kind, itertime, &ni, &Efn, &Efp, &epsilon, &Nb, V_left_B,
			V_right_B, &n, &p, &V);

//    %Calculating the Quasi Fermi Levels
	fermiIteration(N, dx, &ni, &n, &p, &V, V_left_B, V_right_B, v_bias, &mu_p,
			&mu_n, &Efn, &Efp);

	double error = sum(abs(V - V_old));
	return error;
}

void calcJ(const vec* mu_, const vec* v_, const vec* Ef_, double dx, vec* J_) {
	vec f = zeros<vec>(N + 2, 1);
	f.rows(1, N) = *Ef_;
	f.at(0) = Efn_left_bound;
	f.at(N + 1) = Efp_right_bound;

	*J_ = *mu_ % *v_ % (f.rows(2, N + 1) - f.rows(0, N - 1)) / (2 * dx);

}

int parseLine(char* line) {
	int i = strlen(line);
	while (*line < '0' || *line > '9')
		line++;
	line[i - 3] = '\0';
	i = atoi(line);
	return i;
}

int getMemoryUsed() {
	FILE* file = fopen("/proc/self/status", "r");
	int result = -1;
	char line[128];

	while (fgets(line, 128, file) != NULL) {
		if (strncmp(line, "VmRSS:", 6) == 0) {
			result = parseLine(line);
		}
		break;
	}
	fclose(file);
	return result;
}

void process_mem_usage(double& vm_usage, double& resident_set) {
	using std::ios_base;
	using std::ifstream;
	using std::string;

	vm_usage = 0.0;
	resident_set = 0.0;

	// 'file' stat seems to give the most reliable results
	//
	ifstream stat_stream("/proc/self/stat", ios_base::in);

	// dummy vars for leading entries in stat that we don't care about
	//
	string pid, comm, state, ppid, pgrp, session, tty_nr;
	string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	string utime, stime, cutime, cstime, priority, nice;
	string O, itrealvalue, starttime;

	// the two fields we want
	//
	unsigned long vsize;
	long rss;

	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
			>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt >> utime
			>> stime >> cutime >> cstime >> priority >> nice >> O >> itrealvalue
			>> starttime >> vsize >> rss; // don't care about the rest

	stat_stream.close();

	long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	vm_usage = vsize / 1024.0;
	resident_set = rss * page_size_kb;
}
