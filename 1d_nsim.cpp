/*	Author: Matthew J. Metson
 *	Created: 27/9/2018
 *	Last modified: 11/10/2018
 *	Description: Simulates N run-and-tumble random walkers on a one-dimensional
 *	 periodic lattice. The number of particles, number of lattice sites, hop rates,
 *	 tumble rates, and time to run the simulation for are all specified by the
 *	 user.
 */

#include <iostream>
#include <stdlib.h>
#include <random>
#include <time.h>
#include <math.h>
#include <fstream>
#include <algorithm>

class Particle {

	public:
	//Constructor (site is randomly-generated site number;
	//sites is total number of lattice sites) and destructor.
	Particle(int site, int sites, double hopRate, double tumbleRate);
	~Particle();

	//Position and direction of particle.
	int pos;
	int dir;

	//Hop and tumble times (respectively).
	double times[2];
	
	//Says whether a particle should be able to hop or not.
	bool noHop;
	
	//Method to determine value of noHop (left and right refer to the neighbouring particles).
	void findNoHop(Particle* left, Particle* right);

	//Functions that generate hop and tumble times, drawn from exponential distributions.
	void hopTime(double t);
	void tumbleTime(double t);

	//Functions that hop and tumble the particles.
	void hop();
	void tumble();

	//Hopping and tumbling rates.
	double r_h;
	double r_t;

	private:
	//Total sites on lattice.
	int totalSites;

};

Particle::Particle(int site, int sites, double hopRate, double tumbleRate) {

	//Initialising total number of sites.
	totalSites = sites;

	//Assigning hopping and tumbling rates.
	r_h = hopRate;
	r_t = tumbleRate;

	//Assigns position (site number).
	pos = site;

	//Assigns a random direction.
	int temp = rand() % 2;
	dir = temp == 0 ? -1 : 1;

	//Assign initial hop and tumble times.
	hopTime(0);
	tumbleTime(0);

}

Particle::~Particle() {}

void Particle::findNoHop(Particle* left, Particle* right) {
	
	//Useful variables.
	int lPos = left -> pos;
	int rPos = right -> pos;
	
	//Booleans that determine whether the neighbouring particles are next door.
	bool lNextDoor = (lPos + 1) % totalSites == pos;
	bool rNextDoor = (pos + 1) % totalSites == rPos;
	
	//Determining whether particle can hop or not.
	if (lNextDoor && rNextDoor) noHop = true;
	else if (lNextDoor && dir < 0) noHop = true;
	else if (rNextDoor && dir > 0) noHop = true;
	else noHop = false;
	
}

void Particle::hopTime(double t) {

	std::random_device rd;
	std::mt19937 gen(rd());
	std::exponential_distribution<> hRand(r_h);
	times[0] = t + hRand(gen);

}

void Particle::hop() {
	
	//Implements the hop on a periodic lattice.
	pos = (pos + dir + totalSites) % totalSites;

}

void Particle::tumbleTime(double t) {

	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::exponential_distribution<> tRand(r_t);
	times[1] = t + tRand(gen);

}

void Particle::tumble() {

	//Implements the tumble.
	dir = -dir;

}

void generateSites(int* siteArr, int numParts, int sites);
void update(Particle** parts, int pIndex, int htIndex, double t, int numParts);

int main() {

	using std::cout;
	using std::endl;
	using std::cin;
	using std::ofstream;
	using std::sort;

	//Time to run simulation for, number of particles, number of sites on lattice, and rates.
	double tMax_temp = 0.;
	int numParts_temp = 0;
	int sites_temp = 0;
	double hopRate_temp = 0.;
	double tumbleRate_temp = 0.;
	
	cout << "Enter t_max: ";
	cin >> tMax_temp;
	const double tMax = tMax_temp;
	
	cout << "Enter the number of particles: ";
	cin >> numParts_temp;
	const int numParts = numParts_temp;
	
	cout << "Enter the number of sites: ";
	cin >> sites_temp;
	const int sites = sites_temp;

	cout << "Enter the hop rate: ";
	cin >> hopRate_temp;
	const double hopRate = hopRate_temp;
	
	cout << "Enter the tumble rate: ";
	cin >> tumbleRate_temp;
	const double tumbleRate = tumbleRate_temp;

	//Time.
	double t = 0.;

	//Randomly generating distinct ordered site numbers to later be assigned as particle positions.
	//Initialisation of all elements to (sites + 1) is due to the site-generation algorithm.
	int siteArr[numParts];
	for (int i = 0; i < numParts; i++) siteArr[i] = sites + 1;
	srand(time(NULL));
	generateSites(siteArr, numParts, sites);	
	sort(siteArr, siteArr + numParts);	
	
	//Creating particles.
	Particle** parts = new Particle*[numParts];
	for (int n = 0; n < numParts; n++) 
		parts[n] = new Particle(siteArr[n], sites, hopRate, tumbleRate);

	//Determining noHop values.
	for (int n = 0; n < numParts; n++) {
		int leftIndex = n == 0 ? numParts - 1 : n - 1;
		int rightIndex = n == numParts - 1 ? 0 : n + 1;
		parts[n] -> findNoHop(parts[leftIndex], parts[rightIndex]);
	}

	//Table.
	ofstream table;
	table.open("table.txt");
	for (int n = 0; n < numParts; n++) table << "\tp" << n;
	table << "\tt" << endl;

	//Loop that hops/tumbles particles till tMax.
	while (t <= tMax) {

		//Looks for smallest t_t or t_h (for noHop == false).
		int pIndex = 0;
		int htIndex = 1;
		double smallest = parts[0] -> times[1];
		for (int n = 1; n < numParts; n++)
			if (parts[n] -> times[1] < smallest) {
				smallest = parts[n] -> times[1];
				pIndex = n;
				htIndex = 1;
			}
		for (int n = 0; n < numParts; n++)
			if (parts[n] -> noHop == false)
				if (parts[n] -> times[0] < smallest) {
					smallest = parts[n] -> times[0];
					pIndex = n;
					htIndex = 0;
				}
	
		//Updating t and carrying out the hop/tumble on the appropriate particle.
		t = smallest;
		update(parts, pIndex, htIndex, t, numParts);

		//Writing to table.
		for (int n = 0; n < numParts; n++) table << "\t" << parts[n] -> pos;
		table << "\t" << t << endl;

	}
	
	table.close();

	return 0;

}

void generateSites(int* siteArr, int numParts, int sites) {
	
	//The ((RAND_MAX - 1) / RAND_MAX) factor ensures we don't hit site = (sites + 1).
	int randNum = 0;
	bool tryAgain = false;
	for (int i = 0; i < numParts; i++) {
		do {
			tryAgain = false;
			randNum = (rand() / (double)RAND_MAX) * sites
					* ((RAND_MAX - 1) / (double)RAND_MAX);
			for (int j = 0; j < i; j++)
				if (siteArr[j] == randNum) {
					j = i;
					tryAgain = true;
				}
		} while (tryAgain);
		siteArr[i] = randNum;
	}
	
}

void update(Particle** parts, int pIndex, int htIndex, double t, int numParts) {

	//Hops/tumbles the appropriate particle and increases t accordingly.
	if (htIndex == 0) {
		parts[pIndex] -> hopTime(t);
		parts[pIndex] -> hop();
	}
	else {
		parts[pIndex] -> tumbleTime(t);
		parts[pIndex] -> tumble();
	}

	//Updates noHop for all possible affected particles.
	//First, find indices for all appropriate particles.
	int leftLeftIndex = (pIndex - 2 + numParts) % numParts;
	int leftIndex = (pIndex - 1 + numParts) % numParts;
	int rightIndex = (pIndex + 1) % numParts;
	int rightRightIndex = (pIndex + 2) % numParts;

	//Then re-evaluate noHop.
	parts[leftIndex] -> findNoHop(parts[leftLeftIndex], parts[pIndex]);
	parts[pIndex] -> findNoHop(parts[leftIndex], parts[rightIndex]);
	parts[rightIndex] -> findNoHop(parts[pIndex], parts[rightRightIndex]);

	//Checks for particles whose hop times need resetting, and resets where appropriate.
	for (int n = 0; n < numParts; n++)
		if (parts[n] -> times[0] < t) parts[n] -> hopTime(t);

}
