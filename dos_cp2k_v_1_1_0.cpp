/**************************************************************************
 * Program for postprocessing of cp2k pdos files                          *
 * By Antonio Claudio Padilha - claudio.padilha@york.ac.uk                *
 *                                                                        *
 * version 1.1.0, 21th July 2017                                          *
 *                                                                        *
 *  - dependencies: tclap Command line parser (tclap.sourceforge.net)     *
 *                                     -- known to work on version 1.2.1  *
 **************************************************************************/

#include <fstream>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstring>
#include <sstream>
#include <tclap/CmdLine.h>

#define pi 3.14159265
#define au2ev 2.72113838565563E+01
#define version "1.1.0"

using namespace std;

class sDos {
private:
	float * en; // energies
	float * dos_t; // total dos
	float s; // sigma
	float * eig; // original vector w eigenenergies
	int * mos;
	float * occ;
	float ** p_ml; // mls from file
	float ** dos_ml; // adjusted projected dos
	float ** dos_l; // just sums ups px+py+pz=p, etc.
	int _n;
	int _m;
	int n_ml; // # of all projections with m
	int n_l;  // # of all projections (px+py+pz=p, etc. )
	float ef;
	char * at_type;
	string * mls;
	string ls = "spdf";
	float dE;
	void adjust ();
	string fin;
	string fout;
	int count_lines ();
	float find_ef ();
	char * find_kind ();
	int find_ml ();
	float _un;
	sDos ();
	void read_pdos ();
	float g (float x, float x0, float s);
public:
	sDos (const float mult, const float sigma, string f, const float un); // constructor
	~sDos (); // destructor

	void print2DataFile ();
	void print2GnuplotFile ();
};

// sDos::sDos (const float mult, const float sigma, char * filein, char * fileout, const float un)
// : s(sigma), _m(floor(mult)), fin(filein), fout(fileout), _un(un)
sDos::sDos (const float mult, const float sigma, string f, const float un)
: s(sigma), _m(floor(mult)), fin(f), fout(f.substr(0,f.find_last_of('.'))+".out"), _un(un)
{
	_n = count_lines();
	ef = find_ef() * _un;
	at_type = find_kind(); 
	n_ml = find_ml();

	read_pdos();

	en = new float [_n * _m];
	dos_t = new float [_n * _m];

	dos_ml = new float * [_n * _m];
	for (int i = 0; i < _n * _m; i++)
		dos_ml[i] = new float [n_ml];

	float min = eig[0];
	float max = eig[_n - 1];
	dE = (max - min) / (_n * _m);

	for (int i = 0; i < _n * _m; i++)
	{
		en[i] = min + i * dE;
		dos_t[i] = 0.0;
		for (int j = 0; j < n_ml; j++) dos_ml[i][j] = 0.0;
	}

	adjust();

	cout << "Atomic species: " << at_type << endl;
	cout << "Fermi energy = " << ef << endl;
	// take care of projections ml -> l
	if (n_ml == 1)
	{
		puts ("only s");
		n_l = 1;
		dos_l = new float * [_n * _m];
		for (int i = 0; i < _n * _m; i++)
		{
			dos_l[i] = new float [n_l];
			// s is the same
			dos_l[i][0] = dos_ml[i][0];
		}
	}
	else if (n_ml == 4)
	{
		puts ("s + p");
		n_l = 2;
		dos_l = new float * [_n * _m];
		
		// initialize 
		for (int i = 0; i < _n * _m; i++)
		{
			dos_l[i] = new float [n_l];
			// s is the same
			dos_l[i][0] = dos_ml[i][0];
			for (int j = 1; j < n_l; j++)
				dos_l[i][j] = 0.0;
			// deal with p
			for (int j = 1; j < 4; j++)
				dos_l[i][1] += dos_ml[i][j];
		}
	}
	else if (n_ml == 9)
	{
		puts ("s + p + d");

		n_l = 3;
		dos_l = new float * [_n * _m];
		// initialize 

		for (int i = 0; i < _n * _m; i++)
		{
			dos_l[i] = new float [n_l];
			// s is the same
			dos_l[i][0] = dos_ml[i][0];
			for (int j = 1; j < n_l; j++)
				dos_l[i][j] = 0.0;
			// deal with p
			for (int j = 1; j < 4; j++)
				dos_l[i][1] += dos_ml[i][j];
			// deal with d
			for (int j = 4; j < 9; j++)
				dos_l[i][2] += dos_ml[i][j];
		}
	}
	else 
	{
		puts ("s + p + d + f");
		n_l = 4;
		dos_l = new float * [_n * _m];
		// initialize 
		for (int i = 0; i < _n * _m; i++)
		{
			dos_l[i] = new float [n_l];
			// s is the same
			dos_l[i][0] = dos_ml[i][0];
			for (int j = 1; j < n_l; j++)
				dos_l[i][j] = 0.0;
			// deal with p
			for (int j = 1; j < 4; j++)
				dos_l[i][1] += dos_ml[i][j];
			// deal with d
			for (int j = 4; j < 9; j++)
				dos_l[i][2] += dos_ml[i][j];
			// deal with f
			for (int j = 9; j < 16; j++)
				dos_l[i][3] += dos_ml[i][j];
		}
	}
}

sDos::~sDos() {
	delete [] en;
	delete [] dos_t;
	delete [] eig;
	delete [] mos;
	delete [] occ;
	for (int j = 0; j < n_ml; j++) delete [] dos_ml[j];
	delete [] dos_ml;
}

// gaussian function. returns g(x) with mean x0 and sigma s
float sDos::g (float x, float x0, float s)
{
	return exp (-pow(x-x0,2)/(2*pow(s,2))) / sqrt (2 * pi * pow(s,2));
}

void sDos::read_pdos()
{
	ifstream data (fin);
	eig = new float [_n]; // eigenvalue energy
	mos = new int [_n];   // eigenvalue molecular orbital
	occ = new float [_n]; // eigenvalue occupation
	
	// angular momenta projections for each energy
	p_ml = new float * [_n]; 
	for (int i = 0; i < _n; i++) p_ml[i] = new float [n_ml];

	char buffer[1000];
	string temp;

	// read .pdos file
	if (data.is_open()) 
	{
		// ignore first two lines: they are the header
		getline(data,temp);
		getline(data,temp);
		// loop over the eigenvalues
		for (int i = 0; i < _n; i++)
		{
			data >> buffer;
			sscanf (buffer, "%d", &mos[i]);
			data >> buffer;
			sscanf (buffer, "%f", &eig[i]);
			eig[i]*=_un;
			data >> buffer;
			sscanf (buffer, "%f", &occ[i]);
			// loop over the angular momenta read from the header
			for (int j = 0; j < n_ml; j++)
			{
				data >> buffer;
				sscanf (buffer, "%f", &p_ml[i][j]);
			}
		}
	}
}

void sDos::adjust ()
{
	// l is the index of the energy grid point that is 4 * sigma 
	// from the energy in the old grid (the mean value for the gaussian
	// in the new grid) 
	int l = floor(4*s/dE);
	// this is the index of the first energy on the new grid (and also
	// in the old grid). We will use it as starting point for the search
	// for all the gaussian centres ahead
	int i0 = 0;

	// we go over all energies in the old grid, each one is the mean
	// value for the gaussians in the new grid
	for (int i = 0; i < _n; i++)
	{
		// here we are searching for the next energy in the old grid eig
		// by running over the new grid en. Once we are further than the
		// old energy (eig[i]-en[j] becomes < 0) we stop and then our 
		// serach for the next mean value will start from here (i0 = j)
		for (int j = i0; j < _n * _m && eig[i]-en[j] >= 0; j++) i0 = j;
				
		// we check if we are at the borders of the grid. If i0 - l < 0 we are
		// starting before the start of en, and if i0 + l > _n * _m - 1 we are
		// further than the end of en, so we set the values to those limits
		int ii = ((i0 - l) < 0 ? 0 : i0 - l);
		int ie = ((i0 + l) > (_n * _m - 1) ? _n * _m - 1 : i0 + l);
		
		// here we go from mean - 3 * sigma to mean + 4 * sigma to calculate our
		// gaussians and add them to the new grid
		for (int j = ii; j <= ie; j++)
		{
			// we add the gaussians on top of each other on the new grid 
			float temp = g (en[j], en[i0], s);
			dos_t[j] += temp;
			// here we store the projections, which are just (agnular momenta) * gaussian
			for (int k = 0; k < n_ml; k++)
				dos_ml[j][k] += temp * p_ml[i][k];
		}
	}
}

void sDos::print2DataFile ()
{
	ofstream outfile (fout);

	if (outfile.is_open())
	{
		// check if units are eV or a.u. and print the corresponding header
		if (abs(_un-1)<1e-4) outfile << setw(10) << "#Eig[a.u.]" << setw(14) << "Total Dos";
		else  outfile << setw(10) << "#Eig[eV]" << setw(14) << "Total Dos";

		// complete header with angular momenta
		for (int i = 0; i < n_l; i++) outfile << setw(14) << ls[i];
		
		// end of header
		outfile << endl;

		// for each energy on the new grid
		for (int i = 0; i < _n * _m; i++)
		{
			// print the energy and the total dos
			outfile << setw(10) << setprecision(6) << en[i] << setw(14) <<setprecision(6) << dos_t[i];
			// print the projections
			for (int j = 0; j < n_l; j++)
				outfile << setw(14) <<setprecision(6) << dos_l[i][j];
			outfile << endl;
		}
	// tell user where to find the data
	cout << "Data recorded in file " << fout << endl;
	}
	
	// close the data file created
	outfile.close();
}

void sDos::print2GnuplotFile ()
{
	// script file is [inputfile - pdos + gnu]
	string gnuout = fin.substr(0,fin.find_last_of('.'))+".gnu";
	// image file is [inputfile - pdos + eps]
	string epsout = fin.substr(0,fin.find_last_of('.'))+".eps";

	// just dump all that into the gnuplot script file
	ofstream outfile (gnuout);
	if (outfile.is_open())
	{
		outfile << "set terminal postscript eps color enhanced font 24" << endl;
		outfile << "set output \'"<< epsout << "\'" << endl;
		outfile << "set border lw 3" << endl;
		outfile << "set size ratio 0.3" << endl;
		outfile << "#xmin = -6.0" << endl;
		outfile << "#xmax =  6.0" << endl;
		outfile << "#set xrange [xmin:xmax]" << endl;
		outfile << "#set xtics xmin xmax 1.0" << endl;
		outfile << "set format x \"%.1f\"" << endl;
		outfile << "#unset ytics" << endl;
		outfile << "#ymin = -60" << endl;
		outfile << "#ymax =  60" << endl;
		outfile << "#set yrange [ymin:ymax]" << endl;
		outfile << "set style line 1 lt 2 lw 2 lc 0" << endl;
		outfile << "set style arrow 1 nohead ls 1" << endl;
		outfile << "set xlabel \"E-E_F("<< (abs(_un-1)<1e-4?"a.u.":"eV") <<")\"" << endl;
		outfile << "set ylabel \"DoS\"" << endl;
		outfile << "Ef = " << ef << endl;
		outfile << "#set arrow from 0.0, ymin to 0.0, ymax as 1" << endl;
		outfile << "#set label \"" << at_type << "\" at 0.1 * (xmax-xmin) + xmin, ymax - 0.1 * (ymax-ymin) font \",24\"" << endl;
		outfile << "plot \t \"" << fout << "\" u ($1-Ef):2 \t w l lw 3 lt 1 lc 0 \t title \'Total\', \\" << endl;
		for (int i = 0; i < n_l; i++)
			outfile << "\t \"" << fout << "\" u ($1-Ef):" << i + 3 << " \t w l lw 3 lt 1 lc " << i + 2 << " \t title \'" << ls[i] << "\', \\" << endl;

		// tell the user what to do
		cout << "Use script " << gnuout << " to generate the graph using gnuplot" << endl;
	}
}

// count all lines in a file
int sDos::count_lines ()
{
	string line;
	int n = 0;
	ifstream infile (fin);

	// first 2 lines of cp2k.pdos file are comments, do not count
	while (getline(infile, line))
		if (line[0] != '#') n++;

	infile.close();
	return n;
}

// returns fermi energy from header of pdos file
float sDos::find_ef ()
{
	string line;
	char * buff1, * buff2;

	ifstream infile (fin);

	getline(infile,line);

	infile.close();

	buff1 = new char [line.size()+1];
	strcpy(buff1,line.c_str());
	buff2 = strstr (buff1, "Fermi");
	delete buff1;

	buff1 = new char[9];
	strncpy (buff1,&buff2[13],8);
	buff1[9] = '\0';
	
	return strtof (buff1,NULL);
}

// returns atomic species from header of pdos file
char * sDos::find_kind ()
{
	string line;
	char * buff1, * buff2;
	
	// open input file
	ifstream infile (fin);

	// we need the first line 
	getline(infile,line);

	// see ya input file
	infile.close();


	buff1 = new char [line.size()+1];
	strcpy(buff1,line.c_str());
	buff2 = strstr (buff1, "kind");
	// puts (buff2);
	delete buff1;
	buff1 = new char [2];
	strncpy (buff1, &buff2[5], 3);
	buff1[2] = '\0';

	return buff1;
}

// returns angular momenta from header of pdos file
int sDos::find_ml ()
{
	string line;
	int i, n;

	// open input file
	ifstream infile (fin);

	// the stuff we need is in line 2
	getline(infile,line);
	getline(infile,line);
	
	// we don't need the input file anymore
	infile.close();

	// lets parse the second line
	istringstream iss1 (line);

	// we iterate over the stream to find its length
	for (i = 0; iss1; i++)
	{
		string sub;
		iss1 >> sub;
	}

	// n is how many angular momenta we found
	// we subtract 6 since this is the number of
	// strings that appear first in the line before s, p, d, f
	n = i-6;
	mls = new string[n];

	// iterate over the line once more
	istringstream iss2 (line);
	string buff;
	// throw away the first entries
	for (int i = 0; i < 5; i++) iss2 >> buff;
	// that is the part where s, p, d, f are
	for (i = 0; i < n; i++) iss2 >> mls[i];

	// return the number of angular momenta found
	return n;
}

int main (int argc, char ** argv)
{
	float sigma, con, grid;
	
	// adapted from 
	// https://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
	// Wrap everything in a try block.  Do this every time,
    // because exceptions will be thrown for problems.
    try 
    {
	    // Define the command line object, and insert a message
	    // that describes the program. The "Command description message"
	    // is printed last in the help text. The second argument is the
	    // delimiter (usually space) and the last one is the version number.
	    // The CmdLine object parses the argv array based on the Arg objects
	    // that it contains.
	    TCLAP::CmdLine cmd("This program generates cp2k smooth pdos", ' ', version);

	    // Define a value argument and add it to the command line.
	    // A value arg defines a flag and a type of value that it expects,
	    // such as "-n Bishop".
	    TCLAP::ValueArg<std::string> sigmaArg("s","sigma","Gaussian half-width (default 0.001)",false,"0.001","float");
		TCLAP::ValueArg<std::string> unitArg("u","unit","Which unit to use? (eV or a.u. - default a.u.)",false,"a.u.","string");
		TCLAP::ValueArg<std::string> gridArg("g","grid","Multiply grid by how much? (default 1E3)",false,"1E3","int");
		TCLAP::UnlabeledMultiArg<string> filesArg("filenames", "the *.pdos files to be used as input", true, "input files");
	    // Add the argument nameArg to the CmdLine object. The CmdLine object
	    // uses this Arg to parse the command line.
	    cmd.add (sigmaArg);
	    cmd.add (unitArg);
	    cmd.add(gridArg);
	    cmd.add(filesArg);

	    // Define a switch and add it to the command line.
	    // A switch arg is a boolean argument and only defines a flag that
	    // indicates true or false.  In this example the SwitchArg adds itself
	    // to the CmdLine object as part of the constructor.  This eliminates
	    // the need to call the cmd.add() method.  All args have support in
	    // their constructors to add themselves directly to the CmdLine object.
	    // It doesn't matter which idiom you choose, they accomplish the same thing.
	    //TCLAP::SwitchArg reverseSwitch("r","reverse","Print name backwards", cmd, false);

	    // Parse the argv array.
	    cmd.parse( argc, argv );

	    // Get the value parsed by each arg.
	    sigma = strtof (sigmaArg.getValue().c_str(),NULL);
	    string unit (unitArg.getValue());
	    con = (!unit.empty() && strcmp(unit.c_str(), "eV") * strcmp(unit.c_str(), "ev") * strcmp(unit.c_str(), "EV") * strcmp(unit.c_str(), "Ev") == 0 ? au2ev : 1);
	    grid = strtof (gridArg.getValue().c_str(),NULL);
	    vector<string> fileNames = filesArg.getValue();
		
		// loop over filenames 
		for (string c : fileNames)
	    {
			cout << endl << "** input file: " << c << " **" << endl;

			// create a sDos object for each filename
			sDos * d = new sDos (grid, sigma, c, con);

			// creates the outputs for each filename from the sDos objects
			d -> print2DataFile ();
			d -> print2GnuplotFile ();
	    }


	} 
	catch (TCLAP::ArgException &e)  // catch any exceptions
    { 
    	cerr << "error: " << e.error() << " for arg " << e.argId() << endl; 
    }

	
	
	return 0;
}
