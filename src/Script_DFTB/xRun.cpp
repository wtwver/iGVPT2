#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <vector>
#include "Utils.h"
#include "Matrix.h"

#define BFS 1024
using namespace std;

/***********************************************************************************************************************************/
void getNameWorkDir(char* nameWorkDir)
{
        const char* tmpdir = getenv("TMPDIR");
        if(!tmpdir || strstr(tmpdir,"null"))
        sprintf(nameWorkDir,"%s", "/tmp");
        else
        sprintf(nameWorkDir,"%s", getenv("TMPDIR"));
}

/***********************************************************************************************************************************/
void getNameCurrDir(char* nameCurrDir)
{
        sprintf(nameCurrDir,"%s", getenv("PWD"));
}

/***********************************************************************************************************************************/
void addAngularMomentum(fstream& file, vector<string> Symb)
{
        file<<"MaxAngularMomentum = {"<<endl;
        
  	for(int i=0; i<Symb.size(); i++)
  	{
  		if(Symb[i] == string("H")) file<<"H = \"s\""<<endl;
  		if(Symb[i] == string("S")) file<<"S = \"d\""<<endl;
  		if(Symb[i] == string("P")) file<<"P = \"d\""<<endl;
  		if(Symb[i] == string("C")) file<<"C = \"p\""<<endl;
  		if(Symb[i] == string("O")) file<<"O = \"p\""<<endl;
  		if(Symb[i] == string("N")) file<<"N = \"p\""<<endl;
  	}
  	file<<"}"<<endl;

}

/***********************************************************************************************************************************/
void addSpinConstants(fstream& file, vector<string> Symb)
{
	file<<"SpinConstants{"<<endl;
	for(int i = 0; i<Symb.size(); i++)
	{
	if(Symb[i] == string("O")){
	file<<"O = {" << endl;
	file<<"-0.032 -0.028"<<endl;
	file<<"-0.028 -0.027"<<endl;
	file<<"}"<<endl;
	}
	if(Symb[i] == string("H")){
	file<<"H = {" << endl;
	file<<"-0.064 "<<endl;
	file<<"}"<<endl;
	}
	if(Symb[i] == string("C")){
	file<<"C = {" << endl;
	file<<"-0.028 -0.024"<<endl;
	file<<"-0.024 -0.022"<<endl;
	file<<"}"<<endl;
	}
	if(Symb[i] == string("S")){
	file<<"S = {" << endl;
	file<<"-0.019 -0.016 0.000 "<<endl;
	file<<"-0.016 -0.014 0.000"<<endl;
	file<<"0.000 0.000 -0.010"<<endl;
	file<<"}"<<endl;
	}
	if(Symb[i] == string("P")){
	file<<"P = {" << endl;
	file<<"-0.019 -0.016 0.000 "<<endl;
	file<<"-0.016 -0.014 0.000"<<endl;
	file<<"0.000 0.000 -0.010"<<endl;
	file<<"}"<<endl;
	}
	if(Symb[i] == string("N")){
	file<<"N = {" << endl;
	file<<"-0.030 -0.026"<<endl;
	file<<"-0.028 -0.027"<<endl;
	file<<"}"<<endl;
	}
	}
	file<<"}"<<endl;

}
/***********************************************************************************************************************************/
void TypeOfAtoms(vector<string> Atoms, int NbAtoms, vector<string> &Symb)
{
	char buffer[BFS];
        //cout <<"Je cherche a ajouter la geometrie" << endl;
        //vector<string> Symb;
        string test;
        string atom;

        //cout <<"atoms.size() : " << atoms.size() << endl;
        Symb.push_back(Atoms[0].c_str());
        //cout <<"Symb[0] : "<< Symb[0] << endl;
        //cout << "Symb.size() : " << Symb.size() << endl;
        for (int i=0;i<NbAtoms;i++)
        {
                //cout << "i : "<<i<<endl;
                atom = Atoms[i].c_str();
                //cout << "atom : " << atom << endl;
                bool ok = true;
                for(int j=0; j<Symb.size(); j++)
                {
                        //cout << "j : "<<j<<endl;
                        test = Symb[j];
                        //cout <<"test : " <<test << endl;
                        if(Symb.size()> NbAtoms) { ok = false; break;}
                        if( atom == test) {ok = false; break;}
                }//cout<<"atom : " << atom << endl;
                //cout << ok << endl;
                if( ok)
                {
                        Symb.push_back( atom );

                        //cout <<"Symb push_back : "<< Symb[nb-1] << endl;
                }
        }
}

/***********************************************************************************************************************************/
void addXYZGeometryToInputDFTBFile(fstream& file, int NbAtoms, vector<string> Atoms, vector<double> Coords, vector<string> Symb)
{
	char buffer[BFS];
        int index[Symb.size()];
        //cout << "Symb.size() : " << Symb.size() << endl;
        for(int i=0; i<Symb.size(); i++) index[i] = i+1;
        
        for(int i=0; i<Symb.size(); i++) file<<Symb[i].c_str()<<" ";
        //for(int i=0; i<Symb.size(); i++) cout<<"Symb : "<<Symb[i].c_str()<<endl;
        file<<endl;
        int indice = 0;
        for (int i=0;i<NbAtoms;i++)
        {
        	for(int j=0; j<Symb.size(); j++)
        	{
        		if(Atoms[i].c_str() == Symb[j]) indice = j;
                }
                sprintf(buffer,"%d %d %20.14f %20.14f %20.14f\n", i + 1, index[indice], Coords[3*i], Coords[3*i+1], Coords[3*i+2]);
                file<<buffer;
        }
	file << endl << "}" << endl;
        //cout<<"J'ai termine de remplir la geo!" << endl;
}

/***********************************************************************************************************************************/

int main(int argc, char *argv[])
{
	char buffer[BFS];
        char nameWorkDir[BFS];
        char nameCurrDir[BFS];
	int NbOfAtoms = 0;
	int ChargeMol = 0;
	int Multiplicity = 0;
	fstream file;
	fstream fileInput;
	fstream fileOutput;
	vector<string> Atoms;
	
	string Symbol;
	vector<string> MMType;
	vector<string> pdbType;
	vector<string> residueName;
	vector<int> numResidue ;
	vector<double> charge ;
	vector<int> layer ;
	vector<int> optAtom ;
	double XAtom = 0;
	double YAtom = 0;
	double ZAtom = 0;
	vector<int> NbConnect ;
	int Option;
	vector<double> Coordinates;
	vector<string> Connection;
	char* inpName = NULL;
	char* outName = NULL;
	char *optionsName = NULL;
	char* pref = NULL;
	char inpDFTB[BFS];
	char outDFTB[BFS];

	if(argc<3)
	{
		cerr<<" Il faut fournir le nom du fichier d'entree ainsi que le nom du fichier de sortie"<<endl;
		return 1;
	}
	optionsName = argv[1];
	inpName = argv[2];
	outName = argv[3];

	//Lecture du fichier Options.txt
	string Param;
	fstream fileOptions;
	fileOptions.open(optionsName, ios::in);
	if(fileOptions.fail()) cerr << "Impossible de lire le fichier Options.txt " << endl;
	fileOptions >> Param;
	fileOptions.close();
	//Fin de la lecture Options.txt

	file.open(inpName, ios::in);
	if(file.fail()) cerr<< "Impossible d'ouvrir le fichier inp.ici" << endl;
	{
		file.getline(buffer, BFS);	
		vector<string> line = split(buffer);
		Option = atoi(line[0].c_str());

		file.getline(buffer, BFS);	// Le texte Geometry
		file.getline(buffer, BFS);	//Contient le nombre d'atomes, la charge, la multiplicité
		line = split(buffer);
		NbOfAtoms = atoi(line[0].c_str());
		ChargeMol = atoi(line[1].c_str());
		Multiplicity = atoi(line[2].c_str());
		for(int i=0; i< NbOfAtoms; i++)
		{
			file.getline(buffer, BFS);
			line = split(buffer);
			Symbol = line[0];
			MMType.push_back(line[1]);
			pdbType.push_back(line[2]);
			residueName.push_back(line[3]);
			numResidue.push_back(atoi(line[4].c_str()));
			charge.push_back(atof(line[5].c_str()));
			layer.push_back(atoi(line[6].c_str()));
			optAtom.push_back(atoi(line[7].c_str()));
			NbConnect.push_back(atoi(line[11].c_str()));
			Atoms.push_back(Symbol);
			string Connec;
			for(int j=0; j<2*NbConnect[i]; j++)
			{
				Connec += line[12+j];
				Connec += " "; 
			}
			Connection.push_back(Connec);
			XAtom = atof(line[8].c_str());
			YAtom = atof(line[9].c_str());
			ZAtom = atof(line[10].c_str());
			Coordinates.push_back(XAtom);
			Coordinates.push_back(YAtom);
			Coordinates.push_back(ZAtom);
			line = split(buffer);
		}
	}	
	file.close();
			
	//Création du fichier input de dftb
		vector<string> Symb;
		TypeOfAtoms(Atoms, NbOfAtoms, Symb);

		pref = getSuffixNameFile(inpName);
		sprintf(inpDFTB,"%sDFTB.hsd", pref);
		sprintf(outDFTB,"%sDFTB.out", pref);

		fileInput.open(inpDFTB, ios::out);
		if(fileInput.fail()) cerr << "Impossible de créer le fichier dftb_in.hsd" << endl; 
		fileInput << "Geometry = GenFormat {"<<endl;
		fileInput << NbOfAtoms << " C" << endl;
		addXYZGeometryToInputDFTBFile(fileInput, NbOfAtoms, Atoms, Coordinates, Symb);
		fileInput << endl;
		if(Option == 2)
		{
			fileInput << "Driver = ConjugateGradient{" << endl;
			fileInput << "MovedAtoms = 1:-1" << endl;
 			fileInput << "MaxForceComponent = 1.000000000000000E-006  # Extremely small!" << endl;
			fileInput << " MaxSteps = 1000" << endl;
 			fileInput << "OutputPrefix = \"geom.out\""<< endl << "}" << endl;
		}
		fileInput << "Hamiltonian = DFTB {" << endl;
		fileInput<<"Charge = "<<ChargeMol<<endl;
	        fileInput<<"SCC = Yes"<<endl;
        	fileInput<<"SCCTolerance = 1.0E-6 # Extremely small!"<<endl;
	        fileInput<<"MaxSCCIterations = 10000"<<endl;
	        fileInput<<"Mixer = DIIS {"<<endl;
	        fileInput<<"InitMixingParameter = 0.2"<<endl;
	        fileInput<<"Generations = 3"<<endl;
	        fileInput<<"UseFromStart = Yes"<<endl<<"}"<<endl;
        	fileInput<<"SpinPolarisation = Colinear{" << endl;
        	fileInput<<"UnpairedElectrons = "<< Multiplicity - 1  << endl;
        	fileInput<<"}"<<endl;
        	addSpinConstants(fileInput, Symb);
        	addAngularMomentum(fileInput, Symb);
        	fileInput<<"Filling = Fermi {"<<endl;
        	fileInput<<"Temperature [Kelvin] = 300"<<endl;
        	fileInput<<"}"<<endl;
        	fileInput<<"SlaterKosterFiles = Type2FileNames {    # File names from two atom type names"<<endl;
	        fileInput<<"Prefix = \""<<Param<<"\"  # Path as prefix"<<endl;
        	//fileInput << "Prefix = \"/softs/CompChemPackages/DFTBP/parameters/mio-1-1/\"  # Path as prefix" << endl;
        	fileInput<<"Separator = \"-\"         # Dash between type names"<<endl;
        	fileInput<<"Suffix = \".skf\"         # Suffix after second type name"<<endl;
        	fileInput<<"}"<<endl;
        	fileInput<<"}"<<endl;
		if(Option == 1)
		{
			fileInput<<"Options{"<<endl<<"CalculateForces = Yes"<< endl << "}"<<endl;
		}
        	fileInput<<"ParserOptions = {"<<endl;
        	fileInput<<"ParserVersion = 3"<<endl;
        	fileInput<<"}"<<endl;
	fileInput.close();
	//Fin de la création du fichier Input

	//Lancement du calcul
        sprintf(buffer,"runCCDFTB %s",inpDFTB);	
	system(buffer);

	//Creation du nouveau fichier Input
	fileOutput.open(outDFTB, ios::in);
	if(fileOutput.fail()) cerr << "Impossible de lire le fichier "<<outDFTB<< endl;
	double dipole[4];
	double Energy = 0;
	Matrix<double> Coord = Matrix<double> (3, NbOfAtoms);
	Matrix<double> CoordOpt = Matrix<double> (3, NbOfAtoms);
	Matrix<double> Forces = Matrix<double> (NbOfAtoms, 3);
        vector<double> Charge;
	for(int i=0; i<4; i++) dipole[i] = 0;
	while(!fileOutput.eof())
	{
		fileOutput.getline(buffer,BFS);
		if(fileOutput.fail()) break;
		if(strstr(buffer,"Total Energy:"))
		{
			vector<string> line = split(buffer);
			Energy = atof(line[2].c_str());
		}

		if(strstr(buffer,"GEOMETRY") && Option == 2)
		{
			fileOutput.getline(buffer, BFS);
			//cout << "buffer : " << buffer << endl;
			fileOutput.getline(buffer, BFS);
			//cout << "buffer : " << buffer << endl;
			//fileOutput.getline(buffer, BFS);
			//cout << "buffer : " << buffer << endl;

			vector<string> line = split(buffer);
			for(int j=0; j<NbOfAtoms; j++)
			{
				fileOutput.getline(buffer,BFS);
				//cout << "buffer : " << buffer << endl;
				line = split(buffer);
				for(int i = 2; i< line.size(); i++)
				{
					Coord(i-2,j) = atof(line[i].c_str());
					CoordOpt(i-2,j) = atof(line[i].c_str());
				}
			}
		
		}
		if(Option == 0 or Option == 1)
		{	
			for(int i=0; i< NbOfAtoms; i++)
			{
				Coord(0, i) = Coordinates[3*i];
				Coord(1, i) = Coordinates[3*i+1];
				Coord(2, i) = Coordinates[3*i+2];
			}
		}
		if(strstr(buffer," Total Forces") && Option == 1)
		{
			vector<string> line = split(buffer);
			int k=0; 
			while(!line.empty()){
				fileOutput.getline(buffer,BFS);
				line = split(buffer);
				for(int i = 0; i< line.size(); i++)
				{
					Forces(k,i) = atof(line[i].c_str());
				}
				k++;
				}
			k-=1;
		
		}
		if(strstr(buffer,"Net atomic "))
		{
			vector<string> line = split(buffer);
                       	fileOutput.getline(buffer,BFS);
                        for(int i=0; i<Coord.getNumberOfColumns(); i++)
			{
                        	fileOutput.getline(buffer,BFS);
                        	line = split(buffer);
				Charge.push_back(atof(line[1].c_str()));	
                        }
		for(int j = 0; j<Coord.getNumberOfColumns(); j++)
		{
			for(int l = 0; l<Coord.getNumberOfRows(); l++)
			{
				Coord(l,j) = Coord(l,j) * Charge[j];
			}
		}
		for(int i = 0; i<Coord.getNumberOfRows(); i++)
		{
			for(int j = 0; j<Coord.getNumberOfColumns(); j++)
               		{
                        	dipole[i] += Coord(i,j);
                        	
                	}
		} 
		/*for(int i =0; i<3; i++)
		{
		properties.dipole[i] *= 2.54174619;
		}*/
		dipole[3] = sqrt(dipole[0]*dipole[0]+dipole[1]*dipole[1]+dipole[2]*dipole[2]);

		} 
	}
	fileOutput.close();

	fstream Input;
	Input.open(outName, ios::out);
	Input << setprecision(14) << scientific << Energy << endl;
	Input << setprecision(14) << scientific << dipole[0] << "	  " << dipole[1] << "   " << dipole[2] << endl;
	if(Option == 1)
        for(int i=0; i<NbOfAtoms; i++)
	{
		Input << setprecision(14) << scientific << Forces(i, 0) << "  " << Forces(i, 1) << "  " <<  Forces(i, 2) << "  "<< endl;
	}	

	if(Option == 2)
	{
		Input << "Geometry" << endl;
		Input << NbOfAtoms << " " << ChargeMol << " " <<  Multiplicity << endl;
		for(int i=0; i<NbOfAtoms; i++)
		{
			Input << setprecision(14) << scientific  << Atoms[i] << " " << MMType[i] << " " <<  pdbType[i] << " " << residueName[i] << " " << numResidue[i] << " " << charge[i] << " " << layer[i] << " " << optAtom[i] << " " << CoordOpt(0 ,i) << " " << CoordOpt(1, i) << " " << CoordOpt(2, i) << " " << NbConnect[i] << " " << Connection[i] << endl;
		}

	}
	Input.close();

	return 0;
}
