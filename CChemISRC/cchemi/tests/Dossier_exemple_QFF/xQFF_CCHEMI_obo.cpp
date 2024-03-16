#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>


using namespace std;

#define BFS 1024




int main(int argc, char *argv[])
{
	string ListInput = argv[1];
	string ListOutput = argv[2];
	string Output = argv[3];
	vector<double> Freq;
	vector<double> Mass;
	vector<double> Delta;
	vector<double> Energy;
	vector<double> Dipole_X;
	vector<double> Dipole_Y;
	vector<double> Dipole_Z;
	//double ad = 2.54174619;
	double ad = 2.54158059;
	//double hk = 627.50944796;
	double hk = 627.509544796;

	string File ="";
	char buffer[1024];
	int n = 0;
	int k = 0;
	string already = "test";
	
	fstream fileListInput;
	fileListInput.open(ListInput.c_str(), ios::in);
	if(fileListInput.fail()) cerr << "Error while opening  " << ListInput << endl;
	while(!fileListInput.eof())
	{
		fileListInput >> File;
		cout << File << endl;
		if(already == File) break;
		already = File;
		if(n == 6) n = 0;
		fstream file;
	        file.open(File.c_str(), ios::in);
	        if(file.fail()) cerr << "Erreur ouverture  " << File << endl;
	        while(!file.eof())
		{
			file.getline(buffer,1024);
			if(strstr(buffer,"Mode:") && n == 3)
			{
				char dum[1024];
				char freq[1024];
				char mass[1024];
				char delta[1024];
			
				int nb = sscanf(buffer, "%s %s %s %s %s %s %s %s %s %s",&dum, &dum, &freq, &dum, &mass, &dum, &dum, &dum, &delta, &dum);
				cout << nb << endl;
				if(nb == 10) 
				{
					double tmp = atof(freq);
					cout << tmp << endl;
					Freq.push_back(tmp);
					tmp = atof(mass);
					cout << tmp << endl;
					Mass.push_back(tmp);
					tmp = atof(delta);
					cout << tmp << endl;
					Delta.push_back(tmp);
				}
				else
				{
					cerr << "Error : Impossible to read the mode! Pb here : " << File << endl; 
					break;
				}
			}
		}
		n++;
	} 	

	string name = "";
	string unit = "";
	int nb = 0;	
	int nbp = 0;	
	bool Ener = false;
	bool Dipole = false;
	int index_Ener = 0;
	int index_Dip = 0;
	int nb_geo = 0;
	int nb_props = 0;
	int nb_ener = 0;
	int nb_dip = 0;
	double energy = 0;
	double px = 0;
	double py = 0;
	double pz = 0;
	int tempi = 0;

	vector<double> Nombre;

        fstream fileListOutput;
        fileListOutput.open(ListOutput.c_str(), ios::in);
        if(fileListOutput.fail()) cerr << "Erreur ouverture  " << ListOutput << endl;
        while(!fileListOutput.eof())
        {
                fileListOutput >> File;
		cout << File << endl;
                if(already == File) break;
                already = File;
                fstream file;
                file.open(File.c_str(), ios::in);
                if(file.fail()) cerr << "Erreur ouverture  " << File << endl;
                while(!file.eof())
                {
	                file.getline(buffer, BFS);
         	        if(strstr(buffer, "[GEOMS]" ))
	                {
				cout << buffer << endl;
	                	//file.getline(buffer, BFS);
        	                file >> nb_geo  >> nb_props;
				//nb = sscanf(buffer, "%d  %d ",&nb_geo, &nb_props);
				//cout << nb << endl;
                        	cout << nb_geo << "   " << nb_props << endl;
	                        if(nb_geo == 1)
        	                {
                	                for(int i = 0; i < nb_props; i++)
                        	        {
	                			//file.getline(buffer, BFS);
						//cout << buffer << endl;
						//nb = sscanf(buffer, "%s  %s %d",&name, &unit, &nbp);
						//cout << nb << endl;
						//cout << unit << endl;
        	                                file >> name >> unit >> nbp;
						cout << name << " " << unit << " " << nbp << endl;
						Nombre.push_back(nbp);
        	                                if(name == "energy")
                	                        {
                        	                        Ener = true;
                                	                index_Ener = i;
							nb_ener = nbp;
                                        	}
	                                        if(name == "Dipole")
        	                                {
                	                                Dipole = true;
                        	                        index_Dip = i;
							nb_dip = nbp;
                                	        }
	                                }
        	                        for(int i = 0; i < nb_props; i++)
                	                {
						//file.getline(buffer, BFS);
        	                                if(Ener && i == index_Ener)
                	                        {
							cout << "ENERGY" << endl;
							//if(nb_ener == 1) 
							{
								file >> energy;
								//file.getline(buffer, BFS);
								//cout << buffer << endl;
								//nb = sscanf(buffer, "%12.8f ",&energy);
								//cout << nb << endl;
								cout << setprecision(30) << energy << endl;
								energy /= hk;
								cout << setprecision(30) << energy << endl;
								if(fabs(energy) < 1e-08) cout << "WARNING : An energy is equal to zero in : " << File << endl;
								Energy.push_back(energy);
							}
                	                        }
        	                                else if(Dipole && i == index_Dip)
                	                        {
							cout << "DIPOLE" << endl;
							//if(nb_dip == 3) 
							{
								//file.getline(buffer, BFS);
								//cout << buffer << endl;
								file >> px >> py >> pz;
								//nb = sscanf(buffer, "%12.8f %12.8f %12.8f",&px, &py, &pz);
								//cout << nb << endl;
								px /= ad;
								py /= ad;
								pz /= ad;
								cout << px << endl;
								cout << py << endl;
								cout << pz << endl;
								if(fabs(px) < 1e-08 && fabs(py) < 1e-08 && fabs(pz) < 1e-08) cout << "WARNING : A dipole is equal to zero in : " << File << endl;
								Dipole_X.push_back(px);
								Dipole_Y.push_back(py);
								Dipole_Z.push_back(pz);
							}
                	                        }
						else
						{
							cout << "AUTRE" << endl;
							string can;
							for(int l = 0; l < Nombre[i]; l++) file >> can;
						}
						
					}
				}
			}

		}	
		k++;
	} 	
	
	cout << Freq.size() << endl;
	cout << Mass.size() << endl;
	cout << Delta.size() << endl;
	cout << Energy.size() << endl;
	cout << Dipole_X.size() << endl;

	FILE *outputfile;	
	//fstream outputfile;
	outputfile = fopen(Output.c_str(), "w");

	fprintf(outputfile,"Frequencies [cm^-1]\n"); 
	for(int i=0;i<Freq.size();i++) fprintf(outputfile, "%0.12f\n",Freq[i]);

	fprintf(outputfile,"Mass [uma]\n"); 
        for(int i=0;i<Mass.size();i++) fprintf(outputfile,"%0.12f\n",Mass[i]);

	fprintf(outputfile,"Delta [Bohr]\n");
        for(int i=0;i<Delta.size();i++) fprintf(outputfile,"%0.12f\n",Delta[i]);

	fprintf(outputfile,"Energies [Hartree]\n");
	for(int i=0;i<Energy.size();i++) fprintf(outputfile,"%0.12f\n",Energy[i]);

	fprintf(outputfile,"Dipoles [AU]\n");
	for(int i=0;i<Dipole_X.size();i++) fprintf(outputfile,"%0.12f %0.12f %0.12f\n",Dipole_X[i], Dipole_Y[i], Dipole_Z[i]);
	fclose(outputfile);
	
		
	if(Energy.size() != k) cout << "Warning : the number of energies is not equal to the number of files in the list!" << endl;
	if(Dipole_X.size() != k) cout << "Warning : the number of dipoles is not equal to the number of files in the list!" << endl;
	return 0;
}
