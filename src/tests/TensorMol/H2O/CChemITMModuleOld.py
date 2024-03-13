#!/usr/local/bin/python3.6
from __future__ import absolute_import
from __future__ import print_function
from TensorMol import *
import sys

os.environ["CUDA_VISIBLE_DEVICES"]="" # set to use CPU


def printData(a):
	print(a.NAtoms())
	print(a.AtomTypes())
	molName="p"
	a.WriteXYZ(molName)
	print("see datasets/"+molName+".xyz file")
	for mol in a.mols:
		print(mol.properties["name"]+" Energy="+str(mol.properties["energy"]))

def xyzstr(mol):
	lines =""
	natoms = mol.atoms.shape[0]
	lines = lines+(str(natoms)+mol.PropertyString()+"\n")
	for i in range (natoms):
		atom_name =  list(atoi.keys())[list(atoi.values()).index(mol.atoms[i])]
		lines = lines+(atom_name+"   "+str(mol.coords[i][0])+ "  "+str(mol.coords[i][1])+ "  "+str(mol.coords[i][2]))
		for k in range(3):
			lines=lines+"  "+str(mol.properties["gradients"][i,k])
		if (i<natoms-1):
			lines = lines+"\n"

	return lines

def printMols(a):
	try:
		for mol in a.mols:
			lines=xyzstr(mol)
			print(lines+"\n")
	except Exception as Ex:
		print("printMols Failed.", Ex)
		raise Ex
		return


def readCChemI(inputFileName):
	try:
		prefixFN=inputFileName.split(".")[0]
		a = MSet(prefixFN,path_="")
		a.mols.append(Mol())
		mol = a.mols[0]
		f=open(inputFileName,"r")
		lines=f.readlines()
		option=int(lines[0])
		line = lines[2].split()
		natoms=int(line[0])
		mol.natoms=natoms
		charge=int(line[1])
		mult=int(line[2])
		mol.coords = np.zeros(shape=(1,1),dtype=np.float64)
		mol.atoms =  np.zeros(1,dtype=np.uint8)
		mol.atoms.resize((natoms))
		mol.coords.resize((natoms,3))
		for i in range(natoms):
			line = lines[i+3].split()
			#  H H H H 0 0.000000 2 0 1.030675 1.867087 2.816794
			mol.atoms[i]=AtomicNumber(line[0])
			for j in range(3):
				try:
					mol.coords[i,j]=float(line[j+8])
				except:
					mol.coords[i,j]=scitodeci(line[j+8])
		mol.properties["totCharge"] = charge
		mol.properties["spin"] = mult
		mol.properties["gradients"] =  np.zeros(shape=(natoms,3),dtype=np.float64)
		mol.properties["name"] = " "
		mol.properties["energy"] = 0.0
		a.mols[0] = mol
		return a,option
	except Exception as Ex:
		print("Read Failed.", Ex)
		raise Ex
		return


def loadData(inputFileName):
	a, option = readCChemI(inputFileName)
	a.Save() 
	a.Load() 
	return a,option

def setData(name,natoms, charge, mult, symbols, X,Y,Z):
	try:
		prefixFN=name
		a = MSet(prefixFN,path_="")
		a.mols.append(Mol())
		mol = a.mols[0]
		mol.natoms=natoms
		mol.coords = np.zeros(shape=(1,1),dtype=np.float64)
		mol.atoms =  np.zeros(1,dtype=np.uint8)
		mol.atoms.resize((natoms))
		mol.coords.resize((natoms,3))
		for i in range(natoms):
			mol.atoms[i]=AtomicNumber(symbols[i])
			mol.coords[i,0]=float(X[i])
			mol.coords[i,1]=float(Y[i])
			mol.coords[i,2]=float(Z[i])
		mol.properties["totCharge"] = charge
		mol.properties["spin"] = mult
		mol.properties["gradients"] =  np.zeros(shape=(natoms,3),dtype=np.float64)
		mol.properties["name"] = " "
		mol.properties["energy"] = 0.0
		a.mols[0] = mol
		a.Save() 
		a.Load() 
		#printMols(a)
		return a
	except Exception as Ex:
		print("setData Failed.", Ex)
		raise Ex
		return

def setCoordinates(a, X,Y,Z):
	try:
		natoms = a.mols[0].atoms.shape[0]
		for i in range(natoms):
			a.mols[0].coords[i,0]=float(X[i])
			a.mols[0].coords[i,1]=float(Y[i])
			a.mols[0].coords[i,2]=float(Z[i])
		#printMols(a)
		return a
	except Exception as Ex:
		print("setCoordinates Failed.", Ex)
		raise Ex
		return


def usage(progName):
	print("==================================================================")
	print("You must provide 2 arguments : inputFileName outputFileName")
	print("Example "+progName+" mol.ici mol.out")
	print("==================================================================")

def getManager(a):
	if (1): # learning energy, gradient and dipole, energy should be in hartree, gradient should be in hartree/angstrom, dipole should be in a.u.
		random.shuffle(a.mols)
		TreatedAtoms = a.AtomTypes()
		PARAMS["NetNameSuffix"] = "training_H2O"
		PARAMS["learning_rate"] = 0.00001
		PARAMS["momentum"] = 0.95
		PARAMS["max_steps"] = 10000  # Train for 5 epochs in total
		PARAMS["TestRatio"] = 0.5
		PARAMS["batch_size"] =  5
		PARAMS["test_freq"] = 2500 # Test for every epoch
		PARAMS["tf_prec"] = "tf.float64" # double precsion
		PARAMS["EnergyScalar"] = 1.0
		PARAMS["GradScalar"] = 1.0/20.0
		PARAMS["DipoleScaler"] = 1.0
		PARAMS["NeuronType"] = "sigmoid_with_param" # choose activation function
		PARAMS["sigmoid_alpha"] = 100.0  # activation params
		PARAMS["HiddenLayers"] = [10, 10]  # number of neurons in each layer
		PARAMS["EECutoff"] = 15.0
		PARAMS["EECutoffOn"] = 0
		PARAMS["Elu_Width"] = 4.6  # when elu is used EECutoffOn should always equal to 0
		PARAMS["EECutoffOff"] = 15.0
		PARAMS["DSFAlpha"] = 0.18
		PARAMS["AddEcc"] = True
		PARAMS["KeepProb"] = [1.0, 1.0, 1.0, 1.0] # each layer's keep probability for dropout
		PARAMS["learning_rate_dipole"] = 0.0001  # learning rate for dipole learning
		PARAMS["learning_rate_energy"] = 0.00001 # learning rate for energy & grads learning
		PARAMS["SwitchEpoch"] = 5  # Train dipole for 2 epochs, then train energy & grads
		d = MolDigester(TreatedAtoms, name_="ANI1_Sym_Direct", OType_="EnergyAndDipole")
		tset = TensorMolData_BP_Direct_EE_WithEle_Release(a, d, order_=1, num_indis_=1, type_="mol",  WithGrad_ = True)
		manager = TFMolManage("Mol_a0_t_0.01_ANI1_Sym_Direct_fc_sqdiff_BP_Direct_EE_SymFunction_training_H2O", tset,False, "fc_sqdiff_BP_Direct_EE_SymFunction", False, False)
		PARAMS['Profiling']=0
		return manager

def getEnergy(a, manager):
	#print("\nJe suis dans getEnegy\n")
	try:
		#printMols(a)
		#print("\nPA "+ str(PARAMS["AN1_r_Rc"])+"\n")
		#print("\nPA "+ str(PARAMS["AN1_a_Rc"])+"\n")
		#print("\nPA "+ str(PARAMS["EECutoffOff"])+"\n")
		Etotal, Ebp, Ebp_atom, Ecc, Evdw, mol_dipole, atom_charge, gradient = manager.EvalBPDirectEELinearSingle(a.mols[0], PARAMS["AN1_r_Rc"], PARAMS["AN1_a_Rc"], PARAMS["EECutoffOff"], True)
		a.mols[0].properties["energy"] = Etotal[0]
		dipole = mol_dipole[0]
		a.mols[0].properties["dipole"] = np.zeros(3,dtype=np.float64)
		for k in range (3):
			a.mols[0].properties["dipole"][k] = dipole[k]/AUPERDEBYE
		natoms = a.mols[0].atoms.shape[0]

		#print('{:} {:} {:} {:} {:} {:} {:}'.format("\nE = ",Etotal[0]," Dipole = ", a.mols[0].properties["dipole"][0], a.mols[0].properties["dipole"][1], a.mols[0].properties["dipole"][2],"\n"))

		charges = atom_charge[0]
		a.mols[0].properties["atomCharges"] = np.zeros(natoms,dtype=np.float64)
		for i in range (natoms):
			a.mols[0].properties["atomCharges"][i] = charges[i] 

		forces = gradient[0]
		a.mols[0].properties["gradients"] =  np.zeros(shape=(natoms,3),dtype=np.float64)
	#	return a
		listRes = []
		listRes.append(Etotal) 
		for k in range (3):
			listRes.append(a.mols[0].properties["dipole"][k])
		for i in range (natoms):
			for k in range(3):
				listRes.append(a.mols[0].properties["gradients"][i,k])
		return listRes
	except Exception as Ex:
		print("getEnergy Failed.", Ex)
		raise Ex
		return

def getEnergyAndForces(a, manager):
	#print("\nJe suis dans getEnergyAndForces\n")
	try:
		#printMols(a)
		#print("\nPA "+ str(PARAMS["AN1_r_Rc"])+"\n")
		#print("\nPA "+ str(PARAMS["AN1_a_Rc"])+"\n")
		#print("\nPA "+ str(PARAMS["EECutoffOff"])+"\n")
		Etotal, Ebp, Ebp_atom, Ecc, Evdw, mol_dipole, atom_charge, gradient = manager.EvalBPDirectEELinearSingle(a.mols[0], PARAMS["AN1_r_Rc"], PARAMS["AN1_a_Rc"], PARAMS["EECutoffOff"], True)
		a.mols[0].properties["energy"] = Etotal[0]
		dipole = mol_dipole[0]
		a.mols[0].properties["dipole"] = np.zeros(3,dtype=np.float64)
		for k in range (3):
			a.mols[0].properties["dipole"][k] = dipole[k]/AUPERDEBYE
		natoms = a.mols[0].atoms.shape[0]

		#print('{:} {:} {:} {:} {:} {:} {:}'.format("\nE = ",Etotal[0]," Dipole = ", a.mols[0].properties["dipole"][0], a.mols[0].properties["dipole"][1], a.mols[0].properties["dipole"][2],"\n"))

		charges = atom_charge[0]
		a.mols[0].properties["atomCharges"] = np.zeros(natoms,dtype=np.float64)
		for i in range (natoms):
			a.mols[0].properties["atomCharges"][i] = charges[i] 

		#KJPERHARTREE = 2625.499638
		#JOULEPERHARTREE = KJPERHARTREE*1000.0
		forces = gradient[0]
		a.mols[0].properties["gradients"] =  np.zeros(shape=(natoms,3),dtype=np.float64)
		for i in range (natoms):
			for k in range(3):
				a.mols[0].properties["gradients"][i,k] = -forces[i,k]/JOULEPERHARTREE/BOHRPERA

		listRes = []
		listRes.append(Etotal) 
		for k in range (3):
			listRes.append(a.mols[0].properties["dipole"][k])
		for i in range (natoms):
			for k in range(3):
				listRes.append(a.mols[0].properties["gradients"][i,k])
		return listRes
	except Exception as Ex:
		print("getEnergyAndForces Failed.", Ex)
		raise Ex
		return

def optGeom(a, manager):
	def GetEnergyForceForMol(m):
		def EnAndForce(x_, DoForce=True):
			tmpm = Mol(m.atoms,x_)
			Etotal, Ebp, Ebp_atom, Ecc, Evdw, mol_dipole, atom_charge, gradient = manager.EvalBPDirectEELinearSingle(tmpm, PARAMS["AN1_r_Rc"], PARAMS["AN1_a_Rc"], PARAMS["EECutoffOff"], True)
			energy = Etotal[0]
			force = gradient[0]
			if DoForce:
				return energy, force
			else:
				return energy
		return EnAndForce


	print("Optimizing====================> ")
	PARAMS["OptMaxCycles"]=500
	PARAMS["OptThresh"] =0.000001
	#Opt = MetaOptimizer(F,m,Box_=False)
	F = GetEnergyForceForMol(a.mols[0])
	Opt = GeomOptimizer(F)
	m = Opt.Opt(a.mols[0])
	
	a.mols[0].coords = m.coords

	Etotal, Ebp, Ebp_atom, Ecc, Evdw, mol_dipole, atom_charge, gradient = manager.EvalBPDirectEELinearSingle(a.mols[0], PARAMS["AN1_r_Rc"], PARAMS["AN1_a_Rc"], PARAMS["EECutoffOff"], True)
	a.mols[0].properties["energy"] = Etotal[0]
	dipole = mol_dipole[0]
	a.mols[0].properties["dipole"] = np.zeros(3,dtype=np.float64)
	for k in range (3):
		a.mols[0].properties["dipole"][k] = dipole[k]/AUPERDEBYE
	natoms = a.mols[0].atoms.shape[0]

	charges = atom_charge[0]
	a.mols[0].properties["atomCharges"] = np.zeros(natoms,dtype=np.float64)
	for i in range (natoms):
			a.mols[0].properties["atomCharges"][i] = charges[i] 

	#KJPERHARTREE = 2625.499638
	#JOULEPERHARTREE = KJPERHARTREE*1000.0

	forces = gradient[0]
	a.mols[0].properties["gradients"] =  np.zeros(shape=(natoms,3),dtype=np.float64)
	for i in range (natoms):
		for k in range(3):
			a.mols[0].properties["gradients"][i,k] = -forces[i,k]/JOULEPERHARTREE/BOHRPERA

	#printMols(a)

	return a

def computeFreqs(a,manager):
	m = a.mols[0]
	def EnAndForce(x_, DoForce=True):
		tmpm = Mol(m.atoms,x_)
		Etotal, Ebp, Ebp_atom, Ecc, Evdw, mol_dipole, atom_charge, gradient = manager.EvalBPDirectEELinearSingle(tmpm, PARAMS["AN1_r_Rc"], PARAMS["AN1_a_Rc"], PARAMS["EECutoffOff"], True)
		energy = Etotal[0]
		force = gradient[0]
		if DoForce:
			return energy, force
		else:
			return energy

	EnergyForceField =  lambda x: EnAndForce(x)

	def ChargeField(x_):
		mtmp = Mol(m.atoms,x_)
		Etotal, Ebp, Ebp_atom, Ecc, Evdw, mol_dipole, atom_charge, gradient = manager.EvalBPDirectEEUpdateSingle(mtmp, PARAMS["AN1_r_Rc"], PARAMS["AN1_a_Rc"], PARAMS["EECutoffOff"], True)
		energy = Etotal[0]
		force = gradient[0]
		return atom_charge[0]

	def EnergyField(x_):
		return EnAndForce(x_,True)[0]

	def DipoleField(x_):
		q = np.asarray(ChargeField(x_))
		dipole = np.zeros(3)
		for i in  range(0, q.shape[0]):
			dipole += q[i]*x_[i]*BOHRPERA
		return dipole

	#printMols(a)

	#Perform Harmonic frequency analysis
	masses = np.array(list(map(lambda x: ATOMICMASSESAMU[x-1],m.atoms)))
	print("masses size=",str(masses.shape[0])+"\n")
	for i in  range(0, masses.shape[0]):
		print("mass "+str(i)+" = "+str(masses[i])+"\n")

#	w,v = HarmonicSpectra(EnergyField, m.coords, m.atoms, WriteNM_=True, Mu_ = DipoleField)
	w,v = HarmonicSpectra(EnergyField, m.coords, m.atoms)

	print("w size=",str(w.shape[0]))
	for i in  range(0, w.shape[0]):
		print(str(w[i])+" "+str(v[i])+"\n")

	return w,v

def saveOutput(a,outputFileName,option):
	try:
		m = a.mols[0]
		f=open(outputFileName,"w")
		f.write( '{: 0.14f}'.format(m.properties["energy"])+"\n")
		f.write('{: 0.14f} {: 0.14f} {: 0.14f}\n'.format(m.properties["dipole"][0],m.properties["dipole"][1],m.properties["dipole"][2]))
		natoms = m.atoms.shape[0]
		if(option==1):
			for i in range (natoms):
				f.write('{: 0.14f} {: 0.14f} {: 0.14f}'.format(-m.properties["gradients"][i,0],-m.properties["gradients"][i,1],-m.properties["gradients"][i,2]))
				f.write("\n")
			return
		if(option==2):
			f.write("Geometry\n")
			f.write(str(natoms)+" "+str(m.properties["totCharge"])+" "+str(m.properties["spin"])+"\n")
			for i in range (natoms):
				atom_name =  list(atoi.keys())[list(atoi.values()).index(m.atoms[i])]
				f.write('{:5} {:5} {:5} {:5} {:} {:0.14f} {:} {: 0.14f} {: 0.14f} {: 0.14f} {:}'.format(atom_name, atom_name, atom_name,atom_name, " 0 ", m.properties["atomCharges"][i]," 2 1 ", m.coords[i][0],m.coords[i][1],m.coords[i][2]," 0 "))
				f.write('{:10} {: 0.14f} {: 0.14f} {: 0.14f}'.format('GRADIENT',m.properties["gradients"][i,0],m.properties["gradients"][i,1],m.properties["gradients"][i,2]))
				f.write("\n")
			return

	except Exception as Ex:
		print("Write Failed.", Ex)
		raise Ex
		return


def main(argv):
	if len(argv)<3:
		usage(argv[0])
		sys.exit(2) 

	a,option = loadData(argv[1])
	manager = getManager(a)
	if(option==0):
		a = getEnergy(a,manager)

	if(option==1):
		a = getEnergyAndForces(a,manager)

	if(option==2):
		a = optGeom(a,manager)
#		printMols(a)

	if(option==3):
		computeFreqs(a,manager)
	if(option != 3):
		saveOutput(a,argv[2],option)

if __name__ == "__main__":
    main(sys.argv[0:])
