#!/usr/local/bin/python3.6
from __future__ import absolute_import
from __future__ import print_function
from TensorMol import *
import sys

os.environ["CUDA_VISIBLE_DEVICES"]="" # set to use CPU


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

def closeSession():
	try:
		sess = tf.compat.v1.get_default_session()
		sess.close()
	except Exception as Ex:
		print("closeSession Failled.", Ex)
		raise Ex
		return



def getManager(a):
	TreatedAtoms = np.array([1,6,7,8], dtype=np.uint8)
	PARAMS["tf_prec"] = "tf.float64"
	PARAMS["NeuronType"] = "sigmoid_with_param"
	PARAMS["sigmoid_alpha"] = 100.0
	PARAMS["HiddenLayers"] = [2000, 2000, 2000]
	PARAMS["EECutoff"] = 15.0
	PARAMS["EECutoffOn"] = 0
	PARAMS["Elu_Width"] = 4.6  # when elu is used EECutoffOn should always equal to 0
	PARAMS["EECutoffOff"] = 15.0
	PARAMS["AddEcc"] = True
	PARAMS["KeepProb"] = [1.0, 1.0, 1.0, 0.7]
	d = MolDigester(TreatedAtoms, name_="ANI1_Sym_Direct", OType_="EnergyAndDipole")  # Initialize a digester that apply descriptor for the fragme
	tset = TensorMolData_BP_Direct_EE_WithEle(a, d, order_=1, num_indis_=1, type_="mol",  WithGrad_ = True)
	Solvation_ = False
	if Solvation_:
		PARAMS["DSFAlpha"] = 0.18
		manager=TFMolManage("chemspider12_solvation", tset,False,"fc_sqdiff_BP_Direct_EE_ChargeEncode_Update_vdw_DSF_elu_Normalize_Dropout",False,False)
	else:
		PARAMS["DSFAlpha"] = 0.18*BOHRPERA
		manager=TFMolManage("chemspider12_nosolvation", tset,False,"fc_sqdiff_BP_Direct_EE_ChargeEncode_Update_vdw_DSF_elu_Normalize_Dropout",False,False)
	return manager


def getEnergy(a, manager):
	#print("\nJe suis dans getEnegy\n")
	try:
		#printMols(a)
		#print("\nPA "+ str(PARAMS["AN1_r_Rc"])+"\n")
		#print("\nPA "+ str(PARAMS["AN1_a_Rc"])+"\n")
		#print("\nPA "+ str(PARAMS["EECutoffOff"])+"\n")
		Etotal, Ebp, Ebp_atom, Ecc, Evdw, mol_dipole, atom_charge, gradient = manager.EvalBPDirectEEUpdateSingle(a.mols[0], PARAMS["AN1_r_Rc"], PARAMS["AN1_a_Rc"], PARAMS["EECutoffOff"], True)
		a.mols[0].properties["energy"] = Etotal[0]
		dipole = mol_dipole[0]
		a.mols[0].properties["dipole"] = np.zeros(3,dtype=np.float64)
		for k in range (3):
			a.mols[0].properties["dipole"][k] = dipole[k]
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
#		for i in range (natoms):
#			for k in range(3):
#				listRes.append(a.mols[0].properties["gradients"][i,k])
		for i in range (natoms):
				listRes.append(a.mols[0].properties["atomCharges"][i])
	# data must be returned  as an one  list in this order : Energy, Dipole, charges, all in atomic unit
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
		Etotal, Ebp, Ebp_atom, Ecc, Evdw, mol_dipole, atom_charge, gradient = manager.EvalBPDirectEEUpdateSingle(a.mols[0], PARAMS["AN1_r_Rc"], PARAMS["AN1_a_Rc"], PARAMS["EECutoffOff"], True)
		a.mols[0].properties["energy"] = Etotal[0]
		dipole = mol_dipole[0]
		a.mols[0].properties["dipole"] = np.zeros(3,dtype=np.float64)
		for k in range (3):
			a.mols[0].properties["dipole"][k] = dipole[k]
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

		for i in range (natoms):
				listRes.append(a.mols[0].properties["atomCharges"][i])
	# data must be returned  as an one  list in this order : Energy, Dipole, gradients, charges, all in atomic unit
		return listRes
	except Exception as Ex:
		print("getEnergyAndForces Failed.", Ex)
		raise Ex
		return

def optGeom(a, manager):
	def GetEnergyForceForMol(m):
		def EnAndForce(x_, DoForce=True):
			tmpm = Mol(m.atoms,x_)
			Etotal, Ebp, Ebp_atom, Ecc, Evdw, mol_dipole, atom_charge, gradient = manager.EvalBPDirectEEUpdateSingle(tmpm, PARAMS["AN1_r_Rc"], PARAMS["AN1_a_Rc"], PARAMS["EECutoffOff"], True)
			energy = Etotal[0]
			force = gradient[0]
			if DoForce:
				return energy, force
			else:
				return energy
		return EnAndForce


	# print("Optimizing====================> ")
	PARAMS["OptMaxCycles"]=500
	PARAMS["OptThresh"] =0.000001
	#Opt = MetaOptimizer(F,m,Box_=False)
	F = GetEnergyForceForMol(a.mols[0])
	Opt = GeomOptimizer(F)
	m = Opt.Opt(a.mols[0])
	
	a.mols[0].coords = m.coords

	Etotal, Ebp, Ebp_atom, Ecc, Evdw, mol_dipole, atom_charge, gradient = manager.EvalBPDirectEEUpdateSingle(a.mols[0], PARAMS["AN1_r_Rc"], PARAMS["AN1_a_Rc"], PARAMS["EECutoffOff"], True)
	a.mols[0].properties["energy"] = Etotal[0]
	dipole = mol_dipole[0]
	a.mols[0].properties["dipole"] = np.zeros(3,dtype=np.float64)
	for k in range (3):
		a.mols[0].properties["dipole"][k] = dipole[k]
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

	listRes = []
	listRes.append(Etotal) 
	for k in range (3):
		listRes.append(a.mols[0].properties["dipole"][k])
	for i in range (natoms):
		for k in range(3):
			listRes.append(a.mols[0].coords[i,k])
	for i in range (natoms):
		for k in range(3):
			listRes.append(a.mols[0].properties["gradients"][i,k])
		for i in range (natoms):
				listRes.append(a.mols[0].properties["atomCharges"][i])
	# data must be returned  as an one  list in this order : Energy, Dipole, gradients, charges, all in atomic unit, exept coordinates in Ang
	return listRes

