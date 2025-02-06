#pragma once
class CZInChI
{
public' :
	CZInChI(void)
	~CZInChI(void)
	//Input
	CString strInPutInChI


	//Output
	CString strClassL1
	CString strClassL2
	CString strClassL3

	int iClassL1
	int iClassL2
	int IClassL3

	CString strError

	void DoClassification(CString strInput)				//////////////////////////////////
private' :
	//Functions

	//Main Fucnations

	
	void Clearinfo(void)										//////////////////////////////		1										
	void Cutpart(void)											//////////////////////////////		2
	void Countatom(int iSub)									//////////////////////////////		3
	int Findatom(CString str,int iSub)								//////////////////////////		3-1
	int IsOrgSub(int iSub)										//////////////////////////////		4
	void DoL1_1(void)											//////////////////////////////		5
	void DoOrganic(void)											//////////////////////////		5-1
	void Atomclass(void)
	void Oclass1(void)
	void OclassO(void)
	void OclassN(void)
	void OclassH(void)
	void OclassP(void)
	void OclassSi(void)
	void OclassS(void)
	
	void OclassNO(void)
	void OclassOS(void)
	void Subclass(void)

	void DoL1_2(void)											//////////////////////////////		6
	void DoL1_3(void)											//////////////////////////////		7
	void DoL1_4(void)											//////////////////////////////		8
	void DoL1_5(void)											//////////////////////////////		9
	void DoAnion(void)											


	//Sub Functions
	void RingCheck(void)
	void Ringsort(void)
	void Ringdelete(void)
	void Unsatbond(int ithAtom)
	void FIndAtomloc(int ithAtom)
	double BondCheck(int ithAtom)
	double Checkhydrogen(int ithAtom)
	int Countstr(CString string,CString str2)


	CString strAtomclass
	CString strMainclass
	CString strSubclass
	
	int iAtomclass
	int iMainclass
	int iSubclass
	
	//Variables
	CString strTInChI																//InChI string
	CString strTFormula																// Formula part in InChI
	CString strTAtom																// Atom connection part in InChI
	CString strTHydrogen															// Hydrogen atom part in InChI
	CString strProton																// Proton part in InChI
	CString strTCharge																// Chaged part in InChI
	
	CString strInChI																//InChI string
	CString strFormula																// Formula part in InChI
	CString strAtom																// Atom connection part in InChI
	CString strHydrogen															// Hydrogen atom part in InChI
	CString strCharge																// Chaged part in InChI
	
	CString strBond																// Bonding part In InChI
	CString strTetra																// Tetra stereo part in InChI
	CString strStereo																// Stereochemistry part in InChI
	CString strIsotrope															// Isotrope pate in in InChI
	CString strfeature																// Feature of chemical compound.

	
	int iNumSubStruc
	int ithSub
	
	CString strSubInChI[20]
	CString strSubFormula[20]
	CString strSubConnect[20]
	CString strSubHydrogen[20]
	CString strSubCharge[20]
	
	//Numbers of important atoms
	int iSubMetal[20]
	int iNumSubCarbon[20]
	int iNumSubNitrogen[20]
	int iNumSubOxygen[20]
	int iNumSubSilicon[20]
	int iNumSubPhosphine[20]
	int iNumSubSulfur[20]
	int iNumSubHydrogen[20]
	int iNumSubHalogen[20]																// Number of Halogen
	int iNumSubFluorine[20]																// Number of Fluorine
	int iNumSubChlorine[20]																// Number of Chlorine
	int iNumSubBromine[20]																// Number of Bromine
	int iNumSubIodine[20]																	// Number of Iodine	
	int iSubCountAtom[20]
	int iSubAtom[20][255]
	int iOrgSub[20]
	
	//variables for substructure

	
	int iKindKGuest

	int iAtom[255]
	int iNumRing																	// Number of rings
	int iNumUnsatbond																// Number of unsaturated bonds
	int iNumTriple		
	int iNumDouble
	
	int iNumBranch																	// Number of branches
	int iNumAromatic																// Number of aromatics

	int ring[255]
	int ringatom[255][255]
	int iringNUM[255]

	int iloc[20]
	int iBond[20]
	int iUnsatbond[255]
	int iNumbond[20]

	int iOrCheck[10]
	int FGs[255]
	int iMulSub[20]
	int nstdBond
	int iSubAnion

	int iGuestAtom
	int iCountAtom																	// Total atom number in molecule
	int iKindAtom																	// Total kind of atom in molecule
	int iPartNum	

	int iNumCarbon																	// Number of Carbon
	int iNumOxygen																	// Number of Oxygen
	int iNumNitrogen																// Number of Nitrogen
	int iNumSulfur																	// Number of Sulfur
	int iNumSilicon																// Number of Silicon
	int iNumPhosphine																// Number of Phosphorus
	int iNumHalogen																// Number of Halogen
	int iNumFluorine																// Number of Fluorine
	int iNumChlorine																// Number of Chlorine
	int iNumBromine																// Number of Bromine
	int iNumIodine																	// Number of Iodine
	int iNumHydrogen																// Number of Hydrogen
}

