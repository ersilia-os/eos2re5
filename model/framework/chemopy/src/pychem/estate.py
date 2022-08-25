# -*- coding: utf-8 -*-
"""
##############################################################################
This module is to compute the estate fingerprints and values based on Kier 

and Hall's paper. If you have any question please contact me via email.

My email adress is : orientalcds@gmail.com

Created on Tue May 24 14:32:52 2011

@author: Dongsheng Cao
##############################################################################
"""

from rdkit.Chem.EState import Fingerprinter  as ESFP
from rdkit import Chem

import AtomTypes as ATEstate
import numpy

Version=1.0
################################################################

def _CalculateEState(mol,skipH=1):
    """
    #################################################################
    **Internal used only**
    
    Get the EState value of each atom in a molecule
    #################################################################
    """
    mol=Chem.AddHs(mol)
    if skipH==1: 
        mol=Chem.RemoveHs(mol)
    tb1=Chem.GetPeriodicTable()
    nAtoms=mol.GetNumAtoms()
    Is=numpy.zeros(nAtoms,numpy.float)
    for i in range(nAtoms):
        at=mol.GetAtomWithIdx(i)
        atNum=at.GetAtomicNum()
        d=at.GetDegree()
        if d>0:
            h=at.GetTotalNumHs()
            dv=tb1.GetNOuterElecs(atNum)-h         
            #dv=numpy.array(_AtomHKDeltas(at),'d')
            N=_GetPrincipleQuantumNumber(atNum)
            Is[i]=(4.0/(N*N)*dv+1)/d
    dists = Chem.GetDistanceMatrix(mol,useBO=0,useAtomWts=0)
    dists +=1
    accum = numpy.zeros(nAtoms,numpy.float)
    for i in range(nAtoms):
        for j in range(i+1,nAtoms):
            p=dists[i,j]
            if p < 1e6:
                temp=(Is[i]-Is[j])/(p*p)
                accum[i] +=temp
                accum[j] -=temp
    res=accum+Is
    return res


def _GetPrincipleQuantumNumber(atNum):
    """
    #################################################################
    *Internal Use Only*
    
    Get the principle quantum number of atom with atomic
    
    number equal to atNum 
    #################################################################
    """
    if atNum<=2:
        return 1
    elif atNum<=10:
        return 2
    elif atNum<=18:
        return 3
    elif atNum<=36:
        return 4
    elif atNum<=54:
        return 5
    elif atNum<=86:
        return 6
    else:
        return 7



def CalculateHeavyAtomEState(mol):
    
    """
    #################################################################
    The sum of the EState indices over all non-hydrogen atoms
    
    -->Shev
    #################################################################
    """
    
    return round(sum(_CalculateEState(mol)),3)

def _CalculateAtomEState(mol,AtomicNum=6):
    """
    #################################################################
    **Internal used only**
    
    The sum of the EState indices over all atoms 
    #################################################################
    """
    nAtoms=mol.GetNumAtoms()
    Is=numpy.zeros(nAtoms,numpy.float)
    Estate=_CalculateEState(mol)
    
    for i in range(nAtoms):
        at=mol.GetAtomWithIdx(i)
        atNum=at.GetAtomicNum()
        if atNum==AtomicNum:
            Is[i]=Estate[i]
            
    res=sum(Is)
    
    return res
    
def CalculateCAtomEState(mol):
    """
    #################################################################
    The sum of the EState indices over all C atoms
    
    -->Scar
    #################################################################
    """
    return _CalculateAtomEState(mol,AtomicNum=6)



def CalculateHalogenEState(mol):
    
    """
    #################################################################
    The sum of the EState indices over all Halogen atoms
    
    -->Shal
    #################################################################
    """
    
    Nf=_CalculateAtomEState(mol,AtomicNum=9)
    Ncl=_CalculateAtomEState(mol,AtomicNum=17)
    Nbr=_CalculateAtomEState(mol,AtomicNum=35)
    Ni=_CalculateAtomEState(mol,AtomicNum=53)
  
    return round(Nf+Ncl+Nbr+Ni,3)



def CalculateHeteroEState(mol):
    """
    #################################################################
    The sum of the EState indices over all hetero atoms
    
    -->Shet
    #################################################################
    """
    
    Ntotal=sum(_CalculateEState(mol))
    NC=_CalculateAtomEState(mol,AtomicNum=6)
    NH=_CalculateAtomEState(mol,AtomicNum=1)
    
    return round(Ntotal-NC-NH,3)


def CalculateAverageEState(mol):
    """
    #################################################################
    The sum of the EState indices over all non-hydrogen atoms 
    
    divided by the number of non-hydrogen atoms.
    
    -->Save
    #################################################################
    """
    N=mol.GetNumAtoms()
    return round(sum(_CalculateEState(mol))/N,3)

def CalculateMaxEState(mol):

    """
    #################################################################
    Obtain the maximal Estate value in all atoms
    
    -->Smax
    #################################################################
    """
    return round(max(_CalculateEState(mol)),3)


def CalculateMinEState(mol):
    """
    #################################################################
    Obtain the minimal Estate value in all atoms
    
    -->Smin
    #################################################################
    """
    
    return round(min(_CalculateEState(mol)),3)


def CalculateDiffMaxMinEState(mol):
    """
    #################################################################
    The difference between Smax and Smin
    
    -->DS
    #################################################################
    """
    return round(max(_CalculateEState(mol))-min(_CalculateEState(mol)),3)
    

def CalculateAtomTypeEState(mol):
    """
    #################################################################
    Calculation of sum of E-State value of specified atom type
    
    res---->list type
    #################################################################
    """
    AT=ATEstate.GetAtomLabel(mol)
    Estate=_CalculateEState(mol)
    res=[]
    for i in AT:
        if i==[]:
            res.append(0)
        else:
            res.append(sum([Estate[k] for k in i]))
    ESresult={}
    for n,es in enumerate(res):
        ESresult['S'+str(n+1)]=round(es,3)
        
    return ESresult



def CalculateEstateFingerprint(mol):
    """
    #################################################################
    The Calculation of EState Fingerprints.
    
    It is the number of times each possible atom type is hit.
    
    Usage:
        
        result=CalculateEstateFingerprint(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a dict form containing 79 estate fragments.
    #################################################################
    """
    temp=ESFP.FingerprintMol(mol)
    res={}
    for i,j in enumerate(temp[0]):
        res['Sfinger'+str(i+1)]=j
    
    return res


def CalculateEstateValue(mol):
    """
    #################################################################
    The Calculate of EState Values.
    
    It is the sum of the Estate indices for atoms of each type.
    
    Usage:
        
        result=CalculateEstateValue(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a dict form containing 79 estate values.
    #################################################################
    """
    temp=ESFP.FingerprintMol(mol)
    res={}
    for i,j in enumerate(temp[1]):
        res['S'+str(i+1)]=round(j,3)
    
    return res
        

def CalculateMaxAtomTypeEState(mol):
    """
    #################################################################
    Calculation of maximum of E-State value of specified atom type
    
    res---->dict type
    
    Usage:
        
        result=CalculateMaxAtomTypeEState(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a dict form containing 79 max estate values.
    #################################################################
    """
    AT=ATEstate.GetAtomLabel(mol)
    Estate=_CalculateEState(mol)
    res=[]
    for i in AT:
        if i==[]:
            res.append(0)
        else:
            res.append(max([Estate[k] for k in i]))
    ESresult={}
    for n,es in enumerate(res):
        ESresult['Smax'+str(n)]=round(es,3)
        
    return ESresult



def CalculateMinAtomTypeEState(mol):
    """
    #################################################################
    Calculation of minimum of E-State value of specified atom type
    
    res---->dict type
    
    Usage:
        
        result=CalculateMinAtomTypeEState(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a dict form containing 79 min estate values.
    #################################################################
    """
    AT=ATEstate.GetAtomLabel(mol)
    Estate=_CalculateEState(mol)
    res=[]
    for i in AT:
        if i==[]:
            res.append(0)
        else:
            res.append(min([Estate[k] for k in i]))
    ESresult={}
    for n,es in enumerate(res):
        ESresult['Smin'+str(n)]=round(es,3)
        
    return ESresult



def GetEstate(mol):
    """
    #################################################################
    Obtain all descriptors related to Estate.

    Usage:
        
        result=GetEstate(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a dict form containing all estate values.
    #################################################################
    """
    result={}
    #result.update(CalculateEstateFingerprint(mol))
    result.update(CalculateEstateValue(mol))
    result.update(CalculateMaxAtomTypeEState(mol))
    result.update(CalculateMinAtomTypeEState(mol))
    result.update({'Shev':CalculateHeavyAtomEState(mol)})
    result.update({'Scar':CalculateCAtomEState(mol)})
    result.update({'Shal':CalculateHalogenEState(mol)})
    result.update({'Shet':CalculateHeteroEState(mol)})
    result.update({'Save':CalculateAverageEState(mol)})
    result.update({'Smax':CalculateMaxEState(mol)})
    result.update({'Smin':CalculateMinEState(mol)})
    result.update({'DS':CalculateDiffMaxMinEState(mol)})

    
    return result


def _GetHTMLDoc():
    """
    #################################################################
    Write HTML documentation for this module.
    #################################################################
    """
    import pydoc
    pydoc.writedoc('estate')
################################################################
  
if __name__=='__main__':
    

    smi5=['COCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCOCCN','c1ccccc1N']
    smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]']
    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        print index+1
        print smi      
##        print '\t',CalculateEstateFingerprint(m)
##        print '\t',CalculateEstateValue(m)
##        print '\t',CalculateMaxAtomTypeEState(m)
##        print '\t', CalculateMinAtomTypeEState(m)
        
        print GetEstate(m)
        print len(GetEstate(m))

    
