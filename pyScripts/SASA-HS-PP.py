# SASA based hotspot prediction for Protein - Nucleic acid script
# SASA-HS-PNA
# Using VMD and Weka commands
# Authors:
# Cristian R Munteanu, Carlos Fernandez-Lozano (RNASA, Computer Science Faculty, University of A Coruna)
# Irina Moreira (University of Porto)

# import modules
import os,sys
import time
from Bio.PDB import *
from Bio.PDB import PDBParser
from Bio.PDB.PDBParser import PDBConstructionException

# Write PDB ATOM file from original PDB using chain order from Monomers
def WritePDB_ATOM(PDB,Mon1,Mon2):
    # read original PDB
    f = open(PDB,"r")
    li=f.readlines()
    f.close()

    PDB_ATOM=PDB[:-4]+"_ATOM.pdb"
    fo = open(PDB_ATOM, "w+")  # open PDB ATOM to be completed
    for i in range(len(Mon1)): # for each chain in Monomer 1
        for line in li:
            # if line starts with ATOM and chain is from Monomer1
            if (line[:4]=="ATOM" and (line.split())[4]==Mon1[i]): 
                fo.write(line)

    for i in range(len(Mon2)): # for each chain in Monomer 2
        for line in li:
            # if line starts with ATOM and chain is from Monomer1
            if (line[:4]=="ATOM" and (line.split())[4]==Mon2[i]): 
                fo.write(line)

    fo.close()
    return

# Write script_AddH for specific chains
def AddH(sFile,PDB,Monomer1,Monomer2):
    f = open(sFile,'w')
    f.write("mol new "+PDB+".pdb\nmol new "+PDB+".pdb\n")
    f.write(
"""set protein [atomselect top all]
set chains [lsort -unique [$protein get chain]]
foreach chain $chains {
    set sel [atomselect top "protein and chain $chain"]
""")
    f.write("$sel writepdb "+PDB+"_${chain}.pdb\n}\n")
    f.write(
"""

package require psfgen
topology top_na.inp
alias residue HIS HSD
alias residue HOH TIP3
alias residue ZN ZN2
alias atom ILE CD1 CD
alias atom HOH O OH2
pdbalias residue DG GUA
pdbalias residue DC CYT
pdbalias residue DA ADE
pdbalias residue DT THY

foreach bp { GUA CYT ADE THY URA } {
     pdbalias atom $bp "O5\*" O5'
     pdbalias atom $bp "C5\*" C5'
     pdbalias atom $bp "O4\*" O4'
     pdbalias atom $bp "C4\*" C4'
     pdbalias atom $bp "C3\*" C3'
     pdbalias atom $bp "O3\*" O3'
     pdbalias atom $bp "C2\*" C2'
     pdbalias atom $bp "O2\*" O2'
     pdbalias atom $bp "C1\*" C1'
        } 
                                                                                                                                                    
""")

    for i in range(len(Monomer1)):
        f.write("segment "+Monomer1[i]+" {\n     pdb "+PDB+"_"+Monomer1[i]+".pdb\n     first none\n     last none\n}\n")

    for i in range(len(Monomer2)):
        f.write("segment "+Monomer2[i]+" {\n     pdb "+PDB+"_"+Monomer2[i]+".pdb\n     first none\n     last none\n}\n")
    f.write("\n\n")


    for i in range(len(Monomer1)):
        f.write("coordpdb "+PDB+"_"+Monomer1[i]+".pdb "+Monomer1[i]+"\n")

    for i in range(len(Monomer2)):
        f.write("coordpdb "+PDB+"_"+Monomer2[i]+".pdb "+Monomer2[i]+"\n")
    f.write("\n\n")

    f.write("guesscoord\nwritepdb "+PDB+"_Hs.pdb\nquit\nexit")

    f.close()
    return

# --------------------------------------------------------------------
# Write script_create_monomers for specific chains
def Create_monomers(sFile,PDB,Monomer1,Monomer2):
    f = open(sFile,'w')
    f.write("mol new "+PDB+"_Hs.pdb\n")
    print "sFile= ",sFile 
    f.write("set monomer1 [atomselect top \""+"chain "+Monomer1[0])
    for i in range(1,len(Monomer1)):
        f.write(" or chain "+Monomer1[i])
    f.write("\"]\n")
    f.write("$monomer1 writepdb "+PDB[0:4]+"_monomer1.pdb\n")

    f.write("set monomer2 [atomselect top \""+"chain "+Monomer2[0])
    for i in range(1,len(Monomer2)):
        f.write(" or chain "+Monomer2[i])
    f.write("\"]\n")
    f.write("$monomer2 writepdb "+PDB[0:4]+"_monomer2.pdb\n")

    f.write("quit\nexit")

    f.close()
    return

# --------------------------------------------------------------------
# Write script_runSASA
def RunSASA(sFile,PDB):
    f = open(sFile,'w')
    f.write("mol new "+PDB+"_Hs.pdb\nsource _sasa_res.tcl\ngetAllResSASA all 1.4\nquit\nexit")
    return
# --------------------------------------------------------------------
# Write script_runSASA_mon1
def RunSASA_mon1(sFile,PDB):
    f = open(sFile,'w')
    f.write("mol new "+PDB[0:4]+"_monomer1.pdb\nsource _sasa_res.tcl\ngetAllResSASA all 1.4\nquit\nexit")
    return
# --------------------------------------------------------------------
# Write script_runSASA_mon2
def RunSASA_mon2(sFile,PDB):
    f = open(sFile,'w')
    f.write("mol new "+PDB[0:4]+"_monomer2.pdb\nsource _sasa_res.tcl\ngetAllResSASA all 1.4\nquit\nexit")
    return
#-----------------------------------------------------------------------

def unique(s):
    n = len(s)
    if n == 0:
        return []
    u = {}
    try:
        for x in s:
            u[x] = 1
    except TypeError:
        del u  # move on to the next method
    else:
        return u.keys()
    try:
        t = list(s)
        t.sort()
    except TypeError:
        del t  # move on to the next method
    else:
        assert n > 0
        last = t[0]
        lasti = i = 1
        while i < n:
            if t[i] != last:
                t[lasti] = last = t[i]
                lasti += 1
            i += 1
        return t[:lasti]
    # Brute force is all that's left.
    u = []
    for x in s:
        if x not in u:
            u.append(x)
    return u

# --------------------------------------------------------------------
# Write script to generate interface residue file (for monomer 1)
def GetInterfaceResidues1(scriptF,outF,sPDB,Monomer1,Monomer2):
    f = open(scriptF,'w')
    f.write("mol new "+sPDB+"\n")
    f.write("set outfile [open \""+outF+"\" w]\n")
    f.write("set sel1 [atomselect top \"protein and (chain "+Monomer1[0])
    for i in range(1,len(Monomer1)):
        f.write(" or chain "+Monomer1[i])
    f.write(") and within 5 of (chain "+Monomer2[0])
    for i in range(1,len(Monomer2)):
        f.write(" or chain "+Monomer2[i])
    f.write(")\"]\n")
    f.write("$sel1 get {resid resname}\nset sel2 [$sel1 get {resid resname}]\nputs $outfile \"$sel2\"\nclose $outfile\nquit\nexit")

    f.close()
    return

# --------------------------------------------------------------------
# Write script to generate interface residue file (for monomer 2)
def GetInterfaceResidues2(scriptF,outF,sPDB,Monomer1,Monomer2):
    f = open(scriptF,'w')
    f.write("mol new "+sPDB+"\n")
    f.write("set outfile [open \""+outF+"\" w]\n")
    f.write("set sel1 [atomselect top \"protein and (chain "+Monomer2[0])
    for i in range(1,len(Monomer2)):
        f.write(" or chain "+Monomer2[i])
    f.write(") and within 5 of (chain "+Monomer1[0])
    for i in range(1,len(Monomer1)):
        f.write(" or chain "+Monomer1[i])
    f.write(")\"]\n")
    f.write("$sel1 get {resid resname}\nset sel2 [$sel1 get {resid resname}]\nputs $outfile \"$sel2\"\nclose $outfile\nquit\nexit")

    f.close()
    return

# --------------------------------------------------------------------
# Get list of residues from residue files (VMD output)
def GetResidueList(sResidueFile):
    lResidues=unique((open(sResidueFile).read().replace('\n','').replace('} {','#').replace('{','')).split("#"))
    return filter(None, lResidues) # remove null strings

# --------------------------------------------------------------------
# Filter the results.csv => new file
def FilterResults(ResultsFile,ResidueFile1,ResidueFile2,NewResults):
    # get list of residues for P-NA
    lResidues1=GetResidueList(ResidueFile1)
    lResidues2=GetResidueList(ResidueFile2)
    lResidues=unique(lResidues1+lResidues2)

    # read results
    f=open(ResultsFile,"r")
    o=open(NewResults,"w")
    #print header
    o.write("PDBChain,PDBResNo,PDBResName,Res No,Res Name,SASA comp,SASA mon,SASA del,SASA rel,SASA compres,SASA mon/res,SASA del/res,SASA rel/res,SASA comp/ave,SASA mon/ave,SASA del/ave,SASA rel/ave\n")
    lines=f.readlines()
    for i in range(2,len(lines)):
        CurrI=lines[i]
        CurrLine=CurrI.split(",")
        for j in range(len(lResidues)):
            CurrJ=lResidues[j]
            CurrAA=CurrJ.split(" ") # CurrAA[0] = ResNo, CurrAA[1]=ResName
            if (CurrLine[1]==CurrAA[0]) and (CurrLine[2]==CurrAA[1]):
                o.write(CurrI)
    o.close()
    f.close()
    return

# --------------------------------------------------------------------------
# Get PDB residue list fpr protein and nucleic acid residues only as triles
def PDBGetResidueList(PDB,Monomers):
    p=PDBParser(PERMISSIVE=1)
    s=p.get_structure("test",PDB)
    model=s[0]
    AAlist=[] # list with residues for the entire PDB

    for monChain in Monomers:
        for chain in model.get_list():
            for residue in chain.get_iterator():
                if chain.get_id()==monChain:
                    residue_id = residue.get_id()
                    hetfield = residue_id[0]
                    if (hetfield[0]!="H" and hetfield[0]!="W"): # do not print water and hetero residues
                        AAlist.append((chain.get_id(), residue.get_id()[1], residue.get_resname()))
    return AAlist # list of triples (chain, ResNo, ResName)

# ---------------------------------------------------------------------------------------
# Correct the residue numbering for result.csv from PDB for a list of chains (monomers) 
def CorrectResults(sFile,sNewFile,PDB,Monomers):
    PDBResidues = PDBGetResidueList(PDB,Monomers) #print all residues from monomers
    # Open first file to be modified
    f = open(sFile,"r")
    lines=f.readlines()
    print lines
    f.close()

    # Open the new file
    fn = open(sNewFile,"w")
    
    # write first 2 lines
    fn.write("PDBChain,PDBResNo,PDBResName,"+lines[0])
    fn.write(",,,"+lines[1])

    # write the new data
    i=0
    for r in range(2,len(lines)):
        if i<len(PDBResidues):
            (Chain,ResNo,ResName)=PDBResidues[i]
            fn.write(Chain+","+str(ResNo)+","+ResName+","+lines[r])
        i=i+1

    fn.close()
    return

# ---------------------------------------------------------------------------------------
def GetConsurfScores(sFileConsurf):
    # get a list of ResName, ResNo, Chain, Score from consurf files with multiple chain data
    # searching for data with columns: POS,SEQ,3LATOM,SCORE,COLOR,CONFIDENCE,INTERVAL,CONFIDENCE,INTERVAL,COLORS,MSA,DATA,RESIDUE,VARIETY
    ScoreList = [] # final results
    for line in open(sFileConsurf):
        cline = line.split() 
        if (len(cline) == 9 or len(cline) == 10):
            if cline[0] != '-' and cline[0] != 'POS'and cline[2]!='-':
                ResName_No, Chain = cline[2].split(':')
                ScoreList.append([ResName_No[0:3], ResName_No[3:len(ResName_No)], Chain, cline[3]]) # [ResName, ResNo, Chain, Score]  
    return ScoreList # list with ResName, ResNo, Chain, Score

# ---------------------------------------------------------------------------------------
def AddConsurf(Chains,ConsurfList,sResultFilter,sFinalRes):
    # add consurf data to previous filtered results

    # read Filtered results
    inF = open(sResultFilter,"r")
    lines=inF.readlines()
    inF.close()

    # write final results
    outF = open(sFinalRes,'w')
    outF.write((lines[0]).replace('\n', '').replace('\r','')+",Score\n") # write header

    for chain in Chains: # for each chain, check each restult filtered line to match with Consurfs
        for i in range(len(lines)): # for each line in results
            # PDBChain,PDBResNo,PDBResName,Res No,Res Name
            # SASA comp,SASA mon,SASA del,SASA rel,SASA compres,SASA mon/res,SASA del/res,SASA rel/res
            # SASA comp/ave,SASA mon/ave,SASA del/ave,SASA rel/ave

            CurrLine = lines[i]
            cline = (CurrLine.split(','))[0:17] # get only first 17 columns!

            # writing 17 cols as a string
            clineStr = ""
            for item in cline:
                clineStr += str(item)+","
                
            for consurf in ConsurfList: # for each consurf
                # ResName, ResNo, Chain, Score
                if (consurf[2] == cline[0] and chain == cline[0]) and (consurf[0] == cline[2]) and (consurf[1] == cline[1]) : # match chains, residue, res pos
                    outF.write(clineStr+str(consurf[3])+"\n") # write header
                


    outF.close()
    return

# ---------------------------------------------------------------------------------------
def FinalRes(sFile1,sFile2,sFinal):
    # generate the final result file from Weka output and SASA results
    
    f1 = open(sFile1,"r")
    f1lines = f1.readlines()
    f1.close()

    f2 = open(sFile2,"r")
    f2lines = f2.readlines()
    f2.close()

    fres = open(sFinal,"w")
    fres.write("inst#,actual,predicted,error prediction,PDBChain,PDBResNo,PDBResName,SASA comp,SASA mon,SASA del,SASA rel,SASA compres,SASA mon/res,SASA del/res,SASA rel/res,SASA comp/ave,SASA mon/ave,SASA del/ave,SASA rel/ave,Score\n")


    l=1
    for line in f1lines:
        if len(line)>35 and line[1]!='i':
            CurrLine =line[:-1]
            Res = (f2lines[l]).split(",")
            sRes = Res[0]
            for i in range(1,len(Res)):
                if (i!=3) and (i!=4):
                    sRes=sRes+","+Res[i]
            fres.write((CurrLine[0:6]).strip()+","+(CurrLine[7:17]).strip()+","+(CurrLine[18:28]).strip()+","+(CurrLine[29:40]).strip()+","+sRes)

            #print f2lines[l]
            l=l+1
            
    fres.close()
    return

######################################################################################
# MAIN function
######################################################################################

# run a batch of VMD command in order to calculate SASA features for interface residues
# between protein monomers

def RunSASA_VMD(PDBor,Mon1,Mon2,sFileConsurf):
    print (time.strftime("%H:%M:%S"))

    # [0] Get from interface PDB file and monomers
    # PDBor = "1MNM.pdb" # original PDB

    Monomer1= Mon1.split(",") #["A","B","C"] # protein chains
    Monomer2= Mon2.split(",") #["E","F"]     # protein chains

    PDBfile=PDBor[:-4]+"_ATOM.pdb" # PDB ATOM
    PDB=PDBfile[:-4]               # PDB name without ".pdb"
    print "PDB file = ",PDBfile
    print "Monomer 1 = ",Monomer1
    print "Monomer 2 = ",Monomer2
    print "PDBor= ",PDBor
    # [1] Convert PDB in PDB ATOM
    WritePDB_ATOM(PDBor,Monomer1,Monomer2)

    # [2] Add hydrogen
    # Generate spefici AddH script
    AddH("script_AddH",PDB,Monomer1,Monomer2)
    # run AddH script
    cmd='vmd -dispdev text -e script_AddH>script_AddH.out'
    print cmd
    os.system(cmd)

    # [3]Create 2 monomers
    Create_monomers("script_create_monomers",PDB,Monomer1,Monomer2)
    cmd='vmd -dispdev text -e script_create_monomers>script_create_monomers.out'
    print cmd
    os.system(cmd)

    # [4] Calculate SASA features

    # for the complex
    RunSASA("script_runSASA",PDB)
    cmd='vmd -dispdev text -e script_runSASA > script_runSASA.out'
    print cmd
    os.system(cmd)
    os.system('mv res_sasa.dat res_sasa1.dat')

    # for the monomer 1
    RunSASA_mon1("script_runSASA_mon1",PDB)
    cmd='vmd -dispdev text -e script_runSASA_mon1 > script_runSASA_mon1.out'
    print cmd
    os.system(cmd)
    os.system('mv res_sasa.dat res_sasaA1.dat')

    # for the monomer 2
    RunSASA_mon2("script_runSASA_mon2",PDB)
    cmd='vmd -dispdev text -e script_runSASA_mon2 > script_runSASA_mon2.out'
    print cmd
    os.system(cmd)
    os.system('mv res_sasa.dat res_sasaB1.dat')

    # [5]Feature calculation using SASA_features.py => result.csv
    cmd='python SASA_features.py'
    print cmd
    os.system(cmd)

    # [5.x] Correction of result.csv with the residue numbering from PDB
    #       (add 3 new columns at the begining of result.csv)
    CorrectResults("result.csv","result.corrected.csv",PDBor,Monomer1+Monomer2)

    # [6] Interface amino acids
    # Write scripts for P-P interface amino acids
    GetInterfaceResidues1("script_residues1", "residues1",PDBor,Monomer1,Monomer2)
    GetInterfaceResidues2("script_residues2", "residues2",PDBor,Monomer1,Monomer2)

    # Calculate the P-P interface AA
    cmd='vmd -dispdev text -e script_residues1'
    print cmd
    os.system(cmd)

    cmd='vmd -dispdev text -e script_residues2'
    print cmd
    os.system(cmd)

    # filter the corrected results
    FilterResults("result.corrected.csv","residues1","residues2","result.filtered.csv")

    # Add consurf info
    # Get consurf data
    ConsurfList = GetConsurfScores(sFileConsurf)
    
    Chains = (Mon1+","+Mon2).split(",") # all protein residues
    print Chains

    # add consurf data to filterd results ==> final results
    AddConsurf(Chains,ConsurfList,"result.filtered.csv","result.final.csv")

    # Run Weka model
    # --------------------
    #remove from file unused features and convert to arff format
    cmd='java weka.filters.unsupervised.attribute.Remove -R 1,2,3,4,5,6,7,9,10,11,12,14,15,16 -i result.final.csv -o result.final.removed.arff'
    print cmd
    os.system(cmd)

    #create attribute Class, nominal with two levels and add a quotation mark at the end of each row
    cmd='java weka.filters.unsupervised.attribute.Add -T NOM -N Class -L NS,HS -C last -i result.final.removed.arff -o result.final.removed.quotMark.arff'
    print cmd
    os.system(cmd)
    
    #generate output for each value
    cmd='java weka.classifiers.bayes.BayesNet -l ppModel.model -T result.final.removed.quotMark.arff -p 0 > output.txt'
    print cmd
    os.system(cmd)

    # Combine Weka and SASA results into the final file
    sFile1="output.txt"
    sFile2="result.final.csv"
    sFinal="sasa.hs.pp.results.csv"
    FinalRes(sFile1,sFile2,sFinal)

    print (time.strftime("%H:%M:%S"))
    return

#################################################################
# main program
#################################################################

### get the brute Web lists
WebPDBFile=str(sys.argv[1])
WebMon1=str(sys.argv[2])
WebMon2=str(sys.argv[3])
WebConsurf=str(sys.argv[4])
'''import pymol
from pymol import cmd,stored,math
dirname=os.getcwd()+'/'+WebPDBFile.split('.')[0]
print dirname
pymol.finish_launching()
cmd.load(dirname)
cmd.remove('hetatm')
cmd.save(dirname)
cmd.quit()'''

##WebPDBFile="1VFB.pdb"
##WebMon1="A,B"
##WebMon2="C,D"
##WebConsurf="consurf_1VFB_ABC.grades"

# run the main function
RunSASA_VMD(WebPDBFile,WebMon1,WebMon2,WebConsurf)

