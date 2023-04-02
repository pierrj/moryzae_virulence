from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBIO

def cif2pdb(gid, cif):
  #Convert CIF into PDB and then save into the given directory
  parser = MMCIFParser()
  structure = parser.get_structure('A', cif)

  io = PDBIO()
  io.set_structure(structure)
  io.save(f'{gid}.AF2DB.pdb')


#Go through AF2.list and convert the match.cif to pdb 
for line in open('AF2.list', 'r'):
    repSeq = line.split()[2]
    cif = line.split()[1]
    
    try:
        cif2pdb(f'{repSeq}/{repSeq}',
                f'/global/scratch/users/skyungyong/CO_Pierre_MO/Analysis/Structures-DB/AF2/AF-{cif}-F1-model_v3.cif')
    except Exception:
        print(f'{repSeq} match {cif} is missing')
