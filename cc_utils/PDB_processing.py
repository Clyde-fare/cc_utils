
__author__ = 'clyde'


##taken from http://stackoverflow.com/questions/11685716/how-to-extract-chains-from-a-pdb-file
import os
from Bio import PDB


class ChainSplitter:
    def __init__(self, out_dir=None):
        """ Create parsing and writing objects, specify output directory. """
        self.parser = PDB.PDBParser()
        self.writer = PDB.PDBIO()
        if out_dir is None:
            out_dir = os.path.join(os.getcwd(), "chain_PDBs")
        self.out_dir = out_dir

    def make_pdb(self, pdb_path, chain_letters, overwrite=False, struct=None):
        """ Create a new PDB file containing only the specified chains.

        Returns the path to the created file.

        :param pdb_path: full path to the crystal structure
        :param chain_letters: iterable of chain characters (case insensitive)
        :param overwrite: write over the output file if it exists
        """
        chain_letters = [chain.upper() for chain in chain_letters]

        # Input/output files
        (pdb_dir, pdb_fn) = os.path.split(pdb_path)
        pdb_id = pdb_fn[3:7]
        out_name = "pdb%s_%s.ent" % (pdb_id, "".join(chain_letters))
        out_path = os.path.join(self.out_dir, out_name)
        print "OUT PATH:",out_path
        plural = "s" if (len(chain_letters) > 1) else ""  # for printing

        # Skip PDB generation if the file already exists
        if (not overwrite) and (os.path.isfile(out_path)):
            print("Chain%s %s of '%s' already extracted to '%s'." %
                  (plural, ", ".join(chain_letters), pdb_id, out_name))
            return out_path

        print("Extracting chain%s %s from %s..." % (plural,
                                                    ", ".join(chain_letters), pdb_fn))

        # Get structure, write new file with only given chains
        if struct is None:
            struct = self.parser.get_structure(pdb_id, pdb_path)
        self.writer.set_structure(struct)
        self.writer.save(out_path, select=SelectChains(chain_letters))

        return out_path


class SelectChains(PDB.Select):
    """ Only accept the specified chains when saving. """
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)


if __name__ == "__main__":
    """ Parses PDB id's desired chains, and creates new PDB structures. """
    import sys
    if not len(sys.argv) == 2:
        print "Usage: $ python %s 'pdb.txt'" % __file__
        sys.exit()

    pdb_textfn = sys.argv[1]

    pdbList = PDB.PDBList()
    splitter = ChainSplitter("/home/steve/chain_pdbs")  # Change me.

    with open(pdb_textfn) as pdb_textfile:
        for line in pdb_textfile:
            pdb_id = line[:4].lower()
            chain = line[4]
            pdb_fn = pdbList.retrieve_pdb_file(pdb_id)
            splitter.make_pdb(pdb_fn, chain)


from simtk.openmm.app.pdbfile import PDBFile


def get_non_standard_res(file_n, chain=0):
    raw_s = PDBFile(file_n)
    chain = raw_s.topology._chains[chain]
    residues = chain._residues

    non_standard = []
    for r in residues:
        atoms = sorted([a.name for a in r._atoms])
        if expected_atoms(r.name,chain) != atoms:
            non_standard.append(r)
    return non_standard


def expected_atoms(res_name, chain):
    res_bonds = chain.topology._standardBonds[res_name]
    unique_atoms = list(set([b[0] for b in res_bonds] + [b[1] for b in res_bonds]))
    backbone_atoms = [a for a in unique_atoms if 'H' not in a and 'XT' not in a and '-' not in a]
    return sorted(backbone_atoms)


from cc_notebook_utils import unwind
from molmod.io import XYZReader
from chem_utils import mol2_parse
def get_interfacial_atoms(pdb_code, non_std_atom_indexes):
    master_xyz = XYZReader(pdb_code + '_master.xyz').get_first_molecule()
    master_xyz.set_default_graph()
    std_nn_ind = unwind([[n for n in master_xyz.graph.neighbors[i] if n not in non_std_atom_indexes] for i in non_std_atom_indexes])
    std_nnn_ind = [[n for n in master_xyz.graph.neighbors[i] if n not in non_std_atom_indexes] for i in std_nn_ind]
    nstd_nn_ind = [[n for n in master_xyz.graph.neighbors[i] if n in non_std_atom_indexes] for i in std_nn_ind]
    nstd_nnn_ind = [[n for n in master_xyz.graph.neighbors[i] if n in non_std_atom_indexes] for i in unwind(nstd_nn_ind)]

    gaff_atoms, gaff_bonds = mol2_parse(pdb_code +'_master_gaff.mol2')
    amber_atoms, amber_bonds = mol2_parse(pdb_code +'_master_amber.mol2')

    std_nn = [amber_atoms[i] for i in std_nn_ind]
    std_nnn = [[amber_atoms[i] for i in nnn_ind] for nnn_ind in std_nnn_ind]
    nstd_nn = [[gaff_atoms[i] for i in nn_ind] for nn_ind in nstd_nn_ind]
    nstd_nnn = [[gaff_atoms[i] for i in nnn_ind] for nnn_ind in nstd_nnn_ind]

    bonds = []
    for i in range(len(std_nn)):
        for k in range(len(nstd_nn[i])):
            bonds.append(std_nn[i], nstd_nn[i][k])

    dihedrals = []
    for i in range(len(std_nn)):
        for j in range(len(std_nnn[i])):
            for k in range(len(nstd_nn[i])):
                for l in range(len(nstd_nnn[i])):
                    dihedrals.append([std_nnn[i][j], std_nn[i], nstd_nn[i][k], nstd_nnn[i][l]])

    angles = []
    for dihedral in dihedrals:
        ang1, ang2 = dihedral[0:-1], dihedral[1:]
        if ang1 not in angles:
            angles.append(ang1)
        if ang2 not in angles:
            angles.append(ang2)

    return