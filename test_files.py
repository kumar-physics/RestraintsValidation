import sys
from mmcif.io.PdbxReader import PdbxReader
import pynmrstar
import numpy


class TestFiles(object):

    def __init__(self):
        pass

    def test_nef_files(self,nef_file):
        data = pynmrstar.Entry.from_file(nef_file)
        for saveframe in data:
            for loop in saveframe:
                print (saveframe.name,loop.category)

    @staticmethod
    def get_coordinates(cif_file):
        """
        Extract coordinate information from cif file as a dictionary
        {model_id : {(seq_id,chain_id,res_id,atom_id) : array[x,y,x],...},...}
        :param cif_file: Input coordinate file
        :return: dictionary
        """
        cif_data = []
        ifh = open(cif_file, 'r')
        pRd = PdbxReader(ifh)
        pRd.read(cif_data)
        ifh.close()
        c0 = cif_data[0]
        atom_site = c0.getObj('atom_site')
        max_models = int(atom_site.getValue('pdbx_PDB_model_num', -1))
        col_names = atom_site.getAttributeList()
        model_id = col_names.index('pdbx_PDB_model_num')
        x_id = col_names.index('Cartn_x')
        y_id = col_names.index('Cartn_y')
        z_id = col_names.index('Cartn_z')
        atom_id = col_names.index('label_atom_id')
        comp_id = col_names.index('label_comp_id')
        asym_id = col_names.index('label_asym_id')
        entity_id = col_names.index('label_entity_id')
        seq_id = col_names.index('label_seq_id')
        icode_id = col_names.index('pdbx_PDB_ins_code')
        alt_id = col_names.index('label_alt_id')
        aut_seq_id = col_names.index('auth_seq_id')
        aut_aym_id = col_names.index('auth_asym_id')
        pdb_models = {}
        atom_ids = {}
        for model in range(1, max_models + 1):
            pdb = {}
            aid = {}
            for dat in atom_site.getRowList():
                if int(dat[model_id]) == model:
                    aid[(dat[seq_id], dat[asym_id], dat[comp_id], dat[atom_id])] = \
                        (dat[entity_id], dat[asym_id], dat[comp_id], dat[seq_id], dat[aut_seq_id],
                         dat[alt_id], dat[icode_id], dat[aut_aym_id])
                    pdb[(dat[seq_id], dat[asym_id], dat[comp_id], dat[atom_id])] = \
                        numpy.array([float(dat[x_id]), float(dat[y_id]), float(dat[z_id])])
            pdb_models[model] = pdb
            atom_ids[model] = aid
        return pdb_models, atom_ids

    def test_cif_file(self,cif_file):
        pdb,atoms = self.get_coordinates(cif_file)
        print (pdb)
        #print (atoms)



if __name__ == "__main__":
    p = TestFiles()
    p.test_nef_files('vtf_examples/CtR107.nef')
    p.test_cif_file('vtf_examples/CtR107.cif')
