

import json
import pynmrstar
# from monkeypatch import patch_parser
# patch_parser(pynmrstar)
import re

class FixNomenclature(object):

    def __init__(self,star_file):
        (isOk, msg, self.atomDict) = self.load_json_data('NEFTranslator/lib/atomDict.json')
        if not isOk:
            self.logger.error(msg)
        self.chains = self.get_chain_map(star_file)
        self.fix_atom_names(star_file)
    @staticmethod
    def load_json_data(json_file):
        """
        Loads json data files from lib folder
        :param json_file: json file
        :return: dictionay
        """
        try:
            with open(json_file, 'r') as jsonF:
                data_dict = json.loads(jsonF.read())
            is_ok = True
            msg = "{} file is read!".format(json_file)
        except IOError:
            msg = "{} file is missing!".format(json_file)
            is_ok = False
            data_dict = []
        return is_ok, msg, data_dict

    @staticmethod
    def get_chain_map(star_file):
        ent = pynmrstar.Entry.from_file(star_file)
        try:
            entity_assembly = ent.get_loops_by_category('_Entity_assembly')[0]
            col_names = entity_assembly.get_tag_names()
            id_index = col_names.index('_Entity_assembly.ID')
            asym_index = col_names.index('_Entity_assembly.Asym_ID')
            chain_map = {}
            for row in entity_assembly:
                chain_map[row[id_index]] = row[asym_index]
        except IndexError:
            print("Entity_assembly loop no found, No chain mapping created")
            chain_map = {}
        return chain_map

    def get_nmrstar_atom(self, res, nef_atom):
        """
        Returns (atom with out wildcard,[IUPAC atom list],ambiguity code)

        """
        ambiguity_code = 1
        atom_type = None
        try:
            atoms = self.atomDict[res]
            atom_list = []
            try:
                ref_atom = re.findall(r'(\S+)([xyXY])([%*])$|(\S+)([%*])$|(\S+)([xyXY]$)', nef_atom)[0]
                atm_set = [ref_atom.index(i) for i in ref_atom if i != ""]
                pattern = None
                if atm_set == [0, 1, 2]:
                    atom_type = ref_atom[0]
                    pattern = re.compile(r'%s\S\d+' % (ref_atom[0]))
                    alist2 = [i for i in atoms if re.search(pattern, i)]
                    xid = sorted(set([int(i[len(ref_atom[0])]) for i in alist2]))
                    if ref_atom[1] == "x" or ref_atom[1] == "X":
                        atom_list = [i for i in alist2 if int(i[len(ref_atom[0])]) == xid[0]]
                    else:
                        atom_list = [i for i in alist2 if int(i[len(ref_atom[0])]) == xid[1]]
                    ambiguity_code = 2
                elif atm_set == [3, 4]:
                    atom_type = ref_atom[3]
                    if ref_atom[4] == "%":
                        pattern = re.compile(r'%s\d+' % (ref_atom[3]))
                    elif ref_atom[4] == "*":
                        pattern = re.compile(r'%s\S+' % (ref_atom[3]))
                    else:
                        logging.critical("Wrong NEF atom {}".format(nef_atom))
                    atom_list = [i for i in atoms if re.search(pattern, i)]
                    ambiguity_code = 1

                elif atm_set == [5, 6]:
                    atom_type = ref_atom[5]
                    pattern = re.compile(r'%s\S+' % (ref_atom[5]))
                    atom_list = [i for i in atoms if re.search(pattern, i)]
                    if len(atom_list) != 2:
                        atom_list = []
                    elif ref_atom[6] == "y" or ref_atom[6] == "Y":
                        # atom_list.reverse()[]
                        atom_list = atom_list[-1:]
                    elif ref_atom[6] == "x" or ref_atom[6] == "X":
                        atom_list = atom_list[:1]
                    else:
                        logging.critical("Wrong NEF atom {}".format(nef_atom))
                    ambiguity_code = 2

                else:
                    logging.critical("Wrong NEF atom {}".format(nef_atom))
            except IndexError:

                # print nef_atom
                pass
                atom_type = nef_atom
            if len(atom_list) == 0:
                if nef_atom in atoms:
                    atom_list.append(nef_atom)
                else:
                    if nef_atom == "H%":  # To handle terminal protons
                        atom_list = ['H1', 'H2', 'H3']
                        atom_type = "H"
        except KeyError:
            # self.logfile.write("%s\tResidue not found,%s,%s\n"%(self.TimeStamp(time.time()),res,nef_atom))
            # print "Residue not found",res,nef_atom
            if res != ".":
               print ("Non-standard residue found {}".format(res))
            atom_list = []
            atom_type = nef_atom

            if nef_atom == "H%":
                atom_list = ['H1', 'H2', 'H3']
                atom_type = "H"
        return atom_type, atom_list, ambiguity_code

    def fix_atom_names(self,star_file):
        star_data = pynmrstar.Entry.from_file(star_file)
        rest_loops = star_data.get_loops_by_category('_Gen_dist_constraint')
        index = 1

        for dat in rest_loops:
            #dat.add_tag('_Gen_dist_constraint.Index_ID')
            col_names = dat.get_tag_names()
            rest_id = col_names.index('_Gen_dist_constraint.ID')
            logic_id = col_names.index('_Gen_dist_constraint.Member_logic_code')
            seq_id_1 = col_names.index('_Gen_dist_constraint.Comp_index_ID_1')
            asym_id_1 = col_names.index('_Gen_dist_constraint.Auth_asym_ID_1')
            entity_id_1 = col_names.index('_Gen_dist_constraint.Entity_assembly_ID_1')
            strand_id_1 = col_names.index('_Gen_dist_constraint.PDB_strand_ID_1')
            comp_id_1 = col_names.index('_Gen_dist_constraint.Comp_ID_1')
            atom_id_1 = col_names.index('_Gen_dist_constraint.Atom_ID_1')
            auth_atom_id_1 = col_names.index('_Gen_dist_constraint.Auth_atom_ID_1')
            seq_id_2 = col_names.index('_Gen_dist_constraint.Comp_index_ID_2')
            asym_id_2 = col_names.index('_Gen_dist_constraint.Auth_asym_ID_2')
            entity_id_2 = col_names.index('_Gen_dist_constraint.Entity_assembly_ID_2')
            strand_id_2 = col_names.index('_Gen_dist_constraint.PDB_strand_ID_2')
            comp_id_2 = col_names.index('_Gen_dist_constraint.Comp_ID_2')
            atom_id_2 = col_names.index('_Gen_dist_constraint.Atom_ID_2')
            auth_atom_id_2 = col_names.index('_Gen_dist_constraint.Auth_atom_ID_2')
            lb_id = col_names.index('_Gen_dist_constraint.Distance_lower_bound_val')
            ub_id = col_names.index('_Gen_dist_constraint.Distance_upper_bound_val')
            r_dict = {}
            dat2 = []
            for rest2 in dat:
                rest = rest2[:]
                rid = rest[rest_id]
                res1 = rest[comp_id_1]
                atm1 = rest[atom_id_1]
                if 'Q' in atm1:
                    atm1=atm1.replace('Q','H')
                    atm1 += '*'
                elif 'M' in atm1:
                    atm1=atm1.replace('M','H')
                    atm1 += '*'
                else:
                    pass
                (na1,new_atms1,ac1)=self.get_nmrstar_atom(res1,atm1)
                res2 = rest[comp_id_2]
                atm2 = rest[atom_id_2]
                if 'Q' in atm2:
                    atm2 = atm2.replace('Q', 'H')
                    atm2 += '*'
                elif 'M' in atm2:
                    atm2 = atm2.replace('M', 'H')
                    atm2 += '*'
                else:
                    pass
                (na2,new_atms2,ac2) = self.get_nmrstar_atom(res2, atm2)

                rest[asym_id_1] = self.chains[rest[entity_id_1]]
                rest[asym_id_2] = self.chains[rest[entity_id_2]]
                if len(new_atms1)==0 or len(new_atms2)==0:
                    print ("Here",rest)
                if len(new_atms1) > 1 or len(new_atms2) > 1:
                    for a1 in new_atms1:
                        for a2 in new_atms2:
                            rest[rest_id] = rid
                            rest[atom_id_1] = a1
                            rest[atom_id_2] = a2
                            rest[logic_id] = 'OR'
                            dat2.append(rest[:])
                else:
                    dat2.append(rest[:])
            dat.clear_data()
            # print (dat)
            # print (len(dat2),len(dat2[0]))
            for item in dat2:
                if len(item)!=len(dat.tags):
                    print(len(item), len(dat.tags))
            dat.data = dat2
            dat: pynmrstar.Loop
            dat.sort_rows(['_Gen_dist_constraint.ID'])
            dat.add_tag('Index_ID', update_data=True)
            dat.renumber_rows('Index_ID')
            dat.sort_tags()
        with open('test_out.str', 'w') as wstarfile:
            wstarfile.write(str(star_data))




if __name__ == "__main__":
    entry_id = '4ch1'
    url = 'http://www.bmrb.wisc.edu/ftp/pub/bmrb/nmr_pdb_integrated_data/coordinates_restraints_chemshifts/all/' \
          'nmr-star/{}/{}_linked.str'.format(entry_id.lower(), entry_id.lower())
    p = FixNomenclature(url)

