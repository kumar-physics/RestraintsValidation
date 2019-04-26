import sys
import os.path
from numpy import logical_and

# print sys.path
sys.path.append("/home/kumaran/git/py-mmcif")
from mmcif.io.PdbxReader import PdbxReader
from math import sqrt, acos
import numpy
import pynmrstar
from monkeypatch import patch_parser

patch_parser(pynmrstar)
import logging
import ntpath
import operator
import json


class ValidateRestraints:
    """
    NMR restraints validation module
    """
    __version__ = "v0.8"

    def __init__(self, cif_file, star_file):
        pdb = self.get_coordinates(cif_file)
        distance, angle = self.get_restraints(star_file)
        self.validate_distace_restraints(pdb, distance)
        #self.validate_angle_restraints(pdb, angle)

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
        pdb_models = {}
        for model in range(1, max_models + 1):
            pdb = {}
            for dat in atom_site.getRowList():
                if int(dat[model_id]) == model:
                    pdb[(dat[seq_id], dat[asym_id], dat[comp_id], dat[atom_id])] = \
                        numpy.array([float(dat[x_id]), float(dat[y_id]), float(dat[z_id])])
            pdb_models[model] = pdb
        return pdb_models

    @staticmethod
    def get_restraints(star_file):
        lp_flg = 0
        dat_flg = 1
        try:
            star_data = pynmrstar.Entry.from_file(star_file)
            dist_sf = star_data.get_saveframes_by_category('general_distance_constraints')
        except ValueError:
            try:
                sf_data = pynmrstar.Saveframe.from_file(star_file)
                if sf_data.get_tag('Sf_category')[0] == 'general_distance_constraints':
                    dist_sf = [sf_data]
                else:
                    print("Error : file doesn't have restraints data")
                    dat_flg = 0
            except ValueError:
                try:
                    lp_data = pynmrstar.Loop.from_file(star_file)
                    lp_flg = 1
                    if lp_data.category == '_Gen_dist_constraint':
                        dist_sf = [lp_data]
                    else:
                        print("Error : file doesn't have restraints data")
                        dat_flg = 0
                except ValueError:
                    print("ERROR : File contains no valid saveframe or loop")
                    dat_flg = 0
        except IOError:
            print("Error file not found")
            dat_flg = 0
        if dat_flg:
            dist_dict = {}
            for sf in dist_sf:
                if lp_flg:
                    dat = sf
                else:
                    dat = sf.get_loop_by_category('_Gen_dist_constraint')
                col_names = dat.get_tag_names()
                rest_id = col_names.index('_Gen_dist_constraint.ID')
                seq_id_1 = col_names.index('_Gen_dist_constraint.Comp_index_ID_1')
                entity_id_1 = col_names.index('_Gen_dist_constraint.Auth_asym_ID_1')
                # entity_id_1 = col_names.index('Gen_dist_constraint.Entity_assembly_ID_1')
                comp_id_1 = col_names.index('_Gen_dist_constraint.Comp_ID_1')
                atom_id_1 = col_names.index('_Gen_dist_constraint.Atom_ID_1')
                seq_id_2 = col_names.index('_Gen_dist_constraint.Comp_index_ID_2')
                entity_id_2 = col_names.index('_Gen_dist_constraint.Auth_asym_ID_2')
                # entity_id_2 = col_names.index('Gen_dist_constraint.Entity_assembly_ID_2')
                comp_id_2 = col_names.index('_Gen_dist_constraint.Comp_ID_2')
                atom_id_2 = col_names.index('_Gen_dist_constraint.Atom_ID_2')
                lb_id = col_names.index('_Gen_dist_constraint.Distance_lower_bound_val')
                ub_id = col_names.index('_Gen_dist_constraint.Distance_upper_bound_val')
                r_dict = {}
                for rest in dat:
                    if rest[rest_id] not in r_dict.keys():
                        r_dict[rest[rest_id]] = []
                    atom1 = (rest[seq_id_1], rest[entity_id_1], rest[comp_id_1], rest[atom_id_1])
                    atom2 = (rest[seq_id_2], rest[entity_id_2], rest[comp_id_2], rest[atom_id_2])
                    if atom1[1]!=atom2[1]:
                        cat = 'long'
                    elif abs(int(atom1[0])-int(atom2[0]))==0:
                        cat = 'intraresidue'
                    elif abs(int(atom1[0])-int(atom2[0]))==1:
                        cat = 'sequential'
                    elif 1 < abs(int(atom1[0])-int(atom2[0])) < 5:
                        cat = 'medium'
                    elif abs(int(atom1[0])-int(atom2[0]))>=5:
                        cat = 'long'
                    try:
                        lb = float(rest[lb_id])
                    except ValueError:
                        lb = -999.9
                    try:
                        ub = float(rest[ub_id])
                    except ValueError:
                        ub = 999.9
                    if lb == -999.9 and ub == 999.9:
                        print("Error: Distance restraint value not readable for restraint id {}; "
                              "for atoms {},{}".format(rest[rest_id], atom1, atom2))
                    else:
                        r_dict[rest[rest_id]].append([atom1, atom2, lb, ub, cat])
                if lp_flg:
                    dist_dict['distance_restraints'] = r_dict
                else:
                    dist_dict[sf.name] = r_dict
        else:
            dist_dict = None
        dat_flg = 1
        lp_flg = 0
        try:
            star_data = pynmrstar.Entry.from_file(star_file)
            ang_sf = star_data.get_saveframes_by_category('torsion_angle_constraints')
        except ValueError:
            try:
                sf_data = pynmrstar.Saveframe.from_file(star_file)
                if sf_data.get_tag('Sf_category')[0] == 'torsion_angle_constraints':
                    ang_sf = [sf_data]
                else:
                    print("Error : file doesn't have restraints data")
                    dat_flg = 0
            except ValueError:
                try:
                    lp_data = pynmrstar.Loop.from_file(star_file)
                    lp_flg = 1
                    if lp_data.category == '_Torsion_angle_constraint':
                        ang_sf = [lp_data]
                    else:
                        print("Error : file doesn't have restraints data")
                        dat_flg = 0
                except ValueError:
                    print("ERROR : File contains no valid saveframe or loop")
                    dat_flg = 0
        except IOError:
            print("Error file not found")
            dat_flg = 0
        if dat_flg:
            angle_dict = {}
            for sf in ang_sf:
                if lp_flg:
                    dat = sf
                else:
                    dat = sf.get_loop_by_category('_Torsion_angle_constraint')
                col_names = dat.get_tag_names()
                rest_id = col_names.index('_Torsion_angle_constraint.ID')
                rest_name_id = col_names.index("_Torsion_angle_constraint.Torsion_angle_name")
                seq_id_1 = col_names.index('_Torsion_angle_constraint.Comp_index_ID_1')
                entity_id_1 = col_names.index('_Torsion_angle_constraint.Auth_asym_ID_1')
                # entity_id_1 = col_names.index('_Torsion_angle_constraint.Entity_assembly_ID_1')
                comp_id_1 = col_names.index('_Torsion_angle_constraint.Comp_ID_1')
                atom_id_1 = col_names.index('_Torsion_angle_constraint.Atom_ID_1')
                seq_id_2 = col_names.index('_Torsion_angle_constraint.Comp_index_ID_2')
                entity_id_2 = col_names.index('_Torsion_angle_constraint.Auth_asym_ID_2')
                # entity_id_2 = col_names.index('_Torsion_angle_constraint.Entity_assembly_ID_2')
                comp_id_2 = col_names.index('_Torsion_angle_constraint.Comp_ID_2')
                atom_id_2 = col_names.index('_Torsion_angle_constraint.Atom_ID_2')
                seq_id_3 = col_names.index('_Torsion_angle_constraint.Comp_index_ID_3')
                entity_id_3 = col_names.index('_Torsion_angle_constraint.Auth_asym_ID_3')
                # entity_id_3 = col_names.index('_Torsion_angle_constraint.Entity_assembly_ID_3')
                comp_id_3 = col_names.index('_Torsion_angle_constraint.Comp_ID_3')
                atom_id_3 = col_names.index('_Torsion_angle_constraint.Atom_ID_3')
                seq_id_4 = col_names.index('_Torsion_angle_constraint.Comp_index_ID_4')
                entity_id_4 = col_names.index('_Torsion_angle_constraint.Auth_asym_ID_4')
                # entity_id_4 = col_names.index('_Torsion_angle_constraint.Entity_assembly_ID_4')
                comp_id_4 = col_names.index('_Torsion_angle_constraint.Comp_ID_4')
                atom_id_4 = col_names.index('_Torsion_angle_constraint.Atom_ID_4')
                lb_id = col_names.index('_Torsion_angle_constraint.Angle_lower_bound_val')
                ub_id = col_names.index('_Torsion_angle_constraint.Angle_upper_bound_val')
                r_dict = {}
                for rest in dat:
                    if rest[rest_id] not in r_dict.keys():
                        r_dict[rest[rest_id]] = []
                    atom1 = (rest[seq_id_1], rest[entity_id_1], rest[comp_id_1], rest[atom_id_1])
                    atom2 = (rest[seq_id_2], rest[entity_id_2], rest[comp_id_2], rest[atom_id_2])
                    atom3 = (rest[seq_id_3], rest[entity_id_3], rest[comp_id_3], rest[atom_id_3])
                    atom4 = (rest[seq_id_4], rest[entity_id_4], rest[comp_id_4], rest[atom_id_4])
                    rest_name = rest[rest_name_id]
                    try:
                        lb = float(rest[lb_id])
                    except ValueError:
                        lb = -999.9
                    try:
                        ub = float(rest[ub_id])
                    except ValueError:
                        ub = 999.9
                    if lb == -999.9 and ub == 999.9:
                        print("Error: Distance restraint value not readable for restraint id {}; "
                              "for atoms {},{}".format(rest[rest_id], atom1, atom2))
                    else:
                        r_dict[rest[rest_id]].append([atom1, atom2, atom3, atom4, rest_name, lb, ub])
                if lp_flg:
                    angle_dict['angle_restraints'] = r_dict
                else:
                    angle_dict[sf.name] = r_dict
        else:
            angle_dict = None
        return dist_dict, angle_dict

    @staticmethod
    def get_distance(c1, c2):

        """ Calculates the distance between two coordinate values.
        Each coordinate is an array of x,y,z. Returns distance in A assuming the input coordinates are in A"""

        return numpy.linalg.norm(c1 - c2)

    @staticmethod
    def get_dihedral_angle(c1, c2, c3, c4):

        """ Calculates the dihedral angle from the given four coordinate values.
        Each coordinate is an array of x,y,z. Returns angle in degrees"""

        bv12 = c1 - c2
        bv32 = c3 - c2
        bv43 = c4 - c3
        pv13 = numpy.cross(bv12, bv32)
        pv24 = numpy.cross(bv43, bv32)
        pro = numpy.dot(pv13, pv24)
        sqdist13 = numpy.dot(pv13, pv13)
        sqdist24 = numpy.dot(pv24, pv24)
        cosin = pro / sqrt(sqdist13 * sqdist24)
        cosin - min(1.0, max(-1.0, cosin))
        angle = acos(cosin)

        if numpy.dot(pv13, numpy.cross(pv24, bv32)) < 0:
            angle = -angle
        return round(numpy.degrees(angle), 4)

    @staticmethod
    def r6sum(dist_list):
        return (sum([i ** (-6.) for i in dist_list])) ** (-1. / 6.)

    @staticmethod
    def r6average(dist_list):
        return (sum([i ** (-6.) for i in dist_list]) / float(len(dist_list))) ** (-1. / 6.)

    def validate_distace_restraints(self, coordinates, restraints):
        violations = {}
        for k in restraints.keys():
            viol = {}
            for m in restraints[k].keys():
                modl = {}
                for i in coordinates.keys():
                    d = []
                    for n in restraints[k][m]:
                        d.append(self.get_distance(coordinates[i][n[0]], coordinates[i][n[1]]))
                        # print(m, i, n[0], n[1], self.get_distance(coordinates[i][n[0]],coordinates[i][n[1]]))
                    ed = self.r6sum(d)
                    lb = n[2]
                    ub = n[3]
                    if lb <= ed <= ub:
                        err = 0.0
                    elif ed < lb:
                        err = abs(ed - lb)
                    else:
                        err = abs(ed - ub)
                    modl[i] = (lb, ub, ed, err,n[4])
                viol[m] = modl
            violations[k] = viol
        for k in violations.keys():
            for m in violations[k].keys():
                for n in violations[k][m].keys():
                    print(k, m, n, violations[k][m][n])

    def validate_angle_restraints(self, coordinates, restraints):
        violations = {}
        for k in restraints.keys():
            viol = {}
            for m in restraints[k].keys():
                modl = {}
                for i in coordinates.keys():
                    d = []
                    for n in restraints[k][m]:
                        d.append(self.get_dihedral_angle(coordinates[i][n[0]], coordinates[i][n[1]],coordinates[i][n[2]], coordinates[i][n[3]]))
                    ed = numpy.average(d) # No combinatorial angle restraints at present;but for future need
                    lb = n[5]
                    ub = n[6]
                    if lb <= ed <= ub:
                        err = 0.0
                    elif ed < lb:
                        err = abs(ed - lb)
                    else:
                        err = abs(ed - ub)
                    modl[i] = (n[4],lb, ub, ed)
                viol[m] = modl
            violations[k] = viol
        for k in violations.keys():
            for m in violations[k].keys():
                for n in violations[k][m].keys():
                    print(k, m, n, violations[k][m][n])


if __name__ == "__main__":
    p = ValidateRestraints('nef_examples/2l9r.cif', 'nef_examples/2l9r.str')
