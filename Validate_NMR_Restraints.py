import sys
import os.path
from numpy import logical_and

# print sys.path
sys.path.append("~/git/py-mmcif")
from mmcif.io.PdbxReader import PdbxReader
from NEFTranslator import NEFTranslator
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
        chain_dict = self.get_chain_map(star_file)
        nt = NEFTranslator.NEFTranslator()
        print(nt.validate_file(star_file, 'R'))
        sd = pynmrstar.Entry.from_file(star_file)
        psudo_atoms = nt.validate_atom(sd, lp_category='_Gen_dist_constraint', seq_id='Comp_index_ID_1',
                                       res_id='Comp_ID_1', atom_id='Atom_ID_1')
        for i in psudo_atoms:
            res = i[1]
            atm = i[2]
            atm = atm.replace('M', 'H')
            atm = atm.replace('Q', 'H')
            atm = atm + '*'
            print(i, nt.get_nmrstar_atom(res, atm))
        pdb = self.get_coordinates(cif_file)
        distance, angle = self.get_restraints(star_file, chain_dict)
        v=self.calculate_distance_violations(pdb, distance)
        self.calculate_distance_violation_statistics(v)

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
    def get_restraints(star_file, chain_dict):
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
            dist_dict2 = {}
            dist_dict = {}
            for sf in dist_sf:
                if lp_flg:
                    dat = sf
                else:
                    dat = sf.get_loop_by_category('_Gen_dist_constraint')
                col_names = dat.get_tag_names()
                rest_id = col_names.index('_Gen_dist_constraint.ID')
                seq_id_1 = col_names.index('_Gen_dist_constraint.Comp_index_ID_1')
                asym_id_1 = col_names.index('_Gen_dist_constraint.Auth_asym_ID_1')
                entity_id_1 = col_names.index('_Gen_dist_constraint.Entity_assembly_ID_1')
                comp_id_1 = col_names.index('_Gen_dist_constraint.Comp_ID_1')
                atom_id_1 = col_names.index('_Gen_dist_constraint.Atom_ID_1')
                seq_id_2 = col_names.index('_Gen_dist_constraint.Comp_index_ID_2')
                asym_id_2 = col_names.index('_Gen_dist_constraint.Auth_asym_ID_2')
                entity_id_2 = col_names.index('_Gen_dist_constraint.Entity_assembly_ID_2')
                comp_id_2 = col_names.index('_Gen_dist_constraint.Comp_ID_2')
                atom_id_2 = col_names.index('_Gen_dist_constraint.Atom_ID_2')
                lb_id = col_names.index('_Gen_dist_constraint.Distance_lower_bound_val')
                ub_id = col_names.index('_Gen_dist_constraint.Distance_upper_bound_val')
                list_id = col_names.index('_Gen_dist_constraint.Gen_dist_constraint_list_ID')
                r_dict = {}
                for rest in dat:
                    if rest[rest_id] not in r_dict.keys():
                        r_dict[rest[rest_id]] = []
                    if (rest[list_id], rest[rest_id]) not in dist_dict2.keys():
                        dist_dict2[(rest[list_id], rest[rest_id])] = []
                    if rest[asym_id_1] == '.':
                        eid1 = chain_dict[rest[entity_id_1]]
                    else:
                        eid1 = rest[asym_id_1]
                    if rest[asym_id_2] == '.':
                        eid2 = chain_dict[rest[entity_id_2]]
                    else:
                        eid2 = rest[asym_id_2]

                    atom1 = (rest[seq_id_1], eid1, rest[comp_id_1], rest[atom_id_1])
                    atom2 = (rest[seq_id_2], eid2, rest[comp_id_2], rest[atom_id_2])
                    if atom1[1] != atom2[1]:
                        cat = 'long'
                    elif abs(int(atom1[0]) - int(atom2[0])) == 0:
                        cat = 'intraresidue'
                    elif abs(int(atom1[0]) - int(atom2[0])) == 1:
                        cat = 'sequential'
                    elif 1 < abs(int(atom1[0]) - int(atom2[0])) < 5:
                        cat = 'medium'
                    elif abs(int(atom1[0]) - int(atom2[0])) >= 5:
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
                        dist_dict2[(rest[list_id], rest[rest_id])].append([atom1, atom2, cat, lb, ub])
                if lp_flg:
                    dist_dict['distance_restraints'] = r_dict
                else:
                    dist_dict[sf.name] = r_dict
        else:
            dist_dict = None
            dist_dict2 = None
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
            angle_dict2 = {}
            for sf in ang_sf:
                if lp_flg:
                    dat = sf
                else:
                    dat = sf.get_loop_by_category('_Torsion_angle_constraint')
                col_names = dat.get_tag_names()
                rest_id = col_names.index('_Torsion_angle_constraint.ID')
                rest_name_id = col_names.index("_Torsion_angle_constraint.Torsion_angle_name")
                seq_id_1 = col_names.index('_Torsion_angle_constraint.Comp_index_ID_1')
                asym_id_1 = col_names.index('_Torsion_angle_constraint.Auth_asym_ID_1')
                entity_id_1 = col_names.index('_Torsion_angle_constraint.Entity_assembly_ID_1')
                comp_id_1 = col_names.index('_Torsion_angle_constraint.Comp_ID_1')
                atom_id_1 = col_names.index('_Torsion_angle_constraint.Atom_ID_1')
                seq_id_2 = col_names.index('_Torsion_angle_constraint.Comp_index_ID_2')
                asym_id_2 = col_names.index('_Torsion_angle_constraint.Auth_asym_ID_2')
                entity_id_2 = col_names.index('_Torsion_angle_constraint.Entity_assembly_ID_2')
                comp_id_2 = col_names.index('_Torsion_angle_constraint.Comp_ID_2')
                atom_id_2 = col_names.index('_Torsion_angle_constraint.Atom_ID_2')
                seq_id_3 = col_names.index('_Torsion_angle_constraint.Comp_index_ID_3')
                asym_id_3 = col_names.index('_Torsion_angle_constraint.Auth_asym_ID_3')
                entity_id_3 = col_names.index('_Torsion_angle_constraint.Entity_assembly_ID_3')
                comp_id_3 = col_names.index('_Torsion_angle_constraint.Comp_ID_3')
                atom_id_3 = col_names.index('_Torsion_angle_constraint.Atom_ID_3')
                seq_id_4 = col_names.index('_Torsion_angle_constraint.Comp_index_ID_4')
                asym_id_4 = col_names.index('_Torsion_angle_constraint.Auth_asym_ID_4')
                entity_id_4 = col_names.index('_Torsion_angle_constraint.Entity_assembly_ID_4')
                comp_id_4 = col_names.index('_Torsion_angle_constraint.Comp_ID_4')
                atom_id_4 = col_names.index('_Torsion_angle_constraint.Atom_ID_4')
                lb_id = col_names.index('_Torsion_angle_constraint.Angle_lower_bound_val')
                ub_id = col_names.index('_Torsion_angle_constraint.Angle_upper_bound_val')
                list_id = col_names.index('_Torsion_angle_constraint.Torsion_angle_constraint_list_ID')
                r_dict = {}
                for rest in dat:
                    if rest[rest_id] not in r_dict.keys():
                        r_dict[rest[rest_id]] = []
                    if (rest[list_id], rest[rest_id]) not in angle_dict2.keys():
                        angle_dict2[(rest[list_id], rest[rest_id])] = []

                    if rest[asym_id_1] == '.':
                        eid1 = chain_dict[rest[entity_id_1]]
                    else:
                        eid1 = rest[asym_id_1]
                    if rest[asym_id_2] == '.':
                        eid2 = chain_dict[rest[entity_id_2]]
                    else:
                        eid2 = rest[asym_id_2]
                    if rest[asym_id_3] == '.':
                        eid3 = chain_dict[rest[entity_id_3]]
                    else:
                        eid3 = rest[asym_id_3]
                    if rest[asym_id_4] == '.':
                        eid4 = chain_dict[rest[entity_id_4]]
                    else:
                        eid4 = rest[asym_id_4]

                    atom1 = (rest[seq_id_1], eid1, rest[comp_id_1], rest[atom_id_1])
                    atom2 = (rest[seq_id_2], eid2, rest[comp_id_2], rest[atom_id_2])
                    atom3 = (rest[seq_id_3], eid3, rest[comp_id_3], rest[atom_id_3])
                    atom4 = (rest[seq_id_4], eid4, rest[comp_id_4], rest[atom_id_4])
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
                        angle_dict2[(rest[list_id], rest[rest_id])].append([atom1, atom2, atom3, atom4, rest_name, lb, ub])
                if lp_flg:
                    angle_dict['angle_restraints'] = r_dict
                else:
                    angle_dict[sf.name] = r_dict
        else:
            angle_dict = None
            angle_dict2 = None
        return dist_dict2, angle_dict2

    @staticmethod
    def get_distance(c1, c2):

        """ Calculates the distance between two coordinate values.         6167  4830  .  A  478  ARG  HD%   B  2    C    H2'   1  .  .  .  .  4     .

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

    @staticmethod
    def r6sum(dist_list):
        return (sum([i ** (-6.) for i in dist_list])) ** (-1. / 6.)

    @staticmethod
    def r6average(dist_list):
        return (sum([i ** (-6.) for i in dist_list]) / float(len(dist_list))) ** (-1. / 6.)

    def calculate_distance_violations(self, coordinates, restraints):
        violations = {}
        for rest_id in restraints.keys():
            m ={}
            for model in coordinates.keys():
                dist_list=[]
                for rest in restraints[rest_id]:
                    atom_1 = rest[0]
                    atom_2 = rest[1]
                    cat = rest[2]
                    lb = rest[3]
                    ub = rest[4]
                    pos_1 = coordinates[model][atom_1]
                    pos_2 = coordinates[model][atom_2]
                    d = self.get_distance(pos_1,pos_2)
                    dist_list.append(d)
                r6dist = self.r6sum(dist_list)
                if lb <= r6dist <= ub:
                    err = 0.0
                elif r6dist < lb:
                    err = abs(r6dist-lb)
                else:
                    err = abs(r6dist-ub)
                m[model] = err
            violations[rest_id] = m
        return violations

    def calculate_angle_violations(self, coordinates, restraints):
        violations = {}
        for rest_id in restraints.keys():
            m ={}
            for model in coordinates.keys():
                ang_list=[]
                for rest in restraints[rest_id]:
                    atom_1 = rest[0]
                    atom_2 = rest[1]
                    atom_3 = rest[2]
                    atom_4 = rest[3]
                    cat = rest[4]
                    lb = rest[5]
                    ub = rest[6]
                    pos_1 = coordinates[model][atom_1]
                    pos_2 = coordinates[model][atom_2]
                    pos_3 = coordinates[model][atom_3]
                    pos_4 = coordinates[model][atom_4]
                    ang = self.get_dihedral_angle(pos_1,pos_2,pos_3,pos_4)
                    ang_list.append(ang)
                avg_viol = numpy.mean(ang_list)
                if lb <= avg_viol <= ub:
                    err = 0.0
                elif avg_viol < lb:
                    err = abs(avg_viol-lb)
                else:
                    err = abs(avg_viol-ub)
                m[model] = err
            violations[rest_id] = m
        return violations
    def calculate_distance_violation_statistics(self,violations):
        for rest_id in violations.keys():
            v = []
            mkeys = violations[rest_id].keys()
            for model in violations[rest_id].keys():

                if violations[rest_id][model]>0.0:
                    v.append(violations[rest_id][model])
            if len(v)>0:
                print (rest_id,len(v),numpy.mean(v),max(v),min(v))
            else:
                print (rest_id,'Not Violated')
        for model in mkeys:
            v=[]
            for rest_id in violations.keys():
                if violations[rest_id][model] > 0:
                    v.append(violations[rest_id][model])
            if len(v)>0:
                print(model, len(v), numpy.mean(v), max(v), min(v))
            else:
                print(model, 'Not Violated')






    # def validate_distace_restraints(self, coordinates, restraints):
    #     violations = {}
    #     for k in restraints.keys():
    #         viol = {}
    #         for m in restraints[k].keys():
    #             modl = {}
    #             for i in coordinates.keys():
    #                 d = []
    #                 for n in restraints[k][m]:
    #                     d.append(self.get_distance(coordinates[i][n[0]], coordinates[i][n[1]]))
    #                     # print(m, i, n[0], n[1], self.get_distance(coordinates[i][n[0]],coordinates[i][n[1]]))
    #                 ed = self.r6sum(d)
    #                 lb = n[2]
    #                 ub = n[3]
    #                 if lb <= ed <= ub:
    #                     err = 0.0
    #                 elif ed < lb:
    #                     err = abs(ed - lb)
    #                 else:
    #                     err = abs(ed - ub)
    #                 modl[i] = (n[4], lb, ub, ed, err)
    #             viol[m] = modl
    #         violations[k] = viol
    #     return violations
    #

    # def validate_angle_restraints(self, coordinates, restraints):
    #     violations = {}
    #     for k in restraints.keys():
    #         viol = {}
    #         for m in restraints[k].keys():
    #             modl = {}
    #             for i in coordinates.keys():
    #                 d = []
    #                 for n in restraints[k][m]:
    #                     d.append(self.get_dihedral_angle(coordinates[i][n[0]], coordinates[i][n[1]],coordinates[i][n[2]], coordinates[i][n[3]]))
    #                 ed = numpy.average(d) # No combinatorial angle restraints at present;but for future need
    #                 lb = n[5]
    #                 ub = n[6]
    #                 if lb <= ed <= ub:
    #                     err = 0.0
    #                 elif ed < lb:
    #                     err = abs(ed - lb)
    #                 else:
    #                     err = abs(ed - ub)
    #                 modl[i] = (n[4],lb, ub, ed)
    #             viol[m] = modl
    #         violations[k] = viol
    #     return violations
    # def dist_violoation_statistics(self,dist_viol):
    #     types = []
    #     restraint_types={}
    #     for k in dist_viol.keys():
    #         for r in dist_viol[k].keys():
    #             types.append(dist_viol[k][r][1][0])
    #     types_stat = {i:types.count(i) for i in set(types)}
    #     types_stat['total']= len(types)
    #     viol_stat = {}
    #     for k in dist_viol.keys():
    #         v={}
    #         for r in dist_viol[k].keys():
    #             c=0
    #             e=[]
    #             for m in dist_viol[k][r]:
    #                 e.append(dist_viol[k][r][m][-1])
    #                 if dist_viol[k][r][m][-1]>0.0:
    #                     c+=1
    #             v[r]=[dist_viol[k][r][1][0],c,numpy.mean(e),max(e)]
    #         viol_stat[k]=v
    #     for k in viol_stat.keys():
    #         for r in viol_stat[k].keys():
    #             print (k,r,viol_stat[k][r])
    #     viol_stat2={}
    #     for k in dist_viol.keys():
    #         for i in dist_viol[k][list(dist_viol[k].keys())[0]].keys():
    #             c=0
    #             for r in dist_viol[k].keys():
    #                 if dist_viol[k][r][i][-1]>0.0:
    #                     c+=1
    #             print (k,i,r,c)


if __name__ == "__main__":
    # p = ValidateRestraints('nef_examples/2mqq.cif', 'nef_examples/2mqq.str')
    p = ValidateRestraints('nef_examples/2mqq.cif', 'nef_examples/2mqq.str')
