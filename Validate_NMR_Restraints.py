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
import xml.etree.ElementTree as ET
from operator import itemgetter

patch_parser(pynmrstar)

import json


class Validate_NMR_Restraints:
    """
    NMR restraints validation module
    """
    __version__ = "v0.8"

    def __init__(self):
        pass


    def run_validation(self,cif_file,star_file):
        nt = NEFTranslator.NEFTranslator()
        pdb, atom_ids = self.get_coordinates(cif_file)
        max_models = len(pdb.keys())
        distance, angle, chain_dict = self.get_restraints(star_file)
        dv = self.calculate_distance_violations(pdb, distance)
        av = self.calculate_angle_violations(pdb, angle)
        self.write_xml(dv, av, distance, angle, atom_ids)
        # self.write_xml_simple(dv, av)
        dist_viol_stat, dist_viol = self.calculate_violation_statistics(dv)
        ang_viol_stat, ang_viol = self.calculate_violation_statistics(av)
        self.bin_distance_violations(dv)
        sorted_dist_viol_stat = sorted(dist_viol_stat, reverse=True, key=itemgetter(0, 3))
        sorted_dist_viol = sorted(dist_viol, reverse=True, key=itemgetter(0))
        sorted_ang_viol_stat = sorted(ang_viol_stat, reverse=True, key=itemgetter(0, 3))
        sorted_ang_viol = sorted(ang_viol, reverse=True, key=itemgetter(0))
        type_stat_dist = self.restraints_type_statistics(dist_viol_stat, max_models)
        type_stat_ang = self.restraints_type_statistics(ang_viol_stat, max_models)
        json_data = self.generate_json(type_stat_dist, type_stat_ang, sorted_dist_viol_stat, sorted_dist_viol,
                                       sorted_ang_viol_stat, sorted_ang_viol)
        with open('data_json.json', 'w') as write_file:
            json.dump(json_data, write_file)

    @staticmethod
    def generate_json(type_stat_dist, type_stat_ang, dist_viol_stat, dist_viol, ang_viol_stat, ang_viol):
        restraints_validation = {}
        distance = {}
        distance['summary'] = type_stat_dist[0]
        distance['violated'] = type_stat_dist[2]
        distance['consistently_violated'] = type_stat_dist[1]
        distance['sorted_average_violations'] = dist_viol_stat
        distance['sorted_violations'] = dist_viol
        angle = {}
        angle['summary'] = type_stat_ang[0]
        angle['violated'] = type_stat_ang[2]
        angle['consistently_violated'] = type_stat_ang[1]
        angle['sorted_average_violations'] = ang_viol_stat
        angle['sorted_violations'] = ang_viol
        restraints_validation['distance'] = distance
        restraints_validation['angle'] = angle
        return restraints_validation

    @staticmethod
    def restraints_type_statistics(violation_statistics, max_models):
        """
        Counts the number of restraints and violations in each restraints type(distance: intraresidue, sequential,
        medium, long; angle : PHI, PSI, etc.. )
        :param violation_statistics: output from calculate_violation_statisitcs
        :param max_models: Number of models in the ensemble (required to estimate the consistently violated restraitns)
        :return: three statistics in a dictionary format (General summary, Consistently violated, Violated)
        """
        rest_type_v = [i[5] for i in violation_statistics if i[3] != 0]  # Violated at lease in one model
        rest_type_cv = [i[5] for i in violation_statistics if i[3] == max_models]  # Violated in all modes
        rest_type_all = [i[5] for i in violation_statistics]
        uniq_type = list(set(rest_type_all))
        type_count_all = {}
        type_count_v = {}
        type_count_cv = {}
        total_all = len(violation_statistics)
        type_count_all['total'] = (total_all, 100.00)  # Count and percentage
        for i in uniq_type:
            type_count_all[i] = (
                rest_type_all.count(i), round(((float(rest_type_all.count(i)) / float(total_all)) * 100.00), 2))
            type_count_v[i] = (
                rest_type_v.count(i),
                round(((float(rest_type_v.count(i)) / float(rest_type_all.count(i))) * 100.00), 2))
            type_count_cv[i] = (
                rest_type_cv.count(i),
                round(((float(rest_type_cv.count(i)) / float(rest_type_all.count(i))) * 100.00), 2))
        return type_count_all, type_count_cv, type_count_v

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

    @staticmethod
    def get_restraints(star_file):
        """
        Extracts restraints from NMR-STAR file
        :param star_file: NMR-STAR file
        :return: distance restraints in dictionary format, angle restraints in dictionary format, chain map as dictionary
        """
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

        try:
            entity_assembly = star_data.get_loops_by_category('_Entity_assembly')[0]
            col_names = entity_assembly.get_tag_names()
            id_index = col_names.index('_Entity_assembly.ID')
            asym_index = col_names.index('_Entity_assembly.Asym_ID')
            chain_dict = {}
            for row in entity_assembly:
                chain_map[row[id_index]] = row[asym_index]
        except IndexError:
            print("Entity_assembly loop not found, No chain mapping created")
            chain_dict = {}
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
                        angle_dict2[(rest[list_id], rest[rest_id])].append(
                            [atom1, atom2, atom3, atom4, rest_name, lb, ub])
                if lp_flg:
                    angle_dict['angle_restraints'] = r_dict
                else:
                    angle_dict[sf.name] = r_dict
        else:
            angle_dict = None
            angle_dict2 = None
        return dist_dict2, angle_dict2, chain_dict

    @staticmethod
    def get_distance(c1, c2):
        """
        Calculates the distance between two coordinate points
        :param c1: array of x,y,z
        :param c2: array of x,y,z
        :return: distance between two ponts
        """
        return numpy.linalg.norm(c1 - c2)

    @staticmethod
    def get_dihedral_angle(c1, c2, c3, c4):
        """
        Calculates the dihedral angle from the given set of four coordinate values
        :param c1: array of x,y,z
        :param c2: array of x,y,z
        :param c3: array of x,y,z
        :param c4: array of x,y,z
        :return: angle in degrees
        """
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
        """
        Calculates 1/r^6 sum for ambiguous restraints as recommended by NMR VTF
        :param dist_list: list of distances
        :return: r6 sum
        """
        return (sum([i ** (-6.) for i in dist_list])) ** (-1. / 6.)

    def calculate_distance_violations(self, coordinates, restraints):
        """
        Calculates violation for each restraint
        :param coordinates:  output from get_coordinates
        :param restraints:  output from get_restraints
        :return: dictionary { rest identifier : { model no : (value, 'type') }} example {('1', '1210'):
        {1: (0.22874009681108554, 'long'), 2: (0.084672752498432757, 'long')...}}
        """
        violations = {}
        for rest_id in restraints.keys():
            m = {}
            for model in coordinates.keys():
                dist_list = []
                for rest in restraints[rest_id]:
                    atom_1 = rest[0]
                    atom_2 = rest[1]
                    cat = rest[2]
                    lb = rest[3]
                    ub = rest[4]
                    pos_1 = coordinates[model][atom_1]
                    pos_2 = coordinates[model][atom_2]
                    d = self.get_distance(pos_1, pos_2)
                    dist_list.append(d)
                r6dist = self.r6sum(dist_list)
                if lb <= r6dist <= ub:
                    err = 0.0
                elif r6dist < lb:
                    err = abs(r6dist - lb)
                else:
                    err = abs(r6dist - ub)
                m[model] = (err, cat)
            violations[rest_id] = m
        return violations

    def calculate_angle_violations(self, coordinates, restraints):
        """
        Calculates violation for each restraint
        :param coordinates: output from get_coordinates
        :param restraints:  output from get_restraints
        :return: dictionary { rest identifier : { model no : (value, 'type') }} example {('1', '121'):
        {1: (0.22874009681108554, 'PSI'), 2: (0.084672752498432757, 'PHI')...}}
        """
        violations = {}
        for rest_id in restraints.keys():
            m = {}
            for model in coordinates.keys():
                ang_list = []
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
                    ang = self.get_dihedral_angle(pos_1, pos_2, pos_3, pos_4)
                    ang_list.append(ang)
                avg_viol = numpy.mean(ang_list)
                if lb <= avg_viol <= ub:
                    err = 0.0
                elif avg_viol < lb:
                    err = abs(avg_viol - lb)
                else:
                    err = abs(avg_viol - ub)
                m[model] = (err, cat)
            violations[rest_id] = m
        return violations

    @staticmethod
    def write_xml(distance_violations, angle_violations, distance, angle, atom_ids):
        rest_ids = list(distance_violations.keys())
        model_ids = list(distance_violations[rest_ids[0]].keys())
        violations = ET.Element('RestraintsViolations')
        dist_viol = ET.SubElement(violations, 'DistanceViolations')
        for m_id in model_ids[:50]:
            models = ET.SubElement(dist_viol, 'Model')
            models.set('model', str(m_id))
            for r_id in rest_ids[:50]:
                for rest in distance[r_id]:
                    if distance_violations[r_id][m_id][0] > 0.0:
                        model = ET.SubElement(models, 'Violation')
                        model.set('rest_id', str(r_id[1]))
                        model.set('rest_list_id', str(r_id[0]))
                        model.set('model_1', str(m_id))
                        model.set('chain_1', atom_ids[m_id][rest[0]][7])
                        model.set('resnum_1', atom_ids[m_id][rest[0]][4])
                        model.set('ent_1', atom_ids[m_id][rest[0]][0])
                        model.set('altcode_1', atom_ids[m_id][rest[0]][5])
                        model.set('icode_1', atom_ids[m_id][rest[0]][6])
                        model.set('seq_1', atom_ids[m_id][rest[0]][3])
                        model.set('said_1', atom_ids[m_id][rest[0]][1])
                        model.set('resname_1', atom_ids[m_id][rest[0]][2])
                        model.set('model_2', str(m_id))
                        model.set('chain_2', atom_ids[m_id][rest[1]][7])
                        model.set('resnum_2', atom_ids[m_id][rest[1]][4])
                        model.set('ent_2', atom_ids[m_id][rest[1]][0])
                        model.set('altcode_2', atom_ids[m_id][rest[1]][5])
                        model.set('icode_2', atom_ids[m_id][rest[1]][6])
                        model.set('seq_2', atom_ids[m_id][rest[1]][3])
                        model.set('said_2', atom_ids[m_id][rest[1]][1])
                        model.set('resname_2', atom_ids[m_id][rest[1]][2])

                        # model.set('seq_1', str(rest[0][0]))
                        # model.set('chain_1', str(rest[0][1]))
                        # model.set('res_1', str(rest[0][2]))
                        # model.set('atom_1', str(rest[0][3]))
                        # model.set('seq_2', str(rest[1][0]))
                        # model.set('chain_2', str(rest[1][1]))
                        # model.set('res_2', str(rest[1][2]))
                        # model.set('atom_2', str(rest[1][3]))3
                        model.set('violation', str(distance_violations[r_id][m_id][0]))
                        model.set('rest_type', distance_violations[r_id][m_id][1])
        rest_ids = list(angle_violations.keys())
        ang_viol = ET.SubElement(violations, 'AngleViolations')
        for m_id in model_ids[:50]:
            models = ET.SubElement(ang_viol, 'Model')
            models.set('model', str(m_id))
            for r_id in rest_ids[:50]:
                for rest in angle[r_id]:
                    if angle_violations[r_id][m_id][0] > 0.0:
                        model = ET.SubElement(models, 'Violation')
                        model.set('rest_id', str(r_id[1]))
                        model.set('rest_list_id', str(r_id[0]))
                        model.set('model_1', str(m_id))
                        model.set('chain_1', atom_ids[m_id][rest[0]][7])
                        model.set('resnum_1', atom_ids[m_id][rest[0]][4])
                        model.set('ent_1', atom_ids[m_id][rest[0]][0])
                        model.set('altcode_1', atom_ids[m_id][rest[0]][5])
                        model.set('icode_1', atom_ids[m_id][rest[0]][6])
                        model.set('seq_1', atom_ids[m_id][rest[0]][3])
                        model.set('said_1', atom_ids[m_id][rest[0]][1])
                        model.set('resname_1', atom_ids[m_id][rest[0]][2])
                        model.set('model_2', str(m_id))
                        model.set('chain_2', atom_ids[m_id][rest[1]][7])
                        model.set('resnum_2', atom_ids[m_id][rest[1]][4])
                        model.set('ent_2', atom_ids[m_id][rest[1]][0])
                        model.set('altcode_2', atom_ids[m_id][rest[1]][5])
                        model.set('icode_2', atom_ids[m_id][rest[1]][6])
                        model.set('seq_2', atom_ids[m_id][rest[1]][3])
                        model.set('said_2', atom_ids[m_id][rest[1]][1])
                        model.set('resname_2', atom_ids[m_id][rest[1]][2])
                        model.set('model_3', str(m_id))
                        model.set('chain_3', atom_ids[m_id][rest[2]][7])
                        model.set('resnum_3', atom_ids[m_id][rest[2]][4])
                        model.set('ent_3', atom_ids[m_id][rest[2]][0])
                        model.set('altcode_3', atom_ids[m_id][rest[2]][5])
                        model.set('icode_3', atom_ids[m_id][rest[2]][6])
                        model.set('seq_3', atom_ids[m_id][rest[2]][3])
                        model.set('said_3', atom_ids[m_id][rest[2]][1])
                        model.set('resname_4', atom_ids[m_id][rest[2]][2])
                        model.set('model_4', str(m_id))
                        model.set('chain_4', atom_ids[m_id][rest[3]][7])
                        model.set('resnum_4', atom_ids[m_id][rest[3]][4])
                        model.set('ent_4', atom_ids[m_id][rest[3]][0])
                        model.set('altcode_4', atom_ids[m_id][rest[3]][5])
                        model.set('icode_4', atom_ids[m_id][rest[3]][6])
                        model.set('seq_4', atom_ids[m_id][rest[3]][3])
                        model.set('said_4', atom_ids[m_id][rest[3]][1])
                        model.set('resname_4', atom_ids[m_id][rest[3]][2])

                        # model.set('seq_1', str(rest[0][0]))
                        # model.set('chain_1', str(rest[0][1]))
                        # model.set('res_1', str(rest[0][2]))
                        # model.set('atom_1', str(rest[0][3]))
                        # model.set('seq_2', str(rest[1][0]))
                        # model.set('chain_2', str(rest[1][1]))
                        # model.set('res_2', str(rest[1][2]))
                        # model.set('atom_2', str(rest[1][3]))
                        # model.set('seq_3', str(rest[2][0]))
                        # model.set('chain_3', str(rest[2][1]))
                        # model.set('res_3', str(rest[2][2]))
                        # model.set('atom_3', str(rest[2][3]))
                        # model.set('seq_4', str(rest[3][0]))
                        # model.set('chain_4', str(rest[3][1]))
                        # model.set('res_4', str(rest[3][2]))
                        # model.set('atom_4', str(rest[3][3]))
                        model.set('violation', str(angle_violations[r_id][m_id][0]))
                        model.set('rest_type', angle_violations[r_id][m_id][1])
        ET.tostring(violations)
        t = ET.ElementTree(violations)
        t.write('full.xml')


    @staticmethod
    def calculate_violation_statistics(violations):
        """
        Calculates average violation value for each restraint for an ensemble
        :param violations: output from calculate_distance_violations or calculate_angle_violations
        :return: list of average violation for each restraint, list of violations
        example output [0.050748176870277231, 0.00147115297574274, 0.11772253375387987, 7, [2, 5, 7, 11, 14, 18, 19],
         'medium', ('1', '1')],[0.041170250870272262, 2, 'medium', ('1', '1')]
        """
        rest_list = list(violations.keys())
        models = list(violations[rest_list[0]].keys())
        avg_violations = {}
        viol = []
        for rest in rest_list:
            v = []
            m_id = []
            for m in models:
                if violations[rest][m][0] > 0.0:
                    v.append(violations[rest][m][0])
                    viol.append([violations[rest][m][0], m, cat, rest])
                    m_id.append(m)
                cat = violations[rest][m][1]
            if len(v) > 0:
                avg_violations[rest] = [numpy.mean(v), min(v), max(v), len(v), m_id, cat]
            else:
                avg_violations[rest] = [0.0, 0.0, 0.0, 0, m_id, cat]
        avg_viol_list = []
        for rest in rest_list:
            avg_viol_list.append(avg_violations[rest])
            avg_viol_list[-1].append(rest)
        return avg_viol_list, viol

    @staticmethod
    def bin_distance_violations(violations):
        """
        Count the number of violations in different bins (not violated)(0-0.2),(0.2-0.5),(0.5-1.0),(1.0-2.0),(2.0-5.0),(>5)
        :param violations: output from calculate_distance_violations
        :return: disctionary {model id : [count in each bin]} example {1: [1475, 19, 7, 12, 19, 7, 0],
         2: [1469, 30, 7, 9, 15, 9, 0], 3: [1473, 23, 7, 7, 20, 9, 0]...}
        """
        rest_list = list(violations.keys())
        models = list(violations[rest_list[0]].keys())
        stat = {}
        for m in models:
            c = [0, 0, 0, 0, 0, 0, 0]
            for rest in rest_list:
                if violations[rest][m][0] == 0.0:
                    c[0] += 1
                elif 0.0 < violations[rest][m][0] <= 0.2:
                    c[1] += 1
                elif 0.2 < violations[rest][m][0] <= 0.5:
                    c[2] += 1
                elif 0.5 < violations[rest][m][0] <= 1.0:
                    c[3] += 1
                elif 1.0 < violations[rest][m][0] <= 2.0:
                    c[4] += 1
                elif 2.0 < violations[rest][m][0] <= 5.0:
                    c[5] += 1
                elif 5.0 < violations[rest][m][0]:
                    c[6] += 1
                else:
                    print("Error in violation calculation")
            stat[m] = c
        return stat

    @staticmethod
    def bin_angle_violations(violations):
        """
        Count the number of violations in different bins (not violated)(0-5),(5-10),(10-20),(20-40),(40-80),(>80)
        :param violations: output from calculate_angle_violations
        :return: dictionary {model id : [count in each bin]} example {1: [1475, 19, 7, 12, 19, 7, 0],
         2: [1469, 30, 7, 9, 15, 9, 0], 3: [1473, 23, 7, 7, 20, 9, 0]...}
        """
        rest_list = list(violations.keys())
        models = list(violations[rest_list[0]].keys())
        stat = {}
        for m in models:
            c = [0, 0, 0, 0, 0, 0, 0]
            for rest in rest_list:
                if violations[rest][m][0] == 0.0:
                    c[0] += 1
                elif 0.0 < violations[rest][m][0] <= 5.0:
                    c[1] += 1
                elif 5.0 < violations[rest][m][0] <= 10.0:
                    c[2] += 1
                elif 10.0 < violations[rest][m][0] <= 20.0:
                    c[3] += 1
                elif 20.0 < violations[rest][m][0] <= 40.0:
                    c[4] += 1
                elif 40.0 < violations[rest][m][0] <= 80.0:
                    c[5] += 1
                elif 80.0 < violations[rest][m][0]:
                    c[6] += 1
                else:
                    print("Error in violation calculation")
            stat[m] = c
        return stat


if __name__ == "__main__":
    # p = ValidateRestraints('vtf_examples/CtR107.cif', 'vtf_examples/CtR107.nef')
    p = Validate_NMR_Restraints()
    p.run_validation('nef_examples/2m2e.cif', 'nef_examples/2m2e.str')
