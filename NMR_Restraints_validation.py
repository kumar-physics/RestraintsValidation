import sys
import os.path
from numpy import logical_and
#print sys.path
sys.path.append("/home/kumaran/git/py-mmcif")
from mmcif.io.PdbxReader import PdbxReader
from math import sqrt,acos
import numpy
import pynmrstar
from monkeypatch import patch_parser
patch_parser(pynmrstar)
import logging
import ntpath
import operator
import json


class ValidateRestraints(object):
    
    """ NMR Restraint validation module """
    __version__ = "v1.2-17-gef22431"
    
    
    
    def __init__(self,cfile,rfile):
        
        """ Initialize the module with a coordinate( CIF format) and restraint ( NMR-STAR format) files.
        A log file is generated based on input coordinate file name and path""" 
        
        self.cFile = os.path.realpath(cfile)
        self.rFile = os.path.realpath(rfile)
        (filePath,fileName)=ntpath.split(os.path.realpath(self.cFile))
        self.logFile =  filePath+"/"+fileName.split(".")[0]+".log"
        self.outFile =  filePath+"/"+fileName.split(".")[0]+".out"
        self.outFile2 =  filePath+"/"+fileName.split(".")[0]+"_short.out"
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s\t%(levelname)s\t%(message)s')
        self.fh = logging.FileHandler(self.logFile,'w')
        self.fh.setFormatter(formatter)
        self.logger.addHandler(self.fh)
        #logging.basicConfig(filename = self.logFile,level = logging.INFO,
         #                   format='%(asctime)s\t%(levelname)s\t%(message)s',filemode = 'w')
        self.validate()
        self.fh.close()
        self.logger.removeHandler(self.fh)
    @staticmethod   
    def r6sum(dist_list):
        return (sum([i**(-6.) for i in dist_list]))**(-1./6.)
    @staticmethod   
    def r6average(dist_list):
        return (sum([i**(-6.) for i in dist_list])/float(len(dist_list)))**(-1./6.)
    def validate(self):
        
        """ Integrates all steps in a validation process"""
        
        
        if self.checkFiles():
            self.logger.info('Both input files exist')
        else:
            self.logger.critical('Problem with input files! check input files exist ')
            raise IOError("Input files not found")
        self.readCIFFile()
        rest_info= self.readSTARFile()
        if self.rest_flag[0]:
            self.getDistanceRestraints()
            self.distRestraintStat()
            self.calDisanceFromMoels()
            self.distRestAnalysis()
#             for i in range(len(self.model_vs_dist)):
#                 self.logger.info("Violations in bins for Model {} : {}.".format(i+1,self.model_vs_dist[i]))
#             for j in self.dist_vs_model:
#                 for i in j:
#                     self.logger.info('Rest_ID {} violated in {} models;{}'.format(i[0],len(i)-1,i[1:]))
            self.sortDistViolations()
            self.sortAvgDistViolations()
        if self.rest_flag[1]:
            self.getAngleRestraints()
            self.angRestraintStat()
            self.calAngleFromMoels()
            self.angleRestAnalysis()
#             for i in range(len(self.model_vs_ang)):
#                 self.logger.info("Violations in bins for Model {} : {}.".format(i+1,self.model_vs_ang[i]))
#             for j in self.ang_vs_model:
#                 for i in j:
#                     self.logger.info('Rest_ID {} violated in {} models;{}'.format(i[0],len(i)-1,i[1:]))
            self.sortAngleViolations()
            self.sortAvgAngViolations()
        #self.generateReport()
        #self.generateReport2()
        #out=json.loads(self.generateJson())
        #return out
       
    
        
        #self.getAngleRestraints()
        #self.checkAngleRestraints()
        #stat= self.checkDistanceRestraints()
    def checkFiles(self):
        """ Check the existence of input files"""
        
        if os.path.isfile(self.cFile):
            self.logger.info('Coordinate file found {}'.format(self.cFile))
        else:
            self.logger.critical('Cooriantefile NOT found')
            return False
        if os.path.isfile(self.rFile):
            self.logger.info('Restraint file found {}'.format(self.rFile))
        else:
            self.logger.critical('Restraint file NOT found')
            return False
        return True   
    
    def readCIFFile(self):
        """ Reads the input coordinate CIF file using pdbx lightweight parser ref: http://mmcif.wwpdb.org/docs/sw-examples/python/html/index.html """ 
        self.logger.info('Reading {}'.format(self.cFile))
        self.logger.info('Reading coordinate file using {} from {}'.format(PdbxReader.__name__,PdbxReader.__module__))
        self.cifData = []
        ifh = open(self.cFile,'r')
        pRd = PdbxReader(ifh)
        pRd.read(self.cifData)
        ifh.close()
        c0 = self.cifData[0]
        atom_site = c0.getObj('atom_site')
        self.maxModels = int(atom_site.getValue('pdbx_PDB_model_num',-1))
        if self.maxModels == 1:
            self.logger.warn('Coordinate file has only one model')
        elif self.maxModels == 0:
            self.logger.error('Coordinate file has zero models')
        else:
            self.logger.info('Coordinate file has {} models'.format(self.maxModels))    
    def readSTARFile(self):
        
        """ Reads the input NMR-STAR file using pynmrstar library ref : https://github.com/uwbmrb/PyNMRSTAR """
        
        self.logger.info('Reading {}'.format(self.rFile))
        self.logger.info('Reading restraint file using pynmrstar verion {}'.format(pynmrstar.__version__))
        self.starData = pynmrstar.Entry.from_file(self.rFile)  
        cat_list = [saveframe.category for saveframe in self.starData]
        rest_info = [False,False]
        if 'general_distance_constraints' in cat_list:
            self.logger.info('general_distance_constraints saveframe found')
            rest_info[0] = True
        else:
            self.logger.warning('general_distance_constraints saveframe not found')
        if 'torsion_angle_constraints' in cat_list:
            self.logger.info('torsion_angle_constraints saveframe found')
            rest_info[1] = True
        else:
            self.logger.warning('torsion_angle_constraints saveframe not found')
        
        if 'general_distance_constraints' not in cat_list and 'torsion_angle_constraints' not in cat_list:
            self.logger.error('Both general_distance_constraints and torsion_angle_constraints are missing in the STAR file')
            rest_info=[False,False]
        
        self.rest_flag = rest_info
    
    def getSeqlen(self):
        seq = self.starData.get_tag('_Chem_comp_assembly.Comp_index_ID')
        return len(seq)
        
    def checkAngleRestraints(self):
        for modelId in range(1,self.maxModels+1):
            co = self.getCoordinates(modelId)
            for i in self.ang_rest:
                ang = self.getDihedralAngle(co[i[0]], co[i[1]], co[i[2]], co[i[3]]) 
                if ang>=self.ang_rest[i][0] and ang <= self.ang_rest[i][1]:
                    self.logger.info( "Model {} OK".format(modelId))
                else:
                    if ang<self.ang_rest[i][0]:
                        ang_err = self.ang_rest[i][0]-ang
                    elif ang>self.ang_rest[i][1]:
                        ang_err = ang - self.ang_rest[i][1]
                    else:
                        ang_err = -9999    
                    print( "Model {} NOTOK {} {} {} {}".format(modelId,ang_err,i,self.ang_rest[i][0],self.ang_rest[i][1]))
    
    
        
        
    
        
    def getCoordinates(self,modelID):
        
        """Creates a dictionary of coordinates gor a given modelID.
        co[atom] = array([x,y,z])
        atoms are identified by (label_seq_id,label_entity_id,label_comp_id,lable_atom_id)"""
        
        c0 = self.cifData[0]
        atom_site = c0.getObj('atom_site')
        colnames = atom_site.getAttributeList()
        modelid = colnames.index('pdbx_PDB_model_num')
        self.maxModels = int(atom_site.getValue('pdbx_PDB_model_num',-1))
        xid = colnames.index('Cartn_x')
        yid = colnames.index('Cartn_y')
        zid = colnames.index('Cartn_z')
        atomid = colnames.index('label_atom_id')
        compid = colnames.index('label_comp_id')
        asymid = colnames.index('label_asym_id')
        entityid = colnames.index('label_entity_id') # asymid is used instead of entity id because of different definition of cif and star
        seqid = colnames.index('label_seq_id')
        co = {}
        for dat in atom_site.getRowList():
            if int(dat[modelid]) == modelID:
                co[(dat[seqid],dat[asymid],dat[compid],dat[atomid])] = numpy.array([float(dat[xid]),float(dat[yid]),float(dat[zid])])
            
        return co
    
    @staticmethod
    def getDihedralAngle(c1,c2,c3,c4):
        
        """ Calculates the dihedral angle from the given four coordinate values.
        Each coordinate is an array of x,y,z. Returns angle in degrees"""
        
        bv12 = c1-c2
        bv32 = c3-c2
        bv43 = c4-c3
        pv13 = numpy.cross(bv12, bv32)
        pv24 = numpy.cross(bv43, bv32)
        pro = numpy.dot(pv13,pv24)
        sqdist13 = numpy.dot(pv13,pv13)
        sqdist24 = numpy.dot(pv24,pv24)
        cosin = pro / sqrt(sqdist13*sqdist24)
        cosin - min(1.0, max(-1.0, cosin))
        angle = acos(cosin)
        
        if numpy.dot(pv13, numpy.cross(pv24,bv32)) < 0:
            angle = -angle
        return round(numpy.degrees(angle),4)
    @staticmethod
    def getDistance(c1,c2):
        
        """ Calculates the distance between two coordinate values.
        Each coordinate is an array of x,y,z. Returns distance in A assuming the input coordinates are in A""" 
        
        return numpy.linalg.norm(c1-c2)
          #return sqrt(((c1[0]-c2[0])**2)+((c1[1]-c2[1])**2)+((c1[2]-c2[2])**2))  
        
    
    

       
    def getDistanceRestraints(self):
        
        """Creates a dictionary of distance restraints from _Gen_dist_constraint loop.
        dist_rest[atom 1,atom 2] = [lower bound, upper bound]
        atoms are identified by (Comp_index_ID,Entity_assembly_ID,Comp_ID,Atom_ID)"""
        
        dist = self.starData.get_loops_by_category('_Gen_dist_constraint') # list of distance restraint loops
        self.dist_rest = []
        for dl in dist: # for loop in the list
            self.dist_rest.append({}) # each loop info goes into a dictionary 
            rest_id = dl.get_tag_names().index("_Gen_dist_constraint.ID") # rest_id is the key 
            atomid1 = dl.get_tag_names().index("_Gen_dist_constraint.Atom_ID_1") # rest of them are values 
            atomid2 = dl.get_tag_names().index("_Gen_dist_constraint.Atom_ID_2")
            compid1 = dl.get_tag_names().index("_Gen_dist_constraint.Comp_ID_1")
            compid2 = dl.get_tag_names().index("_Gen_dist_constraint.Comp_ID_2")
            #entityid1 = dl.get_tag_names().index("_Gen_dist_constraint.Entity_assembly_ID_1")
            #entityid2 = dl.get_tag_names().index("_Gen_dist_constraint.Entity_assembly_ID_2")
            entityid1 = dl.get_tag_names().index("_Gen_dist_constraint.Auth_asym_ID_1")# asymid is used instead of entity id because of different definition of cif and star
            entityid2 = dl.get_tag_names().index("_Gen_dist_constraint.Auth_asym_ID_2")
            seqid1 = dl.get_tag_names().index("_Gen_dist_constraint.Comp_index_ID_1") 
            seqid2 = dl.get_tag_names().index("_Gen_dist_constraint.Comp_index_ID_2")
            lb = dl.get_tag_names().index("_Gen_dist_constraint.Distance_lower_bound_val")
            ub = dl.get_tag_names().index("_Gen_dist_constraint.Distance_upper_bound_val")
            for i in dl: # for every row in the restraint loop
                if int(i[rest_id]) not in self.dist_rest[-1].keys(): # each rest_id may have multiple rows, because of ambiguity codes and NEF atom nomeclature ; so it will be a list
                    self.dist_rest[-1][int(i[rest_id])]=[]
                try:
                    lbv = float(i[lb])
                except ValueError:
                    lbv = -999
                try:
                    ubv = float(i[ub])
                except ValueError:
                    ubv = 999
                if lbv == -999 and ubv == 999:
                    self.logger.error("Distance restraint value not readable {},{} between ({}{}{}{})-({}{}{}{})".format(i[lb],i[ub],i[seqid1],i[entityid1],i[compid1],i[atomid1],i[seqid2],i[entityid2],i[compid2],i[atomid2]))
                else:
                    self.dist_rest[-1][int(i[rest_id])].append([(i[seqid1],i[entityid1],i[compid1],i[atomid1]),(i[seqid2],i[entityid2],i[compid2],i[atomid2]),(lbv,ubv)])
                #self.dist_rest[(i[seqid1],i[entityid1],i[compid1],i[atomid1]),(i[seqid2],i[entityid2],i[compid2],i[atomid2])]= [lbv,ubv]
        self.logger.info('Number of distance restraints : {}'.format([len(i) for i in self.dist_rest]))
        
    def sortDistViolations(self):
        r_list=[]
        for j in self.dist_vs_model:
                for i in j:
                    if len(i)>1:
                        for k in range(1,len(i)):
                            r_list.append([i[0],i[k][0],i[k][1]])
        self.sorted_dist_rest=sorted(r_list, key = operator.itemgetter(2),reverse = True)
    def sortAvgDistViolations(self):
        r_list=[]
        for j in self.dist_vs_model:
                for i in j:
                    if len(i)>1:
                        rid = i[0]
                        err=[x[1] for x in i[1:]]
                        r_list.append([rid,len(err),sum(err)/len(err)])
                       
        self.sorted_avg_dist_rest=sorted(r_list, key = operator.itemgetter(2),reverse = True)
    def sortAvgAngViolations(self):
        r_list=[]
        for j in self.ang_vs_model:
                for i in j:
                    if len(i)>1:
                        rid = i[0]
                        err=[x[1] for x in i[1:]]
                        r_list.append([rid,len(err),sum(err)/len(err)])
                       
        self.sorted_avg_ang_rest=sorted(r_list, key = operator.itemgetter(2),reverse = True)
        
    def sortAngleViolations(self):
        r_list=[]
        
        for j in self.ang_vs_model:
                for i in j:
                    if len(i)>1:
                        for k in range(1,len(i)):
                            r_list.append([i[0],i[k][0],i[k][1]])
        self.sorted_ang_rest=sorted(r_list, key = operator.itemgetter(2),reverse = True) 
    
    def checkDistRestCateory(self,restid,rlistid):
        cat = -1
        if (restid,rlistid) in self.dist_rest_info[1]:
            cat = 0
        elif (restid,rlistid) in self.dist_rest_info[2]:
            cat = 1
        elif (restid,rlistid) in self.dist_rest_info[3]:
            cat = 2
        elif (restid,rlistid) in self.dist_rest_info[4]:
            cat = 3
        else:
            self.logger.warning("Rest id {} doesn't belong to any category".format((restid,rlistid)))
            cat = -1
        return cat    
    
    def generateReport2(self):
        fo=open(self.outFile2,'w')
        self.logger.info('Generating output file {}'.format(self.outFile2))
        if self.rest_flag[0]:
            fo.write('1. Conformationally restricting restraints\n')
            fo.write('{0:-<45}\n\n'.format(''))
            fo.write('1.1 Distance restraints \n')
            fo.write('{0:-<25}\n\n'.format(''))
            if len(self.dist_vs_model)>1:
                hbonds=min([len(i) for i in self.dist_vs_model])
            else:
                hbonds=0
            self.model_vs_dist2={}
            for i in range(len(self.model_vs_dist)):
                self.model_vs_dist2[i]=[self.model_vs_dist[i][0],self.model_vs_dist[i][1],self.model_vs_dist[i][2],sum(self.model_vs_dist[i][3:4]),self.model_vs_dist[i][7]]
            fo.write('\t{}={} ({}%)\n'.format('{0:<62}'.format('Total '),'{0: >6}'.format(len(self.dist_rest_info[0])),round((float(len(self.dist_rest_info[0]))/float(len(self.dist_rest_info[0])))*100.0,2)))
            fo.write('\t{}={} ({}%)\n'.format('{0:<62}'.format('Intraresidue  |i-j|=0'),'{0: >6}'.format(len(self.dist_rest_info[1])),round((float(len(self.dist_rest_info[1]))/float(len(self.dist_rest_info[0])))*100.0,2)))
            fo.write('\t{}={} ({}%)\n'.format('{0:<62}'.format('Sequential |i-j|=1'),'{0: >6}'.format(len(self.dist_rest_info[2])),round((float(len(self.dist_rest_info[2]))/float(len(self.dist_rest_info[0])))*100.0,2)))
            fo.write('\t{}={} ({}%)\n'.format('{0:<62}'.format('Medium range  |i-j|>1 and |i-j|<5'),'{0: >6}'.format(len(self.dist_rest_info[3])),round((float(len(self.dist_rest_info[3]))/float(len(self.dist_rest_info[0])))*100.0,2)))
            fo.write('\t{}={} ({}%)\n'.format('{0:<62}'.format('Long range |i-j|>=5'),'{0: >6}'.format(len(self.dist_rest_info[4])),round((float(len(self.dist_rest_info[4]))/float(len(self.dist_rest_info[0])))*100.0,2)))
            fo.write('\t{}={} ({}%)\n'.format('{0:<62}'.format('Hydrogen-bond'),'{0: >6}'.format(hbonds),round((float(hbonds)/float(len(self.dist_rest_info[0])))*100.0,2)))
            seq_len=self.getSeqlen()
            if self.rest_flag[1]:
                tot_rest = 0
                for k in self.ang_rest_info.keys():
                    tot_rest+=self.ang_rest_info[k]
                fo.write('\t{}={}\n'.format('{0:<62}'.format('Dihedral angle restraints'),'{0: >6}'.format(tot_rest)))
                
                self.model_vs_ang2={}
                for i in range(len(self.model_vs_ang)):
                    self.model_vs_ang2[i]=[self.model_vs_ang[i][0],sum(self.model_vs_ang[i][1:2]),sum(self.model_vs_ang[i][2:4]),self.model_vs_ang[i][7]]
                ang_per_res = [[float(j)/float(seq_len) for j in self.model_vs_ang2[i]] for i in range(len(self.model_vs_ang2))]
            
            rest_per_res = [[float(j)/float(seq_len) for j in self.model_vs_dist2[i]] for i in range(len(self.model_vs_dist2))]
            fo.write('\t{}={}\n'.format('{0:<62}'.format('No of restraints per residue'),'{0: >6}'.format(round(float(sum([i[0] for i in rest_per_res]))/float(len(rest_per_res)),2))))
            fo.write('\t{}={}\n'.format('{0:<62}'.format('No of long-range restraints per residue'),'{0: >6}'.format(round(float(len(self.dist_rest_info[4]))/float(seq_len),2))))
            fo.write('\n\n')
            fo.write('1.2 Residual restraint violations \n')
            fo.write('{0:-<34}\n\n'.format(''))
            fo.write('\t1.2.1 Average no. of distance violations per structure \n')
            fo.write('\t{0:-<55}\n'.format(''))
            max_val = self.sorted_dist_rest[0][2]
            if max_val<0.2:
                fo.write('\t{}={}(max={})\n'.format('{0:<62}'.format('0.1-0.2A'),'{0: >6}'.format(round(float(sum([i[1] for i in rest_per_res]))/float(len(rest_per_res)),2)),round(max_val,2)))
                fo.write('\t{}={}\n'.format('{0:<62}'.format('0.2-0.5A'),'{0: >6}'.format(round(float(sum([i[2] for i in rest_per_res]))/float(len(rest_per_res)),2))))
                fo.write('\t{}={}\n'.format('{0:<62}'.format('>0.5A'),'{0: >6}'.format(round(float(sum([i[3] for i in rest_per_res]))/float(len(rest_per_res)),2))))
            elif max_val >= 0.2 and max_val < 0.5:
                fo.write('\t{}={}\n'.format('{0:<62}'.format('0.1-0.2A'),'{0: >6}'.format(round(float(sum([i[1] for i in rest_per_res]))/float(len(rest_per_res)),2))))
                fo.write('\t{}={}(max={})\n'.format('{0:<62}'.format('0.2-0.5A'),'{0: >6}'.format(round(float(sum([i[2] for i in rest_per_res]))/float(len(rest_per_res)),2)),round(max_val,2)))
                fo.write('\t{}={}\n'.format('{0:<62}'.format('>0.5A'),'{0: >6}'.format(round(float(sum([i[3] for i in rest_per_res]))/float(len(rest_per_res)),2))))
            else:
                fo.write('\t{}={}\n'.format('{0:<62}'.format('0.1-0.2A'),'{0: >6}'.format(round(float(sum([i[1] for i in rest_per_res]))/float(len(rest_per_res)),2))))
                fo.write('\t{}={}\n'.format('{0:<62}'.format('0.2-0.5A'),'{0: >6}'.format(round(float(sum([i[2] for i in rest_per_res]))/float(len(rest_per_res)),2))))
                fo.write('\t{}={}(max={})\n'.format('{0:<62}'.format('>0.5A'),'{0: >6}'.format(round(float(sum([i[3] for i in rest_per_res]))/float(len(rest_per_res)),2)),round(max_val,2)))
            max_val_ang = self.sorted_ang_rest[0][2]
            fo.write('\t1.2.2 Average no. of dihedral angle violations per structure \n')
            fo.write('\t{0:-<60}\n'.format(''))
            if max_val_ang <=10:
                fo.write('\t{}={}(max={})\n'.format('{0:<62}'.format('0-10'),'{0: >6}'.format(round(float(sum([i[1] for i in ang_per_res]))/float(len(ang_per_res)),2)),round(max_val_ang,2)))
                fo.write('\t{}={}\n'.format('{0:<62}'.format('>10'),'{0: >6}'.format(round(float(sum([i[2] for i in ang_per_res]))/float(len(ang_per_res)),2))))
            else:
                fo.write('\t{}={}\n'.format('{0:<62}'.format('0-10'),'{0: >6}'.format(round(float(sum([i[1] for i in ang_per_res]))/float(len(ang_per_res)),2))))
                fo.write('\t{}={}(max={})\n'.format('{0:<62}'.format('>10'),'{0: >6}'.format(round(float(sum([i[2] for i in ang_per_res]))/float(len(ang_per_res)),2)),round(max_val_ang,2)))
                
            #print round(float(sum([i[3] for i in rest_per_res]))/float(len(rest_per_res)),2)
        fo.close()
        
    def generateJson(self):
       
        if self.rest_flag[0]:
            dist_summary={}
            dist_summary['Total number of restraints']=(len(self.dist_rest_info[0]),round((float(len(self.dist_rest_info[0]))/float(len(self.dist_rest_info[0])))*100.0,2))
            dist_summary['Intraresidue restraints |i-j|=0']=(len(self.dist_rest_info[1]),round((float(len(self.dist_rest_info[1]))/float(len(self.dist_rest_info[0])))*100.0,2))
            dist_summary['Sequential restraints |i-j|=1']=(len(self.dist_rest_info[2]),round((float(len(self.dist_rest_info[2]))/float(len(self.dist_rest_info[0])))*100.0,2))
            dist_summary['Medium range  restraints|i-j|>1 and |i-j|<5']=(len(self.dist_rest_info[3]),round((float(len(self.dist_rest_info[3]))/float(len(self.dist_rest_info[0])))*100.0,2))
            dist_summary['Long range restraints|i-j|>=5']=(len(self.dist_rest_info[4]),round((float(len(self.dist_rest_info[4]))/float(len(self.dist_rest_info[0])))*100.0,2))
            
            dist_violation_each_model={}
            for i in range(len(self.model_vs_dist)):
                dist_violation_each_model[i+1]=[(self.model_vs_dist[i][7],round((float(self.model_vs_dist[i][7])/float(self.model_vs_dist[i][0]))*100.0,2)),
                                                (sum(self.model_vs_dist[i][1:-1]),round((float(sum(self.model_vs_dist[i][1:-1]))/float(self.model_vs_dist[i][0]))*100.0,2)),
                                                (self.model_vs_dist[i][1],round((float(self.model_vs_dist[i][1])/float(self.model_vs_dist[i][0]))*100.0,2)),
                                                (self.model_vs_dist[i][2],round((float(self.model_vs_dist[i][2])/float(self.model_vs_dist[i][0]))*100.0,2)),
                                                (self.model_vs_dist[i][3],round((float(self.model_vs_dist[i][3])/float(self.model_vs_dist[i][0]))*100.0,2)),
                                                (self.model_vs_dist[i][4],round((float(self.model_vs_dist[i][4])/float(self.model_vs_dist[i][0]))*100.0,2)),
                                                (self.model_vs_dist[i][5],round((float(self.model_vs_dist[i][5])/float(self.model_vs_dist[i][0]))*100.0,2)),
                                                (self.model_vs_dist[i][6],round((float(self.model_vs_dist[i][6])/float(self.model_vs_dist[i][0]))*100.0,2)),
                                                ]
            
            
            dist_violation_each_rest={}
            v=0
            nv=0
            ncv=0
            ncv_category=[0,0,0,0]
            v_category=[0,0,0,0]
            nv_category=[0,0,0,0]
            for rlistid in range(len(self.dist_vs_model)):
                j=self.dist_vs_model[rlistid]
                for i in j:
                    cat_n=self.checkDistRestCateory(i[0],rlistid)
                    if len(i)-1 > 0:
                        if len(i)-1 == self.maxModels:
                            ncv+=1
                            if cat_n>=0:
                                ncv_category[cat_n]+=1
                            else:
                                self.logger.warning("Rest id {} doesn't belong to any category".format(i[0]))
                        v+=1
                        if cat_n>=0:
                            v_category[cat_n]+=1
                        else:
                            self.logger.warning("Rest id {} doesn't belong to any category".format(i[0]))
                        dist_violation_each_rest[i[0]]=(len(i)-1,round(sum([float(k[1]) for k in i[1:]])/(len(i)-1),4),','.join([str(k[0]) for k in i[1:]]))
                        #fo.write('\t{}\t{}\t{}\t{}\n'.format('{0:<25}'.format(i[0]),
                        #                     '{0:<25}'.format(len(i)-1),
                         #                    '{0:<25}'.format(round(sum([float(k[1]) for k in i[1:]])/(len(i)-1),4)),
                         #                    '{0:<25}'.format(','.join([str(k[0]) for k in i[1:]]))))
                    else:
                        nv+=1
                        if cat_n>=0:
                            nv_category[cat_n]+=1
                        else:
                            self.logger.warning("Rest id {} doesn't belong to any category".format(i[0]))
           # ncv_p=[ (float(i)/float(sum(ncv_category)))*100.00,2) for i in ncv_category]
            #print ncv_category
            if sum(ncv_category) > 0 :
                ncv_p=['{}%'.format(round((float(i)/float(sum(ncv_category)))*100.00,2)) for i in ncv_category]
            else:
                ncv_p=['{}%'.format(round((float(i)/float(1))*100.00,2)) for i in ncv_category]
            v_p=['{}%'.format(round((float(i)/float(sum(v_category)))*100.00,2)) for i in v_category]
            nv_p=['{}%'.format(round((float(i)/float(sum(nv_category)))*100.00,2)) for i in nv_category]
            ncv_pt=['{}%'.format(round((float(i)/float(float(v+nv)))*100.00,2)) for i in ncv_category]
            v_pt=['{}%'.format(round((float(i)/float(float(v+nv)))*100.00,2)) for i in v_category]
            nv_pt=['{}%'.format(round((float(i)/float(float(v+nv)))*100.00,2)) for i in nv_category]
            dist_viol_summary={}
            dist_viol_summary['Number of consistently violated restraints (violated in all models)']=[(ncv_category[0],ncv_pt[0]),(ncv_category[1],ncv_pt[1]),(ncv_category[2],ncv_pt[2]),(ncv_category[3],ncv_pt[3]),(ncv,round((float(ncv)/float(v+nv))*100,2))]
            dist_viol_summary['Number of violated restraints (at least in one model)']=[(v_category[0],v_pt[0]),(v_category[1],v_pt[1]),(v_category[2],v_pt[2]),(v_category[3],v_pt[3]),(v,round((float(v)/float(v+nv))*100,2))]
            dist_viol_summary['Number of non violated restraints']=[(nv_category[0],nv_pt[0]),(nv_category[1],nv_pt[1]),(nv_category[2],nv_pt[2]),(nv_category[3],nv_pt[3]),(nv,round((float(nv)/float(v+nv))*100,2))]
            dist_viol_summary['Total number of restraints']=[(len(self.dist_rest_info[1]),round((float(len(self.dist_rest_info[1]))/float(len(self.dist_rest_info[0])))*100.0,2)),
                                                             (len(self.dist_rest_info[2]),round((float(len(self.dist_rest_info[2]))/float(len(self.dist_rest_info[0])))*100.0,2)),
                                                             (len(self.dist_rest_info[3]),round((float(len(self.dist_rest_info[3]))/float(len(self.dist_rest_info[0])))*100.0,2)),
                                                             (len(self.dist_rest_info[4]),round((float(len(self.dist_rest_info[4]))/float(len(self.dist_rest_info[0])))*100.0,2)),
                                                             v+nv,round((float(v+nv)/float(v+nv))*100,2)]
            sorted_dist_violations=[(i[0],i[1],round(i[2],4)) for i in self.sorted_dist_rest]
            sorted_avg_dist_violations=[(i[0],i[1],round(i[2],4)) for i in self.sorted_avg_dist_rest]
            
            
       
        if self.rest_flag[1]:
            tot_rest = 0
            for k in self.ang_rest_info.keys():
                tot_rest+=self.ang_rest_info[k]
            ang_summary={}
            ang_summary['Total']=(tot_rest,round((float(tot_rest)/float(tot_rest))*100.0,2))
            for k in self.ang_rest_info.keys():
                ang_summary[k]=(self.ang_rest_info[k],round((float(self.ang_rest_info[k])/float(tot_rest))*100.0,2))
            ang_violation_each_model={}
            for i in range(len(self.model_vs_ang)):
                ang_violation_each_model[i+1]=[(self.model_vs_ang[i][7],round((float(self.model_vs_ang[i][7])/float(self.model_vs_ang[i][0]))*100.0,2)),
                                               (sum(self.model_vs_ang[i][1:-1]),round((float(sum(self.model_vs_ang[i][1:-1]))/float(self.model_vs_ang[i][0]))*100.0,2)),
                                               (self.model_vs_ang[i][1],round((float(self.model_vs_ang[i][1])/float(self.model_vs_ang[i][0]))*100.0,2)),
                                               (self.model_vs_ang[i][2],round((float(self.model_vs_ang[i][2])/float(self.model_vs_ang[i][0]))*100.0,2)),
                                               (self.model_vs_ang[i][3],round((float(self.model_vs_ang[i][3])/float(self.model_vs_ang[i][0]))*100.0,2)),
                                               (self.model_vs_ang[i][4],round((float(self.model_vs_ang[i][4])/float(self.model_vs_ang[i][0]))*100.0,2)),
                                               (self.model_vs_ang[i][5],round((float(self.model_vs_ang[i][5])/float(self.model_vs_ang[i][0]))*100.0,2)),
                                               (self.model_vs_ang[i][6],round((float(self.model_vs_ang[i][6])/float(self.model_vs_ang[i][0]))*100.0,2)),]
                
            ang_violation_each_rest={}
            v=0
            nv=0
            ncv=0
            for j in self.ang_vs_model:
                for i in j:
                    
                    if len(i)-1 > 0:
                        if len(i)-1 == self.maxModels:
                            ncv+=1
                        v+=1
                        ang_violation_each_rest[i[0]]=(len(i)-1,round(sum([float(k[1]) for k in i[1:]])/(len(i)-1),4),','.join([str(k[0]) for k in i[1:]]))
                    else:
                        nv+=1
            ang_viol_summary={}
            ang_viol_summary['Number of consistently violated restraints (violated in all models)']=(ncv,round((float(ncv)/float(v+nv))*100,2))
            ang_viol_summary['Number of violated restraints (at least in one model)']=(v,round((float(v)/float(v+nv))*100,2))
            ang_viol_summary['Number of non violated restraints']=(nv,round((float(nv)/float(v+nv))*100,2))
            ang_viol_summary['Total number of restraints']=(v+nv,round((float(v+nv)/float(v+nv))*100,2))
            
            sorted_ang_violations=[(i[0],i[1],round(i[2],4)) for i in self.sorted_ang_rest]
            sorted_avg_ang_violations=[(i[0],i[1],round(i[2],4)) for i in self.sorted_avg_ang_rest]
            
                
            
                
        jsonlist = json.dumps({'dist_summary':dist_summary,
                                'dist_violation_each_model':dist_violation_each_model,
                                'dist_violation_each_rest':dist_violation_each_rest,
                                'dist_viol_summary':dist_viol_summary,
                                'sorted_dist_violations':sorted_dist_violations,
                                'sorted_avg_dist_violations':sorted_avg_dist_violations,
                                'ang_summary':ang_summary,
                                'ang_violation_each_model':ang_violation_each_model,
                                'ang_violation_each_rest':ang_violation_each_rest,
                                'ang_viol_summary':ang_viol_summary,
                                'sorted_ang_violations':sorted_ang_violations,
                                'sorted_avg_ang_violations':sorted_avg_ang_violations,})
        return jsonlist    
            
            
    def generateReport(self):
        fo = open(self.outFile,'w')
        self.logger.info('Generating output file {}'.format(self.outFile))
        if self.rest_flag[0]:
            fo.write('1. Distance restraints summary\n')
            fo.write('{0:-<30}\n\n'.format(''))
            fo.write('\t{}={} ({}%)\n'.format('{0:<62}'.format('Total number of restraints'),'{0: >6}'.format(len(self.dist_rest_info[0])),round((float(len(self.dist_rest_info[0]))/float(len(self.dist_rest_info[0])))*100.0,2)))
            fo.write('\t{}={} ({}%)\n'.format('{0:<62}'.format('Intraresidue restraints |i-j|=0'),'{0: >6}'.format(len(self.dist_rest_info[1])),round((float(len(self.dist_rest_info[1]))/float(len(self.dist_rest_info[0])))*100.0,2)))
            fo.write('\t{}={} ({}%)\n'.format('{0:<62}'.format('Sequential restraints |i-j|=1'),'{0: >6}'.format(len(self.dist_rest_info[2])),round((float(len(self.dist_rest_info[2]))/float(len(self.dist_rest_info[0])))*100.0,2)))
            fo.write('\t{}={} ({}%)\n'.format('{0:<62}'.format('Medium range  restraints|i-j|>1 and |i-j|<5'),'{0: >6}'.format(len(self.dist_rest_info[3])),round((float(len(self.dist_rest_info[3]))/float(len(self.dist_rest_info[0])))*100.0,2)))
            fo.write('\t{}={} ({}%)\n'.format('{0:<62}'.format('Long range restraints|i-j|>=5'),'{0: >6}'.format(len(self.dist_rest_info[4])),round((float(len(self.dist_rest_info[4]))/float(len(self.dist_rest_info[0])))*100.0,2)))
            fo.write('\n\n')
            fo.write('2. Distance restraint violations\n')
            fo.write('{0:-<32}\n\n'.format(''))
            fo.write('2.1 Number of violated restraints for each model\n')
            fo.write('{0:-<48}\n\n'.format(''))
            fo.write('{}\t{}(%)\t{}(%)\t{}(%)\t{}(%)\t{}(%)\t{}(%)\t{}(%)\t{}(%)\n'.format('{0:>12}'.format('Model'),
                                                                     '{0:>16}'.format('Not Violated'),
                                                                     '{0:>16}'.format('Violated'),
                                                                     '{0:>16}'.format('0A-0.2A'),
                                                                     '{0:>16}'.format('0.2A-0.5A'),
                                                                     '{0:>16}'.format('0.5A-1.0A'),
                                                                     '{0:>16}'.format('1.0A-2.0A'),
                                                                     '{0:>16}'.format('2.0A-5.0A'),
                                                                     '{0:>16}'.format('>5.0A')))
            fo.write('\n')
            for i in range(len(self.model_vs_dist)):
                fo.write('{}\t{}({})\t{}({})\t{}({})\t{}({})\t{}({})\t{}({})\t{}({})\t{}({})\n'.format('{0:>12}'.format(i+1),
                                                                     '{0:>12}'.format(self.model_vs_dist[i][7]),round((float(self.model_vs_dist[i][7])/float(self.model_vs_dist[i][0]))*100.0,2),
                                                                     '{0:>12}'.format(sum(self.model_vs_dist[i][1:-1])),round((float(sum(self.model_vs_dist[i][1:-1]))/float(self.model_vs_dist[i][0]))*100.0,2),
                                                                     '{0:>12}'.format(self.model_vs_dist[i][1]),round((float(self.model_vs_dist[i][1])/float(self.model_vs_dist[i][0]))*100.0,2),
                                                                     '{0:>12}'.format(self.model_vs_dist[i][2]),round((float(self.model_vs_dist[i][2])/float(self.model_vs_dist[i][0]))*100.0,2),
                                                                     '{0:>12}'.format(self.model_vs_dist[i][3]),round((float(self.model_vs_dist[i][3])/float(self.model_vs_dist[i][0]))*100.0,2),
                                                                     '{0:>12}'.format(self.model_vs_dist[i][4]),round((float(self.model_vs_dist[i][4])/float(self.model_vs_dist[i][0]))*100.0,2),
                                                                     '{0:>12}'.format(self.model_vs_dist[i][5]),round((float(self.model_vs_dist[i][5])/float(self.model_vs_dist[i][0]))*100.0,2),
                                                                     '{0:>12}'.format(self.model_vs_dist[i][6]),round((float(self.model_vs_dist[i][6])/float(self.model_vs_dist[i][0]))*100.0,2)))
                
            fo.write('\n\n')
            fo.write('2.2 Number of violated models for each restraint\n')
            fo.write('{0:-<48}\n\n'.format(''))
            fo.write('\t{}\t{}\t{}\t{}\n'.format('{0:<25}'.format('Restraint ID'),
                                             '{0:<25}'.format('No. of Violated Models'),
                                             '{0:<25}'.format('Average Violation(A)'),
                                             '{0:<25}'.format('Violated Models')))
            v=0
            nv=0
            ncv=0
            fo.write("\n")
            ncv_category=[0,0,0,0]
            v_category=[0,0,0,0]
            nv_category=[0,0,0,0]
            for rlistid in range(len(self.dist_vs_model)):
                j=self.dist_vs_model[rlistid]
                for i in j:
                    cat_n=self.checkDistRestCateory(i[0],rlistid)
                    if len(i)-1 > 0:
                        if len(i)-1 == self.maxModels:
                            ncv+=1
                            if cat_n>=0:
                                ncv_category[cat_n]+=1
                            else:
                                self.logger.warning("Rest id {} doesn't belong to any category".format(i[0]))
                        v+=1
                        if cat_n>=0:
                            v_category[cat_n]+=1
                        else:
                            self.logger.warning("Rest id {} doesn't belong to any category".format(i[0]))
                        fo.write('\t{}\t{}\t{}\t{}\n'.format('{0:<25}'.format(i[0]),
                                             '{0:<25}'.format(len(i)-1),
                                             '{0:<25}'.format(round(sum([float(k[1]) for k in i[1:]])/(len(i)-1),4)),
                                             '{0:<25}'.format(','.join([str(k[0]) for k in i[1:]]))))
                    else:
                        nv+=1
                        if cat_n>=0:
                            nv_category[cat_n]+=1
                        else:
                            self.logger.warning("Rest id {} doesn't belong to any category".format(i[0]))
           # ncv_p=[ (float(i)/float(sum(ncv_category)))*100.00,2) for i in ncv_category]
            #print ncv_category
            if sum(ncv_category) > 0 :
                ncv_p=['{}%'.format(round((float(i)/float(sum(ncv_category)))*100.00,2)) for i in ncv_category]
            else:
                ncv_p=['{}%'.format(round((float(i)/float(1))*100.00,2)) for i in ncv_category]
            v_p=['{}%'.format(round((float(i)/float(sum(v_category)))*100.00,2)) for i in v_category]
            nv_p=['{}%'.format(round((float(i)/float(sum(nv_category)))*100.00,2)) for i in nv_category]
            ncv_pt=['{}%'.format(round((float(i)/float(float(v+nv)))*100.00,2)) for i in ncv_category]
            v_pt=['{}%'.format(round((float(i)/float(float(v+nv)))*100.00,2)) for i in v_category]
            nv_pt=['{}%'.format(round((float(i)/float(float(v+nv)))*100.00,2)) for i in nv_category]
            fo.write("\n")
            fo.write('2.3 Distance restraint violations summary\n')
            fo.write('{0:-<41}\n\n'.format(''))
            #fo.write('Number of intraresidue,sequential,medium-range and long-range violations and their percentages are in square bracket\n\n')
            fo.write('\t{}\t{}\t{}\t{}\t{}\t{}\n\n'.format('{0:<65}'.format(' '),
                                                         '{0:<20}'.format('Intrresidue'),
                                                         '{0:<20}'.format('Sequential'),
                                                         '{0:<20}'.format('Medium-range'),
                                                         '{0:<20}'.format('Long-range'),
                                                         '{0:<20}'.format('Total')))
            fo.write('\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('{0:<65}'.format('Number of consistently violated restraints (violated in all models)'),
                                                         '{0:<20}'.format('{} ({})'.format('{0: >6}'.format(ncv_category[0]),ncv_pt[0])),
                                                         '{0:<20}'.format('{} ({})'.format('{0: >6}'.format(ncv_category[1]),ncv_pt[1])),
                                                         '{0:<20}'.format('{} ({})'.format('{0: >6}'.format(ncv_category[2]),ncv_pt[2])),
                                                         '{0:<20}'.format('{} ({})'.format('{0: >6}'.format(ncv_category[3]),ncv_pt[3])),
                                                         '{0:<20}'.format('{} ({}%)'.format('{0: >6}'.format(ncv),round((float(ncv)/float(v+nv))*100,2)))))
            fo.write('\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('{0:<65}'.format('Number of violated restraints (at least in one model)'),
                                                        '{0:<20}'.format('{} ({})'.format('{0: >6}'.format(v_category[0]),v_pt[0])),
                                                         '{0:<20}'.format('{} ({})'.format('{0: >6}'.format(v_category[1]),v_pt[1])),
                                                         '{0:<20}'.format('{} ({})'.format('{0: >6}'.format(v_category[2]),v_pt[2])),
                                                         '{0:<20}'.format('{} ({})'.format('{0: >6}'.format(v_category[3]),v_pt[3])),
                                                         '{0:<20}'.format('{} ({}%)'.format('{0: >6}'.format(v),round((float(v)/float(v+nv))*100,2)))))
            fo.write('\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('{0:<65}'.format('Number of non violated restraints'),
                                                        '{0:<20}'.format('{} ({})'.format('{0: >6}'.format(nv_category[0]),nv_pt[0])),
                                                         '{0:<20}'.format('{} ({})'.format('{0: >6}'.format(nv_category[1]),nv_pt[1])),
                                                         '{0:<20}'.format('{} ({})'.format('{0: >6}'.format(nv_category[2]),nv_pt[2])),
                                                         '{0:<20}'.format('{} ({})'.format('{0: >6}'.format(nv_category[3]),nv_pt[3])),
                                                         '{0:<20}'.format('{} ({}%)'.format('{0: >6}'.format(nv),round((float(nv)/float(v+nv))*100,2)))))
            fo.write('\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('{0:<65}'.format('Total number of restraints'),
                                                         '{0:<20}'.format('{} ({}%)'.format('{0: >6}'.format(len(self.dist_rest_info[1])),round((float(len(self.dist_rest_info[1]))/float(len(self.dist_rest_info[0])))*100.0,2))),
                                                         '{0:<20}'.format('{} ({}%)'.format('{0: >6}'.format(len(self.dist_rest_info[2])),round((float(len(self.dist_rest_info[2]))/float(len(self.dist_rest_info[0])))*100.0,2))),
                                                         '{0:<20}'.format('{} ({}%)'.format('{0: >6}'.format(len(self.dist_rest_info[3])),round((float(len(self.dist_rest_info[3]))/float(len(self.dist_rest_info[0])))*100.0,2))),
                                                         '{0:<20}'.format('{} ({}%)'.format('{0: >6}'.format(len(self.dist_rest_info[4])),round((float(len(self.dist_rest_info[4]))/float(len(self.dist_rest_info[0])))*100.0,2))),
                                                         '{0:<20}'.format('{} ({}%)'.format('{0: >6}'.format(v+nv),round((float(v+nv)/float(v+nv))*100,2)))))
            #fo.write('\t{}={} ({}%) [({})({})({})]\n'.format('{0:<72}'.format('Number of consistently violated restraints (violated in all models) '),'{0: >6}'.format(ncv),round((float(ncv)/float(v+nv))*100,2),",".join([str(i) for i in ncv_category]),",".join(ncv_p),",".join(ncv_pt)))
            #fo.write('\t{}={} ({}%) [({})({})({})]\n'.format('{0:<72}'.format('Number of violated restraints (at least in one model)'),'{0: >6}'.format(v),round((float(v)/float(v+nv))*100,2),",".join([str(i) for i in v_category]),",".join(v_p),",".join(v_pt)))            
            #fo.write('\t{}={} ({}%) [({})({})({})]\n'.format('{0:<72}'.format('Number of non violated restraints'),'{0: >6}'.format(nv),round((float(nv)/float(v+nv))*100,2),",".join([str(i) for i in nv_category]),",".join(nv_p),",".join(nv_p))) 
            #fo.write('\t{}={} ({}%) \n'.format('{0:<72}'.format('Total number of restraints'),'{0: >6}'.format(v+nv),round((float(v+nv)/float(v+nv))*100,2)))      
            fo.write("\n\n")
            fo.write('2.4 List of largest 10 distance violations\n')
            fo.write('{0:-<42}\n\n'.format(''))
            fo.write('\t{}\t{}\t{}\n'.format('{0:<30}'.format('Restraint ID'),
                                             '{0:<30}'.format('Model ID'),
                                             '{0:<30}'.format('|Exp. value - Model value|A')))
            fo.write("\n")
            for i in self.sorted_dist_rest[:10]:
                fo.write('\t{}\t{}\t{}\n'.format('{0:<30}'.format(i[0]),
                                             '{0:<30}'.format(i[1]),
                                             '{0:<30}'.format(round(i[2],4))))
            fo.write("\n\n")
            fo.write('2.5 List of largest 10 average distance violations\n')
            fo.write('{0:-<50}\n\n'.format(''))
            fo.write('\t{}\t{}\t{}\n'.format('{0:<30}'.format('Restraint ID'),
                                             '{0:<30}'.format('Number of Models'),
                                             '{0:<30}'.format('|Exp. value - Model value|A')))
            fo.write("\n")
            for i in self.sorted_avg_dist_rest[:10]:
                fo.write('\t{}\t{}\t{}\n'.format('{0:<30}'.format(i[0]),
                                             '{0:<30}'.format(i[1]),
                                             '{0:<30}'.format(round(i[2],4))))
                
            fo.write('{0:=<120}\n\n'.format(''))
        if self.rest_flag[1]:
            tot_rest = 0
            for k in self.ang_rest_info.keys():
                tot_rest+=self.ang_rest_info[k]
            fo.write('3. Angle restraints summary\n')
            fo.write('{0:-<27}\n\n'.format(''))
            fo.write('\t{}={} ({}%)\n'.format('{0:<62}'.format('Total number of restraints'),'{0: >6}'.format(tot_rest),round((float(tot_rest)/float(tot_rest))*100.0,2)))
            for k in self.ang_rest_info.keys():
                fo.write('\t{}={} ({}%)\n'.format('{0:<62}'.format('Torsion angle Name '+k),'{0: >6}'.format(self.ang_rest_info[k]),round((float(self.ang_rest_info[k])/float(tot_rest))*100.0,2)))
            fo.write('\n\n')
            fo.write('4. Angle restraint violations\n')
            fo.write('{0:-<30}\n\n'.format(''))
            fo.write('4.1 Number of violated restraints for each model\n')
            fo.write('{0:-<48}\n\n'.format(''))
            fo.write('{}\t{} (%)\t{} (%)\t{} (%)\t{} (%)\t{} (%)\t{} (%)\t{} (%)\t{} (%)\n'.format('{0:>12}'.format('Model'),
                                                                     '{0:>12}'.format('Not Violated'),
                                                                     '{0:>12}'.format('Violated'),
                                                                     '{0:>12}'.format('0-5'),
                                                                     '{0:>12}'.format('5-10'),
                                                                     '{0:>12}'.format('10-20'),
                                                                     '{0:>12}'.format('20-40'),
                                                                     '{0:>12}'.format('40-80'),
                                                                     '{0:>12}'.format('>80')))
            fo.write('\n')
            for i in range(len(self.model_vs_ang)):
                fo.write('{}\t{}({})\t{}({})\t{}({})\t{}({})\t{}({})\t{}({})\t{}({})\t{}({})\n'.format('{0:>12}'.format(i+1),
                                                                     '{0:>12}'.format(self.model_vs_ang[i][7]),round((float(self.model_vs_ang[i][7])/float(self.model_vs_ang[i][0]))*100.0,2),
                                                                     '{0:>12}'.format(sum(self.model_vs_ang[i][1:-1])),round((float(sum(self.model_vs_ang[i][1:-1]))/float(self.model_vs_ang[i][0]))*100.0,2),
                                                                     '{0:>12}'.format(self.model_vs_ang[i][1]),round((float(self.model_vs_ang[i][1])/float(self.model_vs_ang[i][0]))*100.0,2),
                                                                     '{0:>12}'.format(self.model_vs_ang[i][2]),round((float(self.model_vs_ang[i][2])/float(self.model_vs_ang[i][0]))*100.0,2),
                                                                     '{0:>12}'.format(self.model_vs_ang[i][3]),round((float(self.model_vs_ang[i][3])/float(self.model_vs_ang[i][0]))*100.0,2),
                                                                     '{0:>12}'.format(self.model_vs_ang[i][4]),round((float(self.model_vs_ang[i][4])/float(self.model_vs_ang[i][0]))*100.0,2),
                                                                     '{0:>12}'.format(self.model_vs_ang[i][5]),round((float(self.model_vs_ang[i][5])/float(self.model_vs_ang[i][0]))*100.0,2),
                                                                     '{0:>12}'.format(self.model_vs_ang[i][6]),round((float(self.model_vs_ang[i][6])/float(self.model_vs_ang[i][0]))*100.0,2)))
            fo.write('\n\n')
            fo.write('4.2 Number of violated models for each restraint\n')
            fo.write('{0:-<48}\n\n'.format(''))
            fo.write('\t{}\t{}\t{}\t{}\n'.format('{0:<25}'.format('Restraint ID'),
                                             '{0:<25}'.format('No. of Violated Models'),
                                             '{0:<25}'.format('Average Violation(A)'),
                                             '{0:<25}'.format('Violated Models')))
            v=0
            nv=0
            ncv=0
            fo.write("\n")
            for j in self.ang_vs_model:
                for i in j:
                    
                    if len(i)-1 > 0:
                        if len(i)-1 == self.maxModels:
                            ncv+=1
                        v+=1
                        fo.write('\t{}\t{}\t{}\t{}\n'.format('{0:<25}'.format(i[0]),
                                             '{0:<25}'.format(len(i)-1),
                                             '{0:<25}'.format(round(sum([float(k[1]) for k in i[1:]])/(len(i)-1),4)),
                                             '{0:<25}'.format(','.join([str(k[0]) for k in i[1:]]))))
                    else:
                        nv+=1
            fo.write("\n")
            fo.write('4.3 Angle restraint violations summary\n')
            fo.write('{0:-<38}\n\n'.format(''))
            fo.write('\t{}={} ({}%)\n'.format('{0:<72}'.format('Number of consistently violated restraints (violated in all models) '),'{0: >6}'.format(ncv),round((float(ncv)/float(v+nv))*100,2)))
            fo.write('\t{}={} ({}%)\n'.format('{0:<72}'.format('Number of violated restraints (at least in one model)'),'{0: >6}'.format(v),round((float(v)/float(v+nv))*100,2)))            
            fo.write('\t{}={} ({}%)\n'.format('{0:<72}'.format('Number of non violated restraints'),'{0: >6}'.format(nv),round((float(nv)/float(v+nv))*100,2))) 
            fo.write('\t{}={} ({}%)\n'.format('{0:<72}'.format('Total number of restraints'),'{0: >6}'.format(v+nv),round((float(v+nv)/float(v+nv))*100,2)))      
            fo.write("\n\n")
            fo.write('4.4 List of largest 10 angle violations\n')
            fo.write('{0:-<39}\n\n'.format(''))
            fo.write('\t{}\t{}\t{}\n'.format('{0:<30}'.format('Restraint ID'),
                                             '{0:<30}'.format('Model ID'),
                                             '{0:<30}'.format('|Exp. value - Model value|A')))
            fo.write("\n")
            for i in self.sorted_ang_rest[:10]:
                fo.write('\t{}\t{}\t{}\n'.format('{0:<30}'.format(i[0]),
                                             '{0:<30}'.format(i[1]),
                                             '{0:<30}'.format(round(i[2],4))))
            fo.write("\n\n")
            fo.write('4.5 List of largest 10 average angle violations\n')
            fo.write('{0:-<47}\n\n'.format(''))
            fo.write('\t{}\t{}\t{}\n'.format('{0:<30}'.format('Restraint ID'),
                                             '{0:<30}'.format('Number of models'),
                                             '{0:<30}'.format('|Exp. value - Model value|A')))
            fo.write("\n")
            for i in self.sorted_avg_ang_rest[:10]:
                fo.write('\t{}\t{}\t{}\n'.format('{0:<30}'.format(i[0]),
                                             '{0:<30}'.format(i[1]),
                                             '{0:<30}'.format(round(i[2],4))))
            
            fo.write('{0:=<120}\n\n'.format(''))
        self.logger.info('Output file {} generated'.format(self.outFile))
        self.logger.info('Report generated successfully')
        
    def distRestAnalysis(self):
        self.model_vs_dist=[]
        for m in range(self.maxModels):
            self.model_vs_dist.append([0,0,0,0,0,0,0,0])
            for rest_list in self.cal_dist:
                for key in sorted(rest_list):
                    self.model_vs_dist[-1][0]+=1
                    dav=[]
                    d_flag=False
                    d_list=[]
                    for rest in rest_list[key]:
                        d_list.append(rest[m+1])
                    
                    r6av=self.r6sum(d_list)  
                    if r6av < rest[0][0] or r6av > rest[0][1]:
                        d_flag = True
                        if r6av < rest[0][0]:
                            dav.append(rest[0][0]-r6av)
                        else:
                            dav.append(r6av-rest[0][1])
                        
#                         d_min = rest[0][0]
#                         d_max = rest[0][1]
#                         d=rest[m+1]
#                         if d>d_max or d<d_min:
#                             d_flag=True
#                             if d<d_max:
#                                 dd=d_max-d
#                             else:
#                                 dd=d-d_min
#                             dav.append(dd)
                    if d_flag:
                        dif=sum(dav)/len(dav)
                        
                        if dif>0 and dif <= 0.2:
                            self.model_vs_dist[-1][1]+=1
                        elif dif > 0.2 and dif <= 0.5:
                            self.model_vs_dist[-1][2]+=1
                        elif dif>0.5 and dif <=1.0:
                            self.model_vs_dist[-1][3]+=1
                        elif dif>1.0 and dif <=2.0:
                            self.model_vs_dist[-1][4]+=1
                        elif dif>2.0 and dif <=5.0:
                            self.model_vs_dist[-1][5]+=1
                        else:
                            self.model_vs_dist[-1][6]+=1
                    else:
                        self.model_vs_dist[-1][7]+=1
        self.dist_vs_model=[]
        for rest_list in self.cal_dist:
            self.dist_vs_model.append([])
            for key in sorted(rest_list):
                mm=[key]
                    #print rest
                for i in range(1,self.maxModels+1):
                    d_flag=False
                    err = []
                    d_list=[]
                    for rest in rest_list[key]:
                        d_list.append(rest[i])
                    r6avg = self.r6sum(d_list)
                    if r6avg<rest[0][0] or r6avg>rest[0][1]:
                        d_flag = True
                        if r6avg < rest[0][0]:
                            err.append(rest[0][0]-r6avg)
                        else:
                            err.append(r6avg-rest[0][1])
#                         if rest[i]<rest[0][0] or rest[i]>rest[0][1]:
#                             d_flag=True
#                             if rest[i]<rest[0][0]:
#                                 err.append(rest[0][0]-rest[i])
#                             else:
#                                 err.append(rest[i]-rest[0][1])
                    if d_flag:
                        err_avg = sum(err)/len(err)
                        mm.append((i,err_avg))
                    
                self.dist_vs_model[-1].append(mm)               
        return [self.model_vs_dist,self.dist_vs_model]           
        #return self.dist_rest
        
    def angleRestAnalysis(self):
        self.model_vs_ang=[]
        for m in range(self.maxModels):
            self.model_vs_ang.append([0,0,0,0,0,0,0,0])
            for rest_list in self.cal_ang:
                for key in sorted(rest_list):
                    self.model_vs_ang[-1][0]+=1
                    dav=[]
                    d_flag=False
                    rest=rest_list[key]
                    d_min = rest[0][0]
                    d_max = rest[0][1]
                    d=rest[m+1]
                    dif=0.0
                    if (d_max > 0 and d_min > 0) and d < 0:
                        d= 360.0+d 
                    if (d_max < 0 and d_min < 0) and d > 0:
                        d= -1*(360.0-d) 
                    if d>d_max :
                        dif = d-d_max
                        d_flag=True
                        self.logger.debug('dmax {}:{},{},{},{},{}'.format(key,d_min,d_max,d,dif,d_flag))
                    elif d<d_min :
                        dif = d_min-d
                        d_flag = True
                        self.logger.debug('dmin {}:{},{},{},{},{}'.format(key,d_min,d_max,d,dif,d_flag))
                    else:
                        d_flag = False
                        self.logger.debug('dmin {}:{},{},{},{},{}'.format(key,d_min,d_max,d,dif,d_flag))
                    self.logger.debug('{}:{},{},{},{}'.format(key,d_min,d_max,d,dif))
                    if d_flag:
                        if dif>0 and dif <= 5:
                            self.model_vs_ang[-1][1]+=1
                        elif dif > 5 and dif <= 10:
                            self.model_vs_ang[-1][2]+=1
                        elif dif>10 and dif <=20:
                            self.model_vs_ang[-1][3]+=1
                        elif dif>20 and dif <=40:
                            self.model_vs_ang[-1][4]+=1
                        elif dif>40 and dif <=80:
                            self.model_vs_ang[-1][5]+=1
                        else:
                            self.model_vs_ang[-1][6]+=1
                    else:
                        self.model_vs_ang[-1][7]+=1
        self.ang_vs_model=[]
        for rest_list in self.cal_ang:
            self.ang_vs_model.append([])
            for key in sorted(rest_list):
                mm=[key]
                    #print rest
                for i in range(1,self.maxModels+1):
                    d_flag=False
                    rest=rest_list[key]
                    rest_ang = rest[i]
                    if ( rest[0][0]> 0 and rest[0][1]> 0 ) and rest[i] <0:
                        rest_ang = 360.0+rest[i]
                    if ( rest[0][0]< 0 and rest[0][1]<0 ) and rest[i] >0:
                        rest_ang = -1*(360.0-rest[i])
                    if rest_ang<rest[0][0] or rest_ang>rest[0][1]:
                        d_flag=True
                        
                        if rest_ang<rest[0][0]:
                            err_ang=rest[0][0]-rest_ang
                        else:
                            err_ang=rest_ang-rest[0][1]
                    if d_flag:
                        mm.append((i,err_ang))
                self.ang_vs_model[-1].append(mm)               
        return [self.model_vs_ang,self.ang_vs_model]           
    def distRestraintStat(self):
        self.dist_rest_info = [[],[],[],[],[]]
        for rlistid  in range(len(self.dist_rest)):
            rest_list=self.dist_rest[rlistid]
            for restid in rest_list.keys():
                for rest in rest_list[restid]: 
                    try:
                        n=abs(int(rest[0][0])-int(rest[1][0]))
                    except ValueError:
                        self.logger.warning("Invalid restraint {}".format(rest))
                self.dist_rest_info[0].append((restid,rlistid))
                if n==0:
                    self.dist_rest_info[1].append((restid,rlistid))
                elif n==1:
                    self.dist_rest_info[2].append((restid,rlistid))
                elif n>1 and n<5:
                    self.dist_rest_info[3].append((restid,rlistid))
                else:
                    self.dist_rest_info[4].append((restid,rlistid))
        #self.logger.debug('Distance restraint statistics {}'.format(self.dist_rest_info))        
        return self.dist_rest_info
    
    def angRestraintStat(self):
        self.ang_rest_info = {}
        
        for rest_list in self.ang_rest:
            for rest in rest_list.keys():
                if rest_list[rest][0] not in self.ang_rest_info.keys():
                    self.ang_rest_info[rest_list[rest][0]]=0
                self.ang_rest_info[rest_list[rest][0]]+=1
        
        self.logger.debug('Angle restraint statistics {}'.format(self.ang_rest_info))
        
    
    
    def getAngleRestraints(self):
        
        """Creates a dictionary of angle restraints from _Torsion_angle_constraint loop.
        ang_rest[atom 1,atom 2, atom 3, atom 4] = [lower bound, upper bound]
        atoms are identified by (Comp_index_ID,Entity_assembly_ID,Comp_ID,Atom_ID)"""
        
        ang = self.starData.get_loops_by_category('_Torsion_angle_constraint')
       
        self.ang_rest = []
        for dl in ang:
            self.ang_rest.append({})
            restid = dl.get_tag_names().index("_Torsion_angle_constraint.ID")
            rest_name = dl.get_tag_names().index("_Torsion_angle_constraint.Torsion_angle_name")
            atomid1 = dl.get_tag_names().index("_Torsion_angle_constraint.Atom_ID_1")
            atomid2 = dl.get_tag_names().index("_Torsion_angle_constraint.Atom_ID_2")
            atomid3 = dl.get_tag_names().index("_Torsion_angle_constraint.Atom_ID_3")
            atomid4 = dl.get_tag_names().index("_Torsion_angle_constraint.Atom_ID_4")
            compid1 = dl.get_tag_names().index("_Torsion_angle_constraint.Comp_ID_1")
            compid2 = dl.get_tag_names().index("_Torsion_angle_constraint.Comp_ID_2")
            compid3 = dl.get_tag_names().index("_Torsion_angle_constraint.Comp_ID_3")
            compid4 = dl.get_tag_names().index("_Torsion_angle_constraint.Comp_ID_4")
            #entityid1 = dl.get_tag_names().index("_Torsion_angle_constraint.Entity_assembly_ID_1")
            #entityid2 = dl.get_tag_names().index("_Torsion_angle_constraint.Entity_assembly_ID_2")
            #entityid3 = dl.get_tag_names().index("_Torsion_angle_constraint.Entity_assembly_ID_3")
            #entityid4 = dl.get_tag_names().index("_Torsion_angle_constraint.Entity_assembly_ID_4")
            entityid1 = dl.get_tag_names().index("_Torsion_angle_constraint.Auth_asym_ID_1")
            entityid2 = dl.get_tag_names().index("_Torsion_angle_constraint.Auth_asym_ID_2")
            entityid3 = dl.get_tag_names().index("_Torsion_angle_constraint.Auth_asym_ID_3")
            entityid4 = dl.get_tag_names().index("_Torsion_angle_constraint.Auth_asym_ID_4")
            seqid1 = dl.get_tag_names().index("_Torsion_angle_constraint.Comp_index_ID_1") 
            seqid2 = dl.get_tag_names().index("_Torsion_angle_constraint.Comp_index_ID_2")
            seqid3 = dl.get_tag_names().index("_Torsion_angle_constraint.Comp_index_ID_3") 
            seqid4 = dl.get_tag_names().index("_Torsion_angle_constraint.Comp_index_ID_4")
            lb = dl.get_tag_names().index("_Torsion_angle_constraint.Angle_lower_bound_val")
            ub = dl.get_tag_names().index("_Torsion_angle_constraint.Angle_upper_bound_val")
            for i in dl:
                try:
                    lbv = float(i[lb])
                except ValueError:
                    lbv = -999
                try:
                    ubv = float(i[ub])
                except ValueError:
                    ubv = -999
                self.ang_rest[-1][int(i[restid])]=[i[rest_name],(i[seqid1],i[entityid1],i[compid1],i[atomid1]),(i[seqid2],i[entityid2],i[compid2],i[atomid2]),(i[seqid3],i[entityid3],i[compid3],i[atomid3]),(i[seqid4],i[entityid4],i[compid4],i[atomid4]),(lbv,ubv)]
        self.logger.info('Number of angle restraints : {}'.format([len(i) for i in self.ang_rest]))
        #return self.ang_rest        
    
    def calDisanceFromMoels(self):
        c=[]
        for modelID in range(1,self.maxModels+1):
            c.append(self.getCoordinates(modelID))
        self.cal_dist=[]
        for rest_list in self.dist_rest:
            self.cal_dist.append({})
            for rest_id in rest_list.keys():
                if rest_id not in self.cal_dist[-1].keys():
                    self.cal_dist[-1][rest_id]=[]
                for rest in rest_list[rest_id]:
                    c_dist = [rest[2]]
                    for m in c:
                        try:
                            d = self.getDistance(m[rest[0]],m[rest[1]])
                        except KeyError:
                            self.logger.warning('One of the atom from the restraint not found in coordinates {},{},{}'.format(rest_id,rest[0],rest[1]))
                            d = '.'
                        c_dist.append(d)
                    self.cal_dist[-1][rest_id].append(c_dist)
    
    def calAngleFromMoels(self):
        c=[]
        for modelID in range(1,self.maxModels+1):
            c.append(self.getCoordinates(modelID))
        self.cal_ang=[]
        for rest_list in self.ang_rest:
            self.cal_ang.append({})
            for rest_id in rest_list.keys():
                rest=rest_list[rest_id]
                
                if rest_id not in self.cal_ang[-1].keys():
                    self.cal_ang[-1][rest_id]=[rest[5]]
                
                for m in c:
                    try:
                        d = self.getDihedralAngle(m[rest[1]],m[rest[2]],m[rest[3]],m[rest[4]])
                    except KeyError:
                        self.logger.warning('One of the atom from the restraint not found in coordinates {},{},{},{},{}'.format(rest_id,rest[1],rest[2],rest[2],rest[2]))
                        d = '.'
                    self.cal_ang[-1][rest_id].append(d)
    
    def checkDistanceRestraints(self):
        mbin=[]
        for modelId in range(1,self.maxModels+1):
            mbin.append([0,0,0,0,0,0,0])
            c = self.getCoordinates(modelId)
            for i in self.dist_rest:
                try:
                    d= self.getDistance(c[i[0]],c[i[1]])
                    if d<self.dist_rest[i][1] and d>self.dist_rest[i][0]:
                        mbin[-1][0]+=1
                        self.logger.info("Rest ok")
                    else:
                        if d>self.dist_rest[i][1]:
                            dif = d - self.dist_rest[i][1]
                        else:
                            dif = self.dist_rest[i][0]-d
                        if dif>0 and dif <= 0.2:
                            mbin[-1][1]+=1
                        elif dif > 0.2 and dif <= 0.5:
                            mbin[-1][2]+=1
                        elif dif>0.5 and dif <=1.0:
                            mbin[-1][3]+=1
                        elif dif>1.0 and dif <=2.0:
                            mbin[-1][4]+=1
                        elif dif>2.0 and dif <=5.0:
                            mbin[-1][5]+=1
                        else:
                            mbin[-1][6]+=1
                    #print (modelId,self.getDistance(c[i[0]],c[i[1]]),self.dist_rest[i][1],diff)
                except KeyError: 
                    self.logger.warning("Unknown atoms found in restraint file {} {}".format(i[0],i[1]))
        return mbin
                

if __name__=="__main__":
    p = ValidateRestraints('pdb_examples/6hvb.cif','pdb_examples/6hvb_linked.str')
    #p = ValidateRestraints('nef_examples/2l9r.cif', 'nef_examples/2l9r.str')

        
    
    
    