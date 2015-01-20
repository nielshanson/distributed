#!/usr/bin/python
# File created on Nov 27 Jan 2012
from __future__ import division

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2010, The metapaths Project"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
     import os
     from os import  makedirs, sys, remove, path
     import re
     from optparse import OptionParser, OptionGroup

     from python_resources.SpeciesComputation import *
     from python_resources.metapaths_utils  import  fprintf, printf
     #from libs.python_modules.sysutil import getstatusoutput
except:
     print """ Could not load some user defined  module functions"""
     print """ """
     sys.exit(3)


objective = re.compile(r'Objective.*=\s*(\d*)')


usage=  sys.argv[0] + " --ncbi-file ncbi_taxfile --pathways-file pathways_list  --enzymes-file enzymes_file"""

glpsol = "/usr/local/bin/glpsol"

parser = OptionParser(usage)
parser.add_option("--ncbi-file", dest="ncbi_file",
                  help='name of the refseq names file')
parser.add_option("--pathways-file", dest="pathways_file",
                  help='name of the pathways file')
parser.add_option("--enzymes-file", dest="enzymes_file",
                  help='name of the enzymes/taxons file')

#parser.add_option("--data-file", dest="data_file",
#                  help='name of the test strings')

def copyList(a, b):
    [ b.append(x) for x in a ]


def get_species(hit):
    species = []
    try:
        m = re.findall(r'\[([^\[]+)\]', hit)
        if m != None:
          copyList(m,species)
    except:
          return None

    if species:
       return species
    else:
       return None


def create_annotation(namefile):

    file = 'ncbi_taxonomy_tree.txt'
    lca = SpeciesComputation(file)

    #taxonomy=lca.getTaxonomy(species)

def read_list_file_to_list(filename,list):
    namefile = open(filename, 'r')
    lines = namefile.readlines()
  
    for line in lines:
       line = line.strip()
       if line:
          list.append(line)


def read_taxonomy_list(filename,namesdict):
    namefile = open(filename, 'r')
    lines = namefile.readlines()
    namefile.close()
  
    for line in lines:
       line = line.strip()
       if line:
          names = get_species(line)
          if names:
            for name in names:
                namesdict[name]=1


def read_enzymes_list(enzymes_file, enzymes):
    enzymesfile = open(enzymes_file, 'r')
    lines = enzymesfile.readlines()
    enzymesfile.close()
  
    for line in lines:
       line = line.strip()
       if line:
          fields = [ x.strip() for x in line.split('\t') ]
          if len(fields[0]) >= 6:
             enzyme_name = fields[0]
             taxonomy = fields[6]
             enzymes[enzyme_name] = taxonomy

def get_buffer():
    return "                                                                                       "

def writeToLine(str1, content, offset):
    i = offset
    strlist = list(str1)
    #print strlist 
    for c in content:
       strlist[i] = c
       i+=1
    str1 = ''.join(strlist)
    return str1
     

def offSize(n):
    if n==0:
       return 0
    if n==1:
       return 1
    if n==2:
       return 4
    if n==3:
       return 14
    if n==4:
       return 24
    if n==5:
       return 39
    if n==5:
       return 49


def field(string, n):
    if n==0:
       return string

    if n==1:
       return " "+ string

    if n==2:
       return "    "+ string

    if n==3:
       return "           "+ string

    if n==4:
       return "                        "+ string

    if n==5:
       return "                                       "+ string

    if n==5:
       return "                                                 "+ string



def read_pathways_list(pathways_file, pathways):
    pathwaysfile = open(pathways_file, 'r')
    lines = pathwaysfile.readlines()
    pathwaysfile.close()
  
    for line in lines:
       line = line.strip()
       if line:
          fields = [ x.strip() for x in line.split('\t') ]
          if fields[0] == 'PATHWAY:':
             pathway_name = fields[1]
             pathway_common_name = fields[2]
             #print pathway_name + ' ' + pathway_common_name
             pathways[ pathway_name] = { 'common_name': pathway_common_name, 'rxns':{}, 'rxns-type':{} }
       

          if fields[0] == 'RXN:':
             rxn_name = fields[1]
             rxn_common_name = fields[2]
             pathways[pathway_name]['rxns'][rxn_name] = {}
             pathways[pathway_name]['rxns-type'][rxn_name] = fields[3].strip()
             if len(fields) > 4:
                 for i in  range(4,len(fields)):
                     enzyme = fields[i]
                     pathways[pathway_name]['rxns'][rxn_name][enzyme] = True
                       

def print_taxonomy_information(pathways, ranked_pathways, enzymes, lca):
        for  key in ranked_pathways:
            print key
            taxons = []
            for rxn in pathways[key]['rxns']:   
               for enzyme in pathways[key]['rxns'][rxn]:   
                    taxons.append(pathways[key]['rxns'][rxn][enzyme])
                    
            lca.build_independent_taxons(taxons)

            for rxn in pathways[key]['rxns']:   
              print '\t' + rxn + '\t' + pathways[key]['rxns-type'][rxn]
              for enzyme in pathways[key]['rxns'][rxn]:   
                 printf("\t\t%s\t%s",  enzyme, pathways[key]['rxns'][rxn][enzyme])
                 printf(" [ ");
                 for child in lca.get_independent_parent(pathways[key]['rxns'][rxn][enzyme]):
                    printf("\t%s", child)
                 printf(" ]\n");
                
             
       #      print '      '+ enzyme + '  ' + pathways[key]['rxns'][rxn][enzyme]
def add_taxonomy_information(pathways, enzymes):
    for  key in pathways:
        for rxn in pathways[key]['rxns']:   
           for enzyme in pathways[key]['rxns'][rxn]:  # get orfs for  
             if enzyme in enzymes: # key to taxonomy
                pathways[key]['rxns'][rxn][enzyme] = enzymes[enzyme] # put taxonomy on reaction
             
       #      print '      '+ enzyme + '  ' + pathways[key]['rxns'][rxn][enzyme]
       # indlist = lca.get_independent_taxons(taxons)
       # for ind in indlist:
       #     print '             ' + ind
         
   
def check_arguments(opts):
    if opts.ncbi_file == None or opts.enzymes_file == None or opts.pathways_file == None:
       print usage
       sys.exit(0)


def compute_min_species(pathways,p, lca):
    taxons = []
    for rxn in pathways[p]['rxns']:   
       for enzyme in pathways[p]['rxns'][rxn]:   
            taxons.append(pathways[p]['rxns'][rxn][enzyme])

    indtaxons = lca.get_independent_taxons(taxons)
    
    #create X variables 
    taxon_variable = {} 
    i = 0
    for x in indtaxons:
      taxon_variable[x] = 'X' + str(i)
      i+=1

    try:
       mpsinput = open('input.mps','w')
    except IOError:
       print """Cannot open \'input.mps\' to write problem"""
       sys.exit(0)
     
    fprintf(mpsinput, "NAME%s\n",field('DISTRIBUTED',3))

    valid_reactions = []
    for rxn in pathways[p]['rxns']:
       if pathways[p]['rxns'][rxn] :  # consider only reactions that have at least one ORF
           valid_reactions.append(rxn)


    # create X variables 
    reaction_names = {} 
    i = 0
    for x in valid_reactions:
      reaction_names[x] = 'RXN' + str(i)
      i+=1
    
    # compute a hash from taxons to reactions
    taxons_reactions = {}
    for rxn in pathways[p]['rxns']:
       if pathways[p]['rxns'][rxn] :  # consider only reactions that have at least one ORF
           for enzyme in  pathways[p]['rxns'][rxn]: 
              for S in lca.get_independent_parent(pathways[p]['rxns'][rxn][enzyme]):
                  if not S in taxons_reactions: 
                     taxons_reactions[S] = {}
                  taxons_reactions[S][rxn] = True 
    
    # write the rows
    fprintf(mpsinput, "ROWS\n")
    str1 = get_buffer()
    str1 = writeToLine(str1,"N",offSize(1))
    str1 = writeToLine(str1,"COST",offSize(2))
    fprintf(mpsinput, "%s\n",str1.rstrip())
    for rxn in valid_reactions:
       str1 = get_buffer()
       str1 = writeToLine(str1,"G", offSize(1))
       str1 = writeToLine(str1, reaction_names[rxn],offSize(2))
       fprintf(mpsinput, "%s\n",str1.rstrip())

    #write the columns
    fprintf(mpsinput, "COLUMNS\n")
    for X in indtaxons:
       str1 = get_buffer()
       str1 = writeToLine(str1, taxon_variable[X],offSize(2))
       str1 = writeToLine(str1,"COST",offSize(3))
       str1 = writeToLine(str1,"1",offSize(4))
       fprintf(mpsinput, "%s\n",str1.rstrip())
       for rxn in taxons_reactions[X]: 
          str1 = get_buffer()
          str1 = writeToLine(str1, taxon_variable[X],offSize(2))
          str1 = writeToLine(str1,reaction_names[rxn],offSize(3))
          str1 = writeToLine(str1,"1",offSize(4))
          fprintf(mpsinput, "%s\n",str1.rstrip())

    #write the columns
    fprintf(mpsinput, "RHS\n")
    for rxn in reaction_names:
        str1 = get_buffer()
        str1 = writeToLine(str1, "RHS1",offSize(2))
        str1 = writeToLine(str1, reaction_names[rxn], offSize(3))
        str1 = writeToLine(str1,"1",offSize(4))
        fprintf(mpsinput, "%s\n",str1.rstrip())

    #write the bounds
    fprintf(mpsinput, "BOUNDS\n")
    for X in indtaxons:
       str1 = get_buffer()
       str1 = writeToLine(str1, "BV", offSize(1))
       str1 = writeToLine(str1,"BND1", offSize(2))
       str1 = writeToLine(str1,taxon_variable[X],offSize(3))
       fprintf(mpsinput, "%s\n",str1.rstrip())

    fprintf(mpsinput, "ENDATA\n")
    mpsinput.close()

    command = glpsol  + " --mps input.mps -o output.sol >> /dev/null"

    os.system(command)

    try:
       glpout = open('output.sol','r')
    except IOError:
       print """Cannot open \'ouptut.sol\' to read solution"""
       sys.exit(0)
    
    solLines = glpout.readlines()
    glpout.close()
    for s in solLines:
        hits = objective.search(s.strip()) 
        if hits:
           value = int(hits.group(1))
           break
     
     
    return value

    #for rxn in pathways[p]['rxns']:
    #   if pathways[p]['rxns'][rxn] :  # consider only reactions that have at least one ORF
    #       for enzyme in  pathways[p]['rxns'][rxn]: 
    #          print rxn +'\t' +  enzyme + '\t' + pathways[p]['rxns'][rxn][enzyme] + '\t' + str(lca.get_independent_parent(pathways[p]['rxns'][rxn][enzyme]))
                 

#+ '\t' + str(lca.get_independent_parent(pathways[key]['rxns'][rxn][enzyme]))
         

# the main function
def main(argv): 

    (opts, args) = parser.parse_args()
    argv = check_arguments(opts)

    namesdict = {} # hash of all the names of the NCBI Taxonomy tree.
    read_taxonomy_list(opts.ncbi_file, namesdict)
    
    # create NCBI Taxonomy tree
    lca = SpeciesComputation(opts.ncbi_file)
    # lca.construct_tree()
    
    enzymes={} # read to taxonomic annotation mapping
    read_enzymes_list(opts.enzymes_file, enzymes)
    
    # Example structure
    # 'GLUTORN-PWY': {'common_name': 'ornithine biosynthesis',
    #                 'rxns-type': {'ACETYLORNDEACET-RXN': '1',
    #                             'N-ACETYLGLUTPREDUCT-RXN': '6',
    #                             'ACETYLGLUTKIN-RXN': '6',
    #                             'ACETYLORNTRANSAM-RXN': '8',
    #                             'N-ACETYLTRANSFER-RXN': '0'},
    #                 'rxns': {'ACETYLORNDEACET-RXN': {},
    #                          'N-ACETYLGLUTPREDUCT-RXN': {'Nmar_1289': True,
    #                                                      'Nlim_1024': True,
    #                                                      'BD31_I1263': True,
    #                                                      'NKOR_07275': True,
    #                                                      'NSED_07155': True},
    #                         'ACETYLGLUTKIN-RXN': {'CSUB_C0425': True,
    #                                               'Nmar_1290': True,
    #                                               'NSED_07160': True,
    #                                               'BD31_I1264': True,
    #                                               'NKOR_07280': True},
    #                         'ACETYLORNTRANSAM-RXN': {'Ngar_c13370': True,
    #                                                  'Nlim_1026': True,
    #                                                  'Nmar_1291': True,
    #                                                  'CSUB_C1307': True,
    #                                                  'NKOR_07285': True,
    #                                                  'MY1_1445': True,
    #                                                  'NSED_07165': True},
    #                         'N-ACETYLTRANSFER-RXN': {}
    #                     }
    #                 }
    
    pathways={} # pathway -> rxns -> orfs
    read_pathways_list(opts.pathways_file, pathways)
    
    # put taxonomy from ORF annotation
    add_taxonomy_information(pathways, enzymes)
    
    distrib_pathways = {}
    for p in pathways:  
       value =  compute_min_species(pathways, p, lca) # perform integer optimization
       distrib_pathways[p] = value

    ranked_pwys = sorted(distrib_pathways, key=distrib_pathways.get, reverse=True)

    for p in ranked_pwys: 
      print p + '\t' + str(distrib_pathways[p]) +  '\t' + pathways[p]['common_name']


    print_taxonomy_information(pathways,ranked_pwys, enzymes, lca)

    #print distrib_pathways 
    sys.exit();

    
    list=[]
    #read_list_file_to_list(opts.data_file,list)
    indlist = lca.get_independent_taxons(list)
    print indlist



# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

