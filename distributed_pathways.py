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

     from SpeciesComputation import *
     from metapaths_utils  import  fprintf, printf, eprintf
     #from libs.python_modules.sysutil import getstatusoutput
except:
     print """ Could not load some user defined  module functions"""
     print """ """
     sys.exit(3)


objective = re.compile(r'Objective.*=\s*(\d*)')


usage=  sys.argv[0] + " --ncbi-file ncbi_taxfile --pathways-file pathways_list  --orf-taxonomy orf-tax-file --equiv-file equivalencefile"""

glpsol = "/usr/local/bin/glpsol"

parser = OptionParser(usage)
parser.add_option("--ncbi-file", dest="ncbi_file",
                  help='name of the refseq names file [OPTIONAL]')

#parser.add_option("--data-file", dest="data_file",
#                  help='name of the test strings')

parser.add_option("--pathways-file", dest="pathways_file",
                  help='name of the pathways file')

parser.add_option("--orf-taxonomy", dest="orf_tax_file",
                  help='map of the enzymes to taxons file')

parser.add_option("--equiv-file", dest="orf_equiv_file",
                  help='orf equivalence by function file [OPTIONAL]')




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

def read_orf_equivalence(orf_equiv_filename, orf_equivalence):
    orfequivfile = open(orf_equiv_filename, 'r')
    lines = orfequivfile.readlines()
    orfequivfile.close()
  
    for line in lines:
       line = line.strip()
       if line:
          fields = [ x.strip() for x in line.split('\t') ]
          if len(fields) > 0:
             orf_equivalence[fields[0]] = []
          for i in range(1, len(fields)):
             orf_equivalence[fields[0]].append(fields[i])
             


def read_enzymes_list(enzymes_file, enzymes, enzymesmap):
    enzymesfile = open(enzymes_file, 'r')
    lines = enzymesfile.readlines()
    enzymesfile.close()
    seen={} 
    count = 0
    for line in lines:
       line = line.strip()
       if line:
          fields = [ x.strip() for x in line.split('\t') ]
          if len(fields) >= 3:
#             enzyme_name = fields[1]
#             taxonomy = fields[2]
# change to accomodate 

             enzyme_name = fields[0]
             taxonomy = fields[3]


             if taxonomy in seen: 
                newtax = seen[taxonomy]
             else:  
                count += 1
                newtax = 'E'+str(count)
                enzymesmap[newtax] = taxonomy
                seen[taxonomy] = newtax

             enzymes[enzyme_name] = newtax
            

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



def read_pathways_list(pathways_file, pathways, orf_equivalence):
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
             pathways[ pathway_name] = { 'common_name': pathway_common_name, 'rxns':{} }
       

          if fields[0] == 'RXN:':
             rxn_name = fields[1]
             rxn_common_name = fields[2]
             pathways[pathway_name]['rxns'][rxn_name] = {}

             if len(fields) > 5:
                for i in  range(5,len(fields)):
                   enzyme = fields[i]
                   pathways[pathway_name]['rxns'][rxn_name][enzyme] = True 
                   if enzyme in orf_equivalence:
                      for equivenzyme in orf_equivalence[enzyme]:
                         pathways[pathway_name]['rxns'][rxn_name][equivenzyme] = True 

                   

def print_taxonomy_information(pathways,ranked_pathways, enzymes, enzymesmap):
    for  key in ranked_pathways:
        taxons = {}
        for rxn in pathways[key]['rxns']:   
           for enzyme in pathways[key]['rxns'][rxn]:   
                taxons[pathways[key]['rxns'][rxn][enzyme]] = True

        #lca.build_independent_taxons(taxons)

        for rxn in pathways[key]['rxns']:   
           print '   ' + rxn
           for enzyme in pathways[key]['rxns'][rxn]:   
               printf("\t%s\t%s\n",  enzyme, enzymesmap[pathways[key]['rxns'][rxn][enzyme]])
             
       #      print '      '+ enzyme + '  ' + pathways[key]['rxns'][rxn][enzyme]
def add_taxonomy_information(pathways, enzymes):
    for  key in pathways:
        for rxn in pathways[key]['rxns']:   
           for enzyme in pathways[key]['rxns'][rxn]:   
             if enzyme in enzymes:
                pathways[key]['rxns'][rxn][enzyme] = enzymes[enzyme]   
             
       #      print '      '+ enzyme + '  ' + pathways[key]['rxns'][rxn][enzyme]
       # indlist = lca.get_independent_taxons(taxons)
       # for ind in indlist:
       #     print '             ' + ind
         
   
def check_arguments(opts):
    if opts.pathways_file == None:
       print usage
       sys.exit(0)


def compute_min_species(pathways,p): #, lca):
    taxons = {}
    for rxn in pathways[p]['rxns']:   
       for enzyme in pathways[p]['rxns'][rxn]:   
            taxons[pathways[p]['rxns'][rxn][enzyme]] = True

    #indtaxons = lca.get_independent_taxons(taxons)
    indtaxons = taxons.keys()
    
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


    #create X variables 
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
              #for S in lca.get_independent_parent(pathways[p]['rxns'][rxn][enzyme]):
              S =pathways[p]['rxns'][rxn][enzyme]
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

    namesdict = {}
    read_taxonomy_list(opts.ncbi_file, namesdict)
#    file = 'ncbi_taxonomy_tree.txt'
#    lca = SpeciesComputation(file)

##    lca.construct_tree()

    orf_equivalence={}
    if opts.orf_equiv_file:
       read_orf_equivalence(opts.orf_equiv_file, orf_equivalence)
       print "read orf equivalence"

    enzymes={}
    enzymesmap={}
    read_enzymes_list(opts.orf_tax_file, enzymes, enzymesmap)
    print "read enzyme list"

    pathways={}
    read_pathways_list(opts.pathways_file, pathways, orf_equivalence)
    print "read pathway  list"

    add_taxonomy_information(pathways, enzymes)
    print "added taxonony"

    distrib_pathways = {}

    count = 0
    for p in pathways:  
       count+=1
#       eprintf("%d\n",str(count))
       value =  compute_min_species(pathways,p)
       #value =  compute_min_species(pathways,p, lca)
       distrib_pathways[p] = value
    print "computed distributed taxonony"

    ranked_pwys = sorted(distrib_pathways, key=distrib_pathways.get, reverse=True)

    for p in ranked_pwys: 
       print p + '\t' + str(distrib_pathways[p]) + '\t' + pathways[p]['common_name']


    print_taxonomy_information(pathways,ranked_pwys, enzymes, enzymesmap)

    #print distrib_pathways 
    sys.exit();

    
    list=[]
    #read_list_file_to_list(opts.data_file,list)
    indlist = lca.get_independent_taxons(list)
    print indlist



# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

