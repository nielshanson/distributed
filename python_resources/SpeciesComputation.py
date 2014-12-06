#!/usr/bin/python

import re
import sys
from metapaths_utils import printf, eprintf

class SpeciesComputation:

    begin_pattern = re.compile("#")

    name_to_id={}
    id_to_name={}
    taxid_to_ptaxid = {}

    pid_to_cid = {}

    tree={}
    leafcount=0

    seen_tids = {}


    Taxons = {}
    Children = {}

    def __init__(self, filename):
       taxonomy_file = open(filename, 'r')
       lines = taxonomy_file.readlines()
       taxonomy_file.close()

       for line in lines:
          if self.begin_pattern.search(line):
              continue
          fields =  [ str(x.strip())  for x in line.rstrip().split('\t')]
          if len(fields) !=3:
              continue
          self.name_to_id[fields[0]] = fields[1]
          self.id_to_name[fields[1]] = fields[0]
          self.taxid_to_ptaxid[fields[1]] = [ fields[2], 0]
         
          if not fields[2] in self.pid_to_cid:
             self.pid_to_cid[fields[2]] = {}
          self.pid_to_cid[fields[2]][fields[1]] = 1
  
#       print self.pid_to_cid
          
       
    def sizeTaxnames(self ):
         return len(self.name_to_id)

    def sizeTaxids(self):
         return len(self.taxid_to_ptaxid)
          
    def get_a_Valid_ID(self, name_group):
        for name in name_group:
           if name in self.name_to_id:
               return  self.name_to_id[name]
        return -1


    def get_lca(self, IDs):
        limit = len(IDs)
        for id in IDs:
           tid = id 
           while( tid in self.taxid_to_ptaxid and tid !='1' ):
               self.taxid_to_ptaxid[tid][1]+=1
               if self.taxid_to_ptaxid[tid][1]==limit:
                  return  self.id_to_name[tid]  
               tid = self.taxid_to_ptaxid[tid][0]
        return ""

    def clear_cells(self, IDs):
        limit = len(IDs)
        for id in IDs:
           tid = id 
           while( tid in self.taxid_to_ptaxid and tid !='1' ):
               if self.taxid_to_ptaxid[tid][1]==0:
                  return  self.id_to_name[tid]  
               self.taxid_to_ptaxid[tid][1]=0
               tid = self.taxid_to_ptaxid[tid][0]
        return ""

    def get_lineage(self, name):
        id = self.get_a_Valid_ID([ name ])
        tid = id
        while( tid in self.taxid_to_ptaxid and tid !='1' ):
             self.taxid_to_ptaxid[tid][1]+=1
             tid = self.taxid_to_ptaxid[tid][0]
             if tid in  self.id_to_name:  
                print '      ' +  self.id_to_name[tid]  
             else:
                return ""
        return ""

    def grow_tree(self, tree, tid, indent):
        if not tid  in self.pid_to_cid:
           return
        for c in self.pid_to_cid[tid]:
           tree[c] = {}
           self.leafcount += 1
#           if self.leafcount %1000 == 0 :
#              eprintf('Count = %d\n',self.leafcount)
 
#           for i in range(indent):
#              printf('\t')
#           if c in self.id_to_name:
#              printf('%s\n',self.id_to_name[c])
#           else:
#              printf('%s\n',c)

           if not c in self.seen_tids:
              self.seen_tids[c] = 1
              self.grow_tree(tree[c],c, indent+1)
           else:
              if c in self.id_to_name:
                   print 'culprit\'s ' + c + ' ' + self.id_to_name[c]

 
    def construct_tree(self):
        self.tree['1'] = {}
        indent = 0
        self.seen_tids = {}
        printf('%s\n',self.id_to_name['1'])
        self.seen_tids['1'] = 1
        self.grow_tree(self.tree['1'],'1', indent+1)
        
        print "Leaf count " + str(self.leafcount)
        #self.print_tree(self.tree, 0)


    def _parent(self, tax):
       if tax in  self.name_to_id:
          id = self.name_to_id[tax]
       else:
          return None

       if id in self.taxid_to_ptaxid:
          pid = self.taxid_to_ptaxid[id][0]
       else:
          return None

       if pid in self.id_to_name:
          return self.id_to_name[pid]

       return None



    def _add_distinct_children(self, tax):
        x = tax
        while x!=None and x != 'all':
          x  = self._parent(x)
          if x:
             #print x
             if x!=tax and x in self.Taxons:
               self.Taxons[x] = False


    def _mark_children(self, tax): 
        x = tax
        while x!=None and x != 'all':
          if x in self.Taxons and x != tax:
             if not x in self.Children:
                self.Children[x]={}
             self.Children[x][tax] = 1
          x  = self._parent(x)

    def build_independent_taxons(self, taxons_list):
        self.Taxons={}
        self.Children={}
        for tax in taxons_list:
           self.Taxons[tax] = True

        for tax in self.Taxons:
           self._add_distinct_children(tax)

        for tax in self.Taxons:
            if self.Taxons[tax]==True:
               self._mark_children(tax)
               if not tax in self.Children:
                 self.Children[tax] ={tax:1}

    def get_independent_parent(self, tax):
        if tax in self.Children:
           return self.Children[tax]
        return None

    def get_independent_taxons(self, taxons_list):
        self.build_independent_taxons(taxons_list)

        inlist = []
        for tax in self.Taxons:
            if self.Taxons[tax]:
               inlist.append(tax)

        return inlist
        


    def print_tree(self, tree, indent):
        for c in tree:
            if c in self.id_to_name:
              for i in range(indent):
                 printf('\t')
              printf('%s\n',self.id_to_name[c])
            self.print_tree(tree[c], indent+1)
        

    def getTaxonomy(self, name_groups):

         IDs = []
         for name_group in name_groups:
            #print name_group
            id = self.get_a_Valid_ID(name_group)
            if id!=-1:
              IDs.append(id)
    
         consensus = self.get_lca(IDs)
         #print "=============="
         self.clear_cells(IDs)
         return consensus

