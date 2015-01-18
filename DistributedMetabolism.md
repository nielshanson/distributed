Distributed Metabolism Notes
----------------------------

* Taking a look at the distributed metabolism folder.

	```  
	python distributed_pathways.v2.py 
	--ncbi-file ncbi_taxonomy_tree.txt 
	--pathways-file BD.pathways2.txt 
	--enzymes-file DISTRIB_CAM/case1/output/BD/results/annotation_table/functional_and_taxonomic_table.txt 
	```

* had to install 'glpk' from homebrew

	```
	brew install homebrew/science/glpk
	```
* got this to work

```
python distributed_pathways.v1.py --ncbi-file ncbi_taxonomy_tree.txt --pathways-file Thaumarchaea_comparison/thamurchaea_pwy.txt --enzymes-file Thaumarchaea_comparison/output/fasta/blast/Thamurchaea_all/results/annotation_table/functional_and_taxonomic_table.txt
```

# Test Cases

* test cases that we have

## DISTRIB_CAM

We have the following samples that we have combined:

* (B) Bacillus subtilis subsp. subtilis str. 168 complete genome (AL009126.3)
* (D) Desulfovibrio vulgaris subsp. vulgaris str. Hildenborough, complete genome (AE017285.1, AE017286.1)
* (M) Methanococcus maripaludis strain S2 (BX950229.1)

and we have combined these in the following ways BD, BM, MD, BDM


* get started tomorrow:

```
cd Dropbox/projects/distributed
python distributed_pathways.v1.py --ncbi-file resources/ncbi_taxonomy_tree.txt --pathways-file data/Thaumarchaea_comparison/thamurchaea_pwy.txt --enzymes-file data/Thaumarchaea_comparison/output/fasta/blast/Thamurchaea_all/results/annotation_table/functional_and_taxonomic_table.txt
```
