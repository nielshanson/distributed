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