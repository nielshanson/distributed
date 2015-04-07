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

# Experiments

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

## mor_trem

Symbionts of Mor and Trem from MetaPathways.

* Extract pathway_reaction information from PGDB using a modified version of `extract_pathways_from_pgdb.py`

```
python extract_pathways_from_pgdb.py -p M_T_PCIT --micro_dex -o exp/mor_trem/M_T_PCIT.rxn.pwy.txt
```

* Run MicroDEX on the distributed pathways

```
python distributed_pathways.v1.py --ncbi-file resources/ncbi_taxonomy_tree.txt --pathways-file exp/mor_trem/M_T_PCIT.rxn.pwy.txt --enzymes-file exp/mor_trem/functional_and_taxonomic_table.txt > exp/mor_trem/mor_trem_initial_results.txt
```

## atcc

Pull ATCC genomes from the NCBI.

* Obtain ATCC genomes by running the utility script `pull_ATCC_NCBI.py`

```
./ATCC_run.sh
```

* Downloaded 507 genomes

```
/Users/nielsh/Downloads/metapathways2/MetaPathways.py -v -i /Users/nielsh/Dropbox/projects/distributed/data/ATCC_genomes -o /Users/nielsh/Dropbox/projects/distributed/data/mp_output/ATCC_genomes -p /Users/nielsh/Downloads/metapathways2/config/template_param.txt -c /Users/nielsh/Downloads/metapathways2/config/template_config.txt -d 8 -s ATCC --runid 14:17:37
```

