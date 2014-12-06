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
     from metapaths_utils  import  fprintf, printf
     #from libs.python_modules.sysutil import getstatusoutput
except:
     print """ Could not load some user defined  module functions"""
     print """ """
     sys.exit(3)


objective = re.compile(r'Objective.*=\s*(\d*)')


usage=  sys.argv[0] + " --ncbi-file ncbi_taxfile --pathways-file pathways_list  --enzymes-file enzymes_file"""

glpsol = "/usr/local/bin/glpsol"

parser = OptionParser(usage = usage, version="v2 (can compute distributedness of composite pathways")
parser.add_option("--ncbi-file", dest="ncbi_file",
                  help='name of the refseq names file')

#parser.add_option("--data-file", dest="data_file",
#                  help='name of the test strings')
parser.add_option("--composite-pathways", dest="composite_pathways_file",
                  help='composite pathways file')

parser.add_option("--pathways-file", dest="pathways_file",
                  help='name of the pathways file as extracted by \"extract_pathway_table_from_pgdb.v2.pl\" ')

parser.add_option("--enzymes-file", dest="enzymes_file",
                  help='name of the enzymes/taxons file. It is expected in the format of the \"functional_taxonomic_table.txt\"')




commonnames = {}
commonnames["PWY-6240"] ="starch plastidic translocation I"
commonnames["PWY-6238"] ="nucleobase ascorbate transport II"
commonnames["PWY-6237"] ="nucleobase ascorbate transport I"
commonnames["PWY-6229"] ="nucleobase derivative transport"
commonnames["PWY-6213"] ="cadmium transport I"
commonnames["PWY-6188"] ="calmodulin regulated calcium transport"
commonnames["PWY-6171"] ="calcium transport II"
commonnames["PWY-6166"] ="calcium transport I"
commonnames["PWY-6162"] ="HAK potassium transport"
commonnames["PWY-6161"] ="sodium transport I"
commonnames["PWY-6159"] ="potassium transport IV"
commonnames["PWY-6156"] ="potassium transport III"
commonnames["PWY-6150"] ="potassium transport II"
commonnames["PWY-6149"] ="potassium transport I"
commonnames["PWY-6137"] ="copper transport II"
commonnames["PWY-6135"] ="copper transport I"
commonnames["PWY-5937"] ="Fe(III)-phytosiderophore transport"
commonnames["PWY-5934"] ="Fe(III)-reduction and Fe(II) transport"
commonnames["TCA"] ="TCA cycle I (prokaryotic)"
commonnames["REDCITCYC"] ="TCA cycle III (helicobacter)"
commonnames["PWY-6969"] ="TCA cycle V (2-oxoglutarate:ferredoxin oxidoreductase)"
commonnames["PWY-5913"] ="TCA cycle VI (obligate autotrophs)"
commonnames["PWY-5690"] ="TCA cycle II (eukaryotic)"
commonnames["P105-PWY"] ="TCA cycle IV (2-oxoglutarate decarboxylase)"
commonnames["GLYOXYLATE-BYPASS"] ="glyoxylate cycle"
commonnames["PWY-6530"] ="perchlorate reduction"
commonnames["PWY-6529"] ="chlorate reduction"
commonnames["PWY-4601"] ="arsenate reduction (respiratory)"
commonnames["ANARESP1-PWY"] ="respiration (anaerobic)"
commonnames["PWY-4302"] ="aerobic respiration (alternative oxidase pathway)"
commonnames["PWY-3781"] ="aerobic respiration (cytochrome c)"
commonnames["PWY-241"] ="C4 photosynthetic carbon assimilation cycle"
commonnames["PWY-181"] ="photorespiration"
commonnames["PWY-101"] ="photosynthesis light reactions"
commonnames["P21-PWY"] ="pentose phosphate pathway (partial)"
commonnames["OXIDATIVEPENT-PWY"] ="pentose phosphate pathway (oxidative branch)"
commonnames["NONOXIPENT-PWY"] ="pentose phosphate pathway (non-oxidative branch)"
commonnames["PWY66-368"] ="ketolysis"
commonnames["PWY66-367"] ="ketogenesis"
commonnames["PWY-6118"] ="glycerol-3-phosphate shuttle"
commonnames["PWY-5261"] ="methanogenesis from tetramethylammonium"
commonnames["PWY-5259"] ="methanogenesis from methanethiol"
commonnames["PWY-5250"] ="methanogenesis from trimethylamine"
commonnames["PWY-5248"] ="methanogenesis from dimethylamine"
commonnames["PWY-5247"] ="methanogenesis from methylamine"
commonnames["PWY-5209"] ="methyl-coenzyme M oxidation to CO<sub>2</sub>"
commonnames["METHFORM-PWY"] ="methyl-coenzyme M reduction to methane"
commonnames["METHANOGENESIS-PWY"] ="methanogenesis from CO<SUB>2</SUB>"
commonnames["METH-ACETATE-PWY"] ="methanogenesis from acetate"
commonnames["PWY-6785"] ="hydrogen production VIII"
commonnames["PWY-6780"] ="hydrogen production VI"
commonnames["PWY-6772"] ="hydrogen production V"
commonnames["PWY-6765"] ="hydrogen production IV"
commonnames["PWY-6759"] ="hydrogen production III"
commonnames["PWY-6758"] ="hydrogen production II"
commonnames["PWY-6744"] ="hydrogen production I"
commonnames["PWY-5484"] ="glycolysis II (from fructose-6P)"
commonnames["PWY-1042"] ="glycolysis IV (plant cytosol)"
commonnames["P341-PWY"] ="glycolysis V (Pyrococcus)"
commonnames["GLYCOLYSIS"] ="glycolysis I"
commonnames["ANAGLYCOLYSIS-PWY"] ="glycolysis III (glucokinase)"
commonnames["PWY-6588"] ="pyruvate fermentation to acetone"
commonnames["PWY-6587"] ="pyruvate fermentation to ethanol III"
commonnames["PWY-5768"] ="pyruvate fermentation to acetate VIII"
commonnames["PWY-5494"] ="pyruvate fermentation to propionate II (acrylate pathway)"
commonnames["PWY-5486"] ="pyruvate fermentation to ethanol II"
commonnames["PWY-5482"] ="pyruvate fermentation to acetate II"
commonnames["PWY-5481"] ="pyruvate fermentation to lactate"
commonnames["PWY-5480"] ="pyruvate fermentation to ethanol I"
commonnames["P108-PWY"] ="pyruvate fermentation to propionate I"
commonnames["CENTFERM-PWY"] ="pyruvate fermentation to butanoate"
commonnames["PWY-6392"] ="<i>meso</i>-butanediol biosynthesis II"
commonnames["PWY-6391"] ="<i>meso</i>-butanediol biosynthesis I"
commonnames["PWY-6390"] ="(<i>S,S</i>)-butanediol biosynthesis"
commonnames["PWY-5951"] ="(<i>R,R</i>)-butanediol biosynthesis"
commonnames["PWY3O-440"] ="acetoin biosynthesis III"
commonnames["PWY-6389"] ="(<i>S</i>)-acetoin biosynthesis"
commonnames["PWY-5939"] ="(<i>R</i>)-acetoin biosynthesis II"
commonnames["PWY-5938"] ="(<i>R</i>)-acetoin biosynthesis I"
commonnames["PWY-6863"] ="pyruvate fermentation to hexanol"
commonnames["PWY-6583"] ="pyruvate fermentation to butanol I"
commonnames["PWY-5677"] ="succinate fermentation to butyrate"
commonnames["P122-PWY"] ="heterolactic fermentation"
commonnames["FERMENTATION-PWY"] ="mixed acid fermentation"
commonnames["PYRUVOX-PWY"] ="pyruvate oxidation pathway"
commonnames["PWY0-1356"] ="formate to dimethyl sulfoxide electron transfer"
commonnames["PWY0-1355"] ="formate to trimethylamine N-oxide electron transfer"
commonnames["PWY0-1353"] ="succinate to cytochrome <i>bd</i> oxidase electron transfer"
commonnames["PWY0-1348"] ="NADH to dimethyl sulfoxide electron transfer"
commonnames["PWY0-1347"] ="NADH to trimethylamine N-oxide electron transfer"
commonnames["PWY0-1336"] ="NADH to fumarate electron transfer"
commonnames["PWY0-1335"] ="NADH to cytochrome <i>bo</i> oxidase electron transfer"
commonnames["PWY0-1334"] ="NADH to cytochrome <i>bd</i> oxidase electron transfer"
commonnames["PWY0-1329"] ="succinate to cytochrome <i>bo</i> oxidase electron transfer"
commonnames["PWY-6512"] ="hydrogen oxidation III (anaerobic, NADP)"
commonnames["PWY-5382"] ="hydrogen oxidation II (aerobic, NAD)"
commonnames["P283-PWY"] ="hydrogen oxidation I (aerobic)"
commonnames["P301-PWY"] ="galena oxidation"
commonnames["P282-PWY"] ="nitrite oxidation"
commonnames["PWY-6970"] ="acetyl-CoA biosynthesis II (NADP-dependent pyruvate dehydrogenase)"
commonnames["PWY-5172"] ="acetyl-CoA biosynthesis III (from citrate)"
commonnames["PWY0-1517"] ="sedoheptulose bisphosphate bypass"
commonnames["PWY-6938"] ="NADH repair"
commonnames["PWY-6883"] ="pyruvate fermentation to butanol II"
commonnames["PWY-6876"] ="isopropanol biosynthesis"
commonnames["PWY-6871"] ="3-methylbutanol biosynthesis"
commonnames["PWY-6859"] ="<i>all-trans</i>-farnesol biosynthesis"
commonnames["PWY-6728"] ="methylaspartate cycle"
commonnames["PWY-5741"] ="ethylmalonyl pathway"
commonnames["PWY-5723"] ="Rubisco shunt"
commonnames["PWY-5197"] ="lactate biosynthesis (archaea)"
commonnames["ENTNER-DOUDOROFF-PWY"] ="Entner-Doudoroff pathway I"
commonnames["PWY-6873"] ="long chain fatty acid ester synthesis for microdiesel production"
commonnames["P641-PWY"] ="phenylmercury acetate degradation"
commonnames["PWY-6726"] ="cyanide detoxification I"
commonnames["PWY-6421"] ="arsenate detoxification III (mycothiol)"
commonnames["PWY-4621"] ="arsenate detoxification II (glutaredoxin)"
commonnames["PWY-4202"] ="arsenate detoxification I (glutaredoxin)"
commonnames["PWY-4521"] ="arsenite oxidation (respiratory)"
commonnames["PWY-6828"] ="linezolid resistance"
commonnames["PWY0-1305"] ="glutamate dependent acid resistance"
commonnames["PWY0-1299"] ="arginine dependent acid resistance"
commonnames["PWY-6900"] ="nectroscordum lachrymator biosynthesis"
commonnames["PWY-6842"] ="glutathione-mediated detoxification II"
commonnames["PWY-6786"] ="detoxification of reactive carbonyls in chloroplasts"
commonnames["PWY-6646"] ="fluoroacetate degradation"
commonnames["PWY-6502"] ="oxidized GTP and dGTP detoxification"
commonnames["PWY-2261"] ="ascorbate glutathione cycle"
commonnames["DETOX1-PWY"] ="superoxide radicals degradation"
commonnames["PWY-6946"] ="cholesterol degradation to androstenedione II (cholesterol dehydrogenase)"
commonnames["PWY-6945"] ="cholesterol degradation to androstenedione I (cholesterol oxidase)"
commonnames["PWY-6948"] ="sitosterol degradation to androstenedione"
commonnames["PWY-6944"] ="androstenedione degradation"
commonnames["PWY-6943"] ="testosterone and androsterone degradation to androstendione"
commonnames["PWY-6526"] ="limonene degradation III (to perillate)"
commonnames["PWY-5924"] ="limonene degradation II (L-limonene)"
commonnames["PWY-5923"] ="limonene degradation I (D-limonene)"
commonnames["PWY-5927"] ="(4<i>S</i>)-carveol and (4<i>S</i>)-dihydrocarveol degradation"
commonnames["PWY-5922"] ="(4<i>R</i>)-carveol and (4<i>R</i>)-dihydrocarveol degradation"
commonnames["PWY-6678"] ="geraniol and nerol degradation"
commonnames["PWY-6672"] ="<i>cis</i>-genanyl-CoA degradation"
commonnames["PWY-6670"] ="citronellol degradation"
commonnames["PWY-6413"] ="ginsenoside degradation III"
commonnames["PWY-6412"] ="ginsenoside degradation II"
commonnames["PWY-6411"] ="ginsenoside degradation I"
commonnames["P601-PWY"] ="D-camphor degradation"
commonnames["6-HYDROXYCINEOLE-DEGRADATION-PWY"] ="1,8-cineole degradation"
commonnames["PWY-6806"] ="carotenoid cleavage dioxygenases"
commonnames["SORBDEG-PWY"] ="sorbitol degradation II"
commonnames["RIBITOLUTIL-PWY"] ="ribitol degradation"
commonnames["PWY-6531"] ="mannitol cycle"
commonnames["PWY-4101"] ="sorbitol degradation I"
commonnames["P562-PWY"] ="<I>myo</I>-inositol degradation"
commonnames["MANNIDEG-PWY"] ="mannitol degradation I"
commonnames["LARABITOLUTIL-PWY"] ="xylitol degradation"
commonnames["GALACTITOLCAT-PWY"] ="galactitol degradation"
commonnames["DARABITOLUTIL-PWY"] ="D-arabitol degradation"
commonnames["PWY-6501"] ="D-glucuronate degradation II"
commonnames["PWY-5525"] ="D-glucuronate degradation I"
commonnames["GLUCUROCAT-PWY"] ="&beta;-D-glucuronide and D-glucuronate degradation"
commonnames["PWY-6499"] ="D-glucarate degradation II"
commonnames["GLUCARDEG-PWY"] ="<i>D</i>-glucarate degradation I"
commonnames["PWY-6491"] ="D-galacturonate degradation III"
commonnames["PWY-6486"] ="D-galacturonate degradation II"
commonnames["GALACTUROCAT-PWY"] ="D-galacturonate degradation I"
commonnames["PWY-6497"] ="D-galactarate degradation II"
commonnames["GALACTARDEG-PWY"] ="D-galactarate degradation I"
commonnames["PWY0-1306"] ="L-galactonate degradation"
commonnames["IDNCAT-PWY"] ="L-idonate degradation"
commonnames["GALACTCAT-PWY"] ="D-galactonate degradation"
commonnames["PWY-4781"] ="phytate degradation II"
commonnames["PWY-4702"] ="phytate degradation I"
commonnames["PWY0-521"] ="fructoselysine and psicoselysine degradation"
commonnames["PWY0-1261"] ="1,6-anhydro-<i>N</i>-acetylmuramic acid recycling"
commonnames["PWY-6507"] ="5-dehydro-4-deoxy-D-glucuronate degradation"
commonnames["PWY-6667"] ="resveratrol degradation"
commonnames["PWY-6684"] ="glucosinolate breakdown (via thiocyanate-forming protein)"
commonnames["PWY-6011"] ="amygdalin and prunasin degradation"
commonnames["PWY-6002"] ="lotaustralin degradation"
commonnames["PWY-5976"] ="dhurrin degradation"
commonnames["PWY-3121"] ="linamarin degradation"
commonnames["PWY-6633"] ="caffeine degradation V (bacteria, to trimethylurate)"
commonnames["PWY-6553"] ="caffeine degradation II"
commonnames["PWY-6552"] ="caffeine degradation I (main, plants)"
commonnames["PWY-6538"] ="caffeine degradation III (bacteria, via demethylation)"
commonnames["PWY-5461"] ="betanidin degradation"
commonnames["PWY-5708"] ="ethiin degradation"
commonnames["PWY-5706"] ="alliin degradation"
commonnames["PWY-31"] ="canavanine degradation"
commonnames["PWY-6848"] ="rutin degradation"
commonnames["PWY0-1527"] ="curcumin degradation"
commonnames["PWY-4441"] ="DIMBOA-glucoside degradation"
commonnames["PWY-6423"] ="hemoglobin degradation"
commonnames["PWY-6018"] ="seed germination protein turnover"
commonnames["PWY-5988"] ="wound-induced proteolysis I"
commonnames["PWY-6341"] ="guaiacylglycerol-&beta;-guaiacyl ether degradation"
commonnames["PWY-743"] ="thiocyanate degradation II"
commonnames["P581-PWY"] ="thiocyanate degradation I"
commonnames["PWY66-221"] ="nicotine degradation III"
commonnames["PWY66-201"] ="nicotine degradation II"
commonnames["P181-PWY"] ="nicotine degradation I"
commonnames["PWY-5336"] ="carbon disulfide oxidation II (aerobic)"
commonnames["PWY-1164"] ="carbon disulfide oxidation I (anaerobic)"
commonnames["PWY66-241"] ="bupropion degradation"
commonnames["PWY1G-1"] ="mycothiol-mediated detoxification"
commonnames["PWY-723"] ="alkylnitronates degradation"
commonnames["PWY-6518"] ="glycocholate metabolism (bacteria)"
commonnames["PWY-6377"] ="&alpha;-tocopherol degradation"
commonnames["PWY-5744"] ="glyoxylate assimilation"
commonnames["PWY-5355"] ="nitroethane degradation"
commonnames["P621-PWY"] ="nylon-6 oligomer degradation"
commonnames["P482-PWY"] ="arsonoacetate degradation"
commonnames["P481-PWY"] ="adamantanone degradation"
commonnames["P401-PWY"] ="cyanide degradation"
commonnames["P344-PWY"] ="acrylonitrile degradation"
commonnames["P221-PWY"] ="octane oxidation"
commonnames["P201-PWY"] ="nitroglycerin degradation"
commonnames["7ALPHADEHYDROX-PWY"] ="cholate degradation (bacteria, anaerobic)"
commonnames["PWY-4061"] ="glutathione-mediated detoxification"
commonnames["P161-PWY"] ="acetylene degradation"
commonnames["THIOSULFOX-PWY"] ="thiosulfate oxidation I (to tetrathionate)"
commonnames["PWY-6677"] ="thiosulfate oxidation IV (multienzyme complex)"
commonnames["PWY-5303"] ="thiosulfate oxidation II (to tetrathionate)"
commonnames["PWY-5296"] ="thiosulfate oxidation III (multienzyme complex)"
commonnames["PWY-5352"] ="thiosulfate disproportionation II (non thiol-dependent)"
commonnames["PWY-5277"] ="thiosulfate disproportionation I (thiol-dependent)"
commonnames["PWY-5350"] ="thiosulfate disproportionation III (rhodanese)"
commonnames["PWY-5359"] ="tetrathionate reductiuon II (to trithionate)"
commonnames["PWY-5358"] ="tetrathionate reduction I (to thiosulfate)"
commonnames["PWY-5364"] ="sulfur reduction II (via polysulfide)"
commonnames["PWY-5332"] ="sulfur reduction I"
commonnames["PWY-6675"] ="sulfur oxidation IV (intracellular sulfur)"
commonnames["SULFUROX-PWY"] ="sulfur oxidation I (aerobic)"
commonnames["FESULFOX-PWY"] ="sulfur oxidation II (Fe<sup>+3</sup>-dependent)"
commonnames["PWY-5302"] ="sulfur disproportionation II (aerobic)"
commonnames["P203-PWY"] ="sulfur disproportionation I (anaerobic)"
commonnames["PWY-6638"] ="sulfolactate degradation III"
commonnames["PWY-6637"] ="sulfolactate degradation II"
commonnames["PWY-6616"] ="sulfolactate degradation I"
commonnames["PWY-6718"] ="sulfoacetaldehyde degradation III"
commonnames["PWY-5982"] ="sulfoacetaldehyde degradation II"
commonnames["PWY-1281"] ="sulfoacetaldehyde degradation I"
commonnames["PWY-5326"] ="sulfite oxidation IV"
commonnames["PWY-5279"] ="sulfite oxidation II"
commonnames["PWY-5278"] ="sulfite oxidation III"
commonnames["PWY-5276"] ="sulfite oxidation I (sulfite oxidoreductase)"
commonnames["PWY-5274"] ="sulfide oxidation II (sulfide dehydrogenase)"
commonnames["PWY-5285"] ="sulfide oxidation III (sulfur dioxygenase)"
commonnames["P222-PWY"] ="sulfide oxidation I (sulfide-quinone reductase)"
commonnames["SULFMETII-PWY"] ="sulfate reduction II (assimilatory)"
commonnames["PWY-6683"] ="sulfate reduction III (assimilatory)"
commonnames["P224-PWY"] ="sulfate reduction V (dissimilatory)"
commonnames["DISSULFRED-PWY"] ="sulfate reduction IV (dissimilatory)"
commonnames["PWY-6048"] ="methylthiopropionate degradation I (cleavage)"
commonnames["PWY-6045"] ="methylthiopropionate degradation II (demethylation)"
commonnames["PWY-6056"] ="dimethylsulfoniopropionate degradation II (cleavage)"
commonnames["PWY-6052"] ="dimethylsulfoniopropionate degradation III (demethylation)"
commonnames["PWY-6046"] ="dimethylsulfoniopropionate degradation I (cleavage)"
commonnames["PWY-6059"] ="dimethyl sulfide degradation II (oxidation)"
commonnames["PWY-6057"] ="dimethyl sulfide degradation III (oxidation)"
commonnames["PWY-6047"] ="dimethyl sulfide degradation I"
commonnames["PWY-6736"] ="sulfur volatiles biosynthesis"
commonnames["PWY-6642"] ="(<i>R</i>)-cysteate degradation"
commonnames["PWY-6634"] ="2,3-dihydroxypropane-1-sulfonate degradation"
commonnames["PWY-6593"] ="sulfoacetate degradation"
commonnames["PWY-6327"] ="tetrathionate oxidation"
commonnames["PWY-6321"] ="homotaurine degradation"
commonnames["PWY-6058"] ="dimethyl sulfone degradation"
commonnames["PWY-6050"] ="dimethyl sulfoxide degradation"
commonnames["PWY-6044"] ="methanesulfonate degradation"
commonnames["PWY-6043"] ="ethanedisulfonate degradation"
commonnames["PWY-5707"] ="isoalliin degradation"
commonnames["PWY-5260"] ="methanogenesis from methylthiopropionate"
commonnames["PWY-5258"] ="methanogenesis from dimethylsulfide"
commonnames["PWY-2601"] ="isethionate degradation"
commonnames["ALKANEMONOX-PWY"] ="two-component alkanesulfonate monooxygenase"
commonnames["PWY-6832"] ="2-aminoethylphosphonate degradation II"
commonnames["PHOSPHONOTASE-PWY"] ="2-aminoethylphosphonate degradation I"
commonnames["PWY0-1533"] ="methylphosphonate degradation"
commonnames["PWY-6357"] ="phosphate utilization in cell wall regeneration"
commonnames["PWY-6348"] ="phosphate acquisition"
commonnames["PWY-5491"] ="diethylphosphate degradation"
commonnames["P483-PWY"] ="phosphonoacetate degradation"
commonnames["PWY490-3"] ="nitrate reduction VI (assimilatory)"
commonnames["PWY-6748"] ="nitrate reduction VII (denitrification)"
commonnames["PWY-5675"] ="nitrate reduction V (assimilatory)"
commonnames["PWY-5674"] ="nitrate reduction IV (dissimilatory)"
commonnames["PWY-381"] ="nitrate reduction II (assimilatory)"
commonnames["DENITRIFICATION-PWY"] ="nitrate reduction I (denitrification)"
commonnames["PWY0-1352"] ="nitrate reduction VIII (dissimilatory)"
commonnames["PWY0-1321"] ="nitrate reduction III (dissimilatory)"
commonnames["PWY-2242"] ="ammonia oxidation III"
commonnames["P303-PWY"] ="ammonia oxidation II (anaerobic)"
commonnames["AMMOXID-PWY"] ="ammonia oxidation I (aerobic)"
commonnames["PWY-6964"] ="ammonia assimilation cycle II"
commonnames["PWY-6963"] ="ammonia assimilation cycle I"
commonnames["PWY-6523"] ="intra-aerobic nitrite reduction"
commonnames["PWY-4984"] ="urea cycle"
commonnames["N2FIX-PWY"] ="nitrogen fixation"
commonnames["PWY-5986"] ="ammonium transport"
commonnames["CYANCAT-PWY"] ="cyanate degradation"
commonnames["PWY-6592"] ="manganese oxidation II"
commonnames["PWY-6591"] ="manganese oxidation I"
commonnames["PWY-6932"] ="selenate reduction"
commonnames["PWY-6692"] ="Fe(II) oxidation"
commonnames["PWY-5949"] ="Fe(II)-nicotianamine transport in phloem"
commonnames["PWY0-1471"] ="uracil degradation III"
commonnames["PWY-6426"] ="uracil degradation I (oxidative)"
commonnames["PWY0-1295"] ="pyrimidine ribonucleosides degradation I"
commonnames["PWY-6556"] ="pyrimidine ribonucleosides degradation II"
commonnames["PWY0-1298"] ="pyrimidine deoxyribonucleosides degradation"
commonnames["PWY-6430"] ="thymine degradation"
commonnames["PWY-6608"] ="guanosine nucleotides degradation III"
commonnames["PWY-6607"] ="guanosine nucleotides degradation I"
commonnames["PWY-6606"] ="guanosine nucleotides degradation II"
commonnames["SALVADEHYPOX-PWY"] ="adenosine nucleotides degradation II"
commonnames["PWY-6617"] ="adenosine nucleotides degradation III"
commonnames["PWY-6596"] ="adenosine nucleotides degradation I"
commonnames["PWY-5532"] ="adenosine nucleotides degradation IV"
commonnames["PWY0-1296"] ="purine ribonucleosides degradation to ribose-1-phosphate"
commonnames["PWY-6019"] ="pseudouridine degradation"
commonnames["PWY-5497"] ="purine nucleotides degradation IV (anaerobic)"
commonnames["P164-PWY"] ="purine nucleotides degradation III (anaerobic)"
commonnames["PWY0-1391"] ="<i>S</i>-methyl-5'-thioadenosine degradation IV"
commonnames["PWY-6756"] ="<i>S</i>-methyl-5'-thioadenosine degradation II"
commonnames["PWY-6754"] ="<i>S</i>-methyl-5'-thioadenosine degradation I"
commonnames["PWY-6753"] ="<i>S</i>-methyl-5'-thioadenosine degradation III"
commonnames["PWY0-1297"] ="purine deoxyribonucleosides degradation"
commonnames["PWY-6227"] ="purine uptake"
commonnames["PWY-6261"] ="thyroid hormone metabolism II (via conjugation and/or degradation)"
commonnames["PWY-6260"] ="thyroid hormone metabolism I (via deiodination)"
commonnames["PWY-2841"] ="cytokinins degradation"
commonnames["PWY-5811"] ="IAA degradation VII"
commonnames["PWY-5797"] ="IAA degradation VI"
commonnames["PWY-5788"] ="IAA degradation V"
commonnames["PWY-5784"] ="IAA conjugate biosynthesis II"
commonnames["PWY-2021"] ="IAA degradation IV"
commonnames["PWY-1981"] ="IAA degradation III"
commonnames["PWY-1962"] ="IAA degradation II"
commonnames["PWY-1961"] ="IAA degradation I"
commonnames["PWY-6400"] ="melatonin degradation III"
commonnames["PWY-6399"] ="melatonin degradation II"
commonnames["PWY-6398"] ="melatonin degradation I"
commonnames["PWY-6688"] ="thyronamine and iodothyronamine metabolism"
commonnames["PWY-6342"] ="noradrenaline and adrenaline degradation"
commonnames["PWY-6313"] ="serotonin degradation"
commonnames["PWY0-1337"] ="oleate &beta;-oxidation"
commonnames["PWY-6837"] ="fatty acid beta-oxidation V (unsaturated, odd number, di-isomerase-dependent)"
commonnames["PWY-5138"] ="fatty acid &beta;-oxidation IV (unsaturated, even number)"
commonnames["PWY-5137"] ="fatty acid &beta;-oxidation III (unsaturated, odd number)"
commonnames["PWY-5136"] ="fatty acid &beta;-oxidation II (core pathway)"
commonnames["PWY-2724"] ="fatty acid &omega;-oxidation"
commonnames["PWY-2501"] ="fatty acid &alpha;-oxidation"
commonnames["FAO-PWY"] ="fatty acid &beta;-oxidation I"
commonnames["PWY-5533"] ="acetone degradation II (to acetoacetate)"
commonnames["PWY-5451"] ="acetone degradation I (to methylglyoxal)"
commonnames["SOPHOROSYLOXYDOCOSANOATE-DEG-PWY"] ="sophorosyloxydocosanoate degradation"
commonnames["PWY6666-1"] ="anandamide degradation"
commonnames["PWY3DJ-11470"] ="sphingosine and sphingosine-1-phosphate metabolism"
commonnames["PWY-6483"] ="ceramide degradation"
commonnames["PWY-6368"] ="3-phosphoinositide degradation"
commonnames["PWY-6111"] ="mitochondrial L-carnitine shuttle pathway"
commonnames["LIPAS-PWY"] ="triacylglycerol degradation"
commonnames["ACETOACETATE-DEG-PWY"] ="acetoacetate degradation (to acetyl CoA)"
commonnames["LIPASYN-PWY"] ="phospholipases"
commonnames["XYLCAT-PWY"] ="xylose degradation I"
commonnames["PWY-6760"] ="xylose degradation III"
commonnames["PWY-5516"] ="xylose degradation II"
commonnames["TREDEGLOW-PWY"] ="trehalose degradation I (low osmolarity)"
commonnames["PWY0-1466"] ="trehalose degradation VI (periplasmic)"
commonnames["PWY0-1182"] ="trehalose degradation II (trehalase)"
commonnames["PWY-2723"] ="trehalose degradation V"
commonnames["PWY-2722"] ="trehalose degradation IV"
commonnames["PWY-2721"] ="trehalose degradation III"
commonnames["SUCUTIL-PWY"] ="sucrose degradation I"
commonnames["SUCROSEUTIL2-PWY"] ="sucrose degradation II"
commonnames["PWY66-373"] ="sucrose degradation V (mammalian)"
commonnames["PWY-621"] ="sucrose degradation III"
commonnames["PWY-5384"] ="sucrose degradation IV"
commonnames["LACTOSEUTIL-PWY"] ="lactose degradation II"
commonnames["BGALACT-PWY"] ="lactose degradation III"
commonnames["RHAMCAT-PWY"] ="L-rhamnose degradation I"
commonnames["PWY-6714"] ="L-rhamnose degradation III"
commonnames["PWY-6713"] ="L-rhamnose degradation II"
commonnames["PWY-5517"] ="L-arabinose degradation III"
commonnames["PWY-5515"] ="L-arabinose degradation II"
commonnames["ARABCAT-PWY"] ="L-arabinose degradation I"
commonnames["PWY-6693"] ="galactose degradation IV"
commonnames["PWY-6317"] ="galactose degradation I (Leloir pathway)"
commonnames["LACTOSECAT-PWY"] ="lactose and galactose degradation I"
commonnames["GALDEG-PWY"] ="galactose degradation II"
commonnames["PWY-5519"] ="D-arabinose degradation III"
commonnames["DARABCATK12-PWY"] ="D-arabinose degradation I"
commonnames["DARABCAT-PWY"] ="D-arabinose degradation II"
commonnames["RIBOKIN-PWY"] ="ribose degradation"
commonnames["PWY0-44"] ="D-allose degradation"
commonnames["PWY0-1314"] ="fructose degradation"
commonnames["PWY0-1301"] ="melibiose degradation"
commonnames["PWY0-1300"] ="2-<I>O</I>-&alpha;-mannosyl-D-glycerate degradation"
commonnames["PWY-3861"] ="mannitol degradation II"
commonnames["PWY-2221"] ="Entner-Doudoroff pathway III (semi-phosphorylative)"
commonnames["PWY-1081"] ="homogalacturonan degradation"
commonnames["P302-PWY"] ="L-sorbose degradation"
commonnames["P124-PWY"] ="Bifidobacterium shunt"
commonnames["NPGLUCAT-PWY"] ="Entner-Doudoroff pathway II (non-phosphorylative)"
commonnames["MANNCAT-PWY"] ="D-mannose degradation"
commonnames["MALTOSECAT-PWY"] ="maltose degradation"
commonnames["LYXMET-PWY"] ="L-lyxose degradation"
commonnames["GLUCOSE1PMETAB-PWY"] ="glucose and glucose-1-phosphate degradation"
commonnames["FUCCAT-PWY"] ="fucose degradation"
commonnames["DHGLUCONATE-PYR-CAT-PWY"] ="glucose degradation (oxidative)"
commonnames["PWY-6812"] ="xyloglucan degradation III (cellobiohydrolase)"
commonnames["PWY-6807"] ="xyloglucan degradation II (exoglucanase)"
commonnames["PWY-6791"] ="xyloglucan degradation I (endoglucanase)"
commonnames["PWY-6789"] ="(1,3)-&beta;-D-xylan degradation"
commonnames["PWY-6717"] ="(1,4)-&beta;-xylan degradation"
commonnames["PWY-6737"] ="starch degradation V"
commonnames["PWY-6735"] ="starch degradation IV"
commonnames["PWY-6731"] ="starch degradation III"
commonnames["PWY-6771"] ="rhamnogalacturonan type I degradation II (bacteria)"
commonnames["PWY-6769"] ="rhamnogalacturonan type I degradation I (fungi)"
commonnames["PWY-6576"] ="dermatan sulfate degradation (metazoa)"
commonnames["PWY-6573"] ="chondroitin sulfate degradation (metazoa)"
commonnames["PWY-6572"] ="chondroitin sulfate and dermatan sulfate degradation I (bacterial)"
commonnames["PWY-5941"] ="glycogen degradation II"
commonnames["PWY-6906"] ="chitin derivatives degradation"
commonnames["PWY-6902"] ="chitin degradation II"
commonnames["PWY-6855"] ="chitin degradation I (archaea)"
commonnames["PWY-6805"] ="cellulose degradation I (cellulosome)"
commonnames["PWY-6788"] ="cellulose degradation II (fungi)"
commonnames["PWY-6822"] ="&iota;-carrageenan degradation"
commonnames["PWY-6821"] ="&kappa;-carrageenan degradation"
commonnames["PWY-6817"] ="&lambda;-carrageenan degradation"
commonnames["PWY-6827"] ="gellan degradation"
commonnames["PWY-6816"] ="agarose degradation"
commonnames["PWY-6815"] ="porphyran degradation"
commonnames["PWY-6813"] ="glucuronoarabinoxylan degradation"
commonnames["PWY-6790"] ="L-arabinan degradation"
commonnames["PWY0-1309"] ="chitobiose degradation"
commonnames["PWY-6527"] ="stachyose degradation"
commonnames["PWY-6028"] ="acetoin degradation"
commonnames["PWY-6927"] ="chlorophyll <i>a</i> degradation II"
commonnames["PWY-5874"] ="heme degradation"
commonnames["PWY-5098"] ="chlorophyll <i>a</i> degradation I"
commonnames["PWY-6370"] ="ascorbate recycling (cytosolic)"
commonnames["PWY-5372"] ="carbon tetrachloride degradation II"
commonnames["PWY-5370"] ="carbon tetrachloride degradation I"
commonnames["PWY-5369"] ="carbon tetrachloride degradation III"
commonnames["PWY-5822"] ="trichloroethylene degradation"
commonnames["PCPDEG-PWY"] ="pentachlorophenol degradation"
commonnames["PCEDEG-PWY"] ="tetrachloroethene degradation"
commonnames["GAMMAHEXCHLORDEG-PWY"] ="&gamma;-hexachlorocyclohexane degradation"
commonnames["12DICHLORETHDEG-PWY"] ="1,2-dichloroethane degradation"
commonnames["PWY0-43"] ="conversion of succinate to propionate"
commonnames["PWY0-42"] ="2-methylcitrate cycle I"
commonnames["PWY-5747"] ="2-methylcitrate cycle II"
commonnames["PROPIONMET-PWY"] ="methylmalonyl pathway"
commonnames["PWY-6698"] ="oxalate degradation V"
commonnames["PWY-6697"] ="oxalate degradation IV"
commonnames["PWY-6696"] ="oxalate degradation III"
commonnames["PWY-6695"] ="oxalate degradation II"
commonnames["PWY-6694"] ="oxalate degradation I"
commonnames["PWY-6060"] ="malonate degradation II (biotin-dependent)"
commonnames["PWY-5794"] ="malonate degradation I (biotin-independent)"
commonnames["PWY-6649"] ="glycolate and glyoxylate degradation III"
commonnames["GLYOXDEG-PWY"] ="glycolate and glyoxylate degradation II"
commonnames["GLYCOLATEMET-PWY"] ="glycolate and glyoxylate degradation I"
commonnames["PWY0-301"] ="L-ascorbate degradation I (bacterial, anaerobic)"
commonnames["PWY-6961"] ="L-ascorbate degradation II (bacterial, aerobic)"
commonnames["PWY-6960"] ="L-ascorbate degradation III"
commonnames["PWY-6959"] ="L-ascorbate degradation V"
commonnames["PWY-6704"] ="L-ascorbate degradation IV"
commonnames["PWY0-1312"] ="acetate formation from acetyl-CoA I"
commonnames["PWY-5536"] ="acetate formation from acetyl-CoA III (succinate)"
commonnames["PWY-5535"] ="acetate formation from acetyl-CoA II"
commonnames["2OXOBUTYRATECAT-PWY"] ="2-oxobutanoate degradation II"
commonnames["PYRUVDEHYD-PWY"] ="acetyl-CoA biosynthesis I (pyruvate dehydrogenase complex)"
commonnames["PWY0-1465"] ="D-malate degradation"
commonnames["PWY0-1313"] ="acetate conversion to acetyl-CoA"
commonnames["PWY-804"] ="glycolate degradation II"
commonnames["PWY-6373"] ="acrylate degradation"
commonnames["PWY-6038"] ="citrate degradation"
commonnames["PWY-6021"] ="nitrilotriacetate degradation"
commonnames["PWY-5749"] ="itaconate degradation"
commonnames["PWY-5654"] ="2-amino-3-carboxymuconate semialdehyde degradation to 2-oxopentenoate"
commonnames["PWY-5652"] ="2-amino-3-carboxymuconate semialdehyde degradation to glutaryl-CoA"
commonnames["PWY-5534"] ="propylene degradation"
commonnames["PWY-5177"] ="glutaryl-CoA degradation"
commonnames["PWY-5162"] ="2-oxopentenoate degradation"
commonnames["PWY-5074"] ="mevalonate degradation"
commonnames["PWY-301"] ="cyclohexane-1-carboxylate degradation (anaerobic)"
commonnames["PWY-2361"] ="3-oxoadipate degradation"
commonnames["GLUCONSUPER-PWY"] ="D-gluconate degradation"
commonnames["PWY-6966"] ="methanol oxidation to formaldehyde I"
commonnames["PWY-6510"] ="methanol oxidation to formaldehyde II"
commonnames["PWY-6509"] ="methanol oxidation to formaldehyde III"
commonnames["PWY-5506"] ="methanol oxidation to formaldehyde IV"
commonnames["PWY-6742"] ="methane oxidation to methanol II"
commonnames["PWY-1641"] ="methane oxidation to methanol I"
commonnames["RUMP-PWY"] ="formaldehyde oxidation I"
commonnames["PWY1G-170"] ="formaldehyde oxidation III (mycothiol-dependent)"
commonnames["PWY-1801"] ="formaldehyde oxidation II (glutathione-dependent)"
commonnames["PWY-1723"] ="formaldehyde oxidation VI (H<sub>4</sub>MPT pathway)"
commonnames["PWY-1722"] ="formaldehyde oxidation V (tetrahydrofolate pathway)"
commonnames["FORMASS-PWY"] ="formaldehyde oxidation IV (thiol-independent)"
commonnames["PWY-1861"] ="formaldehyde assimilation II (RuMP Cycle)"
commonnames["PWY-1622"] ="formaldehyde assimilation I (serine pathway)"
commonnames["P185-PWY"] ="formaldehyde assimilation III (dihydroxyacetone cycle)"
commonnames["PWY-5392"] ="reductive TCA cycle II"
commonnames["P42-PWY"] ="incomplete reductive TCA cycle"
commonnames["P23-PWY"] ="reductive TCA cycle I"
commonnames["PWY-5789"] ="3-hydroxypropionate/4-hydroxybutyrate cycle"
commonnames["PWY-5743"] ="3-hydroxypropionate cycle"
commonnames["CODH-PWY"] ="reductive acetyl coenzyme A pathway"
commonnames["PWYQT-4429"] ="CO<sub>2</sub> fixation into oxaloacetate (anapleurotic)"
commonnames["PWY-5493"] ="reductive monocarboxylic acid cycle"
commonnames["CO2FORM-PWY"] ="methanogenesis from methanol"
commonnames["PWY-1881"] ="formate oxidation to CO<sub>2</sub>"
commonnames["VALDEG-PWY"] ="valine degradation I"
commonnames["PWY-5057"] ="valine degradation II"
commonnames["TYRFUMCAT-PWY"] ="tyrosine degradation I"
commonnames["PWY3O-4108"] ="tyrosine degradation III"
commonnames["PWY-5151"] ="tyrosine degradation II"
commonnames["TRYPDEG-PWY"] ="tryptophan degradation II (via pyruvate)"
commonnames["TRPKYNCAT-PWY"] ="tryptophan degradation IV (via indole-3-lactate)"
commonnames["TRPCAT-PWY"] ="tryptophan degradation I (via anthranilate)"
commonnames["PWY-6307"] ="tryptophan degradation X (mammalian, via tryptamine)"
commonnames["PWY-5651"] ="tryptophan degradation to 2-amino-3-carboxymuconate semialdehyde"
commonnames["PWY-5081"] ="tryptophan degradation VIII (to tryptophol)"
commonnames["PWY-3162"] ="tryptophan degradation V (side chain pathway)"
commonnames["THREONINE-DEG2-PWY"] ="threonine degradation II"
commonnames["THRDLCTCAT-PWY"] ="threonine degradation III (to methylglyoxal)"
commonnames["PWY-5437"] ="threonine degradation I"
commonnames["PWY-5436"] ="threonine degradation IV"
commonnames["SERDEG-PWY"] ="L-serine degradation"
commonnames["PWY0-1535"] ="D-serine degradation"
commonnames["PROUT-PWY"] ="proline degradation"
commonnames["PWY-6318"] ="phenylalanine degradation IV (mammalian, via side chain)"
commonnames["PWY-5079"] ="phenylalanine degradation III"
commonnames["PHENYLALANINE-DEG1-PWY"] ="phenylalanine degradation I (aerobic)"
commonnames["ANAPHENOXI-PWY"] ="phenylalanine degradation II (anaerobic)"
commonnames["TAURINEDEG-PWY"] ="taurine degradation III"
commonnames["PWY0-981"] ="taurine degradation IV"
commonnames["PWY-1264"] ="taurine degradation II"
commonnames["PWY-1263"] ="taurine degradation I"
commonnames["PWY-5029"] ="imidazole-lactate degradation"
commonnames["ORN-AMINOPENTANOATE-CAT-PWY"] ="ornithine degradation I (proline biosynthesis)"
commonnames["OCTOPINEDEG-PWY"] ="octopine degradation"
commonnames["NOPALINEDEG-PWY"] ="nopaline degradation"
commonnames["CITRULLINE-DEG-PWY"] ="citrulline degradation"
commonnames["PWY-701"] ="methionine degradation II"
commonnames["PWY-5082"] ="methionine degradation III"
commonnames["METHIONINE-DEG1-PWY"] ="methionine degradation I (to homocysteine)"
commonnames["PWY0-461"] ="lysine degradation I"
commonnames["PWY-6328"] ="lysine degradation X"
commonnames["PWY-5324"] ="lysine degradation IX"
commonnames["PWY-5314"] ="lysine degradation VIII"
commonnames["PWY-5311"] ="lysine degradation VII"
commonnames["PWY-5298"] ="lysine degradation VI"
commonnames["PWY-5283"] ="lysine degradation V"
commonnames["PWY-5280"] ="lysine degradation IV"
commonnames["LYSINE-DEG1-PWY"] ="lysine degradation II"
commonnames["LYSDEGII-PWY"] ="lysine degradation III"
commonnames["P163-PWY"] ="lysine fermentation to acetate and butyrate"
commonnames["PWY-5076"] ="leucine degradation III"
commonnames["PWY-5075"] ="leucine degradation II"
commonnames["LEU-DEG2-PWY"] ="leucine degradation I"
commonnames["PWY-5078"] ="isoleucine degradation II"
commonnames["ILEUDEG-PWY"] ="isoleucine degradation I"
commonnames["PWY-5031"] ="histidine degradation V"
commonnames["PWY-5030"] ="histidine degradation III"
commonnames["PWY-5028"] ="histidine degradation II"
commonnames["HISTDEG-PWY"] ="histidine degradation IV"
commonnames["HISHP-PWY"] ="histidine degradation VI"
commonnames["HISDEG-PWY"] ="histidine degradation I"
commonnames["GLYCLEAV-PWY"] ="glycine cleavage complex"
commonnames["GLUTAMINEFUM-PWY"] ="glutamine degradation II"
commonnames["GLUTAMINDEG-PWY"] ="glutamine degradation I"
commonnames["PWY-5766"] ="glutamate degradation X"
commonnames["PWY-5087"] ="glutamate degradation VI (to pyruvate)"
commonnames["PWY-4321"] ="glutamate degradation IV"
commonnames["GLUTAMATE-DEG1-PWY"] ="glutamate degradation I"
commonnames["GLUDEG-I-PWY"] ="glutamate degradation III (via 4-aminobutyrate)"
commonnames["P162-PWY"] ="glutamate degradation V (via hydroxyglutarate)"
commonnames["PWY-5329"] ="L-cysteine degradation III"
commonnames["LCYSDEG-PWY"] ="L-cysteine degradation II"
commonnames["CYSTEINE-DEG-PWY"] ="L-cysteine degradation I"
commonnames["PWY-1781"] ="&beta;-alanine degradation II"
commonnames["BETA-ALA-DEGRADATION-I-PWY"] ="&beta;-alanine degradation I"
commonnames["MALATE-ASPARTATE-SHUTTLE-PWY"] ="aspartate degradation II"
commonnames["ASPARTATE-DEG1-PWY"] ="aspartate degradation I"
commonnames["PWY-4002"] ="asparagine degradation II"
commonnames["ASPARAGINE-DEG1-PWY"] ="asparagine degradation I"
commonnames["PWY-6422"] ="D-arginine degradation"
commonnames["PWY-5742"] ="arginine degradation IX (arginine:pyruvate transaminase pathway)"
commonnames["PWY-5024"] ="arginine degradation XI"
commonnames["AST-PWY"] ="arginine degradation II (AST pathway)"
commonnames["ARGDEG-V-PWY"] ="arginine degradation X (arginine monooxygenase pathway)"
commonnames["ARGDEG-IV-PWY"] ="arginine degradation VIII (arginine oxidase pathway)"
commonnames["ARGDEG-III-PWY"] ="arginine degradation IV (arginine decarboxylase/agmatine deiminase pathway)"
commonnames["ARG-GLU-PWY"] ="arginine degradation VII (arginase 3 pathway)"
commonnames["PWY1-2"] ="alanine degradation IV"
commonnames["ALANINE-DEG3-PWY"] ="alanine degradation III"
commonnames["ALADEG-PWY"] ="alanine degradation I"
commonnames["ALACAT2-PWY"] ="alanine degradation II (to D-lactate)"
commonnames["PWY-5159"] ="4-hydroxyproline degradation II"
commonnames["HYDROXYPRODEG-PWY"] ="4-hydroxyproline degradation I"
commonnames["PWY-6344"] ="ornithine degradation II (Stickland reaction)"
commonnames["PWY-6334"] ="L-dopa degradation"
commonnames["PWY-5084"] ="2-ketoglutarate dehydrogenase complex"
commonnames["PWY-5046"] ="branched-chain &alpha;-keto acid dehydrogenase complex"
commonnames["PWY0-1317"] ="L-lactaldehyde degradation (aerobic)"
commonnames["PWY0-1315"] ="L-lactaldehyde degradation (anaerobic)"
commonnames["PWY-901"] ="methylglyoxal degradation II"
commonnames["PWY-5462"] ="methylglyoxal degradation VIII"
commonnames["PWY-5458"] ="methylglyoxal degradation V"
commonnames["PWY-5456"] ="methylglyoxal degradation VII"
commonnames["PWY-5453"] ="methylglyoxal degradation III"
commonnames["PWY-5386"] ="methylglyoxal degradation I"
commonnames["MGLDLCTANA-PWY"] ="methylglyoxal degradation VI"
commonnames["PWY-4261"] ="glycerol degradation I"
commonnames["PWY-6952"] ="glycerophosphodiester degradation"
commonnames["PWY-6131"] ="glycerol degradation II"
commonnames["PWY-6130"] ="glycerol degradation III"
commonnames["GLYCEROLMETAB-PWY"] ="glycerol degradation V"
commonnames["PWY66-21"] ="ethanol degradation II"
commonnames["PWY66-162"] ="ethanol degradation IV"
commonnames["PWY66-161"] ="oxidative ethanol degradation III"
commonnames["ETOH-ACETYLCOA-ANA-PWY"] ="ethanol degradation I"
commonnames["PWY3O-246"] ="(R,R)-butanediol degradation"
commonnames["PWY-6388"] ="(S,S)-butanediol degradation"
commonnames["PWY0-1280"] ="ethylene glycol degradation"
commonnames["PWY-6464"] ="polyvinyl alcohol degradation"
commonnames["CYCLOHEXANOL-OXIDATION-PWY"] ="cyclohexanol degradation"
commonnames["PWY-5731"] ="atrazine degradation III"
commonnames["PWY-5727"] ="atrazine degradation II"
commonnames["P141-PWY"] ="atrazine degradation I (aerobic)"
commonnames["PWY-5726"] ="deethylsimazine degradation"
commonnames["PWY-5171"] ="<i>N</i>-cyclopropylmelamine degradation"
commonnames["PWY-5170"] ="melamine degradation"
commonnames["PWY-5429"] ="<i>p</i>-xylene degradation to <i>p</i>-toluate"
commonnames["PWY-5428"] ="<i>m</i>-xylene degradation to <i>m</i>-toluate"
commonnames["PWY-5691"] ="urate degradation to allantoin"
commonnames["TOLUENE-DEG-DIOL-PWY"] ="toluene degradation to 2-oxopent-4-enoate (<I>via</I> toluene-<I>cis</I>-diol)"
commonnames["TOLUENE-DEG-CATECHOL-PWY"] ="toluene degradation to benzoate"
commonnames["TOLUENE-DEG-4-OH-PWY"] ="toluene degradation to protocatechuate (<I>via</I> <i>p</I>-cresol)"
commonnames["TOLUENE-DEG-3-OH-PWY"] ="toluene degradation to 2-oxopent-4-enoate (<I>via 4-methylcatechol</i>)"
commonnames["TOLUENE-DEG-2-OH-PWY"] ="toluene degradation to 2-oxopent-4-enoate I (<I>via</I> <i>o</I>-cresol)"
commonnames["PWY-81"] ="toluene degradation to benzoyl-CoA (anaerobic)"
commonnames["TOLSULFDEG-PWY"] ="4-toluenesulfonate degradation I"
commonnames["PWY-5165"] ="4-toluenesulfonate degradation II"
commonnames["PWY-6041"] ="4-sulfocatechol degradation"
commonnames["DESULFONATION-PWY"] ="benzenesulfonate degradation"
commonnames["SHIKIMATEDEG-PWY"] ="shikimate degradation I"
commonnames["PWY-6419"] ="shikimate degradation II"
commonnames["PWY-6640"] ="salicylate degradation IV"
commonnames["PWY-6636"] ="salicylate degradation III"
commonnames["PWY-6224"] ="salicylate degradation II"
commonnames["PWY-6183"] ="salicylate degradation I"
commonnames["QUINATEDEG-PWY"] ="quinate degradation I"
commonnames["PWY-6416"] ="quinate degradation II"
commonnames["PWY-6336"] ="protocatechuate degradation III (<I>para</I>-cleavage pathway)"
commonnames["PROTOCATECHUATE-ORTHO-CLEAVAGE-PWY"] ="protocatechuate degradation II (ortho-cleavage pathway)"
commonnames["P184-PWY"] ="protocatechuate degradation I (<I>meta</I>-cleavage pathway)"
commonnames["PWY0-321"] ="phenylacetate degradation I (aerobic)"
commonnames["PWY-1341"] ="phenylacetate degradation II (anaerobic)"
commonnames["PWY-5418"] ="phenol degradation I (aerobic)"
commonnames["PHENOLDEG-PWY"] ="phenol degradation II (anaerobic)"
commonnames["PWY-6690"] ="cinnamate and 3-hydroxycinnamate degradation to 2-oxopent-4-enoate"
commonnames["PWY-6080"] ="4-ethylphenol degradation (anaerobic)"
commonnames["HCAMHPDEG-PWY"] ="3-phenylpropionate and 3-(3-hydroxyphenyl)propionate degradation to 2-oxopent-4-enoate"
commonnames["PWY-5648"] ="2-nitrobenzoate degradation II"
commonnames["PWY-2381"] ="4-nitrobenzoate degradation"
commonnames["PWY-5645"] ="4-chloronitrobenzene degradation"
commonnames["PWY-5640"] ="nitrobenzene degradation II"
commonnames["PWY-5637"] ="nitrobenzene degradation I"
commonnames["PWY-5644"] ="4-nitrotoluene degradation II"
commonnames["P421-PWY"] ="4-nitrotoluene degradation I"
commonnames["PWY-5641"] ="2-nitrotoluene degradation"
commonnames["PWY-5643"] ="2,6-dinitrotoluene degradation"
commonnames["PWY-5642"] ="2,4-dinitrotoluene degradation"
commonnames["PWY-5488"] ="4-nitrophenol degradation II"
commonnames["PWY-5487"] ="4-nitrophenol degradation I"
commonnames["PWY-5636"] ="2-nitrophenol degradation"
commonnames["PWY-6051"] ="2,4,6-trinitrotoluene degradation"
commonnames["PWY-722"] ="nicotinate degradation I"
commonnames["PWY-5055"] ="nicotinate degradation III"
commonnames["PWY-5033"] ="nicotinate degradation II"
commonnames["PWY-1501"] ="mandelate degradation I"
commonnames["P3-PWY"] ="gallate degradation III (anaerobic)"
commonnames["GALLATE-DEGRADATION-II-PWY"] ="gallate degradation I"
commonnames["GALLATE-DEGRADATION-I-PWY"] ="gallate degradation II"
commonnames["PWY-1381"] ="fluorene degradation II"
commonnames["FLUORENE-DEG-9-ONE-PWY"] ="fluorene degradation I"
commonnames["PWY-6192"] ="3,4-dichlorotoluene degradation"
commonnames["PWY-6191"] ="2,5-dichlorotoluene degradation"
commonnames["PWY-6190"] ="2,4-dichlorotoluene degradation"
commonnames["PWY-6104"] ="3-chlorotoluene degradation II"
commonnames["PWY-6103"] ="3-chlorotoluene degradation I"
commonnames["PWY-6214"] ="3-chlorocatechol degradation III (<i>meta</i> pathway)"
commonnames["PWY-6193"] ="3-chlorocatechol degradation II (<i>ortho</i>)"
commonnames["PWY-6089"] ="3-chlorocatechol degradation I (<i>ortho</i>)"
commonnames["PWY-6094"] ="3,4,6-trichlorocatechol degradation"
commonnames["PWY-6093"] ="3,4-dichlorocatechol degradation"
commonnames["PWY-6087"] ="4-chlorocatechol degradation"
commonnames["PWY-6084"] ="3,5-dichlorocatechol degradation"
commonnames["PWY-6228"] ="3-chlorobenzoate degradation III (via gentisate)"
commonnames["PWY-6216"] ="3-chlorobenzoate degradation II (via protocatechuate)"
commonnames["PWY-6088"] ="3-chlorobenzoate degradation I (via chlorocatechol)"
commonnames["PWY-6221"] ="2-chlorobenzoate degradation"
commonnames["PWY-6217"] ="3,4-dichlorobenzoate degradation"
commonnames["PWY-6215"] ="4-chlorobenzoate degradation"
commonnames["PWY-6099"] ="1,2,4,5-tetrachlorobenzene degradation"
commonnames["PWY-6091"] ="1,2,4-trichlorobenzene degradation"
commonnames["PWY-6090"] ="1,2-dichlorobenzene degradation"
commonnames["PWY-6083"] ="chlorobenzene degradation"
commonnames["PWY-6081"] ="1,3-dichlorobenzene degradation"
commonnames["14DICHLORBENZDEG-PWY"] ="1,4-dichlorobenzene degradation"
commonnames["PWY-6200"] ="2,4,5-trichlorophenoxyacetate degradation"
commonnames["PWY-6197"] ="chlorinated phenols degradation"
commonnames["PWY-6178"] ="2,4,6-trichlorophenol degradation"
commonnames["PWY-6107"] ="chlorosalicylate degradation"
commonnames["PWY-6102"] ="5-chloro-3-methyl-catechol degradation"
commonnames["PWY-6086"] ="4-chloro-2-methylphenoxyacetate degradation"
commonnames["PWY-6085"] ="2,4-dichlorophenoxyacetate  degradation"
commonnames["PWY-5419"] ="catechol degradation to 2-oxopent-4-enoate II"
commonnames["P183-PWY"] ="catechol degradation to 2-oxopent-4-enoate I"
commonnames["CATECHOL-ORTHO-CLEAVAGE-PWY"] ="catechol degradation to &beta;-ketoadipate"
commonnames["PWY-1361"] ="benzoyl-CoA degradation I (aerobic)"
commonnames["P321-PWY"] ="benzoyl-CoA degradation III (anaerobic)"
commonnames["CENTBENZCOA-PWY"] ="benzoyl-CoA degradation II (anaerobic)"
commonnames["PWY-283"] ="benzoate degradation II (aerobic and anaerobic)"
commonnames["PWY-2503"] ="benzoate degradation I (aerobic)"
commonnames["PWY-6504"] ="anthranilate degradation IV (aerobic)"
commonnames["PWY-6079"] ="anthranilate degradation I (aerobic)"
commonnames["PWY-6077"] ="anthranilate degradation II (aerobic)"
commonnames["2AMINOBENZDEG-PWY"] ="anthranilate degradation III (anaerobic)"
commonnames["PWY5F9-3233"] ="phthalate degradation"
commonnames["PWY5F9-12"] ="biphenyl degradation"
commonnames["PWY-741"] ="<I>p</I>-cymene degradation to <I>p</I>-cumate"
commonnames["PWY-721"] ="3-methylquinoline degradation"
commonnames["PWY-6941"] ="styrene degradation"
commonnames["PWY-681"] ="dibenzothiophene desulfurization"
commonnames["PWY-6781"] ="chlorogenic acid degradation"
commonnames["PWY-6550"] ="carbazole degradation"
commonnames["PWY-6533"] ="aniline degradation"
commonnames["PWY-6532"] ="diphenylamine degradation"
commonnames["PWY-6343"] ="ferulate degradation"
commonnames["PWY-6340"] ="5,5'-dehydrodivanillate degradation"
commonnames["PWY-6223"] ="gentisate degradation"
commonnames["PWY-6210"] ="2-aminophenol degradation"
commonnames["PWY-6185"] ="4-methylcatechol degradation (<i>ortho</i> cleavage)"
commonnames["PWY-6184"] ="methylsalicylate degradation"
commonnames["PWY-5490"] ="paraoxon degradation"
commonnames["PWY-5489"] ="methyl parathion degradation"
commonnames["PWY-5450"] ="benzene degradation"
commonnames["PWY-5427"] ="naphthalene degradation"
commonnames["PWY-5169"] ="cyanurate degradation"
commonnames["PWY-5163"] ="<I>p</I>-cumate degradation to 2-oxopent-4-enoate"
commonnames["PWY-481"] ="ethylbenzene degradation (anaerobic)"
commonnames["PWY-2421"] ="indole-3-acetate degradation to anthranilate"
commonnames["PWY-142"] ="<i>m</i>-xylene degradation (anaerobic)"
commonnames["PARATHION-DEGRADATION-PWY"] ="parathion degradation"
commonnames["P662-PWY"] ="dibenzofuran degradation"
commonnames["P661-PWY"] ="dibenzo-<I>p</I>-dioxin degradation"
commonnames["P345-PWY"] ="aldoxime degradation"
commonnames["P343-PWY"] ="resorcinol degradation"
commonnames["P342-PWY"] ="orcinol degradation"
commonnames["METHYLGALLATE-DEGRADATION-PWY"] ="methylgallate degradation"
commonnames["M-CRESOL-DEGRADATION-PWY"] ="<I>m</I>-cresol degradation"
commonnames["4TOLCARBDEG-PWY"] ="4-toluenecarboxylate degradation"
commonnames["4-HYDROXYMANDELATE-DEGRADATION-PWY"] ="4-hydroxymandelate degradation"
commonnames["3-HYDROXYPHENYLACETATE-DEGRADATION-PWY"] ="4-hydroxyphenylacetate degradation"
commonnames["2ASDEG-PWY"] ="orthanilate degradation"
commonnames["PWY-5704"] ="urea degradation II"
commonnames["PWY-5703"] ="urea degradation I"
commonnames["PWY-6441"] ="spermine and spermidine degradation III"
commonnames["PWY-6440"] ="spermine and spermidine degradation II"
commonnames["PWY-6117"] ="spermine and spermidine degradation I"
commonnames["PWY0-1221"] ="putrescine degradation II"
commonnames["PWY-3"] ="putrescine degradation V"
commonnames["PWY-2"] ="putrescine degradation IV"
commonnames["PWY-0"] ="putrescine degradation III"
commonnames["PUTDEG-PWY"] ="putrescine degradation I"
commonnames["GLUAMCAT-PWY"] ="<i>N</i>-acetylglucosamine degradation I"
commonnames["PWY-6967"] ="methylamine degradation I"
commonnames["PWY-6965"] ="methylamine degradation II"
commonnames["PWY-4741"] ="creatinine degradation III"
commonnames["PWY-4722"] ="creatinine degradation II"
commonnames["CRNFORCAT-PWY"] ="creatinine degradation I"
commonnames["PWY-3721"] ="choline degradation II"
commonnames["PWY-3641"] ="carnitine degradation III"
commonnames["PWY-3602"] ="carnitine degradation II"
commonnames["CARNMET-PWY"] ="carnitine degradation I"
commonnames["PWY-5698"] ="allantoin degradation to ureidoglycolate II (ammonia producing)"
commonnames["PWY-5697"] ="allantoin degradation to ureidoglycolate I (urea producing)"
commonnames["PWY-6537"] ="4-aminobutyrate degradation II"
commonnames["PWY-6536"] ="4-aminobutyrate degradation III"
commonnames["PWY-6535"] ="4-aminobutyrate degradation I"
commonnames["PWY-6473"] ="4-aminobutyrate degradation IV"
commonnames["PWY-5022"] ="4-aminobutyrate degradation V"
commonnames["PWY6666-2"] ="dopamine degradation"
commonnames["PWY0-1477"] ="ethanolamine utilization"
commonnames["PWY-6534"] ="phenylethylamine degradation II"
commonnames["PWY-6181"] ="histamine degradation"
commonnames["PWY-5736"] ="isopropylamine degradation"
commonnames["PWY-3661"] ="glycine betaine degradation"
commonnames["P561-PWY"] ="stachydrine degradation"
commonnames["2PHENDEG-PWY"] ="phenylethylamine degradation I"
commonnames["PWY0-1324"] ="<i>N</i>-acetylneuraminate and <i>N</i>-acetylmannosamine degradation"
commonnames["PWY-6968"] ="trimethylamine degradation"
commonnames["PWY-6814"] ="acidification and chitin degradation (in carnivorous plants)"
commonnames["PWY-761"] ="rhizobactin 1021 biosynthesis"
commonnames["PWY-6574"] ="achromobactin biosynthesis"
commonnames["PWY-6409"] ="pyoverdine I biosynthesis"
commonnames["PWY-6408"] ="pyochelin biosynthesis"
commonnames["PWY-6407"] ="yersiniabactin biosynthesis"
commonnames["PWY-6381"] ="bisucaberin biosynthesis"
commonnames["PWY-6379"] ="alcaligin biosynthesis"
commonnames["PWY-6378"] ="putrebactin biosynthesis"
commonnames["PWY-6376"] ="desferrioxamine B biosynthesis"
commonnames["PWY-6375"] ="desferrioxamine E biosynthesis"
commonnames["PWY-6289"] ="petrobactin biosynthesis"
commonnames["PWY-5925"] ="hydroxylated mugineic acid phytosiderophore biosynthesis"
commonnames["PWY-5912"] ="2'-deoxymugineic acid phytosiderophore biosynthesis"
commonnames["AEROBACTINSYN-PWY"] ="aerobactin biosynthesis"
commonnames["PWY-5002"] ="tetrahydroxyxanthone biosynthesis (from 3-hydroxybenzoate)"
commonnames["PWY-5001"] ="tetrahydroxyxanthone biosynthesis (from benzoate)"
commonnames["PWY-6115"] ="avenacin biosynthesis"
commonnames["PWY-6109"] ="mangrove triterpenoid biosynthesis"
commonnames["PWY-6105"] ="C30 botryococcene biosynthesis"
commonnames["PWY-6098"] ="diploterol and cycloartenol biosynthesis"
commonnames["PWY-6095"] ="dammara-20,24-diene biosynthesis"
commonnames["PWY-6008"] ="baruol biosynthesis"
commonnames["PWY-6007"] ="arabidiol biosynthesis"
commonnames["PWY-6005"] ="marneral biosynthesis"
commonnames["PWY-5992"] ="thalianol and derivatives biosynthesis"
commonnames["PWY-5774"] ="saponin biosynthesis IV"
commonnames["PWY-5759"] ="saponin biosynthesis III"
commonnames["PWY-5756"] ="saponin biosynthesis II"
commonnames["PWY-5672"] ="ginsenoside biosynthesis"
commonnames["PWY-5377"] ="&alpha;-amyrin biosynthesis"
commonnames["PWY-5203"] ="soybean saponin I biosynthesis"
commonnames["PWY-112"] ="lupeol biosynthesis"
commonnames["PWY-6681"] ="neurosporaxanthin biosynthesis"
commonnames["PWY-5398"] ="crocetin esters biosynthesis"
commonnames["PWY-5397"] ="crocetin biosynthesis"
commonnames["PWY-5305"] ="bixin biosynthesis"
commonnames["PWY-6244"] ="bergamotene biosynthesis II"
commonnames["PWY-6243"] ="bergamotene biosynthesis I"
commonnames["PWY-6836"] ="santalene biosynthesis II"
commonnames["PWY-6669"] ="&delta;-guaiene biosynthesis"
commonnames["PWY-6540"] ="costunolide biosynthesis"
commonnames["PWY-6294"] ="selinene biosynthesis"
commonnames["PWY-6291"] ="valencene and 7-epi-&alpha;-selinene biosynthesis"
commonnames["PWY-6290"] ="&beta;-cubebene biosynthesis"
commonnames["PWY-6278"] ="botrydial biosynthesis"
commonnames["PWY-6275"] ="&beta;-caryophyllene biosynthesis"
commonnames["PWY-6271"] ="eudesmol biosynthesis"
commonnames["PWY-6265"] ="zerumbone biosynthesis"
commonnames["PWY-6258"] ="patchoulol biosynthesis"
commonnames["PWY-6257"] ="curcumene biosynthesis"
commonnames["PWY-6254"] ="santalene biosynthesis I"
commonnames["PWY-6128"] ="cis-calamenene related sesquiterpenoids biosynthesis"
commonnames["PWY-5950"] ="geosmin biosynthesis"
commonnames["PWY-5773"] ="gossypol biosynthesis"
commonnames["PWY-5733"] ="germacrene biosynthesis"
commonnames["PWY-5725"] ="farnesene biosynthesis"
commonnames["PWY-5434"] ="(3<i>E</i>)-4,8-dimethylnona-1,3,7-triene biosynthesis"
commonnames["PWY-5425"] ="oleoresin sesquiterpene volatiles biosynthesis"
commonnames["PWY-5195"] ="artemisinin biosynthesis"
commonnames["PWY-5815"] ="rubber biosynthesis"
commonnames["PWY2OL-4"] ="linalool biosynthesis"
commonnames["PWY-6668"] ="(<i>E,E</i>)-4,8,12-trimethyltrideca-1,3,7,11-tetraene biosynthesis"
commonnames["PWY-6451"] ="3-carene biosynthesis"
commonnames["PWY-6449"] ="fenchone biosynthesis"
commonnames["PWY-6447"] ="trichome monoterpenes biosynthesis"
commonnames["PWY-6445"] ="fenchol biosynthesis II"
commonnames["PWY-6437"] ="fenchol biosynthesis I"
commonnames["PWY-6436"] ="perillyl alcohol biosynthesis"
commonnames["PWY-5928"] ="(4<i>R</i>)-carvone biosynthesis"
commonnames["PWY-5829"] ="geraniol and geranial biosynthesis"
commonnames["PWY-5813"] ="bornyl diphosphate biosynthesis"
commonnames["PWY-5423"] ="oleoresin monoterpene volatiles biosynthesis"
commonnames["PWY-3061"] ="menthol biosynthesis"
commonnames["PWY-3041"] ="monoterpene biosynthesis"
commonnames["PWY-6913"] ="methylbutenol biosynthesis"
commonnames["PWY-5422"] ="isopimaric acid biosynthesis"
commonnames["PWY-5421"] ="dehydroabietic acid biosynthesis"
commonnames["PWY-5414"] ="palustric acid biosynthesis"
commonnames["PWY-5413"] ="neoabietic acid biosynthesis"
commonnames["PWY-5412"] ="levopimaric acid biosynthesis"
commonnames["PWY-5411"] ="abietic acid biosynthesis"
commonnames["PWY-5063"] ="phytyl diphosphate biosynthesis"
commonnames["PWY-6691"] ="plaunotol biosynthesis"
commonnames["PWY-6659"] ="fusicoccins biosynthesis"
commonnames["PWY-6645"] ="labdane-type diterpenes biosynthesis"
commonnames["PWY-6304"] ="casbene biosynthesis"
commonnames["PWY-5660"] ="taxol biosynthesis"
commonnames["PWY-5107"] ="phytol salvage pathway"
commonnames["PWY-6475"] ="<i>trans</i>-lycopene biosynthesis II (plants)"
commonnames["PWY-6809"] ="neoxanthin biosnythesis"
commonnames["PWY-6581"] ="spirilloxanthin and 2,2'-diketo-spirilloxanthin biosynthesis"
commonnames["PWY-6288"] ="zeaxanthin-&beta;-D-diglucoside biosynthesis"
commonnames["PWY-6287"] ="neurosporene biosynthesis"
commonnames["PWY-6286"] ="spheroidene and spheroidenone biosynthesis"
commonnames["PWY-6280"] ="synechoxanthin biosynthesis"
commonnames["PWY-6279"] ="myxol-2' fucoside biosynthesis"
commonnames["PWY-5947"] ="lutein biosynthesis"
commonnames["PWY-5946"] ="&delta;-carotene biosynthesis"
commonnames["PWY-5945"] ="antheraxanthin and violaxanthin biosynthesis"
commonnames["PWY-5944"] ="zeaxanthin biosynthesis"
commonnames["PWY-5943"] ="&beta;-carotene biosynthesis"
commonnames["PWY-5291"] ="canthaxanthin biosynthesis"
commonnames["PWY-5288"] ="astaxanthin biosynthesis"
commonnames["PWY-5175"] ="lactucaxanthin biosynthesis"
commonnames["PWY-5174"] ="capsanthin and capsorubin biosynthesis"
commonnames["PWY-6767"] ="4,4'-diapolycopenedioate biosynthesis"
commonnames["PWY-5875"] ="staphyloxanthin biosynthesis"
commonnames["PWY-5828"] ="lacinilene C biosynthesis"
commonnames["PWY-5827"] ="heliocides biosynthesis"
commonnames["PWY-6295"] ="olivetol biosynthesis"
commonnames["PWY-5808"] ="hyperforin biosynthesis"
commonnames["PWY-5140"] ="cannabinoid biosynthesis"
commonnames["PWY-5133"] ="cohumulone biosynthesis"
commonnames["PWY-5132"] ="humulone biosynthesis"
commonnames["PWY-6372"] ="1D-<i>myo</i>-inositol hexakisphosphate biosynthesis IV (<i>Dictyostelium</i>)"
commonnames["PWY-4661"] ="1D-<i>myo</i>-inositol hexakisphosphate biosynthesis III (<i>Spirodela polyrrhiza</i>)"
commonnames["PWY-6554"] ="1D-<i>myo</i>-inositol hexakisphosphate biosynthesis V (from Ins(1,3,4)P3)"
commonnames["PWY-6362"] ="1D-<i>myo</i>-inositol hexakisphosphate biosynthesis II (mammalian)"
commonnames["PWY-6361"] ="1D-<i>myo</i>-inositol hexakisphosphate biosynthesis I  (from Ins(1,4,5)P3)"
commonnames["PWY-2301"] ="<i>myo</i>-inositol biosynthesis"
commonnames["PWY-6369"] ="inositol pyrophosphates biosynthesis"
commonnames["PWY-6366"] ="D-<i>myo</i>-inositol (1,4,5,6)-tetrakisphosphate biosynthesis"
commonnames["PWY-6365"] ="D-<i>myo</i>-inositol (3,4,5,6)-tetrakisphosphate biosynthesis"
commonnames["PWY-6364"] ="D-<i>myo</i>-inositol (1,3,4)-trisphosphate biosynthesis"
commonnames["PWY-6363"] ="D-<i>myo</i>-inositol (1,4,5)-trisphosphate degradation"
commonnames["PWY-6351"] ="D-<i>myo</i>-inositol (1,4,5)-trisphosphate biosynthesis"
commonnames["PWY-5338"] ="galactosylcyclitol biosynthesis"
commonnames["PWY-5978"] ="kanosamine biosynthesis"
commonnames["PWY-5380"] ="A series fagopyritols biosynthesis"
commonnames["PWY-5379"] ="B series fagopyritols biosynthesis"
commonnames["PWY-6627"] ="salinosporamide A biosynthesis"
commonnames["PWY-6438"] ="phenylphenalenone biosynthesis"
commonnames["PWY-6316"] ="aromatic polyketides biosynthesis"
commonnames["PWY-6314"] ="plumbagin biosynthesis"
commonnames["PWY-6312"] ="barbaloin biosynthesis"
commonnames["PWY-6310"] ="aloesone biosynthesis II"
commonnames["PWY-5810"] ="usnate biosynthesis"
commonnames["PWY-5393"] ="raspberry ketone biosynthesis"
commonnames["PWY-4801"] ="aloesone biosynthesis I"
commonnames["PWY-12"] ="pentaketide chromone biosynthesis"
commonnames["PWY-6888"] ="zealexin biosynthesis"
commonnames["PWY-2961"] ="sesquiterpenoid phytoalexins biosynthesis"
commonnames["PWY-2921"] ="capsidiol biosynthesis"
commonnames["PWY-2981"] ="diterpene phytoalexins precursors biosynthesis"
commonnames["CAMALEXIN-SYN"] ="camalexin biosynthesis"
commonnames["PWY-6418"] ="4-hydroxycoumarin biosynthesis"
commonnames["PWY-5477"] ="gallotannin biosynthesis"
commonnames["PWY-5476"] ="cornusiin E biosynthesis"
commonnames["PWY-5475"] ="pentagalloylglucose biosynthesis"
commonnames["PWY-84"] ="resveratrol biosynthesis"
commonnames["PWY-6939"] ="tetrahydroxystilbenes biosynthesis"
commonnames["PWY-6665"] ="pterostilbene biosynthesis"
commonnames["PWY-5045"] ="pinosylvin metabolism"
commonnames["PWY-5987"] ="sorgoleone biosynthesis"
commonnames["PWY-5802"] ="alizarin biosynthesis"
commonnames["PWY-5801"] ="lawsone biosynthesis"
commonnames["PWY-5795"] ="juglone biosynthesis"
commonnames["PWY-5780"] ="hypericin biosynthesis"
commonnames["PWY-5701"] ="shikonin biosynthesis"
commonnames["PWY-5371"] ="chrysophanol biosynthesis"
commonnames["PWY-116"] ="coniferin metabolism"
commonnames["PWY-83"] ="monolignol glucosides biosynthesis"
commonnames["PWY-6302"] ="dihydroconiferyl alcohol biosynthesis"
commonnames["PWY-6824"] ="justicidin B biosynthesis"
commonnames["PWY-6820"] ="diphyllin biosynthesis"
commonnames["PWY-5479"] ="podophyllotoxin and 6-methoxypodophyllotoxin biosynthesis"
commonnames["PWY-5469"] ="sesamin biosynthesis"
commonnames["PWY-5466"] ="matairesinol biosynthesis"
commonnames["PWY-641"] ="proanthocyanidin biosynthesis from flavanols"
commonnames["PWY-6914"] ="sophoraflavanone G biosynthesis"
commonnames["PWY-5135"] ="xanthohumol biosynthesis"
commonnames["PWY1F-823"] ="leucopelargonidin and leucocyanidin biosynthesis"
commonnames["PWY-5152"] ="leucodelphinidin biosynthesis"
commonnames["PWY-6427"] ="rot-2'-enonate biosynthesis"
commonnames["PWY-6425"] ="rotenoid biosynthesis II"
commonnames["PWY-5775"] ="rotenoid biosynthesis I"
commonnames["PWY-6332"] ="coumestrol biosynthesis"
commonnames["PWY-5825"] ="dalpatein and dalnigrein biosynthesis"
commonnames["PWY-5821"] ="dalcochinin biosynthesis"
commonnames["PWY-5729"] ="vestitol and sativan biosynthesis"
commonnames["PWY-5061"] ="6,7,4'-trihydroxyisoflavone biosynthesis"
commonnames["PWY-4681"] ="kievitone biosynthesis"
commonnames["PWY-4502"] ="wighteone and luteone biosynthesis"
commonnames["PWY-3042"] ="phaseollin biosynthesis"
commonnames["PWY-2762"] ="glyceollin biosynthesis II"
commonnames["PWY-2761"] ="glyceollin biosynthesis I"
commonnames["PWY-2467"] ="pisatin biosynthesis"
commonnames["PWY-2464"] ="maackiain biosynthesis"
commonnames["PWY-2463"] ="medicarpin biosynthesis"
commonnames["PWY-2321"] ="formononetin biosynthesis"
commonnames["PWY-2083"] ="isoflavonoid biosynthesis II"
commonnames["PWY-2002"] ="isoflavonoid biosynthesis I"
commonnames["PWY-6631"] ="<i>O</i> -methylation of tricetin"
commonnames["PWY-6360"] ="flavonol glucosylation I"
commonnames["PWY-6199"] ="quercetinsulphates biosynthesis"
commonnames["PWY-6064"] ="methylquercetin biosynthesis"
commonnames["PWY-6035"] ="2,3-<i>cis</i>-flavanols biosynthesis"
commonnames["PWY-5391"] ="syringetin biosynthesis"
commonnames["PWY-5390"] ="rutin biosynthesis"
commonnames["PWY-5363"] ="chrysin biosynthesis"
commonnames["PWY-5348"] ="kaempferol triglucoside biosynthesis"
commonnames["PWY-5321"] ="quercetin glucoside biosynthesis (Arabidopsis)"
commonnames["PWY-5320"] ="kaempferol glucoside biosynthesis (Arabidopsis)"
commonnames["PWY-5059"] ="pinobanksin biosynthesis"
commonnames["PWY-3101"] ="flavonol biosynthesis"
commonnames["PWY-6239"] ="luteolin glycosides biosynthesis"
commonnames["PWY-6024"] ="isovitexin and isovitexin glycosides biosynthesis"
commonnames["PWY-6015"] ="vitexin and detivatives biosynthesis"
commonnames["PWY-6010"] ="apigenin glycosides biosynthesis"
commonnames["PWY-5793"] ="maysin biosynthesis"
commonnames["PWY-5119"] ="acacetin biosynthesis"
commonnames["PWY-5060"] ="luteolin biosynthesis"
commonnames["PWY-6602"] ="C-glycosylflavone biosynthesis"
commonnames["PWY-5118"] ="ponciretin biosynthesis"
commonnames["PWY-5116"] ="sakuranetin biosynthesis"
commonnames["PWY-5105"] ="hesperitin glycoside biosynthesis"
commonnames["PWY-5094"] ="naringenin glycoside biosynthesis"
commonnames["PWY-6325"] ="echinatin biosynthesis"
commonnames["PWY-5339"] ="chalcone 2'-<i>O</i>-glucoside biosynthesis"
commonnames["PWY-5161"] ="6'-deoxychalcone metabolism"
commonnames["PWY-6401"] ="hispidol biosynthesis"
commonnames["PWY-1901"] ="aurone biosynthesis"
commonnames["PWY-5307"] ="gentiodelphin biosynthesis"
commonnames["PWY-5295"] ="ternatin C5 biosynthesis"
commonnames["PWY-5286"] ="anthocyanidin sophoroside metabolism"
commonnames["PWY-5284"] ="shisonin biosynthesis"
commonnames["PWY-5268"] ="salvianin biosynthesis"
commonnames["PWY-5160"] ="rose anthocyanin biosynthesis"
commonnames["PWY-5153"] ="anthocyanin biosynthesis (delphinidin 3-O-glucoside)"
commonnames["PWY-5139"] ="pelargonidin conjugates biosynthesis"
commonnames["PWY-5125"] ="anthocyanin biosynthesis (pelargonidin 3-O-glucoside, cyanidin 3-O-glucoside)"
commonnames["PWY1F-FLAVSYN"] ="flavonoid biosynthesis"
commonnames["PWY-6787"] ="flavonoid biosynthesis (in equisetum)"
commonnames["PWY-6515"] ="phloridzin biosynthesis"
commonnames["PWY-6232"] ="chrysoeriol biosynthesis"
commonnames["PWY-6029"] ="2,3-<i>trans</i>-flavanols biosynthesis"
commonnames["PWY-5319"] ="coumarin metabolism (to melilotic acid)"
commonnames["PWY-4922"] ="6-methoxymellein biosynthesis"
commonnames["PWY-6792"] ="scopoletin biosynthesis"
commonnames["PWY-5365"] ="linear furanocoumarin biosynthesis"
commonnames["PWY-5349"] ="esculetin biosynthesis"
commonnames["PWY-5176"] ="coumarin biosynthesis (via 2-coumarate)"
commonnames["PWY1F-467"] ="phenylpropanoid biosynthesis, initial reactions"
commonnames["PWY-6040"] ="chlorogenic acid biosynthesis II"
commonnames["PWY-6039"] ="chlorogenic acid biosynthesis I"
commonnames["PWY-5968"] ="cinnamate esters biosynthesis"
commonnames["PWY-5867"] ="<i>t</i>-anethole biosynthesis"
commonnames["PWY-5859"] ="eugenol and isoeugenol biosynthesis"
commonnames["PWY-5168"] ="ferulate and sinapate biosynthesis"
commonnames["PWY-5049"] ="rosmarinic acid biosynthesis II"
commonnames["PWY-5048"] ="rosmarinic acid biosynthesis I"
commonnames["PWY-4421"] ="curcumin glucoside biosynthesis"
commonnames["PWY-4201"] ="volatile cinnamoic ester biosynthesis"
commonnames["PWY-3301"] ="sinapate ester biosynthesis"
commonnames["PWY-2181"] ="free phenylpropanoid acid biosynthesis"
commonnames["PWY-4203"] ="volatile benzenoid biosynthesis I (ester formation)"
commonnames["PWY-6766"] ="salicin biosynthesis"
commonnames["PWY-6763"] ="salicortin biosynthesis"
commonnames["PWY-6444"] ="benzoate biosynthesis II (CoA-independent, non-&beta;-oxidative)"
commonnames["PWY-6835"] ="6-gingerol biosynthesis"
commonnames["PWY-6801"] ="volatile esters biosynthesis (during fruit ripening)"
commonnames["PWY-6673"] ="caffeoylglucarate biosynthesis"
commonnames["PWY-5882"] ="epoxypseudoisoeugenol-2-methylbutyrate biosynthesis"
commonnames["PWY-5665"] ="vanilla biosynthesis"
commonnames["PWY-5021"] ="willardiine and isowillardiine biosynthesis"
commonnames["PWY-5"] ="canavanine biosynthesis"
commonnames["PWY-4985"] ="mimosine biosynthesis"
commonnames["PWY-4961"] ="&beta;-pyrazole-1-ylalanine biosynthesis"
commonnames["PWY-1"] ="lathyrine biosynthesis"
commonnames["PWYQT-4475"] ="glucosinolate biosynthesis from hexahomomethionine"
commonnames["PWYQT-4474"] ="glucosinolate biosynthesis from pentahomomethionine"
commonnames["PWYQT-4473"] ="glucosinolate biosynthesis from tetrahomomethionine"
commonnames["PWYQT-4472"] ="glucosinolate biosynthesis from trihomomethionine"
commonnames["PWYQT-4471"] ="glucosinolate biosynthesis from dihomomethionine"
commonnames["PWYQT-4450"] ="aliphatic glucosinolate biosynthesis, side chain elongation cycle"
commonnames["PWY-601"] ="glucosinolate biosynthesis from tryptophan"
commonnames["PWY-2821"] ="glucosinolate biosynthesis from phenylalanine"
commonnames["PWY-1187"] ="glucosinolate biosynthesis from homomethionine"
commonnames["PWY-861"] ="dhurrin biosynthesis"
commonnames["PWY-5990"] ="lotaustralin biosynthesis"
commonnames["PWY-3022"] ="linamarin biosynthesis"
commonnames["PWY-5388"] ="N-glucosylnicotinate metabolism"
commonnames["PWY-5846"] ="colchicine biosynthesis"
commonnames["PWY-5843"] ="cocaine biosynthesis"
commonnames["PWY-5318"] ="calystegine biosynthesis"
commonnames["PWY-5317"] ="hyoscyamine and scopolamine biosynthesis"
commonnames["PWY-5468"] ="lupanine biosynthesis"
commonnames["PWY-6326"] ="camptothecin biosynthesis"
commonnames["PWY-5848"] ="cinchona alkaloids biosynthesis"
commonnames["PWY-6852"] ="senecionine N-oxide biosynthesis"
commonnames["PWY-5752"] ="piperine biosynthesis"
commonnames["PWY-5316"] ="nicotine biosynthesis"
commonnames["PWY-5110"] ="trigonelline biosynthesis"
commonnames["PWY-5038"] ="caffeine biosynthesis II (via paraxanthine)"
commonnames["PWY-5037"] ="caffeine biosynthesis I"
commonnames["PWY-5040"] ="theobromine biosynthesis II (via xanthine)"
commonnames["PWY-5039"] ="theobromine biosynthesis I"
commonnames["PWY-6337"] ="dehydroscoulerine biosynthesis"
commonnames["PWY-6133"] ="(S)-reticuline biosynthesis II"
commonnames["PWY-5876"] ="magnoflorine biosynthesis"
commonnames["PWY-5472"] ="bisbenzylisoquinoline alkaloid biosynthesis"
commonnames["PWY-5471"] ="laudanine biosynthesis"
commonnames["PWY-5470"] ="palmatine biosynthesis"
commonnames["PWY-5287"] ="sanguinarine and macarpine biosynthesis"
commonnames["PWY-5270"] ="morphine biosynthesis"
commonnames["PWY-3901"] ="berberine biosynthesis"
commonnames["PWY-3581"] ="(S)-reticuline biosynthesis I"
commonnames["PWY-6495"] ="ergotamine biosynthesis"
commonnames["PWY-6493"] ="agroclavine biosynthesis"
commonnames["PWY-5877"] ="beta-carboline biosynthesis"
commonnames["PWY-5467"] ="gramine biosynthesis"
commonnames["PWY-5439"] ="betacyanin biosynthesis (via dopamine)"
commonnames["PWY-5426"] ="betaxanthin biosynthesis"
commonnames["PWY-5404"] ="betaxanthin biosynthesis (via dopaxanthin)"
commonnames["PWY-5403"] ="betaxanthin biosynthesis (via dopamine)"
commonnames["PWY-5400"] ="amaranthin biosynthesis"
commonnames["PWY-5399"] ="betacyanin biosynthesis"
commonnames["PWY-5394"] ="betalamic acid biosynthesis"
commonnames["PWY-5301"] ="ajmaline and sarpagine biosynthesis"
commonnames["PWY-5292"] ="vindoline and vinblastine biosynthesis"
commonnames["PWY-5290"] ="secologanin and strictosidine biosynthesis"
commonnames["PWY-6923"] ="ricinine degradation"
commonnames["PWY-6479"] ="fumigaclavine C biosynthesis"
commonnames["PWY-6027"] ="capsiconiate biosynthesis"
commonnames["PWY-5958"] ="acridone alkaloid biosynthesis"
commonnames["PWY-5883"] ="ephedrine biosynthesis"
commonnames["PWY-5748"] ="&gamma;-coniciene and coniine biosynthesis"
commonnames["PWY-5710"] ="capsaicin biosynthesis"
commonnames["PWY-5666"] ="steroidal glycoalkaloid biosynthesis"
commonnames["PWY-5315"] ="N-methyl-&Delta;<sup>1</sup>-pyrrolinium cation biosynthesis"
commonnames["PWY-6448"] ="hordatine biosynthesis"
commonnames["PWY-6442"] ="spermidine hydroxycinnamic acid conjugates biosynthesis"
commonnames["PWY-5474"] ="hydroxycinnamic acid tyramine amides biosynthesis"
commonnames["PWY-5473"] ="hydroxycinnamic acid serotonin amides biosynthesis"
commonnames["PWY-3181"] ="tryptophan degradation VI (via tryptamine)"
commonnames["PWY-5409"] ="divinyl ether biosynthesis II"
commonnames["PWY-5406"] ="divinyl ether biosynthesis I"
commonnames["PWY-5410"] ="traumatin and (<i>Z</i>)-3-hexen-1-yl acetate biosynthesis"
commonnames["PWY-5408"] ="9-lipoxygenase and 9-hydroperoxide lyase pathway"
commonnames["PWY-5407"] ="9-lipoxygenase and 9-allene oxide synthase pathway"
commonnames["PWY1A0-6325"] ="actinorhodin biosynthesis"
commonnames["PWY-6975"] ="erythromycin biosynthesis"
commonnames["PWY-6971"] ="oleandomycin biosynthesis"
commonnames["PWY-6955"] ="lincomycin biosynthesis"
commonnames["PWY-6919"] ="neopentalenoketolactone and pentalenate biosynthesis"
commonnames["PWY-6915"] ="pentalenolactone biosynthesis"
commonnames["PWY-6722"] ="candicidin biosynthesis"
commonnames["PWY-6721"] ="sangivamycin biosynthesis"
commonnames["PWY-6720"] ="toyocamycin biosynthesis"
commonnames["PWY-6682"] ="dehydrophos biosynthesis"
commonnames["PWY-6679"] ="jadomycin biosynthesis"
commonnames["PWY-6666"] ="pyocyanin biosynthesis"
commonnames["PWY-6511"] ="3-methylarginine biosynthesis"
commonnames["PWY-6346"] ="staurosporine biosynthesis"
commonnames["PWY-6345"] ="K-252 biosynthesis"
commonnames["PWY-6324"] ="rebeccamycin biosynthesis"
commonnames["PWY-6322"] ="phosphinothricin tripeptide biosynthesis"
commonnames["PWY-6003"] ="gramicidin S biosynthesis"
commonnames["PWY-5984"] ="rifamycin B biosynthesis"
commonnames["PWY-5940"] ="streptomycin biosynthesis"
commonnames["PWY-5930"] ="terpentecin biosynthesis"
commonnames["PWY-5929"] ="puromycin biosynthesis"
commonnames["PWY-5887"] ="albaflavenone biosynthesis"
commonnames["PWY-5818"] ="validamycin A biosynthesis"
commonnames["PWY-5776"] ="2-hydroxyphenazine biosynthesis"
commonnames["PWY-5770"] ="phenazine-1-carboxylate biosynthesis"
commonnames["PWY-5757"] ="fosfomycin biosynthesis"
commonnames["PWY-5737"] ="(5<i>R</i>)-carbapenem biosynthesis"
commonnames["PWY-5679"] ="clavulanate biosynthesis"
commonnames["PWY-5633"] ="cephamycin C biosynthesis"
commonnames["PWY-5632"] ="cephalosporin C biosynthesis"
commonnames["PWY-5631"] ="deacetylcephalosporin C biosynthesis"
commonnames["PWY-5630"] ="penicillin K biosynthesis"
commonnames["PWY-5629"] ="isopenicillin N biosynthesis"
commonnames["PWY-6950"] ="DIMBOA-glucoside biosynthesis"
commonnames["PWY-6949"] ="DIBOA-glucoside biosynthesis"
commonnames["PWY-6926"] ="pyrethrin I biosynthesis"
commonnames["PWY-6920"] ="6-gingerol analog biosynthesis"
commonnames["PWY-6839"] ="2-aminoethylphosphonate biosynthesis"
commonnames["PWY-6831"] ="pyrrolnitrin biosynthesis"
commonnames["PWY-6808"] ="dTDP-D-forosamine biosynthesis"
commonnames["PWY-6802"] ="salidroside biosynthesis"
commonnames["PWY-6703"] ="preQ<sub>0</sub> biosynthesis"
commonnames["PWY-6699"] ="oxalate biosynthesis"
commonnames["PWY-6661"] ="4-hydroxy-2(1<i>H</i>)-quinolone biosynthesis"
commonnames["PWY-6660"] ="2-heptyl-3-hydroxy-4(1<i>H</i>)-quinolone biosynthesis"
commonnames["PWY-6644"] ="fluoroacetate and fluorothreonine biosynthesis"
commonnames["PWY-6585"] ="methylketone biosynthesis"
commonnames["PWY-6069"] ="indigo biosynthesis"
commonnames["PWY-6068"] ="indican biosynthesis"
commonnames["PWY-5979"] ="3-amino-5-hydroxybenzoate biosynthesis"
commonnames["PWY-5975"] ="furaneol biosynthesis"
commonnames["PWY-5960"] ="aflatoxins B<sub>2</sub> and G<sub>2</sub> biosynthesis"
commonnames["PWY-5959"] ="aflatoxins B<sub>1</sub> and G<sub>1</sub> biosynthesis"
commonnames["PWY-5956"] ="sterigmatocystin biosynthesis"
commonnames["PWY-5955"] ="versicolorin B biosynthesis"
commonnames["PWY-5954"] ="(1'<i>S</i>,5'<i>S</i>)-averufin biosynthesis"
commonnames["PWY-5935"] ="tuberculosinol biosynthesis"
commonnames["PWY-5826"] ="hypoglycin biosynthesis"
commonnames["PWY-5751"] ="phenylethanol biosynthesis"
commonnames["UDPNACETYLGALSYN-PWY"] ="UDP-<i>N</i>-acetyl-D-glucosamine biosynthesis II"
commonnames["PWY-6834"] ="spermidine biosynthesis III"
commonnames["PWY-6559"] ="spermidine biosynthesis II"
commonnames["BSUBPOLYAMSYN-PWY"] ="spermidine biosynthesis I"
commonnames["PWY-6305"] ="putrescine biosynthesis IV"
commonnames["PWY-46"] ="putrescine biosynthesis III"
commonnames["PWY-43"] ="putrescine biosynthesis II"
commonnames["PWY-40"] ="putrescine biosynthesis I"
commonnames["PWY1F-353"] ="glycine betaine biosynthesis III (plants)"
commonnames["PWY-6004"] ="glycine betaine biosynthesis V (from glycine)"
commonnames["PWY-3722"] ="glycine betaine biosynthesis II (Gram-positive bacteria)"
commonnames["P541-PWY"] ="glycine betaine biosynthesis IV (from glycine)"
commonnames["BETSYN-PWY"] ="glycine betaine biosynthesis I (Gram-negative bacteria)"
commonnames["PWY-4021"] ="&beta;-alanine betaine biosynthesis"
commonnames["CHOLINE-BETAINE-ANA-PWY"] ="choline degradation I"
commonnames["TRYPANOSYN-PWY"] ="trypanothione biosynthesis"
commonnames["PWY0-823"] ="arginine degradation III (arginine decarboxylase/agmatinase pathway)"
commonnames["PWY0-1303"] ="aminopropylcadaverine biosynthesis"
commonnames["PWY-6562"] ="norspermidine biosynthesis"
commonnames["PWY-6456"] ="serinol biosynthesis"
commonnames["PWY-6173"] ="histamine biosynthesis"
commonnames["PWY-5907"] ="homospermidine biosynthesis"
commonnames["PWY-5695"] ="urate biosynthesis/inosine 5'-phosphate degradation"
commonnames["P101-PWY"] ="ectoine biosynthesis"
commonnames["GLYCGREAT-PWY"] ="glycine degradation (creatine biosynthesis)"
commonnames["ARGSPECAT-PWY"] ="spermine biosynthesis"
commonnames["PWY0-662"] ="PRPP biosynthesis I"
commonnames["PWY0-661"] ="PRPP biosynthesis II"
commonnames["PWY-5658"] ="mannosylglycerate biosynthesis II"
commonnames["PWY-5656"] ="mannosylglycerate biosynthesis I"
commonnames["PWY-6687"] ="mannosylglucosylglycerate biosynthesis II"
commonnames["PWY-6686"] ="mannosylglucosylglycerate biosynthesis I"
commonnames["PWY-6685"] ="glucosylglycerate biosynthesis II"
commonnames["PWY-5662"] ="glucosylglycerate biosynthesis I"
commonnames["PWY-6055"] ="dimethylsulfoniopropionate biosynthesis II (Spartina)"
commonnames["PWY-6054"] ="dimethylsulfoniopropionate biosynthesis I (Wollastonia)"
commonnames["PWY-6053"] ="dimethylsulfoniopropionate biosynthesis III (algae)"
commonnames["PWY-6664"] ="di-myo-inositol phosphate biosynthesis"
commonnames["PWY-6157"] ="autoinducer AI-1 biosynthesis"
commonnames["PWY-6154"] ="autoinducer AI-2 biosynthesis II (<i>Vibrio</i>)"
commonnames["PWY-6153"] ="autoinducer AI-2 biosynthesis I"
commonnames["PWY1-3"] ="polyhydroxybutyrate biosynthesis"
commonnames["PWY0-1534"] ="hydrogen sulfide biosynthesis"
commonnames["PWY-6936"] ="seleno-amino acid biosynthesis"
commonnames["PWY-6935"] ="seleno-amino acid detoxification and volatilization II"
commonnames["PWY-6931"] ="seleno-amino acid detoxification and volatilization I"
commonnames["PWY-6841"] ="homophytochelatin biosynthesis"
commonnames["PWY-6840"] ="homoglutathione biosynthesis"
commonnames["PWY-6745"] ="phytochelatins biosynthesis"
commonnames["PWY-6730"] ="methylhalides biosynthesis (plants)"
commonnames["PWY-6657"] ="polyhydroxydecanoate biosynthesis"
commonnames["PWY-6498"] ="eumelanin biosynthesis"
commonnames["PWY-6481"] ="L-dopachrome biosynthesis"
commonnames["PWY-5670"] ="epoxysqualene biosynthesis"
commonnames["PWY-5389"] ="methylthiopropionate biosynthesis"
commonnames["PWY-5109"] ="2-methylbutyrate biosynthesis"
commonnames["PWY-6933"] ="seleno-amino acid detoxification and volatilization III"
commonnames["PWY0-181"] ="salvage pathways of pyrimidine deoxyribonucleotides"
commonnames["PWY0-163"] ="salvage pathways of pyrimidine ribonucleotides"
commonnames["PWY-5686"] ="uridine-5'-phosphate biosynthesis"
commonnames["PWY0-166"] ="pyrimidine deoxyribonucleotides <i>de novo</i> biosynthesis I"
commonnames["PWY-6545"] ="pyrimidine deoxyribonucleotides <i>de novo</i> biosynthesis II"
commonnames["PWY-5687"] ="pyrimidine ribonucleotides interconversion"
commonnames["PWY-6620"] ="guanine and guanosine salvage I"
commonnames["PWY-6618"] ="guanine and guanosine salvage III"
commonnames["PWY-6599"] ="guanine and guanosine salvage II"
commonnames["PWY-6619"] ="adenine and adenosine salvage VI"
commonnames["PWY-6611"] ="adenine and adenosine salvage V"
commonnames["PWY-6610"] ="adenine and adenosine salvage IV"
commonnames["PWY-6609"] ="adenine and adenosine salvage III"
commonnames["PWY-6605"] ="adenine and adenosine salvage II"
commonnames["P121-PWY"] ="adenine and adenosine salvage I"
commonnames["P1-PWY"] ="purine and pyrimidine metabolism"
commonnames["SALVPURINE2-PWY"] ="xanthine and xanthosine salvage"
commonnames["PWY-6124"] ="inosine-5'-phosphate biosynthesis II"
commonnames["PWY-6123"] ="inosine-5'-phosphate biosynthesis I"
commonnames["PWY-6126"] ="adenosine nucleotides <i>de novo</i> biosynthesis"
commonnames["PWY-6125"] ="guanosine nucleotides <i>de novo</i> biosynthesis"
commonnames["PWY-6829"] ="tRNA methylation (yeast)"
commonnames["PWY0-1479"] ="tRNA processing pathway I"
commonnames["PWY-6711"] ="archaeosine biosynthesis"
commonnames["PWY-6700"] ="queuosine biosynthesis"
commonnames["PWY-6689"] ="tRNA splicing"
commonnames["PWY-6122"] ="5-aminoimidazole ribonucleotide biosynthesis II"
commonnames["PWY-6121"] ="5-aminoimidazole ribonucleotide biosynthesis I"
commonnames["PWY-6794"] ="adenosine 5'-phosphoramidate biosynthesis"
commonnames["PWY-6845"] ="nitric oxide biosynthesis I (in plants)"
commonnames["PWY-6405"] ="Rapoport-Luebering glycolytic shunt"
commonnames["PWY-6100"] ="L-carnitine biosynthesis"
commonnames["PPGPPMET-PWY"] ="ppGpp biosynthesis"
commonnames["PWY-6158"] ="creatine-phosphate biosynthesis"
commonnames["PWY-6424"] ="sitosterol biosynthesis"
commonnames["PWY-6132"] ="lanosterol biosynthesis"
commonnames["PWY-6075"] ="ergosterol biosynthesis"
commonnames["PWY-6074"] ="zymosterol biosynthesis"
commonnames["PWY-6061"] ="bile acid biosynthesis, neutral pathway"
commonnames["PWY-6036"] ="cardenolide glucosides biosynthesis"
commonnames["PWY-6032"] ="cardenolide biosynthesis"
commonnames["PWY-6580"] ="L-1-phosphatidyl-inositol biosynthesis (Mycobacteria)"
commonnames["PWY4FS-6"] ="phosphatidylethanolamine biosynthesis II"
commonnames["PWY-6273"] ="phosphatidylethanolamine biosynthesis III"
commonnames["PWY-5669"] ="phosphatidylethanolamine biosynthesis I"
commonnames["PWY4FS-4"] ="phosphatidylcholine biosynthesis IV"
commonnames["PWY4FS-3"] ="phosphatidylcholine biosynthesis III"
commonnames["PWY4FS-2"] ="phosphatidylcholine biosynthesis II"
commonnames["PWY3O-450"] ="phosphatidylcholine biosynthesis I"
commonnames["PWY-6826"] ="phosphatidylcholine biosynthesis VI"
commonnames["PWY-6825"] ="phosphatidylcholine biosynthesis V"
commonnames["PWY-5668"] ="cardiolipin biosynthesis I"
commonnames["PWY-5269"] ="cardiolipin biosynthesis II"
commonnames["PWY0-1319"] ="CDP-diacylglycerol biosynthesis II"
commonnames["PWY-5981"] ="CDP-diacylglycerol biosynthesis III"
commonnames["PWY-5667"] ="CDP-diacylglycerol biosynthesis I"
commonnames["PWY-762"] ="phospholipid desaturation"
commonnames["PWY-6367"] ="D-<i>myo</i>-inositol-5-phosphate metabolism"
commonnames["PWY-6352"] ="3-phosphoinositide biosynthesis"
commonnames["PWY-6349"] ="CDP-archaeol biosynthesis"
commonnames["PWY0-541"] ="cyclopropane fatty acid (CFA) biosynthesis"
commonnames["PWY-4942"] ="cyclopropane and cyclopropene fatty acid biosynthesis"
commonnames["PWY-5375"] ="&alpha;-eleostearic acid biosynthesis"
commonnames["PWY-5374"] ="punicic acid biosynthesis"
commonnames["PWY-5373"] ="calendate biosynthesis"
commonnames["PWY-5368"] ="dimorphecolate biosynthesis"
commonnames["PWY-6958"] ="eicosapentaenoate biosynthesis"
commonnames["PWY-5361"] ="&Delta;5-eicosenoate biosynthesis"
commonnames["PWY-5362"] ="&Delta;6-hexadecenoate biosynthesis"
commonnames["PWY-5367"] ="petroselinate biosynthesis"
commonnames["PWY-6013"] ="crepenynic acid biosynthesis"
commonnames["PWY-6014"] ="vernolic acid biosynthesis"
commonnames["PWY-6433"] ="hydroxylated fatty acid biosynthesis (plants)"
commonnames["PWY-6598"] ="sciadonic acid biosynthesis"
commonnames["PWY-6603"] ="dicranin biosynthesis"
commonnames["PWY-5353"] ="arachidonate biosynthesis"
commonnames["PWY-5989"] ="stearate biosynthesis II (plants)"
commonnames["PWY-5972"] ="stearate biosynthesis I (animals)"
commonnames["PWY-6282"] ="palmitoleate biosynthesis I"
commonnames["PWY-5366"] ="palmitoleate biosynthesis II"
commonnames["PWY-5994"] ="palmitate biosynthesis I (animals)"
commonnames["PWY-5971"] ="palmitate biosynthesis II (bacteria and plants)"
commonnames["PWY-6429"] ="ricinoleate biosynthesis"
commonnames["PWY-5996"] ="oleate biosynthesis II (animals)"
commonnames["PWY-5147"] ="oleate biosynthesis I (plants)"
commonnames["PWY-6001"] ="linoleate biosynthesis II (animals)"
commonnames["PWY-5995"] ="linoleate biosynthesis I (plants)"
commonnames["PWY-6000"] ="&gamma;-linolenate biosynthesis II (animals)"
commonnames["PWY-5998"] ="&gamma;-linolenate biosynthesis I (plants)"
commonnames["PWYG-321"] ="mycolate biosynthesis"
commonnames["PWY0-862"] ="<i>cis</i>-dodecenoyl biosynthesis"
commonnames["PWY-6917"] ="epoxy fatty acid biosynthesis"
commonnames["PWY-6799"] ="fatty acid biosynthesis (plant mitochondria)"
commonnames["PWY-6710"] ="poly-hydroxy fatty acids biosynthesis"
commonnames["PWY-6468"] ="&omega;- hydroxylation of laurate"
commonnames["PWY-6465"] ="&omega;-hydroxylation of caprate and laurate"
commonnames["PWY-5997"] ="&alpha;-linolenate biosynthesis"
commonnames["PWY-5973"] ="<i>cis</i>-vaccenate biosynthesis"
commonnames["PWY-5970"] ="fatty acids biosynthesis (yeast)"
commonnames["PWY-5966"] ="fatty acid biosynthesis initiation II"
commonnames["PWY-5965"] ="fatty acid biosynthesis initiation III"
commonnames["PWY-5080"] ="very long chain fatty acid biosynthesis"
commonnames["PWY-4381"] ="fatty acid biosynthesis initiation I"
commonnames["FASYN-ELONG-PWY"] ="fatty acid elongation -- saturated"
commonnames["PWY-3561"] ="choline biosynthesis III"
commonnames["PWY-3542"] ="choline biosynthesis II"
commonnames["PWY-3385"] ="choline biosynthesis I"
commonnames["TRIGLSYN-PWY"] ="triacylglycerol biosynthesis"
commonnames["SPHINGOLIPID-SYN-PWY"] ="sphingolipid metabolism"
commonnames["SOPHOROSYLOXYDOCOSANOATE-SYN-PWY"] ="sophorosyloxydocosanoate biosynthesis"
commonnames["PWYQT-4427"] ="sulfolipid biosynthesis"
commonnames["PWY3DJ-12"] ="ceramide biosynthesis"
commonnames["PWY3DJ-11281"] ="sphingomyelin metabolism"
commonnames["PWY0-1264"] ="biotin-carboxyl carrier protein assembly"
commonnames["PWY-782"] ="glycolipid desaturation"
commonnames["PWY-6818"] ="ornithine lipid biosynthesis"
commonnames["PWY-6804"] ="diacylglycerol biosynthesis (PUFA enrichment in oilseed)"
commonnames["PWY-6803"] ="phosphatidylcholine acyl editing"
commonnames["PWY-6453"] ="stigma estolide biosynthesis"
commonnames["PWY-6129"] ="dolichol and dolichyl phosphate biosynthesis"
commonnames["PWY-5885"] ="wax esters biosynthesis II"
commonnames["PWY-5884"] ="wax esters biosynthesis I"
commonnames["PWY-5148"] ="acyl-CoA hydrolysis"
commonnames["PWY-5142"] ="acyl-ACP thioesterase pathway"
commonnames["PWY-5129"] ="sphingolipid biosynthesis (plants)"
commonnames["PWY-401"] ="glycolipid biosynthesis"
commonnames["LIPA-CORESYN-PWY"] ="Lipid A-core biosynthesis"
commonnames["KDO-LIPASYN-PWY"] ="(KDO)<sub>2</sub>-lipid A biosynthesis I"
commonnames["PWY-6795"] ="diacylglyceryl-<i>N,N,N</i>-trimethylhomoserine biosynthesis"
commonnames["PWY-6648"] ="rhamnolipid biosynthesis"
commonnames["PWY-735"] ="jasmonic acid biosynthesis"
commonnames["PWY-6297"] ="tuberonate glucoside biosynthesis"
commonnames["PWY-6235"] ="hydroxyjasmonate sulfate biosynthesis"
commonnames["PWY-6233"] ="jasmonoyl-amino acid conjugates biosynthesis II"
commonnames["PWY-6220"] ="jasmonoyl-amino acid conjugates biosynthesis I"
commonnames["PWY-6854"] ="ethylene biosynthesis III (microbes)"
commonnames["PWY-6853"] ="ethylene biosynthesis II (microbes)"
commonnames["ETHYL-PWY"] ="ethylene biosynthesis I (plants)"
commonnames["PWY-5967"] ="lupinate biosynthesis"
commonnames["PWY-2902"] ="cytokinins-<i>O</i>-glucoside biosynthesis"
commonnames["PWY-2901"] ="cytokinins 9-<i>N</i>-glucoside biosynthesis"
commonnames["PWY-2881"] ="cytokinins 7-<i>N</i>-glucoside biosynthesis"
commonnames["PWY-2781"] ="<i>cis</i>-zeatin biosynthesis"
commonnames["PWY-2681"] ="<i>trans</i>-zeatin biosynthesis"
commonnames["PWY-699"] ="brassinosteroid biosynthesis I"
commonnames["PWY-2582"] ="brassinosteroid biosynthesis II"
commonnames["TRPIAACAT-PWY"] ="tryptophan degradation VII (via indole-3-pyruvate)"
commonnames["PWY-6219"] ="indole-3-acetyl-amino acid biosynthesis"
commonnames["PWY-581"] ="IAA biosynthesis I"
commonnames["PWY-5026"] ="IAA biosynthesis V"
commonnames["PWY-5025"] ="IAA biosynthesis IV"
commonnames["PWY-3161"] ="IAA biosynthesis VI (via indole-3-acetamide)"
commonnames["PWY-1921"] ="IAA biosynthesis III"
commonnames["PWY-1822"] ="IAA biosynthesis II"
commonnames["PWY-1741"] ="IAA conjugate biosynthesis I"
commonnames["PWY-5271"] ="phaseic acid biosynthesis"
commonnames["PWY-695"] ="abscisic acid biosynthesis"
commonnames["PWY-6299"] ="aldehyde oxidation I"
commonnames["PWY-6653"] ="<i>ent</i> -kaurene biosynthesis II"
commonnames["PWY-5070"] ="gibberellin biosynthesis I (non C-3, non C-13 hydroxylation)"
commonnames["PWY-5047"] ="gibberellin biosynthesis IV (<i>Gibberella fujikuroi</i>)"
commonnames["PWY-5036"] ="gibberellin biosynthesis II (early C-3 hydroxylation)"
commonnames["PWY-5035"] ="gibberellin biosynthesis III (early C-13 hydroxylation)"
commonnames["PWY-5034"] ="GA<sub>12</sub> biosynthesis"
commonnames["PWY-5032"] ="<i>ent</i>-kaurene biosynthesis I"
commonnames["PWY66-380"] ="estrogen biosynthesis"
commonnames["PWY66-378"] ="androgen biosynthesis"
commonnames["PWY66-377"] ="pregnenolone biosynthesis"
commonnames["PWY66-375"] ="leukotriene biosynthesis"
commonnames["PWY66-301"] ="catecholamine biosynthesis"
commonnames["PWY-6650"] ="juvenile hormone III biosynthesis II"
commonnames["PWY-6575"] ="juvenile hormone III biosynthesis I"
commonnames["PWY-6241"] ="thyroid hormone biosynthesis"
commonnames["PWY-6030"] ="serotonin and melatonin biosynthesis"
commonnames["PYRIDOXSYN-PWY"] ="pyridoxal 5'-phosphate biosynthesis I"
commonnames["PWY-6466"] ="pyridoxal 5'-phosphate biosynthesis II"
commonnames["PLPSAL-PWY"] ="pyridoxal 5'-phosphate salvage pathway"
commonnames["PWY-6875"] ="retinoate biosynthesis II"
commonnames["PWY-6872"] ="retinoate biosynthesis I"
commonnames["PWY-6861"] ="the visual cycle"
commonnames["PWY-6857"] ="retinol biosynthesis"
commonnames["PWY-6910"] ="hydroxymethylpyrimidine salvage"
commonnames["PWY-6899"] ="base-degraded thiamin salvage"
commonnames["PWY-6898"] ="thiamin salvage III"
commonnames["PWY-6896"] ="thiamin salvage I"
commonnames["PWY-6909"] ="thiazole biosynthesis III (eukaryotes)"
commonnames["PWY-6908"] ="thiamin diphosphate biosynthesis IV (eukaryotes)"
commonnames["PWY-6907"] ="thiamin diphosphate biosynthesis III (Staphylococcus)"
commonnames["PWY-6894"] ="thiamin diphosphate biosynthesis I (E. coli)"
commonnames["PWY-6893"] ="thiamin diphosphate biosynthesis II (Bacillus)"
commonnames["PWY-6892"] ="thiazole biosynthesis I (E. coli)"
commonnames["PWY-6890"] ="4-amino-2-methyl-5-diphosphomethylpyrimidine biosynthesis"
commonnames["PWY-6891"] ="thiazole biosynthesis II (Bacillus)"
commonnames["PWY-6654"] ="phosphopantothenate biosynthesis III"
commonnames["PWY-3961"] ="phosphopantothenate biosynthesis II"
commonnames["PANTO-PWY"] ="phosphopantothenate biosynthesis I"
commonnames["PWY-6797"] ="6-hydroxymethyl-dihydropterin diphosphate biosynthesis II (archaea)"
commonnames["PWY-6147"] ="6-hydroxymethyl-dihydropterin diphosphate biosynthesis I"
commonnames["PWY-6614"] ="tetrahydrofolate biosynthesis"
commonnames["PWY-6613"] ="tetrahydrofolate salvage from 5,10-methenyltetrahydrofolate"
commonnames["PWY-6543"] ="p-aminobenzoate biosynthesis"
commonnames["PWY-3841"] ="folate transformations II"
commonnames["PWY-2201"] ="folate transformations I"
commonnames["PWY-2161B"] ="glutamate removal from folates"
commonnames["PWY-2161"] ="folate polyglutamylation"
commonnames["1CMET2-PWY"] ="formylTHF biosynthesis I"
commonnames["RIBOSYN2-PWY"] ="flavin biosynthesis I (bacteria and plants)"
commonnames["PWY66-366"] ="flavin biosynthesis IV (mammalian)"
commonnames["PWY-6168"] ="flavin biosynthesis III (fungi)"
commonnames["PWY-6167"] ="flavin biosynthesis II (archaea)"
commonnames["PWY-6269"] ="adenosylcobalamin salvage from cobinamide II"
commonnames["PWY-6268"] ="adenosylcobalamin salvage from cobalamin"
commonnames["COBALSYN-PWY"] ="adenosylcobalamin salvage from cobinamide I"
commonnames["PWY-5523"] ="5,6-dimethylbenzimidazole biosynthesis"
commonnames["PWY-5509"] ="adenosylcobalamin biosynthesis from cobyrinate <i>a,c</i>-diamide I"
commonnames["PWY-5508"] ="adenosylcobalamin biosynthesis from cobyrinate <i>a,c</i>-diamide II"
commonnames["PWY-5443"] ="aminopropanol phosphate biosynthesis"
commonnames["PWY-6578"] ="7-keto-8-aminopelargonate biosynthesis II"
commonnames["PWY-6519"] ="7-keto-8-aminopelargonate biosynthesis I"
commonnames["PWY0-1507"] ="biotin biosynthesis from 7-keto-8-aminopelargonate"
commonnames["PWY4FS-13"] ="extended VTC2 cycle"
commonnames["PWY4FS-12"] ="VTC2 cycle"
commonnames["PWY4FS-11"] ="ascorbate biosynthesis II (L-gulose pathway)"
commonnames["PWY3DJ-35471"] ="L-ascorbate biosynthesis VI"
commonnames["PWY-882"] ="ascorbate biosynthesis I (L-galactose pathway)"
commonnames["PWY-6415"] ="ascorbate biosynthesis VII"
commonnames["PWY-5521"] ="L-ascorbate biosynthesis V"
commonnames["PWY-6076"] ="1,25-dihydroxyvitamin D<sub>3</sub> biosynthesis"
commonnames["PWY-5189"] ="tetrapyrrole biosynthesis II"
commonnames["PWY-5188"] ="tetrapyrrole biosynthesis I"
commonnames["PWY-5664"] ="tetrahydrobiopterin biosynthesis II"
commonnames["PWY-5663"] ="tetrahydrobiopterin biosynthesis I"
commonnames["PWY-4121"] ="glutathionylspermidine biosynthesis"
commonnames["THIOREDOX-PWY"] ="thioredoxin pathway"
commonnames["PWY8J2-1"] ="bacillithiol biosynthesis"
commonnames["PWY1G-126"] ="mycothiol oxidation"
commonnames["PWY1G-0"] ="mycothiol biosynthesis"
commonnames["PWY-4181"] ="glutathione amide metabolism"
commonnames["PWY-4081"] ="glutathione redox reactions I"
commonnames["GLUTATHIONESYN-PWY"] ="glutathione biosynthesis"
commonnames["GLUT-REDOX-PWY"] ="glutathione redox reactions II"
commonnames["PWY3O-19"] ="ubiquinol-6 biosynthesis (eukaryotic)"
commonnames["PWY-6708"] ="ubiquinol-8 biosynthesis (prokaryotic)"
commonnames["PWY-5873"] ="ubiquinol-7 biosynthesis (eukaryotic)"
commonnames["PWY-5872"] ="ubiquinol-10 biosynthesis (eukaryotic)"
commonnames["PWY-5871"] ="ubiquinol-9 biosynthesis (eukaryotic)"
commonnames["PWY-5870"] ="ubiquinol-8 biosynthesis (eukaryotic)"
commonnames["PWY-5857"] ="ubiquinol-10 biosynthesis (prokaryotic)"
commonnames["PWY-5856"] ="ubiquinol-9 biosynthesis (prokaryotic)"
commonnames["PWY-5855"] ="ubiquinol-7 biosynthesis (prokaryotic)"
commonnames["PWY-5889"] ="rhodoquinone-9 biosynthesis"
commonnames["PWY-5888"] ="rhodoquinone-10 biosynthesis"
commonnames["PWY-6978"] ="plastoquinol-9 biosynthesis II"
commonnames["PWY-1581"] ="plastoquinol-9 biosynthesis I"
commonnames["PWY-5027"] ="phylloquinol biosynthesis"
commonnames["PWY-5895"] ="menaquinol-13 biosynthesis"
commonnames["PWY-5892"] ="menaquinol-12 biosynthesis"
commonnames["PWY-5891"] ="menaquinol-11 biosynthesis"
commonnames["PWY-5890"] ="menaquinol-10 biosynthesis"
commonnames["PWY-5849"] ="menaquinol-6 biosynthesis"
commonnames["PWY-5844"] ="menaquinol-9 biosynthesis"
commonnames["PWY-5839"] ="menaquinol-7 biosynthesis"
commonnames["MENAQUINONESYN-PWY"] ="menaquinol-8 biosynthesis"
commonnames["PWY-6793"] ="demethylmenaquinol-8 biosynthesis III"
commonnames["PWY-6262"] ="demethylmenaquinol-8 biosynthesis II"
commonnames["PWY-5852"] ="demethylmenaquinol-8 biosynthesis I"
commonnames["PWY-5853"] ="demethylmenaquinol-6 biosynthesis"
commonnames["PWY-5851"] ="demethylmenaquinol-9 biosynthesis"
commonnames["PWY-5837"] ="1,4-dihydroxy-2-naphthoate biosynthesis I"
commonnames["PWY-5791"] ="1,4-dihydroxy-2-naphthoate biosynthesis II (plants)"
commonnames["PWY-1422"] ="vitamin E biosynthesis"
commonnames["HEMESYN2-PWY"] ="heme biosynthesis from uroporphyrinogen-III II"
commonnames["HEME-BIOSYNTHESIS-II"] ="heme biosynthesis from uroporphyrinogen-III I"
commonnames["PWY-5531"] ="chlorophyllide <i>a</i> biosynthesis II"
commonnames["CHLOROPHYLL-SYN"] ="chlorophyllide <i>a</i> biosynthesis I"
commonnames["PWY-5086"] ="chlorophyll <i>a</i> biosynthesis I"
commonnames["PWY-5068"] ="chlorophyll cycle"
commonnames["PWY-5064"] ="chlorophyll <i>a</i> biosynthesis II"
commonnames["PWY-5526"] ="bacteriochlorophyll <i>a</i> biosynthesis"
commonnames["PWY-5297"] ="siroheme amide biosynthesis"
commonnames["PWY-5194"] ="siroheme biosynthesis"
commonnames["PWY-5120"] ="geranylgeranyldiphosphate biosynthesis"
commonnames["PWY-5123"] ="<i>trans, trans</i>-farnesyl diphosphate biosynthesis"
commonnames["PWY-6577"] ="farnesylcysteine salvage pathway"
commonnames["PWY-6520"] ="nonaprenyl diphosphate biosynthesis II"
commonnames["PWY-6383"] ="mono-<i>trans</i>, poly-<i>cis</i> decaprenyl phosphate biosynthesis"
commonnames["PWY-5893"] ="tridecaprenyl diphosphate biosynthesis"
commonnames["PWY-5817"] ="dodecaprenyl diphosphate biosynthesis"
commonnames["PWY-5816"] ="all <i>trans</i> undecaprenyl diphosphate biosynthesis"
commonnames["PWY-5807"] ="heptaprenyl diphosphate biosynthesis"
commonnames["PWY-5806"] ="all-<i>trans</i>-decaprenyl diphosphate biosynthesis"
commonnames["PWY-5805"] ="nonaprenyl diphosphate biosynthesis I"
commonnames["PWY-5785"] ="di-<i>trans</i>,poly-<i>cis</i>-undecaprenyl phosphate biosynthesis"
commonnames["PWY-5783"] ="octaprenyl diphosphate biosynthesis"
commonnames["PWY-5122"] ="geranyl diphosphate biosynthesis"
commonnames["NONMEVIPP-PWY"] ="methylerythritol phosphate pathway"
commonnames["HEXPPSYN-PWY"] ="hexaprenyl diphosphate biosynthesis"
commonnames["PYRIDNUCSYN-PWY"] ="NAD biosynthesis I (from aspartate)"
commonnames["PYRIDNUCSAL-PWY"] ="NAD salvage pathway I"
commonnames["PWY3O-4106"] ="NAD salvage pathway III"
commonnames["PWY-5653"] ="NAD biosynthesis from 2-amino-3-carboxymuconate semialdehyde"
commonnames["NAD-BIOSYNTHESIS-III"] ="NAD biosynthesis III"
commonnames["PWY-5381"] ="pyridine nucleotide cycling (plants)"
commonnames["PWY-5083"] ="NAD/NADH phosphorylation and dephosphorylation"
commonnames["NADPHOS-DEPHOS-PWY"] ="NAD phosphorylation and dephosphorylation"
commonnames["PWY-6823"] ="molybdenum cofactor biosynthesis"
commonnames["PWY0-522"] ="lipoate salvage and modification"
commonnames["PWY0-501"] ="lipoate biosynthesis and incorporation I"
commonnames["PWY0-1275"] ="lipoate biosynthesis and incorporation II"
commonnames["PWY-6643"] ="coenzyme M biosynthesis II"
commonnames["P261-PWY"] ="coenzyme M biosynthesis I"
commonnames["PWY-5796"] ="2'-(5'-phosphoribosyl)-3'-dephospho-CoA biosynthesis II (malonate decarboxylase)"
commonnames["P2-PWY"] ="2'-(5'-phosphoribosyl)-3'-dephospho-CoA biosynthesis I (citrate lyase)"
commonnames["COA-PWY"] ="coenzyme A biosynthesis"
commonnames["SAM-PWY"] ="S-adenosyl-L-methionine biosynthesis"
commonnames["PWY0-1433"] ="tetrahydromonapterin biosynthesis"
commonnames["PWY-922"] ="mevalonate pathway I"
commonnames["PWY-6476"] ="cytidylyl molybdenum cofactor biosynthesis"
commonnames["PWY-6420"] ="pyrroloquinoline quinone biosynthesis"
commonnames["PWY-6174"] ="mevalonate pathway II (archaea)"
commonnames["PWY-6148"] ="tetrahydromethanopterin biosynthesis"
commonnames["PWY-6012"] ="acyl carrier protein metabolism"
commonnames["PWY-5964"] ="guanylyl molybdenum cofactor biosynthesis"
commonnames["PWY-5963"] ="thio-molybdenum cofactor biosynthesis"
commonnames["PWY-5917"] ="phycocyanobilin biosynthesis"
commonnames["PWY-5915"] ="phycoerythrobilin biosynthesis"
commonnames["PWY-5499"] ="vitamin B<sub>6</sub> degradation"
commonnames["PWY-5254"] ="methanofuran biosynthesis"
commonnames["PWY-5207"] ="coenzyme B/coenzyme M regeneration"
commonnames["PWY-5199"] ="factor 420 polyglutamylation"
commonnames["PWY-5198"] ="factor 420 biosynthesis"
commonnames["PWY-5196"] ="factor 430 biosynthesis"
commonnames["P241-PWY"] ="coenzyme B biosynthesis"
commonnames["PWY-5800"] ="xylan biosynthesis"
commonnames["PWY-361"] ="phenylpropanoid biosynthesis"
commonnames["PWY-321"] ="cutin biosynthesis"
commonnames["PWY-282"] ="cuticular wax biosynthesis"
commonnames["PWY-6733"] ="sporopollenin precursor biosynthesis"
commonnames["PWY-1121"] ="suberin biosynthesis"
commonnames["UDPNAGSYN-PWY"] ="UDP-<i>N</i>-acetyl-D-glucosamine biosynthesis I"
commonnames["NAGLIPASYN-PWY"] ="lipid IV<sub>A</sub> biosynthesis"
commonnames["PWY0-1338"] ="polymyxin resistance"
commonnames["KDOSYN-PWY"] ="KDO transfer to lipid IV<sub>A</sub> I"
commonnames["PWY-6463"] ="peptidoglycan cross-bridge biosynthesis IV (Weissella viridescens)"
commonnames["PWY-6462"] ="peptidoglycan cross-bridge biosynthesis III (Enterococcus faecalis)"
commonnames["PWY-6461"] ="peptidoglycan cross-bridge biosynthesis II (E. faecium)"
commonnames["PWY-6459"] ="peptidoglycan cross-bridge biosynthesis I (S. aureus)"
commonnames["PWY-6455"] ="vancomycin resistance II"
commonnames["PWY-6454"] ="vancomycin resistance I"
commonnames["PWY-6397"] ="mycolyl-arabinogalactan-peptidoglycan complex biosynthesis"
commonnames["PWY-6387"] ="UDP-<i>N</i>-acetylmuramoyl-pentapeptide biosynthesis III (<i>meso</i>-DAP-containing)"
commonnames["PWY-6386"] ="UDP-<i>N</i>-acetylmuramoyl-pentapeptide biosynthesis II (lysine-containing)"
commonnames["TEICHOICACID-PWY"] ="teichoic acid (poly-glycerol) biosynthesis"
commonnames["PWY-6626"] ="CDP-2-glycerol biosynthesis"
commonnames["PWY-5514"] ="UDP-<i>N</i>-acetyl-D-galactosamine biosynthesis II"
commonnames["PWY-5512"] ="UDP-<i>N</i>-acetyl-D-galactosamine biosynthesis I"
commonnames["ECASYN-PWY"] ="enterobacterial common antigen biosynthesis"
commonnames["PWY-6568"] ="dermatan sulfate biosynthesis (late stages)"
commonnames["PWY-6567"] ="chondroitin sulfate biosynthesis (late stages)"
commonnames["PWY-6566"] ="chondroitin and dermatan biosynthesis"
commonnames["PWY-6558"] ="heparan sulfate biosynthesis (late stages)"
commonnames["PWY-6557"] ="glycoaminoglycan-protein linkage region biosynthesis"
commonnames["PWY-622"] ="starch biosynthesis"
commonnames["PWY-5067"] ="glycogen biosynthesis II (from UDP-D-Glucose)"
commonnames["GLYCOGENSYNTH-PWY"] ="glycogen biosynthesis I (from ADP-D-Glucose)"
commonnames["PWY-6082"] ="alginate biosynthesis II"
commonnames["PWY-6073"] ="alginate biosynthesis I"
commonnames["PWY-6773"] ="callose biosynthesis"
commonnames["PWY-6658"] ="acetan biosynthesis"
commonnames["PWY-6655"] ="xanthan biosynthesis"
commonnames["PWY-6403"] ="carrageenan biosynthesis"
commonnames["PWY-5980"] ="xylogalacturonan biosynthesis"
commonnames["PWY-5936"] ="xyloglucan biosynthesis"
commonnames["PWY-6778"] ="laminaribiose biosynthesis"
commonnames["PWY-6525"] ="stellariose and mediose biosynthesis"
commonnames["PWY-6524"] ="lychnose biosynthesis"
commonnames["PWY-6116"] ="mannosylfructose biosynthesis"
commonnames["PWY-5343"] ="ajugose biosynthesis II (galactinol-independent)"
commonnames["PWY-5342"] ="ajugose biosynthesis I (galactinol-dependent)"
commonnames["PWY-5337"] ="stachyose biosynthesis"
commonnames["MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS"] ="dolichyl-diphosphooligosaccharide biosynthesis"
commonnames["TRESYN-PWY"] ="trehalose biosynthesis I"
commonnames["TREHALOSESYN-PWY"] ="trehalose biosynthesis III"
commonnames["PWY-881"] ="trehalose biosynthesis II"
commonnames["PWY-5985"] ="trehalose biosynthesis VII"
commonnames["PWY-5983"] ="trehalose biosynthesis VI"
commonnames["PWY-2661"] ="trehalose biosynthesis V"
commonnames["PWY-2622"] ="trehalose biosynthesis IV"
commonnames["PWY-5530"] ="sorbitol biosynthesis II"
commonnames["PWY-5054"] ="sorbitol biosynthesis I"
commonnames["PWY-6739"] ="pinitol biosynthesis II"
commonnames["PWY-6738"] ="pinitol biosynthesis I"
commonnames["PWY-3881"] ="mannitol biosynthesis"
commonnames["PWY-3261"] ="UDP-L-rhamnose biosynthesis"
commonnames["PWY-3221"] ="dTDP-L-rhamnose biosynthesis II"
commonnames["GDPRHAMSYN-PWY"] ="GDP-D-rhamnose biosynthesis"
commonnames["DTDPRHAMSYN-PWY"] ="dTDP-L-rhamnose biosynthesis I"
commonnames["PWY-5659"] ="GDP-mannose biosynthesis"
commonnames["PWY-6139"] ="CMP-<i>N</i>-acetylneuraminate biosynthesis II (bacteria)"
commonnames["PWY-6138"] ="CMP-<i>N</i>-acetylneuraminate biosynthesis I (eukaryotes)"
commonnames["PWY0-1241"] ="ADP-L-<i>glycero</i>-&beta;-D-<i>manno</i>-heptose biosynthesis"
commonnames["PWY-82"] ="UDP-L-arabinose biosynthesis II (from L-arabinose)"
commonnames["PWY-6976"] ="dTDP-L-mycarose biosynthesis"
commonnames["PWY-6974"] ="dTDP-L-olivose biosynthesis"
commonnames["PWY-6973"] ="dTDP-D-olivose, dTDP-D-oliose and dTDP-D-mycarose biosynthesis"
commonnames["PWY-6953"] ="dTDP-3-acetamido-3,6-dideoxy-&alpha;-D-galactose biosynthesis"
commonnames["PWY-6942"] ="dTDP-D-desosamine biosynthesis"
commonnames["PWY-6749"] ="CMP-legionaminate biosynthesis"
commonnames["PWY-66"] ="GDP-L-fucose biosynthesis I (from GDP-D-mannose)"
commonnames["PWY-6478"] ="GDP-D-<i>glycero</i>-&alpha;-D-<i>manno</i>-heptose biosynthesis"
commonnames["PWY-63"] ="UDP-L-arabinose biosynthesis I (from UDP-xylose)"
commonnames["PWY-6144"] ="CMP-<i>N</i>-glycoloylneuraminate biosynthesis"
commonnames["PWY-6140"] ="CMP-2-keto-3-deoxy-D-<i>glycero</i>-D-<i>galacto</i>-nononate biosynthesis"
commonnames["PWY-6"] ="GDP-L-fucose biosynthesis II (from L-fucose)"
commonnames["PWY-5834"] ="CDP-tyvelose biosynthesis"
commonnames["PWY-5833"] ="CDP-3,6-dideoxyhexose biosynthesis"
commonnames["PWY-5832"] ="CDP-paratose biosynthesis"
commonnames["PWY-5831"] ="CDP-abequose biosynthesis"
commonnames["PWY-5830"] ="CDP-ascarylose biosynthesis"
commonnames["PWY-5740"] ="GDP-L-colitose biosynthesis"
commonnames["PWY-5739"] ="GDP-&alpha;-D-perosamine biosynthesis"
commonnames["PWY-5738"] ="GDP-6-deoxy-D-talose biosynthesis"
commonnames["PWY-5661"] ="GDP-glucose biosynthesis"
commonnames["PWY-5115"] ="GDP-L-galactose biosynthesis"
commonnames["PWY-3821"] ="galactose degradation III"
commonnames["PWY-5111"] ="CMP-KDO biosynthesis II (from D-arabinose 5-phosphate)"
commonnames["PWY-4861"] ="UDP-D-galacturonate biosynthesis I (from UDP-D-glucuronate)"
commonnames["PWY-4841"] ="UDP-D-glucuronate biosynthesis (from <i>myo</i>-inositol)"
commonnames["PWY-4821"] ="UDP-D-xylose and UDP-D-glucuronate biosynthesis"
commonnames["PWY-4"] ="UDP-D-galacturonate biosynthesis II (from D-galacturonate)"
commonnames["PWY-1269"] ="CMP-KDO biosynthesis I"
commonnames["PWY-5113"] ="UDP-D-apiose biosynthesis (from UDP-D-glucuronate)"
commonnames["GLUCONEO-PWY"] ="gluconeogenesis I"
commonnames["SUCSYN-PWY"] ="sucrose biosynthesis"
commonnames["PWY-842"] ="starch degradation I"
commonnames["PWY-822"] ="fructan biosynthesis"
commonnames["PWY-6724"] ="starch degradation II"
commonnames["PWY-6143"] ="CMP-pseudaminate biosynthesis"
commonnames["PWY-1061"] ="homogalacturonan biosynthesis"
commonnames["PWY-1001"] ="cellulose biosynthesis"
commonnames["GLYCOCAT-PWY"] ="glycogen degradation I"
commonnames["CALVIN-PWY"] ="Calvin-Benson-Bassham cycle"
commonnames["PWY-862"] ="fructan degradation"
commonnames["PWY-6622"] ="heptadecane biosynthesis"
commonnames["PWY-6333"] ="acetaldehyde biosynthesis I"
commonnames["PWY-6330"] ="acetaldehyde biosynthesis II"
commonnames["PWY-5782"] ="2-keto-L-gulonate biosynthesis"
commonnames["PWY-5750"] ="itaconate biosynthesis"
commonnames["TRNA-CHARGING-PWY"] ="tRNA charging"
commonnames["VALSYN-PWY"] ="valine biosynthesis"
commonnames["TYRSYN"] ="tyrosine biosynthesis I"
commonnames["PWY-6134"] ="tyrosine biosynthesis IV"
commonnames["PWY-6120"] ="tyrosine biosynthesis III"
commonnames["PWY-3461"] ="tyrosine biosynthesis II"
commonnames["TRPSYN-PWY"] ="tryptophan biosynthesis"
commonnames["HOMOSER-THRESYN-PWY"] ="threonine biosynthesis from homoserine"
commonnames["PWY0-901"] ="selenocysteine biosynthesis I (bacteria)"
commonnames["PWY-6281"] ="selenocysteine biosynthesis II (archaea and eukaryotes)"
commonnames["SERSYN-PWY"] ="serine biosynthesis"
commonnames["PWY-6196"] ="serine racemization"
commonnames["PWY-4981"] ="proline biosynthesis II (from arginine)"
commonnames["PWY-4281"] ="proline biosynthesis IV"
commonnames["PWY-3341"] ="proline biosynthesis III"
commonnames["PROSYN-PWY"] ="proline biosynthesis I"
commonnames["ARG-PRO-PWY"] ="arginine degradation VI (arginase 2 pathway)"
commonnames["PWY-3462"] ="phenylalanine biosynthesis II"
commonnames["PHESYN"] ="phenylalanine biosynthesis I"
commonnames["PWY-4921"] ="protein citrullination"
commonnames["CITRULBIO-PWY"] ="citrulline biosynthesis"
commonnames["PWY-6922"] ="<i>L-N<sup>&delta;</sup></i>-acetylornithine biosynthesis"
commonnames["PWY-5957"] ="nicotianamine biosynthesis"
commonnames["PWY-5344"] ="homocysteine biosynthesis"
commonnames["PWY-5331"] ="taurine biosynthesis"
commonnames["PWY-1186"] ="homomethionine biosynthesis"
commonnames["HOMOSERSYN-PWY"] ="homoserine biosynthesis"
commonnames["GLUTORN-PWY"] ="ornithine biosynthesis"
commonnames["PWY-5041"] ="<i>S</i>-adenosyl-L-methionine cycle II"
commonnames["PWY-702"] ="methionine biosynthesis II"
commonnames["PWY-5441"] ="S-methylmethionine cycle"
commonnames["HSERMETANA-PWY"] ="methionine biosynthesis III"
commonnames["HOMOSER-METSYN-PWY"] ="methionine biosynthesis I"
commonnames["ADENOSYLHOMOCYSCAT-PWY"] ="methionine salvage II (mammalia)"
commonnames["PWY-6755"] ="<i>S</i>-methyl-5-thio-&alpha;-D-ribose 1-phosphate degradation"
commonnames["PWY-5097"] ="lysine biosynthesis VI"
commonnames["PWY-3081"] ="lysine biosynthesis V"
commonnames["PWY-2942"] ="lysine biosynthesis III"
commonnames["PWY-2941"] ="lysine biosynthesis II"
commonnames["LYSINE-AMINOAD-PWY"] ="lysine biosynthesis IV"
commonnames["DAPLYSINESYN-PWY"] ="lysine biosynthesis I"
commonnames["LEUSYN-PWY"] ="leucine biosynthesis"
commonnames["PWY-5108"] ="isoleucine biosynthesis V"
commonnames["PWY-5104"] ="isoleucine biosynthesis IV"
commonnames["PWY-5103"] ="isoleucine biosynthesis III"
commonnames["PWY-5101"] ="isoleucine biosynthesis II"
commonnames["ILEUSYN-PWY"] ="isoleucine biosynthesis I (from threonine)"
commonnames["HISTSYN-PWY"] ="histidine biosynthesis"
commonnames["GLYSYN-THR-PWY"] ="glycine biosynthesis IV"
commonnames["GLYSYN-PWY"] ="glycine biosynthesis I"
commonnames["GLYSYN-ALA-PWY"] ="glycine biosynthesis III"
commonnames["GLYCINE-SYN2-PWY"] ="glycine biosynthesis II"
commonnames["PWY-6549"] ="glutamine biosynthesis III"
commonnames["PWY-5921"] ="L-glutamine biosynthesis II (tRNA-dependent)"
commonnames["GLNSYN-PWY"] ="glutamine biosynthesis I"
commonnames["PWY-4341"] ="glutamate biosynthesis V"
commonnames["GLUTSYNIII-PWY"] ="glutamate biosynthesis III"
commonnames["GLUTSYN-PWY"] ="glutamate biosynthesis I"
commonnames["GLUTAMATE-SYN2-PWY"] ="glutamate biosynthesis II"
commonnames["GLUGLNSYN-PWY"] ="glutamate biosynthesis IV"
commonnames["ARGASEDEG-PWY"] ="arginine degradation I (arginase pathway)"
commonnames["PWY-6308"] ="cysteine biosynthesis II (RNA-dependent)"
commonnames["HOMOCYSDEGR-PWY"] ="cysteine biosynthesis/homocysteine degradation"
commonnames["CYSTSYN-PWY"] ="cysteine biosynthesis I"
commonnames["PWY-5760"] ="&beta;-alanine biosynthesis IV"
commonnames["PWY-5155"] ="&beta;-alanine biosynthesis III"
commonnames["PWY-3981"] ="&beta;-alanine biosynthesis I"
commonnames["PWY-3941"] ="&beta;-alanine biosynthesis II"
commonnames["PWY-3982"] ="uracil degradation II (reductive)"
commonnames["GLUTDEG-PWY"] ="glutamate degradation II"
commonnames["ASPARTATESYN-PWY"] ="aspartate biosynthesis"
commonnames["ASPSYNII-PWY"] ="cyanide detoxification II"
commonnames["ASPARAGINESYN-PWY"] ="asparagine biosynthesis II"
commonnames["ASPARAGINE-BIOSYNTHESIS"] ="asparagine biosynthesis I"
commonnames["PWY490-4"] ="asparagine biosynthesis III (tRNA-dependent)"
commonnames["PWY-4983"] ="citrulline-nitric oxide cycle"
commonnames["PWY-5154"] ="arginine biosynthesis III"
commonnames["ARGSYNBSUB-PWY"] ="arginine biosynthesis II (acetyl cycle)"
commonnames["ARGININE-SYN4-PWY"] ="arginine biosynthesis IV"
commonnames["PWY0-1021"] ="alanine biosynthesis III"
commonnames["ALANINE-VALINESYN-PWY"] ="alanine biosynthesis I"
commonnames["ALANINE-SYN2-PWY"] ="alanine biosynthesis II"
commonnames["PWY-6482"] ="diphthamide biosynthesis"
commonnames["PWY-5905"] ="hypusine biosynthesis"
commonnames["PWY-6163"] ="chorismate biosynthesis from 3-dehydroquinate"
commonnames["PWY-6435"] ="4-hydroxybenzoate biosynthesis V"
commonnames["PWY-6431"] ="4-hydroxybenzoate biosynthesis IV"
commonnames["PWY-5755"] ="4-hydroxybenzoate biosynthesis II (bacteria and fungi)"
commonnames["PWY-5754"] ="4-hydroxybenzoate biosynthesis I (animals)"
commonnames["PWY-6164"] ="3-dehydroquinate biosynthesis I"
commonnames["PWY-6160"] ="3-dehydroquinate biosynthesis II (archaea)"
commonnames["PWY-981"] ="salicylate biosynthesis II"
commonnames["PWY-6930"] ="phenolic malonylglucosides biosynthesis"
commonnames["PWY-6762"] ="salicylate glucosides biosynthesis IV"
commonnames["PWY-6752"] ="<i>o</i>-diquinones biosynthesis"
commonnames["PWY-6707"] ="gallate biosynthesis"
commonnames["PWY-6624"] ="salicylate glucosides biosynthesis III"
commonnames["PWY-6623"] ="salicylate glucosides biosynthesis II"
commonnames["PWY-6539"] ="petivericin biosynthesis"
commonnames["PWY-6458"] ="benzoyl-CoA biosynthesis"
commonnames["PWY-6457"] ="<i>trans</i>-cinnamoyl-CoA biosynthesis"
commonnames["PWY-6406"] ="salicylate biosynthesis I"
commonnames["PWY-6323"] ="benzoylanthranilate biosynthesis"
commonnames["PWY-6320"] ="phaselate biosynthesis"
commonnames["PWY-5901"] ="2,3-dihydroxybenzoate biosynthesis"
commonnames["PWY-5886"] ="4-hydroxyphenylpyruvate biosynthesis"
commonnames["PWY-5787"] ="oligomeric urushiol biosynthesis"
commonnames["PWY-5765"] ="1,3,5-trimethoxybenzene biosynthesis"
commonnames["PWY-801"] ="homocysteine and cysteine interconversion"
commonnames["PWY-6972"] ="oleandomycin activation/inactivation"
commonnames["PWY-6303"] ="methyl indole-3-acetate interconversion"
commonnames["PWY-5926"] ="afrormosin conjugates interconversion"
commonnames["PWY-2904"] ="formononetin conjugates interconversion"
commonnames["PWY-2861"] ="biochanin A conjugates interconversion"
commonnames["PWY-2701"] ="maackiain conjugates interconversion"
commonnames["PWY-2561"] ="medicarpin conjugates interconversion"
commonnames["PWY-2345"] ="genistein conjugates interconversion"
commonnames["PWY-2343"] ="daidzein conjugates interconversion"
commonnames["PWY-5835"] ="geranyl acetate biosynthesis"
commonnames["PWY-6546"] ="brassinosteroids inactivation"
commonnames["PWY-6494"] ="gibberellin inactivation III (epoxidation)"
commonnames["PWY-6477"] ="gibberellin inactivation II (methylation)"
commonnames["PWY-102"] ="gibberellin inactivation I (2&beta;-hydroxylation)"
commonnames["PWY-5272"] ="abscisic acid glucose ester biosynthesis"
commonnames["PWYQT-4477"] ="indole glucosinolate breakdown (active in intact plant cell)"
commonnames["PWYQT-4476"] ="indole glucosinolate breakdown (insect chewing induced)"
commonnames["PWY-5267"] ="glucosinolate breakdown"
commonnames["PWY-5340"] ="sulfate activation for sulfonation"
commonnames["PWY-5143"] ="fatty acid activation"
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
  
    patt = re.compile(r'BD_1_502')

    for line in lines:
       line = line.strip()
       if line:
          fields = [ x.strip() for x in line.split('\t') ]

          if len(fields) >= 9:
             find = patt.search(fields[0])
             enzyme_name = fields[0]
             # this is the actual taxonomy = fields[8]
             taxonomy = fields[4]
             enzymes[enzyme_name] = taxonomy

def read_composite_pathways_file(composite_pathways_file):
    composite_pathways = []
    try:
       compositepathwaysfile = open(composite_pathways_file, 'r')
       lines = compositepathwaysfile.readlines()
       compositepathwaysfile.close()
    except:
        print "WARNING : Could not open composite pathways definition file " 
        return composite_pathways 

    for line in lines:
       line = line.strip()
       if line:
          fields = [ x.strip() for x in line.split('\t') ]
          composite_pathways.append(fields)  

    return composite_pathways



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
             if len(fields) > 5:
                 for i in  range(5,len(fields)):
                     enzyme = fields[i]
                     pathways[pathway_name]['rxns'][rxn_name][enzyme] = True 
                       

def print_taxonomy_information(pathways,ranked_pathways, enzymes, lca):
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
           for enzyme in pathways[key]['rxns'][rxn]:   
             if enzyme in enzymes:
                pathways[key]['rxns'][rxn][enzyme] = enzymes[enzyme]   
             
                #print '      '+ enzyme + '  ' + pathways[key]['rxns'][rxn][enzyme]
       # indlist = lca.get_independent_taxons(taxons)
       # for ind in indlist:
       #     print '             ' + ind
         
def create_composite_pathways(pathways, enzymes, composite_pathways):
    composite_pathway_list = []
    global commonnames
    for composite in composite_pathways:
        new_name = ''
        common_name = ''
        for pathway in composite:
           new_name = new_name+ "__" +  pathway
           common_name = common_name +  "(" + commonnames[pathway] + ")"

        composite_pathway_list.append(new_name)

        pathways[new_name] = {}
        pathways[new_name]['common_name'] = common_name
        pathways[new_name]['rxns'] = {}
        pathways[new_name]['rxns-type']={}

        incomplete = False
        for pathway in composite:
           if not pathway in pathways:
             incomplete = True
             break;
           else:
             for rxn in pathways[pathway]['rxns']:   
                pathways[new_name]['rxns'][rxn] = {}
                pathways[new_name]['rxns-type'][rxn]=pathways[pathway]['rxns-type'][rxn]
                for enzyme in pathways[pathway]['rxns'][rxn]:   
                   if enzyme in enzymes:
                     pathways[new_name]['rxns'][rxn][enzyme] = enzymes[enzyme]
        if incomplete:
           pathways[new_name]['common_name'] = common_name
           pathways[new_name]['rxns'] = {}
           pathways[new_name]['rxns-type']={}

    return composite_pathway_list 

            
   
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

    composite_pathways = read_composite_pathways_file(opts.composite_pathways_file)

    namesdict = {}
    read_taxonomy_list(opts.ncbi_file, namesdict)
    lca = SpeciesComputation(opts.ncbi_file)
#    lca.construct_tree()


    enzymes={}
    read_enzymes_list(opts.enzymes_file, enzymes)

    pathways={}
    read_pathways_list(opts.pathways_file, pathways)

    add_taxonomy_information(pathways, enzymes)

    composite_pathways_list = create_composite_pathways(pathways, enzymes, composite_pathways)
  
    #for cp in composite_pathways_list:  
    #   print cp + '  ' + pathways[cp]['common_name'] 
    
    distrib_pathways = {}
    for p in pathways:  
    #for p in composite_pathways_list:  
       value =  compute_min_species(pathways,p, lca)
       distrib_pathways[p] = value

    ranked_pwys = sorted(distrib_pathways, key=distrib_pathways.get, reverse=True)

    for p in ranked_pwys: 
      if  p in composite_pathways_list:
         print p + '\t' + str(distrib_pathways[p]) +  '\t' + pathways[p]['common_name']

    sys.exit(0)
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

