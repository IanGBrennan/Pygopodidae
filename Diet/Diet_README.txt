When adding new species to the diet Squamata_Diets_Families.csv
	+ you must also add them to the dietname-to-treename.csv file!
	+ you also must make sure you've spelled the name right and that the taxon is in the tree!
	
Identifying new taxa:
	+ new file: /Users/ianbrennan/Desktop/Pygopodidae_Diet/Squamata_Diets_Families.csv
	+ old file: /Users/ianbrennan/Documents/Science/Title_supplementary_Snake_Singularity/data/3.trait_preparation_diversification/dietModeling/diet.csv
	+ just setdiff the unique predator names
	+ added since 19 August 2025:
		- Amphisbaena_alba
		- Amphisbaena_pretrei
		- Amphisbaena_ibijara
		- Blanus_cinereus
		- Chamaesaura_anguina
		- Platysaurus_intermedius
		- Delma_inornata
		- Delma_fraseri
		- Delma_nasuta
		- Lialis_burtonis
		- Lialis_jicari
		- Pygopus_lepidopodus
		- Pygopus_nigriceps
		- Aprasia_haroldi
		- Aprasia_inaurita
		- Aprasia_parapulchella
		- Aprasia_pulchella
		- Aprasia_repens
		- Aprasia_smithi
		- Aprasia_striolata
		- Paradelma_orientalis
		- Ophidiocephalus_taeniatus
		- Pletholax_gracilis
		- Lerista_desertorum
		- Varanus_mertensi
		- Varanus_panoptes
		- Varanus_mitchelli
		- Varanus_komodoensis
		- Varanus_acanthurus
		- Varanus_glebopalma
		- Varanus_griseus
		- Varanus_rudicollis
		- Varanus_prasinus
		- Varanus_giganteus
		
		TEST:
			Varanus_komodoensis
			Varanus_acanthurus
			Varanus_glebopalma
			Varanus_griseus
			Varanus_niloticus
			Varanus_rudicollis
			Varanus_prasinus
			Varanus_giganteus
		
Run the MCMC chain on the diet data
	+ 1_mcmc-diet.R
	+ before you do this, make sure to check the names of taxa in the Squamata_Diets_Families.csv
		against the tree best_ultrametric_fulltree_ddBD_revision.tre
	+ 
		
Make the diet-proportions-phylo.csv:
	+ 2.make-diet-proportions.R to read in the mcmc files

		