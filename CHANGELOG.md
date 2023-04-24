# CHANGELOG

### Release date 
	* version 2.0.0 on September 30 2022

#### Features
	*	DSL2 based implemetation of tumor WGS
	*	Able to adabpt to both Hematology and Solid WGS questions
	* 	GENS for both normal and tumor
	* 	Able to handle both tumor/normal and tumor only analysis
	* 	COYOTE can handle unfiltered SNV list but tagged SNV from gene panel can be selected
	
	
#### Fixes
	* 	Sharding from 32 to 8 for better I/O specially with bqsr steps in senteion
	* 	Changes with the profile configs
	* 	Changes with the coyote to include sample-wgs
