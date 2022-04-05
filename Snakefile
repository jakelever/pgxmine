
localrules: final_files

assert os.getenv('MODE') in ['full','test'], "Must set environmental variable MODE to full or test"

if os.getenv('MODE') == 'full':
	source_dir = os.getenv('BIOTEXT')
	assert source_dir and os.path.isdir(source_dir), "For full run, must set environmental variable BIOTEXT to directory with BIOTEXT BioC XML files"
	source_dir = source_dir.rstrip('/')
	work_dir = 'working'
elif os.getenv('MODE') == 'test':
	source_dir = 'test_data'
	work_dir = 'test_working'

mesh_files = [ '%s/mesh/%s' % (work_dir,f.replace('.bioc.xml','.json.gz')) for f in os.listdir(source_dir) if f.startswith('pubmed') ]
kb_files = [ '%s/kb/%s' % (work_dir,f.replace('.bioc.xml','.tsv')) for f in os.listdir(source_dir) ]

final_files = [ "%s/%s" % (work_dir,f) for f in ['pgxmine_unfiltered.tsv','pgxmine_collated.tsv','pgxmine_sentences.tsv'] ]

rule final_files:
	input: final_files

rule getRelevantMeSH:
	input: '%s/{f}.bioc.xml' % source_dir
	output: '%s/mesh/{f}.json.gz' % work_dir
	shell: "python getRelevantMeSH.py --inBioc {input} --outJSONGZ {output}"

rule combineRelevantMeSH:
	input: 
		files = mesh_files
	output: '%s/relevantMeSH.json.gz' % work_dir
	shell: "python combineRelevantMeSH.py --inDir %s/mesh --outJSONGZ {output}" % work_dir

rule findPGxSentences:
#	input: '%s/aligned/{f}.bioc.xml' % work_dir
	input: '%s/{f}.bioc.xml' % source_dir
	output: '%s/sentences/{f}.bioc.xml' % work_dir
	shell: "python findPGxSentences.py --inBioc {input} --filterTermsFile pgx_filter_terms.txt --outBioc {output}"

rule createKB:
	input: 
		file = '%s/sentences/{f}.bioc.xml' % work_dir,
		relevant_mesh = '%s/relevantMeSH.json.gz' % work_dir
	output: '%s/kb/{f}.tsv' % work_dir,
	shell: "python createKB.py --trainingFiles data/annotations.variant_star_rs.bioc.xml,data/annotations.variant_other.bioc.xml --inBioC {input.file} --selectedChemicals data/selected_chemicals.json --dbsnp data/dbsnp_selected.tsv --variantStopwords stopword_variants.txt --genes data/gene_names.tsv --relevantMeSH {input.relevant_mesh} --outKB {output}"
   
rule filterAndCollate:
	input: kb_files,
	output: 
		unfiltered = '%s/pgxmine_unfiltered.tsv' % work_dir,
		collated = '%s/pgxmine_collated.tsv' % work_dir,
		sentences = '%s/pgxmine_sentences.tsv' % work_dir
	shell: "python filterAndCollate.py --inData %s/kb --outUnfiltered {output.unfiltered} --outCollated {output.collated} --outSentences {output.sentences}" % work_dir

