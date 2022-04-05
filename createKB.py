import kindred
import argparse
from collections import defaultdict
import sys
import json
import re
import utils
import gzip

def isASCII(s):
	#try:
	#	s.decode('ascii')
	#	return True
	#except UnicodeDecodeError:
	#	return False

	return len(s) == len(s.encode())

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Build classifiers and create the PGmine knowledge base')
	parser.add_argument('--trainingFiles',required=True,type=str,help='3 BioC files (comma-separated) for star_allele, rs and other')
	parser.add_argument('--selectedChemicals',required=True,type=str,help='Which chemicals to filter for')
	parser.add_argument('--dbsnp',required=True,type=str,help='File containing mappings from dbSNP IDs to genes')
	parser.add_argument('--genes',required=True,type=str,help='File containing gene data')
	parser.add_argument('--variantStopwords',required=True,type=str,help='Variants to remove')
	parser.add_argument('--relevantMeSH',required=True,type=str,help='File with MeSH term mappings')
	parser.add_argument('--inBioC',required=True,type=str,help='BioC file to predict things on')
	parser.add_argument('--outKB',required=True,type=str,help='TSV file for KB')
	args = parser.parse_args()

	
	chemMeshIDs = set()
	meshIDsToChemName,meshIDsToPharmGKB,meshIDsToDrugBank = {},{},{}
	cancerDrugMeshIDs = set()
	with open(args.selectedChemicals) as f:
		chemData = json.load(f)
		for chem in chemData:
			for meshID in chem['MeSH']:
				meshID = 'MESH:'+meshID
				chemMeshIDs.add(meshID)

				if chem['isCancerDrug']:
					cancerDrugMeshIDs.add(meshID)

				meshIDsToChemName[meshID] = chem['name']
				meshIDsToDrugBank[meshID] = chem['DrugBank']
				if 'PharmGKB' in chem:
					meshIDsToPharmGKB[meshID] = chem['PharmGKB']
					

	dbsnp = {}
	with open(args.dbsnp) as f:
		for line in f:
			rsid,geneInfos = line.strip('\n').split('\t')
			geneInfos = [ tuple(geneInfo.split(':')) for geneInfo in geneInfos.split('|') ]
			geneNames = [ geneName for geneName,geneID in geneInfos ]
			geneIDs = [ geneID for geneName,geneID in geneInfos ]

			dbsnp[rsid] = (geneNames,geneIDs)

	with open(args.variantStopwords) as f:
		variantStopwords = set( line.strip().lower() for line in f )

	geneID2Name = {}
	proteinCodingGenes = set()
	with open(args.genes) as f:
		for line in f:
			entrezID,name,geneType = line.strip('\n').split('\t')
			geneID2Name[entrezID] = name
			if geneType == 'protein-coding':
				proteinCodingGenes.add(entrezID)

	print("Loaded chemical, gene and variant data")

	pediatricTerms = set(['Pediatrics','Infant','Infant, Newborn','Child','Child, Preschool','Adolescent','Birth Cohort'])
	adultTerms = set(['Adult','Aged','Middle Aged','Young Adult'])

	with gzip.open(args.relevantMeSH,'rt') as f:
		relevantMeSH = json.load(f)

	pediatricPMIDs = set( int(pmid) for pmid,terms in relevantMeSH.items() if any( t in pediatricTerms for t in terms ) )
	adultPMIDs = set( int(pmid) for pmid,terms in relevantMeSH.items() if any( t in adultTerms for t in terms ) )
	print("Loaded mesh PMIDs for pediatric/adult terms")

	# Fix mapping of some popular variants to the correct SNP
	variantFixes = {
		'rs121913377':'rs113488022' # BRAF V600E
	}

	modes = ["star_rs","other"]
	trainingFiles = args.trainingFiles.split(',')
	assert len(trainingFiles) == 2, "Must provide 2 files (comma-separated) for star_rs and other"

	hlaGeneIDs = {"3105","3106","3107","3108","3109","3111","3112","3113","3115","3117","3118","3119","3120","3121","3122","3123","3125","3126","3127","3133","3134","3135"}

	obviousMistakes = {('Abacavir','HLA-B*15:02'),('Allopurinol','HLA-B*15:02'),('Carbamazepine','HLA-B*57:01'),('Allopurinol','HLA-B*57:01'),('Carbamazepine','HLA-B*58:01'),('Abacavir','HLA-B*58:01')}
	chemicalExclusions = {'cc and tc', 'cc+tc', 'cc + tc','whitehall ii','rvoto','lev-pae','oxaipn','luminal b','oxc-mpe','rapid stemi','vp40e'}

	headers = ['pmid','title','journal','journal_short','year','month','day','is_pediatric_paper','is_adult_paper','section','subsection','chemical_mesh_id','chemical_pharmgkb_id','chemical_drugbank_id','chemical_text','chemical_normalized','chemical_position','variant_id','variant_type','variant_text','variant_normalized','variant_position','gene_ids','gene_names','score','sentence','formatted_sentence']
	with open(args.outKB,'w') as outF:
		outF.write("\t".join(headers) + "\n")

		for mode,trainingData in zip(modes,trainingFiles):
			print("Creating classifier for %s" % mode)
			predictedCount = 0

			trainCorpus = kindred.load('biocxml',trainingData)
			corpus = kindred.load('biocxml',args.inBioC)

			for doc in trainCorpus.documents:
				for relation in doc.relations:
					relation.relationType = 'ChemicalMutation'

			for doc in corpus.documents:
				if mode == 'star_allele':
					doc.entities = [ e for e in doc.entities if not (e.entityType == 'Mutation' and not e.text.strip().startswith('*')) ]
				elif mode == 'rs':
					doc.entities = [ e for e in doc.entities if not (e.entityType == 'Mutation' and not e.text.strip().startswith('rs')) ]
				elif mode == 'star_rs':
					doc.entities = [ e for e in doc.entities if not (e.entityType == 'Mutation' and not (e.text.strip().startswith('rs') or e.text.strip().startswith('*'))) ]
				else:
					doc.entities = [ e for e in doc.entities if not (e.entityType == 'Mutation' and (e.text.strip().startswith('*') or e.text.strip().startswith('rs'))) ]

				doc.entities = [ e for e in doc.entities if e.position[0][0] >= 0 ]
				doc.entities = [ e for e in doc.entities if e.position[0][1] <= len(doc.text) ]
				doc.entities = [ e for e in doc.entities if not (e.entityType == 'Chemical' and not e.metadata['conceptid'] in chemMeshIDs) ]
				doc.entities = [ e for e in doc.entities if not (e.entityType == 'Chemical' and len(e.text) <= 4) ]
				doc.entities = [ e for e in doc.entities if not (e.entityType == 'Mutation' and utils.normalizeMutation(e.text) is None) ]
				doc.entities = [ e for e in doc.entities if not (e.entityType == 'Mutation' and e.text.lower() in variantStopwords) ]

				doc.entities = [ e for e in doc.entities if not (e.entityType == 'Chemical' and re.match("^[A-Z][\[\]\-\(\)\d]*[A-Z]$",e.text) ) ]
				doc.entities = [ e for e in doc.entities if not (e.entityType == 'Chemical' and e.text.lower() in chemicalExclusions) ]


			parser = kindred.Parser(model='en_core_sci_sm')
			parser.parse(trainCorpus)
			parser.parse(corpus)

			chemicalVariantClassifier = kindred.RelationClassifier(classifierType='LogisticRegression',threshold=0.5,acceptedEntityTypes=[('Chemical','Mutation')])
			chemicalVariantClassifier.train(trainCorpus)
			chemicalVariantClassifier.predict(corpus)

			for doc in corpus.documents:

				pmid = doc.metadata['pmid']
				title = doc.metadata['title']
				journal = doc.metadata['journal']
				year = doc.metadata['year']
				month = doc.metadata['month']
				day = doc.metadata['day']
				section = doc.metadata['section']

				is_pediatric_paper = pmid and int(pmid) in pediatricPMIDs
				is_adult_paper = pmid and int(pmid) in adultPMIDs

				journal_short = journal
				if len(journal_short) > 50:
					journal_short = journal_short[:50] + '...'

				if 'subsection' in doc.metadata:
					subsection = doc.metadata['subsection'].replace('’',"'")
				elif doc.metadata['section'] == 'abstract':
					subsection = 'abstract'
				elif doc.metadata['section'] == 'title':
					subsection = 'title'

				if subsection == 'None':
					subsection = 'unknown'

				groups = defaultdict(set)
				scores = defaultdict(lambda : 0)

				# We want to group pairs of chemical/variants together so if we don't create redundant relations explaining the same relation where there are multiples of the same chemical/variant across the sentence
				chemicalVariantRelations = [ r for r in doc.relations if r.relationType == 'ChemicalMutation' ]
				for chemicalVariant in chemicalVariantRelations:
					chemical,variant = chemicalVariant.entities

					chemical_mesh_id = chemical.metadata['conceptid']
					variant_concept_id = variant.metadata['conceptid']

					if ';' in chemical_mesh_id:
						continue
					elif ';' in variant_concept_id:
						continue

					key = (chemical_mesh_id,variant_concept_id)
					groups[key].add(chemical)
					groups[key].add(variant)
					scores[key] = max(scores[key],chemicalVariant.probability)


				for key,chemicalVariants in groups.items():

					score = scores[key]

					# Sort by location in sentence
					chemicalVariants = sorted(chemicalVariants, key = lambda x: x.position[0] )

					chemicals = [ e for e in chemicalVariants if e.entityType == 'Chemical' ]
					variants = [ e for e in chemicalVariants if e.entityType == 'Mutation' ]

					chemical,variant = chemicals[0],variants[0]

					chemical_text = chemical.text
					chemical_mesh_id = chemical.metadata['conceptid']
					chemical_pharmgkb_id = meshIDsToPharmGKB[chemical_mesh_id] if chemical_mesh_id in meshIDsToPharmGKB else 'NA'
					chemical_normalized = meshIDsToChemName[chemical_mesh_id]
					chemical_drugbank_id = meshIDsToDrugBank[chemical_mesh_id]

					# Remap statins
					chemical_text_lower = chemical_text.lower()
					if chemical_text_lower in ['statin','statins']:
						chemical_pharmgkb_id = 'PA133950441'
						chemical_normalized = 'HMG-CoA reductase inhibitors'
						chemical_drugbank_id = ''
					elif chemical_text_lower == 'citalopram':
						chemical_pharmgkb_id = 'PA449015'
						chemical_normalized = 'Citalopram'
						chemical_drugbank_id = 'DB00215'
					elif chemical_text_lower == 'levomilnacipran':
						chemical_pharmgkb_id = 'PA166182150'
						chemical_normalized = 'Levomilnacipran'
						chemical_drugbank_id = 'DB08918'

					variant_text = variant.text
					variant_normalized = utils.normalizeMutation(variant_text)
					if variant_normalized is None:
						continue

					variant_metadata = variant.metadata['conceptid'].split(';')
					corresponding_rsids = [ x for x in variant_metadata if re.match(r'rs\d+',x) ]
					corresponding_genes = [ x for x in variant_metadata if re.match(r'CorrespondingGene:(?P<id>\d+)',x) ]

					variant_id = ''
					genes,gene_names,gene_ids = [],'',''
					if len(corresponding_rsids) == 1:
						variant_id = corresponding_rsids[0]
						if variant_id in dbsnp:
							gene_names,gene_ids = dbsnp[variant_id]

							proteinCoding = [ (gene_id,gene_name) for gene_id,gene_name in zip(gene_ids,gene_names) if gene_id in proteinCodingGenes ]
							if len(proteinCoding) > 0:
								# Only include the protein coding if there are any
								gene_ids = [ gene_id for gene_id,gene_name in proteinCoding ]
								gene_names = [ gene_name for gene_id,gene_name in proteinCoding ]

							genes = [ e for e in doc.entities if e.entityType == 'Gene' and e.metadata['conceptid'] in gene_ids ]

							gene_names = ",".join(gene_names)
							gene_ids = ",".join(gene_ids)

					if len(corresponding_genes) == 1:
						tmp_gene_id = corresponding_genes[0].split(':')[1]
						if tmp_gene_id in geneID2Name:
							gene_names = geneID2Name[tmp_gene_id]
							gene_ids = tmp_gene_id

					if variant_id in variantFixes:
						variant_id = variantFixes[variant_id]

					chemical_position = ";".join( "%s,%s" % c.position[0] for c in chemicals )
					variant_position = ";".join( "%s,%s" % v.position[0] for v in variants )

					if variant_text.startswith('rs') and variant_text == variant_id:
						variant_normalized = variant_id

					# Skip variants that start with asterisks but don't have metadata for a star allele - likely a mistake
					if variant_text.strip().startswith('*') and not 'associated_gene' in variant.metadata:
						continue

					variant_type = 'unclear'
					if variant_normalized.startswith('*'):
						variant_type = 'star_allele'
					elif variant_normalized.startswith('p.'):
						variant_type = 'protein'
					elif variant_normalized.startswith('c.') or variant_normalized.startswith('g.') or variant_normalized.startswith('m.'):
						variant_type = 'dna'
					elif variant_normalized.startswith('rs'):
						variant_type = 'rs_id'



					if variant_type == 'star_allele':
						variant_normalized = variant.metadata['conceptid']

						associated_gene = variant.metadata['associated_gene']
						gene_ids,gene_names = None,None
						gene_ids = [ gene_id for gene_id in associated_gene.split(';') if gene_id in geneID2Name ]
						gene_names = [ geneID2Name[gene_id] for gene_id in gene_ids ]

						if len(gene_ids) != 1:
							continue

						gene_ids = gene_ids[0]
						gene_names = gene_names[0]

						genes = [ e for e in doc.entities if e.entityType == 'Gene' and e.metadata['conceptid'] == gene_ids ]

						isHLAGene = gene_ids in hlaGeneIDs
						if isHLAGene:
							variant_normalized = variant_normalized[1:].lstrip('0').replace(':','')
							if len(variant_normalized) % 2 == 1: # Pad the variant name with a zero to make it an even length
								variant_normalized = "0" + variant_normalized

							variant_normalized = re.sub("(\d)(?=(\d{2})+(?!\d))", r"\1:", variant_normalized) # Put in ':' every two digits
							variant_normalized = '*' + variant_normalized
							
						variant_id = gene_names + variant_normalized

					# Skip cancer drugs that are associated with a DNA/protein variant as likely somatic
					if chemical_mesh_id in cancerDrugMeshIDs and variant_type in ['dna','protein']:
						continue

					# Remove some very frequent mismatches
					if (chemical_normalized,variant_id) in obviousMistakes:
						continue

					sentence = doc.text.replace('’',"'")
					formatted_sentence = utils.getFormattedDoc(doc, chemicals + variants + genes)

					outData = [ pmid, title, journal, journal_short, year, month, day, is_pediatric_paper, is_adult_paper, section, subsection, chemical_mesh_id, chemical_pharmgkb_id, chemical_drugbank_id, chemical_text, chemical_normalized, chemical_position, variant_id, variant_type, variant_text, variant_normalized, variant_position, gene_ids, gene_names, score, sentence, formatted_sentence ]

					allowedUnicode = {'title','journal','journal_short','chemical_text','variant_text','sentence','formatted_sentence'}

					assert len(outData) == len(headers)
					for h,v in zip(headers,outData):
						if not (h in allowedUnicode or isASCII(str(v))):
							print('WARNING: Found non-ASCII "%s" in column "%s" (PMID=%s)' % (str(v),h,pmid))



					outF.write("\t".join(map(str,outData)) + "\n")
					predictedCount += 1

			print("Predicted %d association(s) for %s variants" % (predictedCount, mode))
			

