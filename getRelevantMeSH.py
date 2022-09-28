import argparse
import bioc
import json
import gzip


if __name__ == '__main__':
	parser = argparse.ArgumentParser('Check if a document has relevant MeSH tags')
	parser.add_argument('--inBioc',required=True,type=str,help='Input BioC file')
	parser.add_argument('--outJSONGZ',required=True,type=str,help='JSON file with PMIDs mapped to relevant MeSH terms')
	args = parser.parse_args()

	pmids = set()

	print("Loaded PMIDs from corpus file...")

	relevantTerms = ['Pediatrics','Infant','Infant, Newborn','Child','Child, Preschool','Adolescent','Birth Cohort']
	relevantTerms += ['Adult','Aged','Middle Aged','Young Adult']
	relevantTerms = set(relevantTerms)

	print("Searching for MeSH terms in: ", sorted(relevantTerms))
	print()

	seenPMIDs = set()	
	pmidToMesh = {}
	with bioc.biocxml.iterparse(open(args.inBioc, 'rb')) as parser:
		for i,doc in enumerate(parser):
			if not ('pmid' in doc.infons and doc.infons['pmid'] and doc.infons['pmid'] != 'None'):
				continue
			if not ('meshHeadings' in doc.infons and doc.infons['meshHeadings'] and doc.infons['meshHeadings'] != 'None'):
				continue

			pmid = int(doc.infons['pmid'])
			if pmid in seenPMIDs:
				continue
			seenPMIDs.add(pmid)

			meshHeadings = [ heading.split('~') for heading in doc.infons['meshHeadings'].split('\t') ]

			descriptorNames = []
			for descriptorQualifiers in meshHeadings:
				descriptor = descriptorQualifiers[0].split('|')
				assert len(descriptor) == 4, "Expected four pipe-delimited columns. Got: %s" % descriptorQualifiers[0]
				_, meshID, isMajorYN, name = descriptor
				descriptorNames.append(name)

			descriptorNames = [ x for x in descriptorNames if x in relevantTerms ]
			if descriptorNames:
				pmidToMesh[pmid] = descriptorNames

	print("Found %d PubMed ID(s) with relevant MeSH terms" % len(pmidToMesh))

	with gzip.open(args.outJSONGZ,'wt') as f:
		json.dump(pmidToMesh,f)

