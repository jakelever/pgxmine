import argparse
import bioc


if __name__ == '__main__':
	parser = argparse.ArgumentParser('Check if a document has pediatric MeSH terms')
	parser.add_argument('--inBioc',required=True,type=str,help='Input BioC file')
	parser.add_argument('--outTxt',required=True,type=str,help='Text file with PMIDs for pediatric documents')
	args = parser.parse_args()

	pmids = set()

	print("Loaded PMIDs from corpus file...")

	pediatricTerms = set(['Pediatrics','Infant','Infant, Newborn','Child','Child, Preschool','Adolescent'])
	pediatricPMIDs = []

	#with bioc.BioCXMLDocumentReader(args.inBioc) as parser:
	with open(args.inBioc,'rb') as f:
		parser = bioc.BioCXMLDocumentReader(f)
		for i,doc in enumerate(parser):
			if not ('pmid' in doc.infons and doc.infons['pmid'] and doc.infons['pmid'] != 'None'):
				continue
			if not ('meshHeadings' in doc.infons and doc.infons['meshHeadings'] and doc.infons['meshHeadings'] != 'None'):
				continue

			pmid = int(doc.infons['pmid'])

			meshHeadings = [ heading.split('~') for heading in doc.infons['meshHeadings'].split('\t') ]

			for descriptorQualifiers in meshHeadings:
				descriptor = descriptorQualifiers[0].split('|')
				assert len(descriptor) == 4, "Expected four pipe-delimited columns. Got: %s" % descriptorQualifiers[0]
				_, meshID, isMajorYN, name = descriptor
				if name in pediatricTerms:
					pediatricPMIDs.append(pmid)

				#print(descriptor)

	pediatricPMIDs = sorted(set(pediatricPMIDs))
	print("Found %d PubMed ID(s) with pediatric MeSH terms" % len(pediatricPMIDs))

	with open(args.outTxt,'w') as f:
		for pmid in pediatricPMIDs:
			f.write("%d\n" % pmid)

