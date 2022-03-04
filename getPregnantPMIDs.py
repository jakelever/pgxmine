import argparse
import bioc


if __name__ == '__main__':
	parser = argparse.ArgumentParser('Check if a document has pregnant MeSH terms')
	parser.add_argument('--inBioc',required=True,type=str,help='Input BioC file')
	parser.add_argument('--outTxt',required=True,type=str,help='Text file with PMIDs for pregnant documents')
	args = parser.parse_args()

	pmids = set()

	print("Loaded PMIDs from corpus file...")

	#pregnantTerms = set(['Pregnant Women'])
	pregnantTerms = set(['Pregnant Women','Pregnancy'])
	pregnantPMIDs = []

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
				if name in pregnantTerms:
					pregnantPMIDs.append(pmid)

				#print(descriptor)

	pregnantPMIDs = sorted(set(pregnantPMIDs))
	print("Found %d PubMed ID(s) with pregnant MeSH terms" % len(pregnantPMIDs))

	with open(args.outTxt,'w') as f:
		for pmid in pregnantPMIDs:
			f.write("%d\n" % pmid)

