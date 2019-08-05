import argparse
import bioc
import pickle
from collections import defaultdict,Counter
import re
import sys
import string

import datetime
def now():
	return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def createRegex(mention):
	mention = re.sub('\s+',' ',mention.strip())
	m = re.escape(mention)
	m = m.replace('alpha','(alpha|α)')
	m = m.replace('beta','(beta|β)')
	m = m.replace('gamma','(gamma|γ)')
	m = m.replace('delta','(delta|δ)')
	m = m.replace('\ ','\s+')
	return '\\b%s\\b' % m

if __name__ == '__main__':
	parser = argparse.ArgumentParser('Text align PubTator annotations against a BioC file')
	parser.add_argument('--inBioc',required=True,type=str,help='Input BioC file')
	parser.add_argument('--annotations',required=True,type=str,help='Pre-pickled annotations')
	parser.add_argument('--outBioc',required=True,type=str,help='Output BioC file')
	args = parser.parse_args()

	pmids = set()

	with bioc.iterparse(args.inBioc) as parser:
		for i,doc in enumerate(parser):
			if 'pmid' in doc.infons and doc.infons['pmid'] != 'None':
				pmid = int(doc.infons['pmid'])
				pmids.add(pmid)
	
	pmidToAnnotations = defaultdict(list)
	with open(args.annotations) as f:
		for line in f:
			split = line.strip('\n').split('\t')
			pmid,annotationType,conceptid,mentions,database = split
			mentions = mentions.strip()
			pmid = int(pmid)
			if len(mentions) > 0 and pmid in pmids:
				pmidToAnnotations[pmid].append((annotationType,conceptid,mentions))




	currentID = 1
	writer = bioc.iterwrite(args.outBioc)
	with bioc.iterparse(args.inBioc) as parser:
		for i,doc in enumerate(parser):
			for passage in doc.passages:
				passage.annotations = []

			if 'pmid' in doc.infons and doc.infons['pmid'] != 'None':
				pmid = int(doc.infons['pmid'])

				print(now(),i,pmid)
				sys.stdout.flush()

				for passage in doc.passages:
					candidates = defaultdict(lambda : defaultdict(list))

					for annotationType,conceptid,mentions in pmidToAnnotations[pmid]:
						mentions = mentions.split('|')
						regexs = [ createRegex(m) for m in mentions ]
						for mention,regex in zip(mentions,regexs):
							for match in re.finditer(regex, passage.text):
								start,end = match.span()

								candidates[(start,end)][annotationType].append(conceptid)

					locations = sorted([ (end-start,start,end) for start,end in candidates.keys() ], reverse=True)
					nonoverlapping = []
					for _,start1,end1 in locations:
						overlapping = False
						for start2,end2 in nonoverlapping:
							if start1 > end2:
								pass
							elif start2 > end1:
								pass
							else:
								overlapping = True
								break

						if not overlapping:
							nonoverlapping.append ((start1,end1))

					for start,end in nonoverlapping:
						for annotationType,conceptids in candidates[(start,end)].items():
							conceptid = conceptids = ";".join(sorted(list(set(conceptids))))

							a = bioc.BioCAnnotation()
							a.text = passage.text[start:end]
							a.infons = {'type':annotationType, 'conceptid': conceptid}
							a.id = 'T%d' % currentID
							currentID += 1
							a.locations.append(bioc.BioCLocation(offset=start, length=(end-start)))
							passage.annotations.append(a)

			writer.writedocument(doc)

	print ('Done!')

