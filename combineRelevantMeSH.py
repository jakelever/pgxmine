import argparse
import json
from tqdm import tqdm
import os
import gzip


if __name__ == '__main__':
	parser = argparse.ArgumentParser('Combine all the MeSH data into a single file')
	parser.add_argument('--inDir',required=True,type=str,help='Input directory of MeSH terms')
	parser.add_argument('--outJSONGZ',required=True,type=str,help='JSON GZ file with PMIDs mapped to relevant MeSH terms')
	args = parser.parse_args()

	pmids = set()

	print("Loaded PMIDs from corpus file...")

	filenames = sorted( f for f in os.listdir(args.inDir) if f.endswith('.json.gz') )

	pmidToMesh = {}
	for filename in tqdm(filenames):
		with gzip.open(os.path.join(args.inDir,filename),'rt') as f:
			tmp = json.load(f)
			pmidToMesh.update(tmp)

	print("Combined %d PubMed ID(s) with relevant MeSH terms" % len(pmidToMesh))

	with gzip.open(args.outJSONGZ,'wt') as f:
		json.dump(pmidToMesh,f)

