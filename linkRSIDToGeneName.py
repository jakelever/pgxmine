import argparse
import re

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Extract mappings from rsID to genes')
	parser.add_argument('--dbsnp',required=True,type=str,help='dbSNP VCF file')
	parser.add_argument('--pubtator',required=True,type=str,help='PubTator file to know which RS ids to check')
	parser.add_argument('--outFile',required=True,type=str,help='Output file')
	args = parser.parse_args()

	#regex = re.compile("GENEINFO=[^;]*")

	wanted = []
	with open(args.pubtator) as inF:
		for line in inF:
			split = line.strip().split('\t')
			if split[1] == 'Mutation' and split[2].startswith('rs'):
				wanted.append(split[2])
	wanted = set(wanted)

	with open(args.dbsnp) as inF, open(args.outFile,'w') as outF:
		for line in inF:
			if line[0] == '#':
				continue

			split = line.strip('\n').split('\t')
			rsid = split[2]
			if not rsid in wanted:
				continue

			details = split[7]

			#match = re.search(regex,details)
			#if match:
			#	geneInfo = match.group().split('=')[1]
			#	outF.write("%s\t%s\n" % (rsid, geneInfo))

			geneInfos = [ d for d in details.split(';') if d.startswith('GENEINFO') ]
			if geneInfos:
				geneInfo = geneInfos[0].split('=')[1]
				outF.write("%s\t%s\n" % (rsid, geneInfo))


