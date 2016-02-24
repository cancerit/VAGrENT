#!/usr/bin/env python
import os
import sys
import optparse
import warnings
import traceback
import re

"""
Get SV annotation map for the genome
Make sure input file is sorted in chr coordinate order
"""
#chr1 start1 end1 chr2 start2 end2 name quality strand1 strand2 filter info : total 12 fields

BEDPE_LINE='{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'
ID_FORMAT='{}:{}:{}:{}'



def createAnnotationMap(out_dir, annotation_bed, index_file,annotation_class):
	print "Dir:%s,Bed:%s,index:%s,ann_class:%s" %(out_dir,annotation_bed,index_file,annotation_class)
	#get chr length
	dict_chrLen = _getChrLen(index_file)
	fh_bed = open(annotation_bed, 'r')
	# create another file handler to seek data
	fh_bed_reset = open(annotation_bed, 'r')
	fh_out = open('annotation_map.bed', 'w')
	fh_out.write(BEDPE_LINE.format('chrL','startL','endL','chrH','startH','endH','id','score','strandL','strandH','Type','featureL','featureH'))
	svclass_counter=0
	while 1:
		line=fh_bed.readline()
		if not line:
			break
		fields = line.split("\t")
		#21	10156200	10906200	-g1	.	-	ID=ENSG00000166157.12;gene_id=ENSG00000166157.12;transcript_id=ENSG00000166157.12;gene_type=protein_coding;gene_status=KNOWN;gene_name=TPTE;transcript_type=protein_coding;transcript_status=KNOWN;transcript_name=TPTE;level=2;havana_gene=OTTHUMG00000074127.5
		#get fields to variables
		chr=fields[0]; start=fields[1]; end=fields[2]; strand=fields[5]
		#seek file handler to current bed file location
		fh_bed_reset.seek(fh_bed.tell())
		#match only with genes
		if re.match('^g',fields[3]):
			svclass_counter+=1
			gene=re.split(';|=',fields[6])[1]
			(chrL,startL,endL,chrH,startH,endH)=_getCopyLoss(chr,start,end,dict_chrLen[fields[0]])
			fh_out.write(BEDPE_LINE.format(chrL,startL,endL,chrH,startH,endH,ID_FORMAT.format(chrL,'CL',svclass_counter,chrH),'.',strand,strand,'CL',gene,gene))
			(chrL,startL,endL,chrH,startH,endH)=_getInternalDel(chr,start,end)
			fh_out.write(BEDPE_LINE.format(chrL,startL,endL,chrH,startH,endH,ID_FORMAT.format(chrL,'ID',svclass_counter,chrH),'.',strand,strand,'ID',gene,gene))
			(chrL,startL,endL,chrH,startH,endH)=_get5pOverhang(chr,start,end,strand,dict_chrLen[fields[0]])
			fh_out.write(BEDPE_LINE.format(chrL,startL,endL,chrH,startH,endH,ID_FORMAT.format(chrL,'OH',svclass_counter,chrH),'.',strand,strand,'OH',gene,gene))
			svclass_counter=_getGeneFusion(chr,start,end,strand,gene,svclass_counter,fh_bed_reset,fh_out)
			# set file handler back to start
			fh_bed_reset.seek(0)
		#enhancer apposition : pair downstream coordinates(+e) of enhancer with +ve strand genes upstream coordinates(-g) 
		if re.match('^\+e',fields[3]):
			fh_bed_reset.seek(fh_bed.tell())
			svclass_counter=_getEnhancerApposiotionToGenes(chr,start,end,strand,'.',svclass_counter,fh_bed_reset,fh_out)
			fh_bed_reset.seek(0)
		#pair  downstream cooridnates(+g) of -ve strand gene with enhancers upstream coordinates(-e) 
		if (re.match('^\+g',fields[3]) and strand == '-'):
			gene=re.split(';|=',fields[6])[1]
			fh_bed_reset.seek(fh_bed.tell())
			svclass_counter=_getGeneApposiotionToEnhancers(chr,start,end,strand,gene,svclass_counter,fh_bed_reset,fh_out)
			fh_bed_reset.seek(0)
		# set file handler back to start				
		fh_bed_reset.seek(0)	
	print "Dict %s" %(dict_chrLen['10'])
	fh_bed.close()
	fh_bed_reset.close()
	return 0

#get chromosome length
def _getChrLen(index_file):
	dict_chrLen = {}
	fh_input = open(index_file, 'r')
	for line in fh_input.readlines():
		fields = line.split("\t")
		#fields = [f.strip() for f in fields]
		dict_chrLen[fields[0]] = fields[1]
	fh_input.close()
	return dict_chrLen

def _getCopyLoss(chr,start,strand,chrLen):
	if strand == '-':
		return (chr,'0',end,chr,(int(end)-1),chrLen)
	else:
		return (chr,'0',(int(start)+1),chr,start,chrLen)
	
def _getInternalDel(chr,start,end):
	return (chr,start,end,chr,start,end)
	
def _get5pOverhang(chr,start,end,strand,chrLen):
	if strand == '-':
		return (chr,'0',start,chr,start,end)
	else:
		return (chr,start,end,chr,end,chrLen)			
		
def _getGeneFusion(chrL,startL,endL,strandL,geneL,svclass_counter,fh_bed_reset,fh_out):
	for line in fh_bed_reset:
		fields = line.split("\t")
		chrH=fields[0];
		if(re.match('^g',fields[3]) and chrL == chrH):
			startH=fields[1]; endH=fields[2];strandH=fields[5]; geneH=re.split(';|=',fields[6])[1]
			if strandL == strandH:
				svclass_counter+=1
				fh_out.write(BEDPE_LINE.format(chrL,startL,endL,chrH,startH,endH,ID_FORMAT.format(chrL,'GF',svclass_counter,chrH),'.',strandL,strandH,'GF',geneL,geneH))
	return svclass_counter
	
def _getEnhancerApposiotionToGenes(chrL,startL,endL,strandL,geneL,svclass_counter,fh_bed_reset,fh_out):
	for line in fh_bed_reset:
		fields = line.split("\t")
		chrH=fields[0];startH=fields[1]; endH=fields[2];strandH=fields[5]; 
		if(re.match('^\-g',fields[3]) and chrL == chrH and strandH == '+' ):
			geneH=re.split(';|=',fields[6])[1]
			svclass_counter+=1
			fh_out.write(BEDPE_LINE.format(chrL,startL,endL,chrH,startH,endH,ID_FORMAT.format(chrL,'EA',svclass_counter,chrH),'.','.',strandH,'EA',geneL,geneH))
	return svclass_counter

def _getGeneApposiotionToEnhancers(chrL,startL,endL,strandL,geneL,svclass_counter,fh_bed_reset,fh_out):
	for line in fh_bed_reset:
		fields = line.split("\t")
		chrH=fields[0];startH=fields[1]; endH=fields[2]
		if(re.match('^\-e',fields[3]) and chrL == chrH):
			svclass_counter+=1
			fh_out.write(BEDPE_LINE.format(chrL,startL,endL,chrH,startH,endH,ID_FORMAT.format(chrL,'EA',svclass_counter,chrH),'.',strandL,'.','EA',geneL,'.'))
	return svclass_counter
	
#####################################################################################################################
## -- FUNCTION MAIN: 
## -- To produce SV annotation map for a genome of interest
#####################################################################################################################
def main():
   
   optParser = optparse.OptionParser( 
      
      usage = "%prog [options] -o outdir -b input_annotation -i genome_index_file -c annotation_class ",
      
      description=
         "This script takes as input an output directory (outdir), " +
         "input annotation bed file (input_annotation), " + 
				 "genome index file (genome_index_file), " +
				 "annotation class (annotation_class), " +
         "set of coordinates were reported for give annotation class")

   optParser.add_option( "-o", "--outputDirectory", type="string", dest="outdir",
      default = "", help = "Output directory " )

   optParser.add_option( "-b", "--annotationBed", type="string", dest="input_annotation",
      default = "", help = "Input bed file containing various annotations types and associated gene(s) " )

   optParser.add_option( "-i", "--indexFile", type="string", dest="genome_index_file",
      default = "", help = "genome index file containing chr and length " )

   optParser.add_option( "-c", "--annotationClass", type="string", dest="annotation_class",
      default = "", help = "Annotation calss [del,td,ins]" )
   optParser.add_option( "-v", "--verbose", action="store_true", dest="verbose",
      help = "suppress progress report and warnings" )

   if len( sys.argv ) == 1:
      print "pass here"
      optParser.print_help()
      sys.exit( 1 )

   (opts, args) = optParser.parse_args()
   
   if (not opts.outdir) or (not opts.input_annotation) or (not opts.genome_index_file) or (not opts.annotation_class):
      sys.stderr.write( sys.argv[0] + ": Error - Please see the list of required parameters to provide on the command line \n" )
      sys.stderr.write( "Call with '-h' to get usage information.\n" )
      sys.exit( 1 )
      
   try:
			sys.stderr.write("Analysing data.\n")
			createAnnotationMap(opts.outdir,opts.input_annotation,opts.genome_index_file,opts.annotation_class)  
    #r2g.peakrescue_probabilistic_assignment(opts.outdir, opts.peak_filename, opts.mappings_reads2genes_filename, opts.gene_length_filename, opts.readtype)
   except:
      sys.stderr.write( "Error: %s\n" % str( sys.exc_info()[1] ) )
      sys.stderr.write( "[Exception type: %s, raised in %s:%d]\n" % 
         ( sys.exc_info()[1].__class__.__name__, 
           os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]), 
           traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
      sys.exit( 1 )


if __name__ == '__main__':
	main()

#python createAnnotationMap.py -o test_ou -b gene_enhancer.bed -i genome_h37.fa.fai -c del
