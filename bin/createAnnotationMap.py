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
svclass_counter=0
enhancer_separation_flank=500000
ehnancer_distance=1000001



annotationClass_dict={'del':'interstitial_deletion','td':'tandem_duplication','ct':'chromothripsis','ut':'unbalanced_translocation','bt':'balanced_translocation','inv':'inversion','fb':'foldback'}

EnhancerApposition='EA'
EnhancerAppositionOverlap='EO'
EnhancerSeparartion='ES'
Overhang='OH'
GeneFusion='GF'
CopyLoss='CL'
CopyGain='CG'
Complex='CX'
InternalDeletion='IDe'
InternalDeletion='IDu'

# compile regulare expressions
ug=re.compile('^\-g', flags=0) # upstream gene coordinates irrespective of strand
ag=re.compile('^g', flags=0) # actual gene coordinates
dg=re.compile('^\-g', flags=0) # downstream gene coordinates irrespective of strand
ue=re.compile('^\-e', flags=0) # upstream enhancer coordinates
ae=re.compile('^e', flags=0) # actual enhancer coordinates
de=re.compile('^\+e', flags=0)	# downstream enhancer coordinates
attr=re.compile(';|=', flags=0) # gtf attributes
allg=re.compile('g', flags=0) # any gene location

def createAnnotationMap(out_dir,annotation_bed, index_file,annotation_class,input_chr):
	print "Dir:%s,Bed:%s,index:%s,ann_class:%s" %(out_dir,annotation_bed,index_file,annotation_class)
	input_chr=re.sub(r'(?i)chr','',input_chr)
	#get chr length
	global svclass_counter;chr_tmp=0
	
	
	
	dict_chrLen = _getChrLen(index_file)
	fh_bed = open(annotation_bed, 'r')
	# create another file handler to seek data
	fh_bed_reset = open(annotation_bed, 'r')
	#write header information only for first chromosome
	fh_out = open('annotation_map_header.bed', 'w')
	fh_out.write(BEDPE_LINE.format('#chrL','startL','endL','chrH','startH','endH','id','score','strandL','strandH','Type','featureL','featureH'))
	seek_to_chr=0
		
		
	while 1:
		line=fh_bed.readline()
		if not line:
			break
		fields = line.split("\t")
		#21	10156200	10906200	-g1	.	-	ID=ENSG00000166157.12;gene_id=ENSG00000166157.12;transcript_id=ENSG00000166157.12;gene_type=protein_coding;gene_status=KNOWN;gene_name=TPTE;transcript_type=protein_coding;transcript_status=KNOWN;transcript_name=TPTE;level=2;havana_gene=OTTHUMG00000074127.5
		#get fields to variables
		chr=fields[0];
		if input_chr!=chr:
			continue
		if chr_tmp!=chr:
			#open per chromosome out file 
			fh_out.close()
			fh_out = open('{}{}{}.bed'.format(annotationClass_dict[annotation_class],'_',chr), 'w')
			#reset counter to zero for every chromosome
			chr_tmp=chr
			svclass_counter=0
		#seek file handler to current bed file location
		fh_bed_reset.seek(fh_bed.tell())
		if svclass_counter == 0:
			seek_to_chr=fh_bed.tell()
		svclass_counter+=1	
		start=fields[1]; end=fields[2]; strand=fields[5]; chr_arm=fields[7].strip();centromere=fields[8];
		#enhancer apposition : pair upstream coordinates(-e) of enhancer with downstream coordinates of gene(+g) 
		
		list_tmp_output=[]
		print (fields)
		if ue.match(fields[3]):
			if annotation_class == 'td':
				list_tmp_output=_getEnhancerApposiotionToGenesTd(chr,start,end,strand,'.',fh_bed_reset,list_tmp_output,EnhancerApposition,EnhancerAppositionOverlap)
				fh_bed_reset.seek(0)
						
		#enhancer apposition : pair downstream coordinates(+e) of enhancer with upstream coordinates of gene(-g) 
		elif de.match(fields[3]):
			if annotation_class == 'del':
				list_tmp_output=_getEnhancerApposiotionToGenesDel(chr,start,end,strand,'.',fh_bed_reset,list_tmp_output,EnhancerApposition,EnhancerAppositionOverlap)
				fh_bed_reset.seek(0)
					  			
		#match only with genes
		elif ag.match(fields[3]):
			gene=attr.split(fields[6])[1]
			#deletion...
			if annotation_class == 'del':
				list_tmp_output=_getCopyLossDel(chr,start,end,strand,gene,dict_chrLen[chr],list_tmp_output,CopyLoss)
				list_tmp_output=_getInternalDel(chr,start,end,strand,gene,list_tmp_output,InternalDeletion)
				list_tmp_output=_get5pOverhangDel(chr,start,end,strand,gene,dict_chrLen[chr],list_tmp_output,Overhang)
				list_tmp_output=_getGeneFusionDel(chr,start,end,strand,gene,fh_bed_reset,list_tmp_output,list_tmp_output,GeneFusion)
			elif annotation_class == 'td':
				list_tmp_output=_getCopyGainTd(chr,start,end,strand,gene,dict_chrLen[chr],list_tmp_output,CopyGain)
				list_tmp_output=_getInternalTd(chr,start,end,strand,gene,list_tmp_output,InternalDeletion)
				list_tmp_output=_getGeneFusionTd(chr,start,end,strand,gene,fh_bed_reset,list_tmp_output,GeneFusion)
				fh_bed_reset.seek(0)
				list_tmp_output=_get5pOverhangTd(chr,start,end,strand,gene,dict_chrLen[chr],list_tmp_output,Overhang)
				fh_bed_reset.seek(fh_bed.tell())
				list_tmp_output=_getEnhancerSeparationTd(chr,start,end,strand,gene,fh_bed_reset,list_tmp_output,EnhancerSeparartion)
				fh_bed_reset.seek(0)
			elif annotation_class == 'bt':
				list_tmp_output=_getCopyLossBt(chr,start,end,strand,gene,list_tmp_output,dict_chrLen,CopyLoss)
				list_tmp_output=_get5pOverhangBt(chr,start,end,strand,gene,list_tmp_output,dict_chrLen,Overhang)
				list_tmp_output=_getGeneFusionBt(chr,start,end,strand,gene,chr_arm,fh_bed_reset,list_tmp_output,GeneFusion)
				fh_bed_reset.seek(0)
			elif annotation_class == 'ut':
				list_tmp_output=_getCopyLossUt(chr,start,end,strand,gene,chr_arm,centromere,list_tmp_output,dict_chrLen,CopyLoss)
			#Chromothripsis...
			elif annotation_class == 'ct':
				list_tmp_output=_getChromothripsis(chr,start,end,strand,gene,list_tmp_output,dict_chrLen[chr],Complex)	
			# set file handler back to start
			fh_bed_reset.seek(0)
		#pair uptream cooridnates of gene(-g) with downstream coordinates(+e) of enhancers 
		elif ug.match(fields[3]):
			gene=attr.split(fields[6])[1]
			if annotation_class == 'td':
				list_tmp_output=_getGeneApposiotionToEnhancersTd(chr,start,end,strand,gene,fh_bed_reset,list_tmp_output,EnhancerApposition)
				fh_bed_reset.seek(0)
		#pair downstream cooridnates of gene(+g) with upstream coordinates(-e) of enhancers 
		elif dg.match(fields[3]):
			gene=attr.split(fields[6])[1]
			if annotation_class == 'del':
				list_tmp_output=_getGeneApposiotionToEnhancersDel(chr,start,end,strand,gene,fh_bed_reset,list_tmp_output,EnhancerApposition)
				fh_bed_reset.seek(0)
		#enhancer coordiantes
		elif ae.match(fields[3]):
			if annotation_class == 'del':
				list_tmp_output=_getEnhancerSeparationDel(chr,start,end,strand,'.',list_tmp_output,EnhancerSeparartion)
			elif annotation_class == 'td':
				list_tmp_output=_getEnhancerSeparationTd(chr,start,end,strand,'.',fh_bed_reset,list_tmp_output,EnhancerSeparartion)
			elif annotation_class == 'bt':
				fh_bed_reset.seek(0)	
				list_tmp_output=_getEnhancerAppositionBt(chr,start,end,strand,'.',chr_arm,fh_bed_reset,list_tmp_output,dict_chrLen,EnhancerApposition)
				fh_bed_reset.seek(0)
				if fields[6] == '.':
					continue
				#ENSG00000162576_MXRA8 ENSG00000242590_NA ENSG00000107404_DVL1 ENSG00000187608_ISG15
				# seek to the beginning of the chromosome
				fh_bed_reset.seek(seek_to_chr)
				list_tmp_output=_getEnhancerSeparationBt(chr,start,end,strand,'.',fields[6],fh_bed_reset,list_tmp_output,dict_chrLen,EnhancerSeparartion)
				fh_bed_reset.seek(0)	
		# set file handler back to start				
		fh_bed_reset.seek(0)
		fh_out.writelines(list_tmp_output)
	print "Complete annotation map for:{}".format(input_chr)
	fh_bed.close()
	fh_bed_reset.close()

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

def _getEnhancerGeneDict(enhancer_genes):
	dict_egene={}
	list_egene=re.split(' ',enhancer_genes)
	for eg in list_egene:
		egene=re.split('_',eg)[0]	
		dict_egene[egene]=egene
	return dict_egene
#Del----------
def _getCopyLossDel(chrL,startL,endL,strandL,geneL,chrLen,list_tmp_output,sv_class):
	global svclass_counter
	svclass_counter+=1
	if strandL == '-':
			list_tmp_output+=[BEDPE_LINE.format(chrL,'0',endL,chrL,int(startL)-1,chrLen,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrL),'.',strandL,strandL,sv_class,geneL,geneL)]
	else:
			list_tmp_output+=[BEDPE_LINE.format(chrL,'0',int(startL)+1,chrL,startL,chrLen,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrL),'.',strandL,strandL,sv_class,geneL,geneL)]
	return list_tmp_output

def _getInternalDel(chrL,startL,endL,strandL,geneL,list_tmp_output,sv_class):
	global svclass_counter
	svclass_counter+=1
	list_tmp_output+=[BEDPE_LINE.format(chrL,startL,endL,chrL,startL,endL,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrL),'.',strandL,strandL,sv_class,geneL,geneL)]
	return list_tmp_output
	
def _get5pOverhangDel(chrL,startL,endL,strandL,geneL,chrLen,list_tmp_output,sv_class):
	global svclass_counter
	svclass_counter+=1
	if strand == '-':
		list_tmp_output+=[BEDPE_LINE.format(chrL,'0',startL,chrL,startL,endL,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrL),'.',strandL,strandL,sv_class,geneL,geneL)]
	else:
		list_tmp_output+=[BEDPE_LINE.format(chrL,startL,endL,chrL,endL,chrLen,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrL),'.',strandL,strandL,sv_class,geneL,geneL)]
	return list_tmp_output	
		
def _getGeneFusionDel(chrL,startL,endL,strandL,geneL,fh_bed_reset,list_tmp_output,sv_class):
	global svclass_counter
	for line in fh_bed_reset:
		fields = line.split("\t")
		chrH=fields[0];
		if not(ag.match(fields[3])):
			continue
		elif(chrL == chrH):
			startH=fields[1]; endH=fields[2];strandH=fields[5]; geneH=attr.split(fields[6])[1]
			if strandL == strandH:
				svclass_counter+=1
				list_tmp_output+=[BEDPE_LINE.format(chrL,startL,endL,chrH,startH,endH,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrH),'.',strandL,strandH,sv_class,geneL,geneH)]
		else:
			return list_tmp_output
	return list_tmp_output
def _getEnhancerApposiotionToGenesDel(chrL,startL,endL,strandL,geneL,fh_bed_reset,list_tmp_output,sv_class,sv_class2):
	global svclass_counter
	flag=0
	for line in fh_bed_reset:
		fields = line.split("\t")
		chrH=fields[0];startH=fields[1]; endH=fields[2];strandH=fields[5]; 
		if not(ug.match(fields[3])):
			continue
		elif( chrL == chrH ):
			geneH=attr.split(fields[6])[1]
			svclass_counter+=1
			if flag:
				list_tmp_output+=[BEDPE_LINE.format(chrL,startL,endL,chrH,startH,endH,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrH),'.','.',strandH,sv_class,geneL,geneH)]
			elif((int(startH) - int(endL)) < 0):
				flag=0
				#enhancer Overlap
				list_tmp_output+=[BEDPE_LINE.format(chrL,startL,endH,chrH,startL,endH,ID_FORMAT.format(chrL,sv_class2,svclass_counter,chrH),'.','.',strandH,sv_class2,geneL,geneH)]
			else:
				flag=1
				list_tmp_output+=[BEDPE_LINE.format(chrL,startL,endL,chrH,startH,endH,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrH),'.','.',strandH,sv_class,geneL,geneH)]
	  	#return if done with the underlyig chromosome
		else:
			return list_tmp_output
	return list_tmp_output
	
def _getGeneApposiotionToEnhancerDel(chrL,startL,endL,strandL,geneL,fh_bed_reset,list_tmp_output,sv_class,sv_class2):
	global svclass_counter
	flag=0
	for line in fh_bed_reset:
		fields = line.split("\t")
		chrH=fields[0];startH=fields[1]; endH=fields[2]
		if not(ue.match(fields[3])):
			continue
		elif(chrL == chrH):
			svclass_counter+=1
			if flag:
				list_tmp_output+=[BEDPE_LINE.format(chrL,startL,endL,chrH,startH,endH,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrH),'.',strandL,'.',sv_class,geneL,'.')]
			elif((int(startH) - int(endL)) < 0):
				flag=0
				#enhancer Overlap
				list_tmp_output+=[BEDPE_LINE.format(chrL,startL,endH,chrH,startL,endH,ID_FORMAT.format(chrL,sv_class2,svclass_counter,chrH),'.',strandL,'.',sv_class2,geneL,'.')]
			else:
				flag=1
				list_tmp_output+=[BEDPE_LINE.format(chrL,startL,endL,chrH,startH,endH,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrH),'.',strandL,'.',sv_class,geneL,'.')]
		#return if done with the underlyig chromosome
		elif chrL!=chrH:
			return list_tmp_output
	return list_tmp_output

def _getEnhancerSeparationDel(chrL,startL,endL,strandL,geneL,list_tmp_output,sv_class):
	global svclass_counter
	svclass_counter+=1
	list_tmp_output+=[BEDPE_LINE.format(chrL,startL,endL,chrL,startL,endL,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrL),'.',strandL,strandL,sv_class,geneL,geneL)]

#Ct----------
def _getChromothripsis(chr,start,end,strand,gene,list_tmp_output,chrLen,sv_class):
	global svclass_counter
	#Lchromo
	svclass_counter+=1
	list_tmp_output+=[BEDPE_LINE.format(chr,'0',start,chr,start,end,ID_FORMAT.format(chr,sv_class,svclass_counter,chr),'.',strand,strand,sv_class,gene,gene)]
	#Mchromo
	svclass_counter+=1
	list_tmp_output+=[BEDPE_LINE.format(chr,start,end,chr,start,end,ID_FORMAT.format(chr,sv_class,svclass_counter,chr),'.',strand,strand,sv_class,gene,gene)]
	#Rchromo
	svclass_counter+=1
	list_tmp_output+=[BEDPE_LINE.format(chr,start,end,chr,end,chrLen,ID_FORMAT.format(chr,sv_class,svclass_counter,chr),'.',strand,strand,sv_class,gene,gene)]
	return list_tmp_output

#Td----------
def _getCopyGainTd(chrL,startL,endL,strandL,geneL,chrLen,list_tmp_output,sv_class):
	global svclass_counter
	svclass_counter+=1
	list_tmp_output+=[BEDPE_LINE.format(chrL,'0',startL,chrL,endL,chrLen,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrL),'.',strandL,strandL,sv_class,geneL,geneL)]
	return list_tmp_output

def _getInternalTd(chrL,startL,endL,strandL,geneL,list_tmp_output,sv_class):
	global svclass_counter
	svclass_counter+=1
	list_tmp_output+=[BEDPE_LINE.format(chrL,startL,endL,chrL,startL,endL,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrL),'.',strandL,strandL,sv_class,geneL,geneL)]
	return list_tmp_output

def _getGeneFusionTd(chrL,startL,endL,strandL,geneL,fh_bed_reset,list_tmp_output,sv_class):
	global svclass_counter
	for line in fh_bed_reset:
		fields = line.split("\t")
		chrH=fields[0];
		if not(ag.match(fields[3])):
			continue
		elif(chrL == chrH):
			startH=fields[1]; endH=fields[2];strandH=fields[5]; geneH=re.split(';|=',fields[6])[1]
			if strandL == strandH:
				svclass_counter+=1
				list_tmp_output+=[BEDPE_LINE.format(chrL,startL,endL,chrH,startH,endH,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrH),'.',strandL,strandH,sv_class,geneL,geneH)]
		else:
			return list_tmp_output
	return list_tmp_output
			
def _get5pOverhangTd(chrL,startL,endL,strandL,geneL,chrLen,list_tmp_output,sv_class):
	global svclass_counter
	svclass_counter+=1
	if strand == '-':
		list_tmp_output+=[BEDPE_LINE.format(chrL,startL,endL,chrL,endL,chrLen,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrL),'.',strandL,strandL,sv_class,geneL,geneL)]
	else:
		list_tmp_output+=[BEDPE_LINE.format(chrL,'0',startL,chrL,startL,endL,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrL),'.',strandL,strandL,sv_class,geneL,geneL)]
	return list_tmp_output

def _getEnhancerApposiotionToGenesTd(chrL,startL,endL,strandL,geneL,fh_bed_reset,list_tmp_output,sv_class):
	global svclass_counter
	for line in fh_bed_reset:
		fields = line.split("\t")
		chrH=fields[0];startH=fields[1]; endH=fields[2];strandH=fields[5];
		if not(dg.match(fields[3])):
				continue
		elif(chrL == chrH ):
			geneH=re.split(';|=',fields[6])[1]
			svclass_counter+=1
			list_tmp_output+=[BEDPE_LINE.format(chrL,startL,endL,chrH,startH,endH,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrH),'.','.',strandH,sv_class,geneL,geneH)]
	 	#return if done with the underlyig chromosome
		# valid condition as fh_bed_reset will always seek to current bed line
		else:
			return list_tmp_output
	return list_tmp_output

def _getGeneApposiotionToEnhancersTd(chrL,startL,endL,strandL,geneL,fh_bed_reset,list_tmp_output,sv_class):
	global svclass_counter
	for line in fh_bed_reset:
		fields = line.split("\t")
		chrH=fields[0];startH=fields[1]; endH=fields[2]
		if not(de.match(fields[3])):
			continue
		elif(chrL == chrH):
			svclass_counter+=1
			list_tmp_output+=[BEDPE_LINE.format(chrL,startL,endL,chrH,startH,endH,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrH),'.',strandL,'.',sv_class,geneL,'.')]
		else:
			list_tmp_output
	return list_tmp_output
	
def _getEnhancerSeparationTd(chrL,startL,endL,strandL,geneL,fh_bed_reset,list_tmp_output,sv_class):
	global svclass_counter
	counter=0
	for line in fh_bed_reset:
		flag=0;counter+=1
		fields = line.split("\t")
		chrH=fields[0];startH=fields[1];endH=fields[2];strandH=fields[5]; 
		if( ag.match(fields[3]) and chrL == chrH and (int(startH) - int(startL)) > enhancer_separation_flank and (int(startH) - int(startL)) < ehnancer_distance ):
			geneH=attr.split(fields[6])[1];svclass_counter+=1;flag=1
			list_tmp_output+=[BEDPE_LINE.format(chrL,endL,(int(startH)- enhancer_separation_flank),chrH,(int(endL)+ enhancer_separation_flank),startH,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrH),'.',strandL,strandH,sv_class,'.',geneH)]
		elif(ae.match(fields[3]) and chrL == chrH and (int(startH) - int(endL)) > enhancer_separation_flank and (int(startH) - int(endL)) < ehnancer_distance):
			svclass_counter+=1;flag=1
			list_tmp_output+=[BEDPE_LINE.format(chrL,endL,(int(startH)- enhancer_separation_flank),chrH,(int(endL)+ enhancer_separation_flank),startH,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrH),'.',strandL,strandH,sv_class,geneL,'.')]
		elif (not flag and counter > 100): # continue for few lines to capture the nearest gene if any
			return list_tmp_output
	return list_tmp_output

#Bt------------		
def _getCopyLossBt(chrH,startH,endH,strandH,geneH,list_tmp_output,dict_chrLen,sv_class):
	global svclass_counter
	for chr in dict_chrLen:
		svclass_counter+=1
		list_tmp_output+=[BEDPE_LINE.format(chr,'0',dict_chrLen[chr],chrH,startH,endH,ID_FORMAT.format(chr,sv_class,svclass_counter,chrH),'.','.',strandH,sv_class,'.',geneH)]
	return list_tmp_output
	
def _get5pOverhangBt(chrH,startH,endH,strandH,geneH,list_tmp_output,dict_chrLen,sv_class):
	global svclass_counter
	for chr in dict_chrLen:
		svclass_counter+=1
		list_tmp_output+=[BEDPE_LINE.format(chr,'0',dict_chrLen[chr],chrH,(int(startH)-1),endH,ID_FORMAT.format(chr,sv_class,svclass_counter,chrH),'.','.',strandH,sv_class,'.',geneH)]
	return list_tmp_output
 
def _getGeneFusionBt(chrL,startL,endL,strandL,geneL,chr_arm,fh_bed_reset,list_tmp_output,sv_class):
	global svclass_counter
	for line in fh_bed_reset:
		fields = line.split("\t")
		chrH=fields[0];startH=fields[1]; endH=fields[2];strandH=fields[5]; 
		if not(ag.match(fields[3])):
			continue
		elif(fields[7] == chr_arm and strandL == strandH):
			geneH=attr.split(fields[6])[1];svclass_counter+=1
			list_tmp_output+=[BEDPE_LINE.format(chrL,startL,endL,chrH,startH,endH,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrH),'.',strandL,strandH,sv_class,geneL,geneH)]
		elif(fields[7] != chr_arm and strandL != strandH):
			geneH=attr.split(fields[6])[1];svclass_counter+=1
			list_tmp_output+=[BEDPE_LINE.format(chrL,startL,endL,chrH,startH,endH,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrH),'.',strandL,strandH,sv_class,geneL,geneH)]
	return list_tmp_output
	
def _getEnhancerSeparationBt(chrL,startL,endL,strandL,geneL,egenes,fh_bed_reset,list_tmp_output,dict_chrLen,sv_class):
	global svclass_counter
	egene_dict=_getEnhancerGeneDict(egenes)
	for egene in egene_dict:
		for line in fh_bed_reset:
			fields = line.split("\t")
			chrH=fields[0]
			if not(ag.match(fields[3])):
				continue 
			elif(chrL == chrH):
				geneH=attr.split(fields[6])[1]
				startH=fields[1]; endH=fields[2];strandH=fields[5]
				geneTmp=re.split('\.',geneH)[0]
				if (geneTmp == egene) and (endL < startH): # enhancer upstream of a gene
					svclass_counter+=1
					for chr in dict_chrLen:
						list_tmp_output+=[BEDPE_LINE.format(chr,'0',dict_chrLen[chr],chrL,endL,(int(startH)-1),ID_FORMAT.format(chr,sv_class,svclass_counter,chrH),'.','.',strandH,sv_class,geneL,geneH)]
					break
				if (geneTmp == egene) and (startL > endH): # enhancer downstream of gene
					svclass_counter+=1
					for chr in dict_chrLen:
						list_tmp_output+=[BEDPE_LINE.format(chr,'0',dict_chrLen[chr],chrL,endH,startL,ID_FORMAT.format(chr,sv_class,svclass_counter,chrH),'.','.',strandH,sv_class,geneL,geneH)]
					break
			else: # if different chromosome
				break
	return list_tmp_output
			
def _getEnhancerAppositionBt(chrL,startL,endL,strandL,geneL,chr_armL,fh_bed_reset,list_tmp_output,dict_chrLen,sv_class):
	global svclass_counter
	(up_flank,down_flank)=getEnhancerFLank(chrL,startL,endL,dict_chrLen)
	for line in fh_bed_reset:
		fields = line.split("\t")
		chrH=fields[0]
		if not(allg.search(fields[3])):
			continue
		elif(chrL==chrH):
			continue
		elif(dg.match(fields[3]) ):
			startH=fields[1]; endH=fields[2];strandH=fields[5];chr_armH=fields[7].strip();geneH=attr.split(fields[6])[1]
			svclass_counter+=1
			if (chr_armL == chr_armH): #flank before e and flank after g
				list_tmp_output+=[BEDPE_LINE.format(chrL,up_flank,startL,chrH,startH,endH,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrH),'.','.',strandH,sv_class,'.',geneH)]
			else: #flank after e and flank after g
				list_tmp_output+=[BEDPE_LINE.format(chrL,startL,down_flank,chrH,startH,endH,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrH),'.','.',strandH,sv_class,'.',geneH)]
		elif(ag.match(fields[3])):
			startH=fields[1]; endH=fields[2];strandH=fields[5];chr_armH=fields[7].strip();geneH=attr.split(fields[6])[1]
			svclass_counter+=1
			if (chr_armL == chr_armH): #flank after e and flank before g
				list_tmp_output+=[BEDPE_LINE.format(chrL,startL,down_flank,chrH,startH,endH,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrH),'.','.',strandH,sv_class,'.',geneH)]
			else: #flank before e and flank before g
				list_tmp_output+=[BEDPE_LINE.format(chrL,up_flank,startL,chrH,startH,endH,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrH),'.','.',strandH,sv_class,'.',geneH)]
	return list_tmp_output
	
def getEnhancerFLank(chr,start,end,dict_chrLen):
	enhancer_flank=750000
	enhancer_upstream=(int(start)-enhancer_flank)
	enhancer_downstream=(int(end)+enhancer_flank)
	if (enhancer_upstream < 0):
		enhancer_upstream=0
	elif (enhancer_downstream > int(dict_chrLen[chr])):
		enhancer_downstream=dict_chrLen[chr]
	return (enhancer_upstream,enhancer_downstream)

#Ut-----------
def _getCopyLossUt(chrH,startH,endH,strandH,geneH,chr_arm,centromere,list_tmp_output,dict_chrLen,sv_class):
	global svclass_counter
	# explaination if gene on -ve strand ??	
	if chr_arm == 'q': # p arm retained
		for chr in dict_chrLen:
			svclass_counter+=1
			list_tmp_output+=[BEDPE_LINE.format(chr,'0',dict_chrLen[chr],chrH,centromere,endH,ID_FORMAT.format(chr,sv_class,svclass_counter,chrH),'.','.',strandH,sv_class,'.',geneH)]
	else: # p arm lost
			for chr in dict_chrLen:
			svclass_counter+=1
			list_tmp_output+=[BEDPE_LINE.format(chr,'0',dict_chrLen[chr],chrH,startH,centromere,ID_FORMAT.format(chr,sv_class,svclass_counter,chrH),'.','.',strandH,sv_class,'.',geneH)]

	return list_tmp_output

def _get5pOverhangUt(chrH,startH,endH,strandH,geneH,chr_arm,centromere,list_tmp_output,dict_chrLen,sv_class):
	global svclass_counter
		# explaination if gene on -ve strand ??	
	if (chr_arm == 'q' and strandH == '+'): # p arm retained 
		for chr in dict_chrLen:
			svclass_counter+=1
			list_tmp_output+=[BEDPE_LINE.format(chr,'0',dict_chrLen[chr],chrH,centromere,int(endH)-1,ID_FORMAT.format(chr,sv_class,svclass_counter,chrH),'.','.',strandH,sv_class,'.',geneH)]
	elif(chr_arm == 'q'): # p arm lost
			for chr in dict_chrLen:
			svclass_counter+=1
			list_tmp_output+=[BEDPE_LINE.format(chr,'0',dict_chrLen[chr],chrH,int(startH)-1,dict_chrLen[chr],ID_FORMAT.format(chr,sv_class,svclass_counter,chrH),'.','.',strandH,sv_class,'.',geneH)]
	elif (chr_arm == 'p' and strandH == '-'):
			for chr in dict_chrLen:
			svclass_counter+=1
			list_tmp_output+=[BEDPE_LINE.format(chr,'0',dict_chrLen[chr],chrH,int(startH)-1,centromere,ID_FORMAT.format(chr,sv_class,svclass_counter,chrH),'.','.',strandH,sv_class,'.',geneH)]
	else:
			for chr in dict_chrLen:
			svclass_counter+=1
			list_tmp_output+=[BEDPE_LINE.format(chr,'0',dict_chrLen[chr],chrH,'0',int(end)-1,ID_FORMAT.format(chr,sv_class,svclass_counter,chrH),'.','.',strandH,sv_class,'.',geneH)]
	
	return list_tmp_output

def _getGeneFusionUt(chrL,startL,endL,strandL,geneL,chr_arm,fh_bed_reset,list_tmp_output,sv_class):
	global svclass_counter
	for line in fh_bed_reset:
		fields = line.split("\t")
		chrH=fields[0];startH=fields[1]; endH=fields[2];strandH=fields[5]; 
		if not(ag.match(fields[3])):
			continue
		elif(fields[7] == chr_arm and strandL == strandH):
			geneH=attr.split(fields[6])[1];svclass_counter+=1
			list_tmp_output+=[BEDPE_LINE.format(chrL,startL,endL,chrH,startH,endH,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrH),'.',strandL,strandH,sv_class,geneL,geneH)]
		elif(fields[7] != chr_arm and strandL != strandH):
			geneH=attr.split(fields[6])[1];svclass_counter+=1
			list_tmp_output+=[BEDPE_LINE.format(chrL,startL,endL,chrH,startH,endH,ID_FORMAT.format(chrL,sv_class,svclass_counter,chrH),'.',strandL,strandH,sv_class,geneL,geneH)]
	return list_tmp_output

def _getEnhancerSeparationUt(chrL,startL,endL,strandL,geneL,egenes,fh_bed_reset,list_tmp_output,dict_chrLen,sv_class):
	global svclass_counter
	egene_dict=_getEnhancerGeneDict(egenes)
	for egene in egene_dict:
		for line in fh_bed_reset:
			fields = line.split("\t")
			chrH=fields[0]
			if not(ag.match(fields[3])):
				continue 
			elif(chrL == chrH):
				geneH=attr.split(fields[6])[1]
				startH=fields[1]; endH=fields[2];strandH=fields[5]
				geneTmp=re.split('\.',geneH)[0]
				if (geneTmp == egene) and (endL < startH): # enhancer upstream of a gene
					svclass_counter+=1
					for chr in dict_chrLen:
						list_tmp_output+=[BEDPE_LINE.format(chr,'0',dict_chrLen[chr],chrL,endL,(int(startH)-1),ID_FORMAT.format(chr,sv_class,svclass_counter,chrH),'.','.',strandH,sv_class,geneL,geneH)]
					break
				if (geneTmp == egene) and (startL > endH): # enhancer downstream of gene
					svclass_counter+=1
					for chr in dict_chrLen:
						list_tmp_output+=[BEDPE_LINE.format(chr,'0',dict_chrLen[chr],chrL,endH,startL,ID_FORMAT.format(chr,sv_class,svclass_counter,chrH),'.','.',strandH,sv_class,geneL,geneH)]
					break
			else: # if different chromosome
				break
	return list_tmp_output

#####################################################################################################################
## -- FUNCTION MAIN: 
## -- To produce SV annotation map for a genome of interest
#####################################################################################################################
def main():
   
   optParser = optparse.OptionParser( 
      
      usage = "%prog [options] -o outdir -b input_annotation -i genome_index_file -a annotation_class -c restrict_chr ",
      
      description=
         "This script takes as input an output directory (outdir), " +
         "input annotation bed file (input_annotation), " + 
				 "genome index file (genome_index_file), " +
				 "annotation class (annotation_class), " +
				 "restrict to a chromosome  (restrict_chr), " +
         "set of coordinates were reported for give annotation class")

   optParser.add_option( "-o", "--outputDirectory", type="string", dest="outdir",
      default = "", help = "Output directory " )

   optParser.add_option( "-b", "--annotationBed", type="string", dest="input_annotation",
      default = "", help = "Input bed file containing various annotations types and associated gene(s) " )

   optParser.add_option( "-i", "--indexFile", type="string", dest="genome_index_file",
      default = "", help = "genome index file containing chr and length " )

   optParser.add_option( "-a", "--annotationClass", type="string", dest="annotation_class",
      default = "", help = "Annotation calss [del,td,ins]" )
   optParser.add_option( "-c", "--restrictChr", type="str", dest="restrict_chr",
      default = "", help = "restrict to specific Chromosome" )
   optParser.add_option( "-v", "--verbose", action="store_true", dest="verbose",
      help = "suppress progress report and warnings" )

   if len( sys.argv ) == 1:
      print "pass here"
      optParser.print_help()
      sys.exit( 1 )

   (opts, args) = optParser.parse_args()
   
   if (not opts.outdir) or (not opts.input_annotation) or (not opts.genome_index_file) or (not opts.annotation_class) or (not opts.restrict_chr):
      sys.stderr.write( sys.argv[0] + ": Error - Please see the list of required parameters to provide on the command line \n" )
      sys.stderr.write( "Call with '-h' to get usage information.\n" )
      sys.exit( 1 )
      
   try:
			sys.stderr.write("Analysing data.\n")
			createAnnotationMap(opts.outdir,opts.input_annotation,opts.genome_index_file,opts.annotation_class,opts.restrict_chr)  
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
