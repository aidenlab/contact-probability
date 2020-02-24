#!/usr/bin/env python
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import sys
from bisect import bisect_left, bisect_right
from scipy import stats
import random

"""
1st argument is name of experiment
2nd argument is location of file containing contacts
3rd argument is genome (hg19 or mm10)

contacts should be prefiltered using awk
example with medium format merged_nodups:
<readname> <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <mapq2>
awk '{if($3 == $7 && $10 > 29 && $11 > 29){print $2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8}}' merged_nodups.txt > filtered_merged_nodups.txt

example with long format merged_nodups:
<str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <cigar1> <sequence1> <mapq2> <cigar2> <sequence2> <readname1> <readname2>
awk '{if($2 == $6 && $9 > 29 && $12 > 29){print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7}}' merged_nodups.txt > filtered_merged_nodups.txt

"""

description = sys.argv[1]
file_location = sys.argv[2]
genome = str(sys.argv[3])

# chrom.sizes
chrom_dot_sizes = [249904550,243199373,198022430,191535534,180915260,\
	171115067,159321559,146440111,141696573,135534747,\
	135046619,133851895,115169878,107349540,102531392,\
	90354753,81529607,78081510,59380841,63025520,48157577,\
	51304566,155270560,59373566] 

if genome == "mm10":
	chrom_dot_sizes = [195471971,182113224,160039680,156508116,151834684,\
		149736546,145441459,129401213,124595110,130694993,\
		122082543,120129022,120421639,124902244,104043685,\
		98207768,94987271,90702639,61431566,171031299]


def get_chrom_length(number):
	if number=="X":
		number=23
		if genome == "mm10":
			number=20
	return chrom_dot_sizes[int(number)-1]

def totalpairs(chrlength,d1,d2):
	if d2<=d1:	poss=0
	else: poss = (d2-d1-1)*chrlength-(d2*(d2-1)-d1*(d1+1))/2
	return poss


def run_code(exp_name, hicfile):
	figname=exp_name+'_genome_wide_cp.png'
	figname_split=exp_name+'_iolr_genome_wide_cp.png'
	histname=exp_name+'_hist_genome_wide_cp.txt'

	actualhist=numpy.zeros(1000)
	actualhist_split={}

	for type in ['inner','outer','left','right']:
		actualhist_split[type]=numpy.zeros(1000)

	possiblehist=numpy.zeros(1000)
	bin_edges=numpy.logspace(0.0,(numpy.log(get_chrom_length(1))/numpy.log(10)),1001)

	# compute actual contacts

	f=open(hicfile, 'r')
	next=f.readline()
	while next!="":
		s1,chr1,x1,s2,chr2,x2=next.split()[:6]
		chr1 = chr1.replace("chr","")
		chr2 = chr2.replace("chr","")
		if chr1 == chr2 and chr1 != "MT" and chr1 != "M" and chr1 != "Y":
			actualhist[bisect_left(bin_edges,abs(int(x2)-int(x1)))-1]+=1
			# check to make sure x2 corresponds to a further position
			if int(x2)<int(x1):
				s1, s2 = s2, s1

			if s1=='0'  and s2=='0':  
				type='right'
			elif s1=='16' and s2=='16': 
				type='left'
			elif s1=='0'  and s2=='16': 
				type='inner'
			elif s1=='16' and s2=='0':  
				type='outer'
			actualhist_split[type][bisect_left(bin_edges,abs(int(x2)-int(x1)))-1]+=1
		next=f.readline()
	f.close()

	# compute possible contacts
	num_max_chroms = len(chrom_dot_sizes)
	for i in range(1,num_max_chroms):
		chr_num=str(i)
		if i==num_max_chroms:
			chr_num='X'
		chrom_length=get_chrom_length(chr_num)
		for j in range(1000):
			possiblehist[j]+=totalpairs(chrom_length,int(bin_edges[j]),int(bin_edges[j+1])+1)
	#g=open('poss_contacts.txt','w')
	#for i in range(1000):
	#	g.write(str(possiblehist[i])+'\t'+str(bin_edges[i])+'\n')
	#g.close()
	# compute contact probability to make the final histogram
	histfile=open(histname,"w")
	for i in range(1000):
		actual=actualhist[i]
		if possiblehist[i]==0:
			actualhist[i]=0
		else:
			actualhist[i]=actualhist[i]/float(possiblehist[i])
		histfile.write(str(actualhist[i])+' '+str(bin_edges[i])+' '+str(actual)+' '+str(possiblehist[i])+'\n')
	histfile.close()

	###
	# Plot histogram
	###
	plotxleft = 100
	plotxright= 100000000
	fit1left =  30000
	fit1right = 300000
	fit2left =  300000
	fit2right = 3000000

	plt.loglog(bin_edges[0:1000],actualhist)
	plt.axis([plotxleft,plotxright,actualhist[-1]/100,actualhist[bisect_left(bin_edges,plotxleft)]*100])
	plt.xlabel('Distance(bp)')
	plt.ylabel('Contact Probability')	

	if plotxright > fit1right:
		x1=bin_edges[bisect_left(bin_edges,fit1left):bisect_right(bin_edges,fit1right)]
		x1c=actualhist[bisect_left(bin_edges,fit1left):bisect_right(bin_edges,fit1right)]
		scal1=numpy.polyfit(numpy.log(x1),numpy.log(x1c),1)
		y1=numpy.exp(scal1[0]*numpy.log(x1)+scal1[1]+.5)
		r1=stats.pearsonr(x1c,y1)
		plt.loglog(x1,y1,'r',linewidth=2.0)
		plt.legend((exp_name,str(scal1[0])[:8]+" "+str(r1[0])[:8]))
	else:
		plt.legend((exp_name))
	if plotxright > fit2right:
		x2=bin_edges[bisect_left(bin_edges,fit2left):bisect_right(bin_edges,fit2right)]
		x2c=actualhist[bisect_left(bin_edges,fit2left):bisect_right(bin_edges,fit2right)]
		scal2=numpy.polyfit(numpy.log(x2),numpy.log(x2c),1)
		y2=numpy.exp(scal2[0]*numpy.log(x2)+scal2[1]+.5)
		r2=stats.pearsonr(x2c,y2)
		plt.loglog(x2,y2,'g',linewidth=2.0)
		plt.legend((exp_name, str(scal1[0])[:8]+" "+str(r1[0])[:8], str(scal2[0])[:8]+" "+str(r2[0])[:8]))
	plt.savefig(figname)
	plt.close('all')

	# make plot for left/right/inner/outer-CP
	for type,hist in actualhist_split.items():
		for i in range(1000):
			if possiblehist[i]==0:
				hist[i]=0
			else:
				hist[i]=hist[i]/float(possiblehist[i])
		plt.loglog(bin_edges[0:1000],hist,label=type)
	plt.axis([plotxleft,plotxright,10e-10,10e-3])
	plt.xlabel('Distance(bp)')
	plt.ylabel('Contact Probability')
	plt.title(exp_name)
	plt.legend()
	plt.savefig(figname_split)


if __name__=="__main__":
	run_code(description, file_location)
