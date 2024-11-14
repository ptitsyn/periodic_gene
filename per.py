# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 11:59:24 2018

@author: Andrey Ptitsyn
"""
import numpy as np
import math
from scipy import stats
from scipy import signal
from scipy.special import comb
from scipy.stats import chi2

#Savitzky-Golay Polynomial fit reduced for practicality and performance
#low sampling rates in gene expression data make most options obsolete
def polysmooth3(v):
    sv=np.array(v)
    sv[0]=(39*sv[0]+8*sv[1]-4*(sv[2]+sv[3]+sv[4])+sv[5]-2*sv[6])/42      
    sv[1]=(8*sv[0]+19*sv[1]+16*sv[2]+6*sv[3]-4*sv[4]-7*sv[5]+4*sv[7])/42
    sv[2]=(-4*sv[0]+16*sv[1]+19*sv[2]+12*sv[3]+2*sv[4]-4*sv[5]+sv[6])/42
    for i in range(3,sv.size-3): 
            sv[i]=(7*sv[i]+6*(sv[i+1]+sv[i-1])+3*(sv[i+2]+sv[i-2])-2*(sv[i+3]+sv[i-3]))/21;
    sv[-3]=(sv[-7]-4*sv[-6]+2*sv[-5]+12*sv[-4]+19*sv[-3]+16*sv[-2]-4*sv[-1])/42;
    sv[-2]=(4*sv[-7]-7*sv[-6]-4*sv[-5]+6*sv[-4]+16*sv[-3]+19*sv[-2]+8*sv[-1])/42;
    sv[-1]=(-2*sv[-7]+4*sv[-6]+sv[-5]-4*sv[-4]-4*sv[-3]+8*sv[-2]+39*sv[-1])/42;
    return sv

#Savitzky-Golay Polynomial fit reduced for practicality and performance
#low sampling rates in gene expression data make most options obsolete
#Here the integer number of periods is assumed in the timeline, so the last few
#values are used to smooth the first and vice versa     
def cyclicpolysmooth3(v):
    sv=np.array(v)
    sv[0]=(7*sv[0]+6*(sv[1]+sv[-1])+3*(sv[2]+sv[-2])-2*(sv[3]+sv[-3]))/21;
    sv[1]=(7*sv[1]+6*(sv[2]+sv[0])+3*(sv[3]+sv[-1])-2*(sv[4]+sv[-2]))/21;
    sv[2]=(7*sv[2]+6*(sv[3]+sv[1])+3*(sv[4]+sv[0])-2*(sv[5]+sv[-1]))/21;
    for i in range(3,sv.size-3): 
            sv[i]=(7*sv[i]+6*(sv[i+1]+sv[i-1])+3*(sv[i+2]+sv[i-2])-2*(sv[i+3]+sv[i-3]))/21;
    sv[-3]=(7*sv[-3]+6*(sv[-2]+sv[-4])+3*(sv[-1]+sv[-5])-2*(sv[0]+sv[-6]))/21;
    sv[-2]=(7*sv[-2]+6*(sv[-1]+sv[-3])+3*(sv[0]+sv[-4])-2*(sv[1]+sv[-5]))/21;
    sv[-1]=(7*sv[-1]+6*(sv[0]+sv[-2])+3*(sv[1]+sv[-3])-2*(sv[2]+sv[-4]))/21;
    return sv

#the function to estimate the most likely phase of oscillation


#Autocorrelation test
#returns the correleation corefficient of a time series to itself with a cyclic shift by one preiod
#and the p-value for that correlation
def autocorr(timeline, period):
    x = timeline.size // period - 1
    corr=np.zeros(x)
    ruler=timeline.values
    p_value=np.zeros(x)
    for i in range(x):
        ruler = np.roll(ruler, period)
        corr[i], p_value[i] = stats.pearsonr(timeline.values, ruler)
    ph=corr.argmax()
    best_corr=corr[ph]
    best_p=p_value[ph]
    return(best_corr,best_p)


#fisher's g-test for periodicity
#timeline is a series of measurements, period is the vave length to be tested, 
#measured in timepoints, i.g. 6 timepoints per day in a timeline of 12 points over two days 
def Fg_test(timeline,period):
    frequency=int(timeline.size/period)
    N = len(timeline)
    f, px_den = signal.periodogram(timeline)
    px_den=px_den[1:N//2]
    g=px_den[frequency]/px_den.sum()
    p_value = 1 - chi2.cdf(g * (N - 2), 2)
    return(p_value);

#permuted timeline test
#maxp is the number of recombinations in-phase, but not more then some reasonable number
#here this number is reduced to 100 for speed and simplicity: 6 bad out of 100 is like p=0.06,
#failure to reject the null hypothesis, etc.
def Pt_test(timeline, period):
    frequency=int(timeline.size/period)
    test=np.array(timeline)
    f, px_den = signal.periodogram(timeline)
    original=px_den[frequency]
    maxp=1000
    tally=0
    for i in range(maxp):
        np.random.shuffle(test)
        f, px_den = signal.periodogram(test)
        if(px_den[frequency] >= original):
            tally+=1
    return(tally/maxp)
        

def assign_phase(timeline, period):
    corr=np.zeros(period)
    ruler=np.zeros(timeline.size)
    t = np.arange(0.0, timeline.size, 1.0)
    ruler = np.cos(2 * np.pi * t / period)
    for i in range(period):
        corr[i], p_value = stats.pearsonr(timeline.values, ruler)
        ruler=np.roll(ruler,1)
    ph=corr.argmax()
    return(ph)


# tests periodicity in a phase continuum, i.e. in the continuous profile where all genes oscillating
# in the same phase are concatenated in order of decreasing signal/noise ratio
# returns the array of p-values in the same order as concatenated expression profiles
def slide_test_periodicity(timeline, num_genes, period, framelength, method):
    points_per_gene = timeline.size // num_genes
    results = np.zeros(num_genes)
    #frame = np.array(timeline[:framelength])
    if (method == 1):
        # useFisher's g-test
        #slide the frame by one gene at a time, although it might be reasonable to slide by a point or a period in some cases
        for i in range(num_genes-framelength//points_per_gene):
            framestart=i*points_per_gene
            results[i] = Fg_test(timeline[ framestart: framestart + framelength], period)  # test
        # last genes in the last frame receive the same test results for all
        results[i:] = results[i]
    elif (method == 2):
        # use Pt-test
        for i in range(num_genes-framelength//points_per_gene):
            framestart = i * points_per_gene
            results[i] = Pt_test(timeline[ framestart: framestart + framelength], period)  # test
        # last genes in the last frame receive the same test results for all
        results[i:]=results[i]
    return (results)


# testing all timelines in a phase continuum
# parameters: timelines - the matrix of expression profiles, period- number of timepoints per complete cycle
# framelength - how many genes to test in a single frame, i.e. multiple of the number of timepoints per gene expression profile
# cutoff is the p-value for the test and
# method defines the test, 1 for Fisher's g-test 2 for Pt-test
def phase_continuum_test(timelines, period, framelength, cutoff, method):
    num_genes = timelines.shape[0]
    timepoints = timelines.columns
    num_timepoints = timelines.shape[1]
    frequency = num_timepoints / period
    qual = np.zeros(num_genes)
    prob = np.zeros(num_genes)
    phases = np.zeros(num_genes, dtype=int)
    phase_sizes = np.zeros(period, dtype=int)
    # first test to get signal/noise ratio estimated by autocorrelation with a phase shift
    i=0
    for gene in timelines.index:
        qual[i], prob[i] = autocorr(timelines.loc[gene], period)
        phases[i] = assign_phase(timelines.loc[gene], period)
        phase_sizes[phases[i]] += 1
        i+=1
    # allocate space for phase continuums
    timelines['phase']=phases
    timelines['signal/noise']=np.abs(qual) #take best autocorrelation for esimation of signal to noise ratio
    timelines=timelines.sort_values(by=['phase','signal/noise'])
    pvals=[]
    #take one phase at a time
    for phase in range(period):
        phase_df=timelines[timelines['phase']==phase]
        #concat profiles in the order of decreasing signal/noise
        phase_cont=phase_df[timepoints].values.reshape(-1)
        #test for periodicity in a sliding frame, a few profiles at a time
        pvals.append( slide_test_periodicity(phase_cont, phase_df.shape[0], period, framelength, method) )
    pvalues=np.concatenate(pvals)
    #add the results column to the dataframe
    timelines['PhC p-values']=pvalues
    return (timelines)

