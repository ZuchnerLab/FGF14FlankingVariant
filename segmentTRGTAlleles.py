
import pandas
import time
import sys, getopt, os

def main ( argv ):
    SAMPLE=""
    try:
        opts, args = getopt.getopt(argv,"h",["sample=","help"])
    except getopt.GetoptError:
        print('segmentTRGTAlleles.py --sample=<SampleID>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('--sample'):
            SAMPLE=arg
        elif opt in ('h','--help'):
            print('segmentTRGTAlleles.py --sample=<SampleID>')
            sys.exit()
        else:
            print('segmentTRGTAlleles.py --sample=<SampleID>')
            sys.exit()

    # define functions

    # identify all motifs that account for at least 10% of the repeat region length with motifs either equal to those specified, or double that length, remove redundant entries
    def collectRepeats(row):
        from itertools import groupby
        def find_kgrams(repeatRegion, k):
            string=str(repeatRegion)
            kgrams = sorted(
                string[j:j+k]
                    for i in range(k)
                        for j in range(i, (len(string) - i) // k * k, k)
            )
            groups = [(k, len(list(g))) for k, g in groupby(kgrams)]
            return sorted(groups, key=lambda i: i[1], reverse=True)
        repeatLength=row.repeatLength
        motifLengths=set()
        approvedMotifs={}
        approvedMotifCounts={}
        for i in row.MotifLengths:
            motifLengths.add(i)
            motifLengths.add(2*i)
        multiplierFactor=round(max(motifLengths)/min(motifLengths))
        for k in motifLengths:
            motifCountsThisK=find_kgrams(row.Allele[19:-25],k)
            # keep only entries that account for at least 10% of the repeat region span and have a count of at least 3
            for motif, count in motifCountsThisK:
                if (((k*count)>=(repeatLength/10)) & (count>2)):
                    approvedMotifCounts[motif]=count
        return sorted(approvedMotifCounts.items(), key=lambda item: item[1], reverse=True)


    # fuzzy matching to look for repetitive stretches using the motifs identified above
    # edit: removed reverse complementation from the simplest_motif function because that operation only applies when comparing with EHDn, not when segmenting
    def defineFuzzySegments(row):
        def simplest_motif(motif):
            i = (motif+motif).find(motif, 1, -1)
            if i != -1:
                motif = motif[:i]
            motifs = []
            b = len(motif)
            for i in range(b):
                c = motif[i:]+motif[:i]
                motifs.append(c)
            motifs = [x.upper() for x in motifs]
            return(min(motifs))
        import regex as re
        def recursive_merge(inter, start_index = 0):
            for i in range(start_index, len(inter) - 1):
                if inter[i][1] >= inter[i+1][0]:
                    new_start = inter[i][0]
                    new_end = inter[i+1][1]
                    inter[i] = [new_start, new_end]
                    del inter[i+1]
                    return recursive_merge(inter.copy(), start_index=i)
            return inter    
        motifCounts=sorted(row.motifCounts, key=lambda item: len(item[0]), reverse=True) # sort by largest motif first, then by most counts
        motifSegments={}
        bestMotifVersionSpans={}
        bestMergedIntervals={}
        mismatchTolerance=4 # allow one variant every N motif copies with this value, but minimum of 12bp 
        for i in range(len(motifCounts)):
            thisMotif=motifCounts[i][0]
            if len(thisMotif)==2:
                mismatchTolerance=int(mismatchTolerance*1.5)
            elif len(thisMotif)==1:
                mismatchTolerance=mismatchTolerance*3
            intervals=[]
            for matchObj in re.finditer('(' + (thisMotif*mismatchTolerance) + '){e<=1}',row.Allele):
                intervals.append(matchObj.span())
            sorted_on_start = sorted(intervals)
            mergedIntervals = recursive_merge(sorted_on_start.copy())
            thisMotif=simplest_motif(thisMotif)
            spanCovered=0
            for l in range(len(mergedIntervals)):
                spanCovered=spanCovered+(mergedIntervals[l][1]-mergedIntervals[l][0])
            if thisMotif in bestMotifVersionSpans.keys():
                if spanCovered>bestMotifVersionSpans[thisMotif]:
                    bestMotifVersionSpans[thisMotif]=spanCovered
                    bestMergedIntervals[thisMotif]=mergedIntervals
            else:
                bestMotifVersionSpans[thisMotif]=spanCovered
                bestMergedIntervals[thisMotif]=mergedIntervals
        for motif in bestMergedIntervals.keys():
            if bestMergedIntervals[motif]!=[]:
                motifSegments[motif]=bestMergedIntervals[motif]

        # prune overlapping intervals: switch from dictionary to list of lists
        motifSegmentsDict=motifSegments
        motifSegments=[]
        for motif in motifSegmentsDict.keys():
            for j in range(len(motifSegmentsDict[motif])):
                motifSegments.append([motif,motifSegmentsDict[motif][j][0],motifSegmentsDict[motif][j][1],motifSegmentsDict[motif][j][1]-motifSegmentsDict[motif][j][0]])

        # prune overlapping intervals: sort by motif length, then by span length
        motifSegments=sorted(motifSegments, key=lambda item: item[3], reverse=True)
        motifSegments=sorted(motifSegments, key=lambda item: len(item[0]), reverse=True)
        i=0
        while i<len(motifSegments):
            # check if this interval overlaps any before it, and if so, trim it down
            for j in range(0,i):
                if ((motifSegments[i][1]<=motifSegments[j][2]) & (motifSegments[i][1]>motifSegments[j][1]) & (motifSegments[i][2]>motifSegments[j][2]+1)):
                    # trim left edge
                    motifSegments[i][1]=motifSegments[j][2]+1
                elif ((motifSegments[i][2]>=motifSegments[j][1]) & (motifSegments[i][1]<motifSegments[j][1]-1) & (motifSegments[i][2]<motifSegments[j][2])):
                    # trim right edge
                    motifSegments[i][2]=motifSegments[j][1]-1
                elif ((motifSegments[i][1]>=motifSegments[j][1]-1) & (motifSegments[i][2]<=motifSegments[j][2]+1)):
                    # discard this segment, it is subsumed by an earlier one
                    motifSegments.pop(i)
                    i=i-1
                    break
            # if the trimming of this segment has made it shorter than 12bp, drop it
            if ((motifSegments[i][2]-motifSegments[i][1])<12):
                motifSegments.pop(i)
                i=i-1
            i=i+1
        motifSegments=sorted(motifSegments, key=lambda item: item[1], reverse=False)

        # swtich back to dictionary
        motifSegmentsDict={}
        for i in range(len(motifSegments)):
            motif = motifSegments[i][0]
            if motif in motifSegmentsDict.keys():
                motifSegmentsDict[motif].append([motifSegments[i][1],motifSegments[i][2]])
            else:
                motifSegmentsDict[motif]=[[motifSegments[i][1],motifSegments[i][2]]]
        return motifSegmentsDict

    def sizeMajorMotifs(row):
        coordinates=row.motifSegments
        newDict = {}
        for x, y in coordinates.items():
            total = 0
            for i in y:
                diff = i[1] - i[0]
                total += diff
            newDict[x] = total
        return newDict

    def encodeMotifSegments(row):
        from Levenshtein import distance as levenshtein_distance
        motifSegmentsDict=row.motifSegments
        motifSegments=[]
        for motif in motifSegmentsDict.keys():
            for j in range(len(motifSegmentsDict[motif])):
                motifSegments.append([motif,motifSegmentsDict[motif][j][0],motifSegmentsDict[motif][j][1],motifSegmentsDict[motif][j][1]-motifSegmentsDict[motif][j][0]])
        # sort by motif length, then by span length
        motifSegments=sorted(motifSegments, key=lambda item: item[3], reverse=True)
        motifSegments=sorted(motifSegments, key=lambda item: len(item[0]), reverse=True)
        i=0
        while i<len(motifSegments):
            # check if this interval overlaps any before it, and if so, trim it down
            for j in range(0,i):
                if ((motifSegments[i][1]<=motifSegments[j][2]) & (motifSegments[i][1]>motifSegments[j][1]) & (motifSegments[i][2]>motifSegments[j][2]+1)):
                    # trim left edge
                    motifSegments[i][1]=motifSegments[j][2]+1
                elif ((motifSegments[i][2]>=motifSegments[j][1]) & (motifSegments[i][1]<motifSegments[j][1]-1) & (motifSegments[i][2]<motifSegments[j][2])):
                    # trim right edge
                    motifSegments[i][2]=motifSegments[j][1]-1
                elif ((motifSegments[i][1]>=motifSegments[j][1]-1) & (motifSegments[i][2]<=motifSegments[j][2]+1)):
                    # discard this segment, it is subsumed by an earlier one
                    motifSegments.pop(i)
                    i=i-1
                    break
            i=i+1
        motifSegments=sorted(motifSegments, key=lambda item: item[1], reverse=False)
        startPos=motifSegments[0][1]
        endPos=motifSegments[0][2]
        # determine whether this sequence has reference, variant of interest, or something else at the 5' end
        referenceSeq='TCTATGCAACCAACTTTCTGTGA'
        referenceSeq2='TCTATGCAACCAACTTTCTGAGA'
        referenceSeq3='TCTATGCAACCAACTTTCTGTGGA'
        C4AA='TCTATGCAACCAACTTTCTGTTAGTCATAGTACCCCAA'
        C4A='TCTATGCAACCAACTTTCTGTTAGTCATAGTACCCCA'
        C5AA='TCTATGCAACCAACTTTCTGTTAGTCATAGTACCCCCAA'
        C5A='TCTATGCAACCAACTTTCTGTTAGTCATAGTACCCCCA'
        C3AA='TCTATGCAACCAACTTTCTGTTAGTCATAGTACCCAA'
        C3A='TCTATGCAACCAACTTTCTGTTAGTCATAGTACCCA'
        C2AA='TCTATGCAACCAACTTTCTGTTAGTCATAGTACCAA'
        C2A='TCTATGCAACCAACTTTCTGTTAGTCATAGTACCA'
        C1AA='TCTATGCAACCAACTTTCTGTTAGTCATAGTACAA'
        C1A='TCTATGCAACCAACTTTCTGTTAGTCATAGTACA'
        shortSeq='TCTATGCAACCAACTTTCTGAA'
        
        refDist=levenshtein_distance(referenceSeq,row.Allele[:23])
        refDist2=levenshtein_distance(referenceSeq2,row.Allele[:23])
        refDist3=levenshtein_distance(referenceSeq3,row.Allele[:24])
        altDistC4AA=levenshtein_distance(C4AA,row.Allele[:38])
        altDistC4A=levenshtein_distance(C4A,row.Allele[:37])
        altDistC5AA=levenshtein_distance(C5AA,row.Allele[:39])
        altDistC5A=levenshtein_distance(C5A,row.Allele[:38])
        altDistC3AA=levenshtein_distance(C3AA,row.Allele[:37])
        altDistC3A=levenshtein_distance(C3A,row.Allele[:36])
        altDistC2AA=levenshtein_distance(C2AA,row.Allele[:36])
        altDistC2A=levenshtein_distance(C2A,row.Allele[:35])
        altDistC1AA=levenshtein_distance(C1AA,row.Allele[:35])
        altDistC1A=levenshtein_distance(C1A,row.Allele[:34])
        shortDist=levenshtein_distance(shortSeq,row.Allele[:22])
    
        encodedSequence=''
        if (altDistC4AA==0):
            encodedSequence='(C4AA)'
        elif (altDistC4A==0):
            encodedSequence='(C4A)'
        elif (altDistC5AA==0):
            encodedSequence='(C5AA)'
        elif (altDistC5A==0):
            encodedSequence='(C5A)'
        elif (altDistC3AA==0):
            encodedSequence='(C3AA)'
        elif (altDistC3A==0):
            encodedSequence='(C3A)'
        elif (altDistC2AA==0):
            encodedSequence='(C2AA)'
        elif (altDistC2A==0):
            encodedSequence='(C2A)'
        elif (altDistC1AA==0):
            encodedSequence='(C1AA)'
        elif (altDistC1A==0):
            encodedSequence='(C1A)'
        elif (refDist==0):
            encodedSequence='(5\'-RFS_' + str(refDist) + 'edits)'
        elif (refDist2==0):
            encodedSequence='(5\'-RFS_' + str(refDist2) + 'edits)'
        elif (refDist3==0):
            encodedSequence='(5\'-RFS_' + str(refDist3) + 'edits)'
        elif (shortDist==0):
            encodedSequence='(shortFlank)'
        else:
            encodedSequence='(otherFlank_' + str(refDist) + 'editsFromReference)'
        # write (motif)copies and any intervening sequence
        i=1
        fullSequence=row.Allele
        GAARepeatUnits=0
        encodedSequence=encodedSequence + "(" + motifSegments[0][0] + ")" + str(round(motifSegments[0][3]/len(motifSegments[0][0])))
        if (motifSegments[0][0]) in ('GAA','AGA','AAG'):
            GAARepeatUnits=round(motifSegments[0][3]/len(motifSegments[0][0]))
        while i<len(motifSegments):
            if (motifSegments[i][1]>(motifSegments[i-1][2]+1)):
                encodedSequence=encodedSequence + fullSequence[motifSegments[i-1][2]+1:motifSegments[i][1]]
            encodedSequence=encodedSequence + "(" + motifSegments[i][0] + ")" + str(round(motifSegments[i][3]/len(motifSegments[i][0])))
            if (motifSegments[i][0]) in ('GAA','AGA','AAG'):
                GAARepeatUnits=GAARepeatUnits+round(motifSegments[i][3]/len(motifSegments[i][0]))
            if (motifSegments[i][2]>endPos):
                endPos=motifSegments[i][2]
            i=i+1
        # calculate distance between remaining sequence and the regular 3' flank
        dist3P=levenshtein_distance('GAAATGTGTTTAAGAATTCCTCAATAAG',fullSequence[motifSegments[i-1][2]:])
        encodedSequence=encodedSequence + "(3Prime_" + str(dist3P) + "edits)"
        repeatLength=endPos-startPos
        # re-calculate motif purity
        purity=3*GAARepeatUnits/repeatLength
        if purity>1:
            purity=1
        elif purity<0:
            purity=0
        return pandas.Series([encodedSequence, repeatLength, purity],index=['encodedSequence','repeatLength','motifPurity'])
    


    starttime=time.time()
    sample=pandas.read_csv(SAMPLE + '.vcf.gz',sep='\t',low_memory=True,header=None,comment='#',
                           names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE'],
                           dtype={'CHROM':str,'POS':int,'ID':str,'REF':str,'ALT':str,'QUAL':str,'FILTER':str,'INFO':str,'FORMAT':str,'SAMPLE':str})
    sample[['TRID','END','MOTIFS','STRUC']]=sample.loc[:,'INFO'].str.split(';',expand=True)
    sample['TRID']=sample.loc[:,'TRID'].str[5:]
    # only evaluate the FGF14-SCA27B locus
    sample=sample.loc[sample['TRID']=='chr13_102161544_102161756',:].reset_index(drop=True)
    sample['END']=sample.loc[:,'END'].str[4:]
    sample['MOTIFS']=sample.loc[:,'MOTIFS'].str[7:]
    sample['MOTIFS']=sample.loc[:,'MOTIFS'].str.split(',',expand=False) # make MOTIFS a list within each cell
    sample['STRUC']=sample.loc[:,'STRUC'].str[6:]
    sample['MotifLengths']=sample.loc[:,'MOTIFS'].apply(lambda x: [len(entry) for entry in x]) # list of lengths of motifs
    sample[['GT','AL','ALCI','SD','MC','MS','AP','AM']]=sample.loc[:,'SAMPLE'].str.split(':',expand=True)
    sample['refLeftFlank']=sample.loc[:,'REF'].str[:25]
    sample['refRightFlank']=sample.loc[:,'REF'].str[-25:]
    sample[['Allele1','Allele2']]=sample.loc[:,'ALT'].str.split(',',expand=True)
    # fill in empty calls
    sample.loc[sample['Allele1'].isnull(),'Allele1']=sample.loc[sample['Allele1'].isnull(),'REF']
    sample.loc[sample['Allele1']=='.','Allele1']=sample.loc[sample['Allele1']=='.','REF']
    sample.loc[sample['Allele2'].isnull(),'Allele2']=sample.loc[sample['Allele2'].isnull(),'REF']
    sample.loc[sample['Allele2']=='.','Allele2']=sample.loc[sample['Allele2']=='.','REF']
    sample[['Allele1RepeatLength','Allele2RepeatLength']]=sample.loc[:,'AL'].str.split(',',expand=True)
    sample['Allele1RepeatLength']=sample.loc[:,'Allele1RepeatLength'].astype(int)-50
    sample['Allele2RepeatLength']=sample.loc[:,'Allele2RepeatLength'].astype(int)-50

    # convert to a one allele per row format. You can adjust this to keep whatever info you need going forward. 
    time1=time.time()
    print('checkpoint 1: ' + str(round(time1-starttime)) + ' seconds elapsed')
    melted=pandas.concat([sample.loc[:,['Allele1','Allele1RepeatLength','TRID','MOTIFS','STRUC','refLeftFlank','refRightFlank']].rename(columns={'Allele1':'Allele','Allele1RepeatLength':'repeatLength'}),sample.loc[:,['Allele2','Allele2RepeatLength','TRID','MOTIFS','STRUC','refLeftFlank','refRightFlank']].rename(columns={'Allele2':'Allele','Allele2RepeatLength':'repeatLength'})],ignore_index=True)
    starttime=time.time()
    melted['MotifLengths']=melted.apply(lambda row: [3],axis=1)
    melted['motifCounts']=melted.apply(collectRepeats,axis=1)
    endtime=time.time()
    print('Collect motifs: ' + str(endtime-starttime) + ' seconds')
    # identify segments
    starttime=time.time()
    melted['motifSegments']=melted.apply(defineFuzzySegments,axis=1)
    endtime=time.time()
    print('Identify segments: ' + str(endtime-starttime) + ' seconds')
    # get major motifs
    starttime=time.time()
    melted['majorMotifs']=melted.loc[:,'motifSegments'].apply(lambda x: list(x.keys()))
    endtime=time.time()
    print('Get major motifs: ' + str(endtime-starttime) + ' seconds')
    # sum up bases covered by each major motif
    starttime=time.time()
    melted['majorMotifsSizes']=melted.apply(sizeMajorMotifs,axis=1)
    endtime=time.time()
    print('Sum up bases: ' + str(endtime-starttime) + ' seconds')
    # encode motif segmentation
    starttime=time.time()
    melted[['encodedSequence','dist5P','dist3P']]=melted.apply(encodeMotifSegments,axis=1)
    endtime=time.time()
    print('Encode sequence: ' + str(endtime-starttime) + ' seconds')

    melted.to_csv(SAMPLE + '_FGF14_meltedWithMotifInfo.txt',sep='\t',index=False)
    return

if __name__ == "__main__":
    main(sys.argv[1:])
