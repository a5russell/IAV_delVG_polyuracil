
import pandas as pd
import os

def reverseComplement(sequence):
    rc = ''
    sequence = sequence[::-1]
    for letter in sequence:
        if letter == 'A':
            rc += 'T'
        elif letter == 'T':
            rc += 'A'
        elif letter == 'C':
            rc += 'G'
        elif letter == 'G':
            rc += 'C'
    return rc

#get R1 and R2 sequiences and open relevant file.
def rowAssign(row, basename, fastaFiles, matchLength):
    templateFasta = os.path.join(fastaFiles, row.template) + '.fa'
    with open(templateFasta, 'r') as input:
        input.readline()
        template = input.readline().upper() #shsould have no newline characters, but just in case...
        while not template[-1:] in ['A','G','C','T']:
            template = template[:-1]
    if row.forward_primer == 'R1':
        R1 = template[:matchLength]
        R2 = reverseComplement(template[-matchLength:])
    else:
        R2 = template[:matchLength]
        R1 = reverseComplement(template[-matchLength:])
    filenameR1 = '_'.join([basename, row.template + '_R1.fastq'])
    filenameR2 = '_'.join([basename, row.template + '_R2.fastq'])
    if os.path.exists(filenameR1):
        return [R1, R2, open(filenameR1, 'a'), open(filenameR2, 'a')]
    else:
        return [R1, R2, open(filenameR1, 'w'), open(filenameR2, 'w')]

def rowAssignV2(row, fastaFiles, matchLength):
    templateFasta = os.path.join(fastaFiles, row.template) + '.fa'
    with open(templateFasta, 'r') as input:
        input.readline()
        template = input.readline().upper() #shsould have no newline characters, but just in case...
        while not template[-1:] in ['A','G','C','T']:
            template = template[:-1]
    if row.forward_primer == 'R1':
        R1 = template[:matchLength]
        R2 = reverseComplement(template[-matchLength:])
    else:
        R2 = template[:matchLength]
        R1 = reverseComplement(template[-matchLength:])
    return [R1, R2, [], []]

def rowPrintV2(row, basename):
    filenameR1 = '_'.join([basename, row.template + '_R1.fastq'])
    filenameR2 = '_'.join([basename, row.template + '_R2.fastq'])
    if os.path.exists(filenameR1):
        with open(filenameR1, 'a') as R1file:
            for element in row.openFileR1:
                R1file.write(element)
        with open(filenameR2, 'a') as R2file:
            for element in row.openFileR2:
                R2file.write(element)
    else:
        with open(filenameR1, 'w') as R1file:
            for element in row.openFileR1:
                R1file.write(element)
        with open(filenameR2, 'w') as R2file:
            for element in row.openFileR2:
                R2file.write(element)


def closeFiles(row):
    row.openFileR1.close()
    row.openFileR2.close()
    
#sorts into the different fastq files. RACE will need a second sort but too much of a pain to get it on first go.
def sortFastq(*, R1, R2, segment, database, basename, fastaFiles, matchLength, plasmid):
    #bit modified to not lose NA, annoying but whatever
    database = pd.read_csv(database, keep_default_na=False)
    database = database.loc[database.segment == segment,]
    if plasmid:
        database = database.loc[database.template.str[1:4] != 'RNA']
    else:
        database = database.loc[database.template.str[:4] != 'plas']
    files = {}        
    #get all the relevant sequences to aassign files, as well as open relevant files
    database[['R1','R2','openFileR1', 'openFileR2']] = database.apply(lambda row: rowAssign(row, basename, fastaFiles,matchLength), axis=1, result_type='expand')
    with open(R1, 'r') as R1, open(R2, 'r') as R2:
        R1name = R1.readline().split(' ')[0] + '\n'
        R2name = R2.readline().split(' ')[0] + '\n'
        while R1name[0] == '@':
            try:
                if R1name != R2name:
                    raise ValueError
            except:                
                print('dischordant fastq names')
                print(R1name)
                print(R2name)
                
            seqR1 = R1.readline()
            seqR2 = R2.readline()
            plusR1 = R1.readline()
            plusR2 = R2.readline()
            scoreR1 = R1.readline()
            scoreR2 = R2.readline()
            matchFrame = pd.merge(database, pd.DataFrame({'R1':[seqR1[:matchLength]], 'R2':[seqR2[:matchLength]]}), how='inner')
            try:
                if len(matchFrame) > 1:
                    raise ValueError
            except:
                print('ambiguous assignments, primers match? Maybe extend length or reconsider this pipeline')
                print(matchFrame)
            if len(matchFrame) == 1:
                file1 = matchFrame.openFileR1.iloc[0]
                file2 = matchFrame.openFileR2.iloc[0]
                file1.write(R1name)
                file1.write(seqR1)
                file1.write(plusR1)
                file1.write(scoreR1)
                file2.write(R2name)
                file2.write(seqR2)
                file2.write(plusR2)
                file2.write(scoreR2)
            R1name = R1.readline().split(' ')[0] + '\n'
            R2name = R2.readline().split(' ')[0] + '\n'
    database.apply(closeFiles, axis=1)
            

#sorts into the different fastq files. RACE will need a second sort but too much of a pain to get it on first go.
#v2 loads everything into a list and then writes to try and improve runtime at the cost of RAM usage.
def sortFastqV2(*, R1, R2, segment, database, basename, fastaFiles, matchLength, plasmid):
    database = pd.read_csv(database)
    database = database.loc[database.segment == segment,]
    count = 0
    if plasmid:
        database = database.loc[database.template.str[1:4] != 'RNA']
    else:
        database = database.loc[database.template.str[:4] != 'plas']
    files = {}        
    #get all the relevant sequences to aassign files, as well as open relevant files
    database[['R1','R2','openFileR1', 'openFileR2']] = database.apply(lambda row: rowAssignV2(row, fastaFiles,matchLength), axis=1, result_type='expand')
    with open(R1, 'r') as R1, open(R2, 'r') as R2:
        R1name = R1.readline().split(' ')[0] + '\n'
        R2name = R2.readline().split(' ')[0] + '\n'
        while R1name != '':
            try:
                if R1name != R2name:
                    raise ValueError
            except:                
                print('dischordant fastq names')
                print(R1name)
                print(R2name)
                
            seqR1 = R1.readline()
            seqR2 = R2.readline()
            plusR1 = R1.readline()
            plusR2 = R2.readline()
            scoreR1 = R1.readline()
            scoreR2 = R2.readline()
            matchFrame = pd.merge(database, pd.DataFrame({'R1':[seqR1[:matchLength]], 'R2':[seqR2[:matchLength]]}), how='inner')
            try:
                if len(matchFrame) > 1:
                    raise ValueError
            except:
                print('ambiguous assignments, primers match? Maybe extend length or reconsider this pipeline')
                print(matchFrame)
            if len(matchFrame) == 1:
                file1 = matchFrame.openFileR1.iloc[0]
                file2 = matchFrame.openFileR2.iloc[0]
                file1 += [''.join([R1name,seqR1,plusR1,scoreR1])]
                file2 += [''.join([R2name,seqR2,plusR2,scoreR2])]
            R1name = R1.readline().split(' ')[0] + '\n'
            R2name = R2.readline().split(' ')[0] + '\n'
            count += 1
            if count % 1000 == 0:
                print(count)
    database.apply(lambda row: rowPrintV2(row, basename), axis=1, result_type='expand')            



#list of bases recorded at each position that meet a minimum qscore. can then compare against genome later. 1-index rather than 0
def positionalData(*, sortedBamfile, outfile, qvalCutoff):
    tempFile = sortedBamfile.split('.')[0] + '_temp.sam'
    output = {}
    command  = ' '.join(['samtools view -f 2 -F 524', sortedBamfile,">", tempFile])
    os.system(command)
    with open(tempFile, 'r') as infile:
        while True:
            line1 = infile.readline()
            line2 = infile.readline()
            if not line2:
                break
            firstMate = line1.split('\t')
            secondMate = line2.split('\t')
            if firstMate[0] != secondMate[0]:
                print('Odd number of matepairs in file...fix')
                break
            #for each, store mapped positions (and small deletions, not counting them) in one list, and deletions in another. Then compare.
            #ignore chimeras
            if firstMate[2] == secondMate[2]:
                currentSeq = {}
                insertions ={}
                segment = firstMate[2]
                for read in [firstMate, secondMate]:
                    CIGAR = read[5]
                    currentPosSeq = int(read[3])
                    score = read[10]
                    sequence = read[9]
                    #start positions and mapping. Just go with lowest score at position rather than handling in a Bayesian framework. There is a lot of discussion about that.
                    CIGARval = []
                    curstring = ''
                    for character in CIGAR:
                        if not character.isalpha():
                            curstring += character
                        else:
                            CIGARval.append((int(curstring), character))
                            curstring = ''
                    currentPosRead = 1
                    #For D, N, advance position in sequence but not read. For I, S, advance position in read but not sequence
                    #For M, add in
                    for value in CIGARval:
                        if value[1] == 'M':
                            #add all mapped positions to a list
                            basesAdded = 0
                            while basesAdded < value[0]:
                                if currentPosSeq not in currentSeq or ord(score[currentPosRead-1]) - 33 > currentSeq[currentPosSeq][1] or currentSeq[currentPosSeq] == 'D':
                                    currentSeq[currentPosSeq] = (sequence[currentPosRead-1], ord(score[currentPosRead-1]) - 33)
                                basesAdded += 1
                                currentPosSeq += 1
                                currentPosRead +=1
                        elif value[1] == 'D':
                            basesAdded = 0
                            #only count if quality score to either side meets standard. Only called if base not already called at position.
                            if currentPosSeq not in currentSeq and ord(score[currentPosRead-2]) - 33 >= qvalCutoff and ord(score[currentPosRead]) - 33 >= qvalCutoff:
                                while basesAdded < value[0]:
                                    currentSeq[currentPosSeq] = ('D', qvalCutoff)
                                    basesAdded += 1
                                    currentPosSeq += 1
                                    #do not advance place in read
                            else:
                                currentPosSeq += value[0]
                        elif value[1] == 'N':
                            currentPosSeq += value[0]
                        #inserted sequence gets added at the position if bases before and after are a good score
                        elif value[1] == 'I':
                            if ord(score[currentPosRead-2]) - 33 >= qvalCutoff and ord(score[currentPosRead + value[0] - 1] ) - 33 >= qvalCutoff:
                                insertions[currentPosSeq] = value[0]
                                currentPosRead += value[0]
                        elif value[1] == 'S':
                            currentPosRead += value[0]

                for position in currentSeq:
                    if segment not in output:
                        output[segment] = {}
                    if position not in output[segment]:
                        output[segment][position] = {'A':0, 'G':0,'C':0,'T':0, 'D':0, 'I':0}
                    if currentSeq[position][1] >= qvalCutoff:
                        #for any reasonable qValCutoff none should be N's, but just to be anal retentive
                        base = currentSeq[position][0]
                        if base != 'N':
                            output[segment][position][base] += 1
                for position in insertions:
                    output[segment][position]['I'] += 1
        with open(outfile, 'w') as outfile:
            #might as well put it all in order in a tsv
            outfile.write('\t'.join(['segment','position','A','G','C','T', 'D','I']) + '\n')
            for segment in sorted(output.keys()):
                for position in sorted(output[segment]):
                    writeline = '\t'.join([segment, str(position)])
                    for base in ['A','G','C','T','D','I']:
                        writeline += '\t' + str(output[segment][position][base])
                    writeline += '\n'
                    outfile.write(writeline)
        arguments = f"rm -f {tempFile}"
        os.system(arguments)

#ignore positions covered by primers. Use known subamplicon coordinates and ignore given positions. Ignorepos given as a dictionary.
def positionalDataTrim(*, sortedBamfile, outfile, qvalCutoff, ignorePos):
    tempFile = sortedBamfile.split('.')[0] + '_temp.sam'
    output = {}
    command  = ' '.join(['samtools view -f 2 -F 524', sortedBamfile,">", tempFile])
    os.system(command)
    with open(tempFile, 'r') as infile:
        while True:
            line1 = infile.readline()
            line2 = infile.readline()
            if not line2:
                break
            firstMate = line1.split('\t')
            secondMate = line2.split('\t')
            if firstMate[0] != secondMate[0]:
                print('Odd number of matepairs in file...fix')
                break
            #for each, store mapped positions (and small deletions, not counting them) in one list, and deletions in another. Then compare.
            #ignore chimeras
            if firstMate[2] == secondMate[2]:
                currentSeq = {}
                insertions ={}
                segment = firstMate[2]
                for read in [firstMate, secondMate]:
                    CIGAR = read[5]
                    currentPosSeq = int(read[3])
                    score = read[10]
                    sequence = read[9]
                    #start positions and mapping. Just go with lowest score at position rather than handling in a Bayesian framework. There is a lot of discussion about that.
                    CIGARval = []
                    curstring = ''
                    for character in CIGAR:
                        if not character.isalpha():
                            curstring += character
                        else:
                            CIGARval.append((int(curstring), character))
                            curstring = ''
                    currentPosRead = 1
                    #For D, N, advance position in sequence but not read. For I, S, advance position in read but not sequence
                    #For M, add in
                    for value in CIGARval:
                        if value[1] == 'M':
                            #add all mapped positions to a list
                            basesAdded = 0
                            while basesAdded < value[0]:
                                if currentPosSeq not in currentSeq or ord(score[currentPosRead-1]) - 33 > currentSeq[currentPosSeq][1] or currentSeq[currentPosSeq] == 'D':
                                    currentSeq[currentPosSeq] = (sequence[currentPosRead-1], ord(score[currentPosRead-1]) - 33)
                                basesAdded += 1
                                currentPosSeq += 1
                                currentPosRead +=1
                        elif value[1] == 'D':
                            basesAdded = 0
                            #only count if quality score to either side meets standard. Only called if base not already called at position.
                            if currentPosSeq not in currentSeq and ord(score[currentPosRead-2]) - 33 >= qvalCutoff and ord(score[currentPosRead]) - 33 >= qvalCutoff:
                                while basesAdded < value[0]:
                                    currentSeq[currentPosSeq] = ('D', qvalCutoff)
                                    basesAdded += 1
                                    currentPosSeq += 1
                                    #do not advance place in read
                            else:
                                currentPosSeq += value[0]
                        elif value[1] == 'N':
                            currentPosSeq += value[0]
                        #inserted sequence gets added at the position if bases before and after are a good score
                        elif value[1] == 'I':
                            if ord(score[currentPosRead-2]) - 33 >= qvalCutoff and ord(score[currentPosRead + value[0] - 1] ) - 33 >= qvalCutoff:
                                insertions[currentPosSeq] = value[0]
                                currentPosRead += value[0]
                        elif value[1] == 'S':
                            currentPosRead += value[0]
                ignoreList = []
                #now "trim" Use a tuple-keyed list of positions unique to each amplicon and then ignore positions defined in thing
                for amplicon in ignorePos.keys():
                    if currentSeq.keys() & amplicon:
                        ignoreList += ignorePos[amplicon]
                if len(ignoreList) != 0:
                    for position in currentSeq:
                        if position not in ignoreList:
                            if segment not in output:
                                output[segment] = {}
                            if position not in output[segment]:
                                output[segment][position] = {'A':0, 'G':0,'C':0,'T':0, 'D':0, 'I':0}
                            if currentSeq[position][1] >= qvalCutoff:
                                #for any reasonable qValCutoff none should be N's, but just to be anal retentive
                                base = currentSeq[position][0]
                                if base != 'N':
                                    output[segment][position][base] += 1
                    for position in insertions:
                        if position not in ignoreList:
                            output[segment][position]['I'] += 1
        with open(outfile, 'w') as outfile:
            #might as well put it all in order in a tsv
            outfile.write('\t'.join(['segment','position','A','G','C','T', 'D','I']) + '\n')
            for segment in sorted(output.keys()):
                for position in sorted(output[segment]):
                    writeline = '\t'.join([segment, str(position)])
                    for base in ['A','G','C','T','D','I']:
                        writeline += '\t' + str(output[segment][position][base])
                    writeline += '\n'
                    outfile.write(writeline)
        arguments = f"rm -f {tempFile}"


#quick and dirty to handle caps before I write something better that takes it pre-alignemnt
def positionalDataIgnoreInsert(*, sortedBamfile, outfile, qvalCutoff, insertionCutoff):
    tempFile = sortedBamfile.split('.')[0] + '_temp.sam'
    output = {}
    command  = ' '.join(['samtools view -f 2 -F 524', sortedBamfile,">", tempFile])
    os.system(command)
    with open(tempFile, 'r') as infile:
        while True:
            line1 = infile.readline()
            line2 = infile.readline()
            if not line2:
                break
            firstMate = line1.split('\t')
            secondMate = line2.split('\t')
            if firstMate[0] != secondMate[0]:
                print('Odd number of matepairs in file...fix')
                break
            #for each, store mapped positions (and small deletions, not counting them) in one list, and deletions in another. Then compare.
            #ignore chimeras
            if firstMate[2] == secondMate[2]:
                currentSeq = {}
                insertions ={}
                segment = firstMate[2]
                insertionTooLarge = False
                for read in [firstMate, secondMate]:
                    CIGAR = read[5]
                    currentPosSeq = int(read[3])
                    score = read[10]
                    sequence = read[9]
                    #start positions and mapping. Just go with lowest score at position rather than handling in a Bayesian framework. There is a lot of discussion about that.
                    CIGARval = []
                    curstring = ''
                    for character in CIGAR:
                        if not character.isalpha():
                            curstring += character
                        else:
                            if character == 'I':
                                if int(curstring) >= insertionCutoff:
                                    insertionTooLarge = True
                            CIGARval.append((int(curstring), character))
                            curstring = ''
                    currentPosRead = 1
                    #For D, N, advance position in sequence but not read. For I, S, advance position in read but not sequence
                    #For M, add in
                    for value in CIGARval:
                        if value[1] == 'M':
                            #add all mapped positions to a list
                            basesAdded = 0
                            while basesAdded < value[0]:
                                if currentPosSeq not in currentSeq or ord(score[currentPosRead-1]) - 33 > currentSeq[currentPosSeq][1] or currentSeq[currentPosSeq] == 'D':
                                    currentSeq[currentPosSeq] = (sequence[currentPosRead-1], ord(score[currentPosRead-1]) - 33)
                                basesAdded += 1
                                currentPosSeq += 1
                                currentPosRead +=1
                        elif value[1] == 'D':
                            basesAdded = 0
                            #only count if quality score to either side meets standard. Only called if base not already called at position.
                            if currentPosSeq not in currentSeq and ord(score[currentPosRead-2]) - 33 >= qvalCutoff and ord(score[currentPosRead]) - 33 >= qvalCutoff:
                                while basesAdded < value[0]:
                                    currentSeq[currentPosSeq] = ('D', qvalCutoff)
                                    basesAdded += 1
                                    currentPosSeq += 1
                                    #do not advance place in read
                            else:
                                currentPosSeq += value[0]
                        elif value[1] == 'N':
                            currentPosSeq += value[0]
                        #inserted sequence gets added at the position if bases before and after are a good score
                        elif value[1] == 'I':
                            if ord(score[currentPosRead-2]) - 33 >= qvalCutoff and ord(score[currentPosRead + value[0] - 1] ) - 33 >= qvalCutoff:
                                insertions[currentPosSeq] = value[0]
                                currentPosRead += value[0]
                        elif value[1] == 'S':
                            currentPosRead += value[0]
                ignoreList = []
                #now "trim" Use a tuple-keyed list of positions unique to each amplicon and then ignore positions defined in thing
                if not insertionTooLarge:
                    for position in currentSeq:
                        if position not in ignoreList:
                            if segment not in output:
                                output[segment] = {}
                            if position not in output[segment]:
                                output[segment][position] = {'A':0, 'G':0,'C':0,'T':0, 'D':0, 'I':0}
                            if currentSeq[position][1] >= qvalCutoff:
                                #for any reasonable qValCutoff none should be N's, but just to be anal retentive
                                base = currentSeq[position][0]
                                if base != 'N':
                                    output[segment][position][base] += 1
                    for position in insertions:
                        if position not in ignoreList:
                            output[segment][position]['I'] += 1
        with open(outfile, 'w') as outfile:
            #might as well put it all in order in a tsv
            outfile.write('\t'.join(['segment','position','A','G','C','T', 'D','I']) + '\n')
            for segment in sorted(output.keys()):
                for position in sorted(output[segment]):
                    writeline = '\t'.join([segment, str(position)])
                    for base in ['A','G','C','T','D','I']:
                        writeline += '\t' + str(output[segment][position][base])
                    writeline += '\n'
                    outfile.write(writeline)
        arguments = f"rm -f {tempFile}"



#seperates presumptive cap-containing sequences, and their mate-pairs from uncapped.
#Due to a simple fuzzy-matching rather than a true alignment, this technically will have issues with indels.
#We accept that issue.
#Do prior to any trimming, expect all reads have a mate.
#trims to presumptive 5' end. Keep caps in seperate file, keep cap lengths 
def capSortFastq(*, capRead, pairedRead, outfile_capped_prefix, outfile_uncapped_prefix, capLength, vRNAseq, cRNAseq, adapterSeq, toleratedfivePrime, toleratedAdapter, outfile_caps_capped_prefix, outfile_caps_uncapped_prefix, outStatsPrefix, primerAdapter, primerAdapterLimit, uncapLimit):
    #append prefix to existing base file names for outfiles
    outCapped = outfile_capped_prefix + os.path.basename(capRead)
    outCappedpaired = outfile_capped_prefix + os.path.basename(pairedRead)
    outUncapped = outfile_uncapped_prefix + os.path.basename(capRead)
    outUncappedpaired =outfile_uncapped_prefix + os.path.basename(pairedRead)
    outCappedcaps = outfile_caps_capped_prefix + 'capSequences_' + os.path.basename(capRead)
    outUncappedcaps = outfile_caps_uncapped_prefix + 'uncappedCapSequences_' + os.path.basename(capRead)
    capLengthscRNA = []
    capLengthsvRNA = []
    uncapLengths = []

    
    #open all files, will write as a stream to each. Throw away reads that do not meet the indicated requirements for now.
    #keep track of numbers with outStats just so we have some kind of QC if we need to alter code or at least a heads-up on what is going on
    with open(capRead, 'r') as capRead, open(pairedRead, 'r') as pairedRead, open(outCapped, 'w') as outCapped, open(outCappedpaired, 'w') as outCappedpaired, open(outUncapped, 'w') as outUncapped, open(outUncappedpaired, 'w') as outUncappedpaired, open(outCappedcaps, 'w') as outCappedcaps, open(outUncappedcaps, 'w') as outUncappedcaps:
        numPairs = 0
        numCapped = 0
        numUncapped = 0
        #order is name, sequence, throw away line ,and score. For this ignore score (but output it)
        nameCap = capRead.readline()
        namePair = pairedRead.readline()
        statistics = {'Total_sequences':0, 'retained_uncapped':0,'retained_capped':0, 'uncapped_vRNA':0, 'capped_vRNA':0,'uncapped_cRNA': 0, 'capped_cRNA': 0}
        while nameCap != '':
            #make sure read names match. Sometimes stupid things like R1 and R2 get appended, so just grab the chip coordinates
            assert nameCap.split(' ')[0] == namePair.split(' ')[0]
            seqCap = capRead.readline()
            seqPair = pairedRead.readline()
            toss = capRead.readline()
            toss = pairedRead.readline()
            scoreCap = capRead.readline()
            scorePair = pairedRead.readline()
            # identify if adapter present in read. If not, chuck.
            numDif = len(adapterSeq)
            adapterPos = 0
            currPos = 0
            #identify best match to adapter sequence, does not tolerate indels
            while currPos <= len(seqCap) - len(adapterSeq):
                comparisonSeq = seqCap[currPos:currPos + len(adapterSeq)]
                newDif  = sum(1 for a, b in zip(adapterSeq, comparisonSeq) if a != b)
                if newDif < numDif:
                    adapterPos = currPos
                    numDif = newDif
                currPos += 1
            if numDif <= toleratedAdapter:
                currPos = adapterPos + len(adapterSeq)
                fiveEnd = 0
                numDif = len(vRNAseq)
                fivePrimepos = 0
                while currPos <= len(seqCap) - len(vRNAseq):
                    #test v and c go with best match
                    comparisonSeq = seqCap[currPos:currPos + len(vRNAseq)]
                    newDif  = sum(1 for a, b in zip(vRNAseq, comparisonSeq) if a != b)
                    if newDif < numDif:
                        fiveEnd = currPos
                        numDif = newDif
                    currPos += 1
                currPos = adapterPos + len(adapterSeq)
                #flag for statistics. For distributing in cap vs non actually irrelevant
                cRNA= False

                while currPos <= len(seqCap) - len(cRNAseq):
                    #test v and c go with best match
                    comparisonSeq = seqCap[currPos:currPos + len(cRNAseq)]
                    newDif  = sum(1 for a, b in zip(cRNAseq, comparisonSeq) if a != b)
                    if newDif < numDif:
                        fiveEnd = currPos
                        numDif = newDif
                        cRNA = True
                    currPos += 1
                #to pass this gate we need to meet requirements for sequencing adapter (TSO) and the five prime of either vRNA or cRNA. Might as well report that too
                if numDif <= toleratedfivePrime:
                    #if insertion meets criteria store in capped
                    insertionLength = fiveEnd - adapterPos - len(adapterSeq)
                    if insertionLength >= capLength:
                        #check if insertion is an erroneous primer
                        capSeq = seqCap[adapterPos + len(adapterSeq):fiveEnd]
                        numDif = len(primerAdapter)
                        currPos = 0
                        while currPos <= len(capSeq) - len(primerAdapter):
                            comparisonSeq = capSeq[currPos:currPos + len(primerAdapter)]
                            newDif  = sum(1 for a, b in zip(primerAdapter, comparisonSeq) if a != b)
                            if newDif < numDif:
                                numDif = newDif
                            currPos += 1
                        #only write if the cap does not just look like a primer
                        if numDif > primerAdapterLimit:
                            if cRNA:
                                capLengthscRNA += [insertionLength] 
                            else:
                                capLengthsvRNA += [insertionLength] 
                            outCapped.write(nameCap)
                            outCapped.write(seqCap[fiveEnd:])
                            outCapped.write(toss)
                            outCapped.write(scoreCap[fiveEnd:])
                            outCappedpaired.write(namePair)
                            outCappedpaired.write(seqPair)
                            outCappedpaired.write(toss)
                            outCappedpaired.write(scorePair)
                            if cRNA:
                                outCappedcaps.write(nameCap[:-1] + ' cRNA\n')
                            else:
                                outCappedcaps.write(nameCap[:-1] + ' vRNA\n')
    
                            #write just cap sequence, exclude adapter. This will include the G's from TSO, so keep that in mind for any processing
                            outCappedcaps.write(seqCap[adapterPos + len(adapterSeq):fiveEnd] + '\n' )
                            outCappedcaps.write(toss)
                            outCappedcaps.write(scoreCap[adapterPos + len(adapterSeq):fiveEnd] + '\n' )
                            statistics['retained_capped'] += 1
                            if cRNA:
                                statistics['capped_cRNA'] += 1
                            else:
                                statistics['capped_vRNA'] += 1
                    elif insertionLength == uncapLimit:
                        uncapLengths += [insertionLength] 
                        outUncapped.write(nameCap)
                        outUncapped.write(seqCap[fiveEnd:])
                        outUncapped.write(toss)
                        outUncapped.write(scoreCap[fiveEnd:])
                        outUncappedpaired.write(namePair)
                        outUncappedpaired.write(seqPair)
                        outUncappedpaired.write(toss)
                        outUncappedpaired.write(scorePair)
                        statistics['retained_uncapped'] += 1
                        outUncappedcaps.write(nameCap)
                        #write just cap sequence, exclude adapter. This will include the G's from TSO, so keep that in mind for any processing
                        outUncappedcaps.write(seqCap[adapterPos + len(adapterSeq):fiveEnd] + '\n' )
                        outUncappedcaps.write(toss)
                        outUncappedcaps.write(scoreCap[adapterPos + len(adapterSeq):fiveEnd] + '\n' )
                        if cRNA:
                            statistics['uncapped_cRNA'] += 1
                        else:
                            statistics['uncapped_vRNA'] += 1
                        outUncappedcaps
            statistics['Total_sequences'] += 1
            nameCap = capRead.readline()
            namePair = pairedRead.readline() 
    pd.DataFrame(statistics.items(), columns=['Category', 'number']).to_csv(outStatsPrefix + '_totalStats.csv')
    pd.DataFrame({'Lengths':capLengthscRNA}).to_csv(outStatsPrefix + '_cappedLengthscRNA.csv')  
    pd.DataFrame({'Lengths':capLengthsvRNA}).to_csv(outStatsPrefix + '_cappedLengthsvRNA.csv')  

    pd.DataFrame({'Lengths':uncapLengths}).to_csv( outStatsPrefix + '_uncappedLengths.csv')  




#seperates presumptive cap-containing sequences, and their mate-pairs from uncapped.
#Due to a simple fuzzy-matching rather than a true alignment, this technically will have issues with indels.
#We accept that issue.
#Do prior to any trimming, expect all reads have a mate.
#trims to presumptive 5' end. Keep caps in seperate file, keep cap lengths 
#added a line of code to remove from 3' end if run-on
def capSortFastqOneEnd(*,capReadNumber, capRead, pairedRead, outfileCappedPrefix, outfileUnCappedPrefix, capLength, cRNAseq, adapterSeq, toleratedfivePrime, toleratedAdapter, outStatsPrefix, primerAdapter, primerAdapterLimit, uncapLimit, trim3primeCapRead, trim3primepairedRead):
    #append prefix to existing base file names for outfiles
    if capReadNumber == 'R2':
        capPaired = 'R1'
    else:
        capPaired = 'R2'
    outCapped = outfileCappedPrefix + '_' + capReadNumber + '.fastq'
    outCappedpaired = outfileCappedPrefix + '_' + capPaired + '.fastq'
    outUncapped = outfileUnCappedPrefix + '_' + capReadNumber + '.fastq'
    outUncappedpaired =outfileUnCappedPrefix + '_' + capPaired + '.fastq'
    capLengthscRNA = []
    capLengthsvRNA = []
    uncapLengths = []

    
    #open all files, will write as a stream to each. Throw away reads that do not meet the indicated requirements for now.
    #keep track of numbers with outStats just so we have some kind of QC if we need to alter code or at least a heads-up on what is going on
    with open(capRead, 'r') as capRead, open(pairedRead, 'r') as pairedRead, open(outCapped, 'w') as outCapped, open(outCappedpaired, 'w') as outCappedpaired, open(outUncapped, 'w') as outUncapped, open(outUncappedpaired, 'w') as outUncappedpaired:
        numPairs = 0
        numCapped = 0
        numUncapped = 0
        #order is name, sequence, throw away line ,and score. For this ignore score (but output it)
        nameCap = capRead.readline()
        namePair = pairedRead.readline()
        statistics = {'Total_sequences':0, 'retained_uncapped':0,'retained_capped':0, 'uncapped_cRNA': 0, 'capped_cRNA': 0}
        while nameCap != '':
            #make sure read names match. Sometimes stupid things like R1 and R2 get appended, so just grab the chip coordinates
            assert nameCap.split(' ')[0] == namePair.split(' ')[0]
            seqCap = capRead.readline()[:-trim3primeCapRead-1] + '\n'
            seqPair = pairedRead.readline()[:-trim3primepairedRead-1] + '\n'
            toss = capRead.readline()
            toss = pairedRead.readline()
            scoreCap = capRead.readline()[:-trim3primeCapRead-1] + '\n'
            scorePair = pairedRead.readline()[:-trim3primepairedRead-1] + '\n'
            # identify if adapter present in read. If not, chuck.
            numDif = len(adapterSeq)
            adapterPos = 0
            currPos = 0
            #identify best match to adapter sequence, does not tolerate indels
            while currPos <= len(seqCap) - len(adapterSeq):
                comparisonSeq = seqCap[currPos:currPos + len(adapterSeq)]
                newDif  = sum(1 for a, b in zip(adapterSeq, comparisonSeq) if a != b)
                if newDif < numDif:
                    adapterPos = currPos
                    numDif = newDif
                currPos += 1
            if numDif <= toleratedAdapter:
                currPos = adapterPos + len(adapterSeq)
                fiveEnd = 0
                numDif = len(cRNAseq)
                fivePrimepos = 0
                while currPos <= len(seqCap) - len(cRNAseq):
                    comparisonSeq = seqCap[currPos:currPos + len(cRNAseq)]
                    newDif  = sum(1 for a, b in zip(cRNAseq, comparisonSeq) if a != b)
                    if newDif < numDif:
                        fiveEnd = currPos
                        numDif = newDif
                    currPos += 1
                currPos = adapterPos + len(adapterSeq)
                #flag for statistics. For distributing in cap vs non actually irrelevant
                cRNA= True
                #to pass this gate we need to meet requirements for sequencing adapter (TSO) and the five prime of either vRNA or cRNA. Might as well report that too
                if numDif <= toleratedfivePrime:
                    #if insertion meets criteria store in capped
                    insertionLength = fiveEnd - adapterPos - len(adapterSeq)
                    if insertionLength >= capLength:
                        #check if insertion is an erroneous primer
                        capSeq = seqCap[adapterPos + len(adapterSeq):fiveEnd]
                        numDif = len(primerAdapter)
                        currPos = 0
                        while currPos <= len(capSeq) - len(primerAdapter):
                            comparisonSeq = capSeq[currPos:currPos + len(primerAdapter)]
                            newDif  = sum(1 for a, b in zip(primerAdapter, comparisonSeq) if a != b)
                            if newDif < numDif:
                                numDif = newDif
                            currPos += 1
                        #only write if the cap does not just look like a primer
                        if numDif > primerAdapterLimit:
                            if cRNA:
                                capLengthscRNA += [insertionLength] 
                            else:
                                capLengthsvRNA += [insertionLength] 
                            outCapped.write(nameCap)
                            outCapped.write(seqCap[fiveEnd:])
                            outCapped.write(toss)
                            outCapped.write(scoreCap[fiveEnd:])
                            outCappedpaired.write(namePair)
                            outCappedpaired.write(seqPair)
                            outCappedpaired.write(toss)
                            outCappedpaired.write(scorePair)

                            statistics['retained_capped'] += 1
                            if cRNA:
                                statistics['capped_cRNA'] += 1
                            else:
                                statistics['capped_vRNA'] += 1
                    elif insertionLength == uncapLimit:
                        uncapLengths += [insertionLength] 
                        outUncapped.write(nameCap)
                        outUncapped.write(seqCap[fiveEnd:])
                        outUncapped.write(toss)
                        outUncapped.write(scoreCap[fiveEnd:])
                        outUncappedpaired.write(namePair)
                        outUncappedpaired.write(seqPair)
                        outUncappedpaired.write(toss)
                        outUncappedpaired.write(scorePair)
                        statistics['retained_uncapped'] += 1
                        if cRNA:
                            statistics['uncapped_cRNA'] += 1
                        else:
                            statistics['uncapped_vRNA'] += 1
            statistics['Total_sequences'] += 1
            nameCap = capRead.readline()
            namePair = pairedRead.readline() 
    pd.DataFrame(statistics.items(), columns=['Category', 'number']).to_csv(outStatsPrefix + '_totalStats.csv')
    pd.DataFrame({'Lengths':capLengthscRNA}).to_csv(outStatsPrefix + '_cappedLengthscRNA.csv')  
    pd.DataFrame({'Lengths':uncapLengths}).to_csv( outStatsPrefix + '_uncappedLengths.csv')  


def capSortFastqOtherEnd(*,capReadNumber, capRead, pairedRead, outfileCappedPrefix, outfileUnCappedPrefix, capLength, vRNAseq, adapterSeq, toleratedfivePrime, toleratedAdapter, outStatsPrefix, primerAdapter, primerAdapterLimit, uncapLimit, trim3primeCapRead, trim3primepairedRead):
    #append prefix to existing base file names for outfiles
    if capReadNumber == 'R2':
        capPaired = 'R1'
    else:
        capPaired = 'R2'
    outCapped = outfileCappedPrefix + '_' + capReadNumber + '.fastq'
    outCappedpaired = outfileCappedPrefix + '_' + capPaired + '.fastq'
    outUncapped = outfileUnCappedPrefix + '_' + capReadNumber + '.fastq'
    outUncappedpaired =outfileUnCappedPrefix + '_' + capPaired + '.fastq'
    capLengthscRNA = []
    capLengthsvRNA = []
    uncapLengths = []

    
    #open all files, will write as a stream to each. Throw away reads that do not meet the indicated requirements for now.
    #keep track of numbers with outStats just so we have some kind of QC if we need to alter code or at least a heads-up on what is going on
    with open(capRead, 'r') as capRead, open(pairedRead, 'r') as pairedRead, open(outCapped, 'w') as outCapped, open(outCappedpaired, 'w') as outCappedpaired, open(outUncapped, 'w') as outUncapped, open(outUncappedpaired, 'w') as outUncappedpaired:
        numPairs = 0
        numCapped = 0
        numUncapped = 0
        #order is name, sequence, throw away line ,and score. For this ignore score (but output it)
        nameCap = capRead.readline()
        namePair = pairedRead.readline()
        statistics = {'Total_sequences':0, 'retained_uncapped':0,'retained_capped':0, 'uncapped_vRNA': 0, 'capped_vRNA': 0}
        while nameCap != '':
            #make sure read names match. Sometimes stupid things like R1 and R2 get appended, so just grab the chip coordinates
            assert nameCap.split(' ')[0] == namePair.split(' ')[0]
            seqCap = capRead.readline()[:-trim3primeCapRead-1] + '\n'
            seqPair = pairedRead.readline()[:-trim3primepairedRead-1] + '\n'
            toss = capRead.readline()
            toss = pairedRead.readline()
            scoreCap = capRead.readline()[:-trim3primeCapRead-1] + '\n'
            scorePair = pairedRead.readline()[:-trim3primepairedRead-1] + '\n'
            # identify if adapter present in read. If not, chuck.
            numDif = len(adapterSeq)
            adapterPos = 0
            currPos = 0
            #identify best match to adapter sequence, does not tolerate indels
            while currPos <= len(seqCap) - len(adapterSeq):
                comparisonSeq = seqCap[currPos:currPos + len(adapterSeq)]
                newDif  = sum(1 for a, b in zip(adapterSeq, comparisonSeq) if a != b)
                if newDif < numDif:
                    adapterPos = currPos
                    numDif = newDif
                currPos += 1
            if numDif <= toleratedAdapter:
                currPos = adapterPos + len(adapterSeq)
                fiveEnd = 0
                numDif = len(vRNAseq)
                fivePrimepos = 0
                while currPos <= len(seqCap) - len(vRNAseq):
                    comparisonSeq = seqCap[currPos:currPos + len(vRNAseq)]
                    newDif  = sum(1 for a, b in zip(vRNAseq, comparisonSeq) if a != b)
                    if newDif < numDif:
                        fiveEnd = currPos
                        numDif = newDif
                    currPos += 1
                currPos = adapterPos + len(adapterSeq)
                #flag for statistics. For distributing in cap vs non actually irrelevant
                cRNA= False
                #to pass this gate we need to meet requirements for sequencing adapter (TSO) and the five prime of either vRNA or cRNA. Might as well report that too
                if numDif <= toleratedfivePrime:
                    #if insertion meets criteria store in capped
                    insertionLength = fiveEnd - adapterPos - len(adapterSeq)
                    if insertionLength >= capLength:
                        #check if insertion is an erroneous primer
                        capSeq = seqCap[adapterPos + len(adapterSeq):fiveEnd]
                        numDif = len(primerAdapter)
                        currPos = 0
                        while currPos <= len(capSeq) - len(primerAdapter):
                            comparisonSeq = capSeq[currPos:currPos + len(primerAdapter)]
                            newDif  = sum(1 for a, b in zip(primerAdapter, comparisonSeq) if a != b)
                            if newDif < numDif:
                                numDif = newDif
                            currPos += 1
                        #only write if the cap does not just look like a primer
                        if numDif > primerAdapterLimit:
                            if cRNA:
                                capLengthscRNA += [insertionLength] 
                            else:
                                capLengthsvRNA += [insertionLength] 
                            outCapped.write(nameCap)
                            outCapped.write(seqCap[fiveEnd:])
                            outCapped.write(toss)
                            outCapped.write(scoreCap[fiveEnd:])
                            outCappedpaired.write(namePair)
                            outCappedpaired.write(seqPair)
                            outCappedpaired.write(toss)
                            outCappedpaired.write(scorePair)

                            statistics['retained_capped'] += 1
                            if cRNA:
                                statistics['capped_cRNA'] += 1
                            else:
                                statistics['capped_vRNA'] += 1
                    elif insertionLength == uncapLimit:
                        uncapLengths += [insertionLength] 
                        outUncapped.write(nameCap)
                        outUncapped.write(seqCap[fiveEnd:])
                        outUncapped.write(toss)
                        outUncapped.write(scoreCap[fiveEnd:])
                        outUncappedpaired.write(namePair)
                        outUncappedpaired.write(seqPair)
                        outUncappedpaired.write(toss)
                        outUncappedpaired.write(scorePair)
                        statistics['retained_uncapped'] += 1
                        if cRNA:
                            statistics['uncapped_cRNA'] += 1
                        else:
                            statistics['uncapped_vRNA'] += 1
            statistics['Total_sequences'] += 1
            nameCap = capRead.readline()
            namePair = pairedRead.readline() 
    pd.DataFrame(statistics.items(), columns=['Category', 'number']).to_csv(outStatsPrefix + '_totalStats.csv')
    pd.DataFrame({'Lengths':capLengthsvRNA}).to_csv(outStatsPrefix + '_cappedLengthsvRNA.csv')  
    pd.DataFrame({'Lengths':uncapLengths}).to_csv( outStatsPrefix + '_uncappedLengths.csv')  



#Similar to above for fairness to plasmid samples, but with the modification that you don't bother with an adapter or capsequences. 
def capSortFastqPlas(*, capRead, pairedRead, file_prefix,  vRNAseq, cRNAseq, toleratedfivePrime, outStatsPrefix):
    #append prefix to existing base file names for outfiles
    outfile = file_prefix + os.path.basename(capRead)
    outPaired = file_prefix + os.path.basename(pairedRead)


    
    #open all files, will write as a stream to each. Throw away reads that do not meet the indicated requirements for now.
    #keep track of numbers with outStats just so we have some kind of QC if we need to alter code or at least a heads-up on what is going on
    with open(capRead, 'r') as capRead, open(pairedRead, 'r') as pairedRead, open(outfile, 'w') as outcap, open(outPaired, 'w') as outpair:
        numPairs = 0
        numCapped = 0
        numUncapped = 0
        #order is name, sequence, throw away line ,and score. For this ignore score (but output it)
        nameCap = capRead.readline()
        namePair = pairedRead.readline()
        statistics = {'Total_sequences':0, 'retained':0,'retained_capped':0, 'vRNA':0, 'cRNA': 0}
        while nameCap != '':
            #make sure read names match. Sometimes stupid things like R1 and R2 get appended, so just grab the chip coordinates
            assert nameCap.split(' ')[0] == namePair.split(' ')[0]
            seqCap = capRead.readline()
            seqPair = pairedRead.readline()
            toss = capRead.readline()
            toss = pairedRead.readline()
            scoreCap = capRead.readline()
            scorePair = pairedRead.readline()
            currPos = 0
            fiveEnd = 0
            numDif = len(vRNAseq)
            fivePrimepos = 0
            while currPos <= len(seqCap) - len(vRNAseq):
                #test v and c go with best match
                comparisonSeq = seqCap[currPos:currPos + len(vRNAseq)]
                newDif  = sum(1 for a, b in zip(vRNAseq, comparisonSeq) if a != b)
                if newDif < numDif:
                    fiveEnd = currPos
                    numDif = newDif
                currPos += 1
            currPos = 0
            #flag for statistics. For distributing in cap vs non actually irrelevant
            cRNA= False
            while currPos <= len(seqCap) - len(cRNAseq):
                #test v and c go with best match
                comparisonSeq = seqCap[currPos:currPos + len(cRNAseq)]
                newDif  = sum(1 for a, b in zip(cRNAseq, comparisonSeq) if a != b)
                if newDif < numDif:
                    fiveEnd = currPos
                    numDif = newDif
                    cRNA = True
                currPos += 1
            #just testing for if meets criteria for five prime
            if numDif <= toleratedfivePrime:
                #just write to file the sequence that matches 5prime end same as above
                outcap.write(nameCap)
                outcap.write(seqCap[fiveEnd:])
                outcap.write(toss)
                outcap.write(scoreCap[fiveEnd:])
                outpair.write(namePair)
                outpair.write(seqPair)
                outpair.write(toss)
                outpair.write(scorePair)                
                statistics['retained'] += 1
                if cRNA:
                    statistics['cRNA'] += 1
                else:
                    statistics['vRNA'] += 1
            statistics['Total_sequences'] += 1
            nameCap = capRead.readline()
            namePair = pairedRead.readline() 
    pd.DataFrame(statistics.items(), columns=['Category', 'number']).to_csv(outStatsPrefix + '_totalStats.csv')




#due to issues with potential phagfe contaimination in an enzyme, used single-end reads for RACE from native samples, so need to modify to simply call from single.
    
def positionalDataSingleEnd(*, sortedBamfile, outfile, qvalCutoff):
    tempFile = sortedBamfile.split('.')[0] + '_temp.sam'
    output = {}
    command  = ' '.join(['samtools view -F 260', sortedBamfile,">", tempFile])
    os.system(command)
    with open(tempFile, 'r') as infile:
        for line in infile:
            currentSeq = {}
            insertions ={}
            read = line.split('\t')

            segment = read[2]
            CIGAR = read[5]
            currentPosSeq = int(read[3])
            score = read[10]
            sequence = read[9]
            CIGARval = []
            curstring = ''
            for character in CIGAR:
                if not character.isalpha():
                    curstring += character
                else:
                    CIGARval.append((int(curstring), character))
                    curstring = ''
                currentPosRead = 1
                #For D, N, advance position in sequence but not read. For I, S, advance position in read but not sequence
                #For M, add in
            for value in CIGARval:
                if value[1] == 'M':
                    #add all mapped positions to a list
                    basesAdded = 0
                    while basesAdded < value[0]:
                        if currentPosSeq not in currentSeq or ord(score[currentPosRead-1]) - 33 > currentSeq[currentPosSeq][1] or currentSeq[currentPosSeq] == 'D':
                            currentSeq[currentPosSeq] = (sequence[currentPosRead-1], ord(score[currentPosRead-1]) - 33)
                        basesAdded += 1
                        currentPosSeq += 1
                        currentPosRead +=1
                elif value[1] == 'D':
                    basesAdded = 0
                    #only count if quality score to either side meets standard. Only called if base not already called at position.
                    if currentPosSeq not in currentSeq and ord(score[currentPosRead-2]) - 33 >= qvalCutoff and ord(score[currentPosRead]) - 33 >= qvalCutoff:
                        while basesAdded < value[0]:
                            currentSeq[currentPosSeq] = ('D', qvalCutoff)
                            basesAdded += 1
                            currentPosSeq += 1
                            #do not advance place in read
                    else:
                        currentPosSeq += value[0]
                elif value[1] == 'N':
                    currentPosSeq += value[0]
                #inserted sequence gets added at the position if bases before and after are a good score
                elif value[1] == 'I':
                    if ord(score[currentPosRead-2]) - 33 >= qvalCutoff and ord(score[currentPosRead + value[0] - 1] ) - 33 >= qvalCutoff:
                        insertions[currentPosSeq] = value[0]
                        currentPosRead += value[0]
                elif value[1] == 'S':
                    currentPosRead += value[0]

            for position in currentSeq:
                if segment not in output:
                    output[segment] = {}
                if position not in output[segment]:
                    output[segment][position] = {'A':0, 'G':0,'C':0,'T':0, 'D':0, 'I':0}
                if currentSeq[position][1] >= qvalCutoff:
                    #for any reasonable qValCutoff none should be N's, but just to be anal retentive
                    base = currentSeq[position][0]
                    if base != 'N':
                        output[segment][position][base] += 1
            for position in insertions:
                output[segment][position]['I'] += 1
    with open(outfile, 'w') as outfile:
        #might as well put it all in order in a tsv
        outfile.write('\t'.join(['segment','position','A','G','C','T', 'D','I']) + '\n')
        for segment in sorted(output.keys()):
            for position in sorted(output[segment]):
                writeline = '\t'.join([segment, str(position)])
                for base in ['A','G','C','T','D','I']:
                    writeline += '\t' + str(output[segment][position][base])
                writeline += '\n'
                outfile.write(writeline)
    arguments = f"rm -f {tempFile}"
