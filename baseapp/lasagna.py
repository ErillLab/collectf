"""
LASAGNA Python implementation

http://biogrid.engr.uconn.edu/lasagna_search/

"""


import gzip
import string
import numpy

complement = string.maketrans('acgtACGT-', 'tgcaTGCA-')
def ReverseComplement(s):
    return s[::-1].translate(complement)

def ComputeCounts(sites, SCOPE, counts_2mer=None, alphabet='acgt-'):
    l = len(sites[0])
    bgCounts = {}
    counts = numpy.array([{} for i in xrange(l)]) # a dict at each position
    scopeStart = 1
    if SCOPE > 0:
        if counts_2mer is None:
            counts_2mer = numpy.array( [[{} for i in xrange(l - scope)]\
                                        for scope in xrange(1, SCOPE + 1)] )
        else:
            scopeStart = len(counts_2mer) + 1
            tmp = numpy.array( [[{} for i in xrange(l - scope)]\
                                for scope in xrange(len(counts_2mer) + 1, SCOPE + 1)] )
            counts_2mer = numpy.concatenate((counts_2mer, tmp))
    nSites = len(sites)
    for i in xrange(nSites): # for each site
        assert len(sites[i]) == l, '\n'.join(sites)
        for k in xrange(l):
            c = sites[i][k].lower()
            if not c in alphabet: continue
            counts[k][c] = counts[k].get(c, 0.0) + 1.0
            bgCounts[c] = bgCounts.get(c, 0.0) + 1.0
        for scope in xrange(scopeStart, SCOPE + 1):
            for k in xrange(l - scope):
                c1 = sites[i][k].lower(); c2 = sites[i][k+scope].lower()
                if not (c1 in alphabet and c2 in alphabet):
                    continue
                p = (c1, c2)
                counts_2mer[scope - 1][k][p] = counts_2mer[scope - 1][k].get(p, 0.0) + 1.0
    return bgCounts, counts, counts_2mer

def ComputeFreqs(sites, SCOPE, alphabet='acgt-', bgCounts=None, counts=None,
                 counts_2mer=None):
    if bgCounts is None or counts is None or (SCOPE > 0 and counts_2mer is None):
        bgCounts, counts, counts_2mer = ComputeCounts(sites, SCOPE, alphabet=alphabet)
    l = len(counts) if counts is not None else len(sites[0])
    aSize = len(alphabet)
    bgFreqs = {'a': 0.0, 'c': 0.0, 'g': 0.0, 't': 0.0, '-': 0.0}
    n = float( sum( bgCounts.values() ) )
    for c in alphabet:
        # Smooth bgFreqs by a pseudo count of 1.0
        bgFreqs[c] = (bgCounts.get(c, 0.0) + 1.0 / float(aSize)) / (n + 1)
    freqs = numpy.array([{'a': 0.0, 'c': 0.0, 'g': 0.0, 't': 0.0, '-': 0.0}\
                         for i in xrange(l)]) # a dict at each position
    for i in xrange(l):
        n = float( sum( counts[i].values() ) )
        for c in alphabet:
            freqs[i][c] = (counts[i].get(c, 0.0) + bgFreqs.get(c, 0.0)) / (n + 1.0)
    freqs_2mer = None
    if SCOPE > 0:
        freqs_2mer = numpy.array( [[{} for i in xrange(l - scope)]\
                                   for scope in xrange(1, SCOPE + 1)] )
    for scope in xrange(1, SCOPE + 1):
        for i in xrange(l - scope):
            n = float( sum( counts_2mer[scope - 1][i].values() ) )
            for c1 in alphabet:
                for c2 in alphabet:
                    p = (c1, c2)
                    freqs_2mer[scope - 1][i][p] = (counts_2mer[scope - 1][i].get(p, 0.0) + bgFreqs[c1] * bgFreqs[c2]) / (n + 1)
    return bgFreqs, freqs, freqs_2mer

def e(s, n):
    return (s - 1.0) / (2.0 * numpy.log(2) * n)

## With small sample correction
def ComputeIC(sites, SCOPE, counts=None, counts_2mer=None):
    alphabet = 'acgt'
    aSize = len(alphabet)
    if counts is None or (SCOPE > 0 and counts_2mer is None):
        bgCounts, counts, counts_2mer = ComputeCounts(sites, SCOPE, alphabet=alphabet)
    l = len(counts) if counts is not None else len(sites[0])
    IC_2mer = None
    IC = numpy.ones(l) * (-numpy.log2(1.0 / float(aSize))) # information content
    for k in xrange(l):
        n = float( sum( [ counts[k].get(c, 0.0) for c in alphabet ]) )
        if n == 0.0:
            IC[k] = 0.0
            continue
        for c in alphabet:
            if not (c in counts[k]) or counts[k][c] == 0.0:
                continue
            freq = counts[k][c] / n
            IC[k] += freq * numpy.log2(freq)
        IC[k] = max(0.0, IC[k] - e(aSize, n))
    if SCOPE > 0:
        IC_2mer = [numpy.ones(l-scope) * (-numpy.log2(1.0 / float(aSize ** 2)))\
                   for scope in xrange(1, SCOPE + 1)]    # information content for pairs of positions
    for scope in xrange(1, SCOPE + 1):
        for k in xrange(l-scope):
            n = float( sum( [ counts_2mer[scope - 1][k].get((c1, c2), 0.0)\
                              for c1 in alphabet for c2 in alphabet] ) )
            if n == 0.0:
                IC_2mer[scope - 1][k] = 0.0
                continue
            for c1 in alphabet:
                for c2 in alphabet:
                    p = (c1, c2)
                    if not (p in counts_2mer[scope - 1][k]) or counts_2mer[scope - 1][k][p] == 0.0:
                        continue
                    freq = counts_2mer[scope - 1][k][p] / n
                    IC_2mer[scope - 1][k] += freq * numpy.log2(freq)
            IC_2mer[scope - 1][k] = max(0.0, IC_2mer[scope - 1][k] - e(aSize**2, n))
    return IC, IC_2mer

def ComputeCoverage(counts, nSites, alphabet='acgt'):
    coverages = []
    for cnt in counts:
        n = float( sum([cnt.get(c, 0) for c in alphabet]) )
        coverages.append(n / nSites)
    return coverages

def SuggestSize(IC, coverages, ICThres=0.0, covThres=0.4):
    maxL = len(IC)
    if isinstance(ICThres, str) and ICThres.startswith('mean'):
        tmp = ICThres[4:]
        if tmp.isdigit(): maxL = int(tmp)
        ICThres = numpy.mean(IC)
    if covThres == "mean":
        covThres = numpy.mean(coverages)
    l = len(IC)
    trimmedL = l
    cntL = cntR = 0
    for i in xrange(l):
        if IC[i] > ICThres and coverages[i] > covThres:
            break
        trimmedL -= 1
        cntL += 1
    for j in xrange(l - 1, i, -1):
        if IC[j] > ICThres and coverages[j] > covThres:
            break
        trimmedL -= 1
        cntR += 1
    j = l-cntR-1
    while trimmedL > 0 and (trimmedL > maxL or IC[cntL] == 0 or IC[j] == 0):
        if IC[cntL] < IC[j]:
            cntL += 1
        else:
            cntR += 1; j -= 1
        trimmedL -= 1
    return trimmedL, cntL, cntR

def TrimPSSM(model, cntL, cntR):
    model.l = model.l - cntL - cntR
    start = cntL
    end = cntL + model.l
    model.scores = model.scores[start:end]
    if hasattr(model, 'scores_2mer') and model.scores_2mer != None:
        for scope in xrange(model.SCOPE, 0, -1):
            end = cntL + model.l - scope
            if scope + 1 > model.l:
                del model.scores_2mer[scope - 1]
                continue
            model.scores_2mer[scope - 1] = model.scores_2mer[scope - 1][start:end]
        model.SCOPE = min(model.SCOPE, model.l - 1)
    return model

class modelPSSM:
    pass

def PSSM(sites, SCOPE, alphabet='acgt-', bgCounts=None, counts=None,
         counts_2mer=None, trim=False, ICThres=0.0, covThres=0.4,
         A=True):
    aSize = len(alphabet)
    l = len(sites[0])
    assert SCOPE >=0, 'SCOPE must be non-negative!'
    SCOPE = min(SCOPE, l-1)
    score_2mer = None
    if SCOPE > 0:
        scores_2mer = [numpy.array( [{} for i in xrange(l - scope)] )\
                       for scope in xrange(1, SCOPE + 1)]
    if bgCounts is None or counts is None:
        bgCounts, counts, counts_2mer = ComputeCounts(sites, SCOPE,
                                                      alphabet=alphabet)
    bgFreqs, freqs, freqs_2mer = ComputeFreqs(sites, SCOPE, alphabet=alphabet,
                                              bgCounts=bgCounts, counts=counts,
                                              counts_2mer=counts_2mer)
    scores = numpy.array( [{} for i in xrange(l)] )
    for i in xrange(l):
        for c in alphabet:
            scores[i][c] = numpy.log2(freqs[i].get(c, 0.0) / bgFreqs.get(c, 0.0))
        if A: scores[i]['-'] = min(scores[i].values())
    scores_2mer = [numpy.array( [{} for i in xrange(l - scope)] )\
                   for scope in xrange(1, SCOPE + 1)] if SCOPE > 0 else None
    for scope in xrange(1, SCOPE + 1):
        for i in xrange(l - scope):
            for c1 in alphabet:
                for c2 in alphabet:
                    p = (c1, c2)
                    scores_2mer[scope - 1][i][p] = numpy.log2(freqs_2mer[scope - 1][i][p] / (bgFreqs[c1] * bgFreqs[c2]))
            if A:
                m = min(scores_2mer[scope - 1][i].values())
                scores_2mer[scope - 1][i][('-', '-')] = m
                for c in alphabet:
                    scores_2mer[scope - 1][i][('-', c)] = m
                    scores_2mer[scope - 1][i][(c, '-')] = m
    model = modelPSSM()
    model.l = l
    model.scores = scores
    model.scores_2mer = scores_2mer
    model.SCOPE = SCOPE
    model.alphabet = alphabet
    model.cntL = model.cntR = 0
    if trim:
        IC, IC_2mer = ComputeIC(sites, 0, counts=counts, counts_2mer=counts_2mer)
        if IC != None:
            coverages = ComputeCoverage(counts, len(sites), alphabet='acgt')
            model.sugL, model.cntL, model.cntR = SuggestSize(IC, coverages,
                                                             ICThres=ICThres,
                                                             covThres=covThres)
            if model.sugL == 0:
                model.cntL = model.cntR = 0
                model.sugL = len(sites[0])
        else:
            (model.sugL, model.cntL, model.cntR) = (model.l, 0, 0)
        model = TrimPSSM(model, model.cntL, model.cntR)
        model.trimmedIC = IC[model.cntL:(model.cntL+model.sugL)]
    return model

def ScoreByPSSM(seq, model):
    l = len(seq)
    score = 0.0
    for i in xrange(l):
        if not (seq[i] in model.scores[i]):
            s = min( model.scores[i].values() )
            model.scores[i][seq[i]] = s
        else:
            s = model.scores[i][seq[i]]
        score += s
    for scope in xrange(1, model.SCOPE + 1):
        for i in xrange(l - scope):
            p = (seq[i], seq[i + scope])
            if not (p in model.scores_2mer[scope - 1][i]):
                s = min( model.scores_2mer[scope - 1][i].values() )
                model.scores_2mer[scope - 1][i][p] = s
            else:
                s = model.scores_2mer[scope - 1][i][p]
            score += s
    return score

def SlideScoreByPSSM(seq, model):
    seq = seq.lower()
    seqL = len(seq)
    bestScore = None
    jStart = 0 # first position of the PSSM
    for iStart in xrange(seqL - 1, -1, -1): # last letter of the sequence
        subseqL = min(seqL - iStart, model.l)
        subseq = seq[iStart:(iStart + subseqL)]
        subseq = subseq + '-' * (model.l - subseqL)
        score = ScoreByPSSM(subseq, model)
        if bestScore == None or score > bestScore:
            bestScore = score
            bestAlign = (iStart, jStart)
    for jStart in xrange(1, model.l):
        subseqL = min(model.l - jStart, seqL)
        subseq = seq[:subseqL]
        subseq = '-' * jStart + subseq + '-' * (model.l - subseqL - jStart)
        score = ScoreByPSSM(subseq, model)
        if bestScore == None or score > bestScore:
            bestScore = score
            bestAlign = (iStart, jStart)
    return bestScore, bestAlign

def NextSite(model, sites, lens, indices, maxIt=5):
    bestScore = None
    k = 0
    while k < maxIt and\
          (k < len(indices) and lens[indices[k]] == lens[indices[0]]):
        i = indices[k]
        # Score sites[i] with model
        for site, strand in [(sites[i], '+'),
                             (ReverseComplement(sites[i]), '-')]:
            score, (jSeq, jPSSM) = SlideScoreByPSSM(site, model)
            if bestScore == None or score > bestScore:
                bestScore = score
                bestK = k
                bestJSeq = jSeq
                bestJPSSM = jPSSM
                bestSite = site
                bestStrand = strand
        k += 1
    return bestSite, bestK, bestStrand, bestJSeq, bestJPSSM, k

def NMatches(s1, s2, i1, i2):
    l1 = len(s1); l2 = len(s2)
    s1 = s1.lower(); s2 = s2.lower()
    cnt = 0
    while i1 < l1 and i2 < l2:
        if s1[i1] == s2[i2]:
            cnt += 1
        i1 += 1; i2 += 1
    return cnt

def AddGaps(aligned, jSeq, jPSSM,
            bgCounts=None, counts=None, counts_2mer=None,
            SCOPE=0, alphabet='acgt-'):
    nAligned = len(aligned)
    PSSML = len(aligned[0])
    seqL = len(aligned[-1])
    leftDiff = jPSSM - jSeq
    rightDiff = (PSSML - jPSSM) - (seqL - jSeq)
    leftGaps = '-' * abs(leftDiff)
    rightGaps = '-' * abs(rightDiff)
    if leftDiff < 0 or rightDiff < 0:
        leftGaps = '-' * abs(min(0, leftDiff))
        rightGaps = '-' * abs(min(0, rightDiff))
        for i in xrange(len(aligned) - 1):
            aligned[i] = leftGaps + aligned[i] + rightGaps
        ## Update counts
        nLeft = len(leftGaps)
        nRight = len(rightGaps)
        if isinstance(bgCounts, dict):
            bgCounts['-'] = bgCounts.get('-', 0) + (nLeft + nRight) * (nAligned - 1)
        if isinstance(counts, numpy.ndarray):
            counts = numpy.concatenate((numpy.array([{'-': float(nAligned - 1)} for i in xrange(nLeft)]), counts,
                                        numpy.array([{'-': float(nAligned - 1)} for i in xrange(nRight)])))
        l = len(aligned[0])
        if counts_2mer != None and SCOPE > 0:
            if type(counts_2mer) == numpy.ndarray:
                counts_2mer = counts_2mer.tolist()
            for scope in xrange(1, SCOPE + 1):
                if len(counts_2mer[scope - 1]) == 0:
                    counts_2mer[scope - 1] = [{} for k in xrange(l-scope)]
                else:
                    counts_2mer[scope - 1] = [{} for k in xrange(nLeft)] + counts_2mer[scope - 1] + [{} for k in xrange(nRight)]
            for i in xrange(nAligned - 1):
                for scope in xrange(1, SCOPE + 1):
                    for k in range(min(nLeft, l - scope)) + range(max(nLeft, nLeft + PSSML - scope), l-scope):
                        c1 = aligned[i][k].lower(); c2 = aligned[i][k+scope].lower()
                        if not (c1 in alphabet and c2 in alphabet):
                            continue
                        p = (c1, c2)
                        counts_2mer[scope - 1][k][p] = counts_2mer[scope - 1][k].get(p, 0.0) + 1.0
    if leftDiff > 0 or rightDiff > 0:
        leftGaps = '-' * max(0, leftDiff)
        rightGaps = '-' * max(0, rightDiff)
        aligned[-1] = leftGaps + aligned[-1] + rightGaps
    if isinstance(bgCounts, dict):
        assert type(counts) == numpy.ndarray, 'type(counts) != numpy.ndarray'
        l = len(aligned[-1])
        for k in xrange(l):
            c = aligned[-1][k].lower()
            if not (c in alphabet):
                continue
            counts[k][c] = counts[k].get(c, 0.0) + 1.0
            bgCounts[c] = bgCounts.get(c, 0.0) + 1.0
        for scope in xrange(1, SCOPE + 1):
            for k in xrange(l - scope):
                c1 = aligned[-1][k].lower(); c2 = aligned[-1][k+scope].lower()
                if not (c1 in alphabet and c2 in alphabet):
                    continue
                p = (c1, c2)
                counts_2mer[scope - 1][k][p] = counts_2mer[scope - 1][k].get(p, 0.0) + 1.0
    if not (bgCounts == None and counts == None and counts_2mer == None):
        return bgCounts, counts, counts_2mer

def FormAlignment(info):
    minOffset = None
    for (i, strand, offset, site) in info:
        minOffset = offset if minOffset == None else min(minOffset, offset)
    info_new = []
    maxL = 0
    for (i, strand, offset, site) in info:
        offset -= minOffset
        l = offset + len(site)
        maxL = max(l, maxL)
        info_new.append( (i, strand, offset, l, site) )
    aligned = []
    idxAligned = []
    strands = []
    for (i, strand, offset, l, site) in info_new:
        aligned.append('-' * offset + site + '-' * (maxL - l) )
        idxAligned.append(i)
        strands.append( strand )
    return aligned, idxAligned, strands        

def LASAGNA(sites, SCOPE, seedIdx=-1, trim=False, ICThres=0.0, covThres=0.4):
    sizeThres = 1
    nSites = len(sites)
    lens = numpy.array([len(site) for site in sites])
    indices = lens.argsort().tolist()
    aligned = []; idxAligned = []; strands = []
    if seedIdx != -1:
        aligned.append(sites[seedIdx])
        idxAligned.append( seedIdx )
        strands.append('+')
        indices.remove(seedIdx)
    switched = False
    while indices:
        if len(aligned) == 1:
            bgCounts, counts, counts_2mer = ComputeCounts(aligned, SCOPE,
                                                          alphabet='acgt-')
        i = indices[0]
        if not aligned:
            aligned.append(sites[i])
            idxAligned.append(i)
            strands.append('+')
            del indices[0]
        else:
            model = PSSM(aligned, SCOPE, bgCounts=bgCounts, counts=counts,
                         counts_2mer=counts_2mer, trim=trim and len(aligned) > 2,
                         ICThres=ICThres, covThres=covThres)
            bestSite, bestK, bestStrand, bestJSeq, bestJPSSM, k = NextSite(model, sites, lens, indices)
            if len(aligned) == 1:
                if (not switched) and (NMatches(aligned[0], bestSite, bestJPSSM, bestJSeq) == 0):
                    if k < len(indices):
                        indices = [indices[k]] + indices[:k] + indices[(k+1):]
                        switched = True
                        continue
            aligned.append(bestSite)
            idxAligned.append( indices[bestK] )
            strands.append(bestStrand)
            bgCounts, counts, counts_2mer = AddGaps(aligned, bestJSeq, bestJPSSM + model.cntL, SCOPE = SCOPE,
                                                    bgCounts = bgCounts, counts = counts, counts_2mer = counts_2mer)
            del indices[bestK]
    indices = numpy.argsort(idxAligned)
    aligned = numpy.array(aligned)[indices].tolist()
    idxAligned = numpy.array(idxAligned)[indices].tolist()
    strands = numpy.array(strands)[indices].tolist()
    return aligned, idxAligned, strands

def LASAGNA_ChIP(sites, SCOPE, seedIdx=None, maxIt1=6, maxIt2=50,
                 ICThres='mean15', covThres='mean'):
    if seedIdx is None:
        lens = numpy.array([len(site) for site in sites])
        i = numpy.argmin(lens)
        indices = [i] + range(i) + range(i+1, len(sites))
        best = None, None, None, 0
        for seedIdx in indices[:maxIt1]:
            attempt = LASAGNA_ChIP(sites, SCOPE, seedIdx=seedIdx,
                                   ICThres=ICThres, covThres=covThres)
            if attempt[-1] > best[-1]:
                best = attempt
        aligned, idxAligned, strands, sumIC = best
        return aligned, idxAligned, strands
    else:
        aligned, idxAligned, strands = LASAGNA(sites, SCOPE, seedIdx=seedIdx,
                                               trim=True, ICThres=ICThres,
                                               covThres=covThres)
        nAligned = len(aligned)
        it = 0; noChange = 0; oldICSum = 0
        while it < maxIt2:
            model = PSSM(aligned, SCOPE, trim=True, ICThres=ICThres,
                         covThres=covThres)
            ICSum = sum(model.trimmedIC)
            if ICSum == oldICSum:
                noChange += 1
            else:
                oldICSum = ICSum; noChange = 0
            if noChange == 2: break
            alignment = []
            for i in xrange(nAligned):
                j = idxAligned[i]
                site = sites[j]
                strand = '+'
                score, (jSeq, jPSSM) = SlideScoreByPSSM(site, model)
                site_rc = ReverseComplement(site)
                score_rc, (jSeq_rc, jPSSM_rc) = SlideScoreByPSSM(site_rc, model)
                if score_rc > score:
                    (jSeq, jPSSM) = (jSeq_rc, jPSSM_rc)
                    site = site_rc
                    strand = '-'
                offset = jPSSM - jSeq
                alignment.append( (j, strand, offset, site) )
            aligned, idxAligned, strands = FormAlignment(alignment)
            it += 1
        return aligned, idxAligned, strands, ICSum
        
def TrimAlignment(aligned, ICThres=0.0, covThres=0.4):
    nAligned = len(aligned)
    counts = ComputeCounts(aligned, 0, alphabet = 'acgt')[1]
    IC = ComputeIC(aligned, 0, counts=counts, counts_2mer=None)[0]
    coverages = ComputeCoverage(counts, nAligned, alphabet = 'acgt')
    sugL, cntL, cntR = SuggestSize(IC, coverages, ICThres=ICThres,
                                   covThres=covThres)
    trimmed = [ aSite[cntL:(cntL+sugL)] for aSite in aligned]
    return trimmed

def ReadSeq(path):
    seqs = []; headers = []
    inp = gzip.open(path, 'r') if path.endswith('.gz') else\
          open(path, 'r')
    tmp = inp.read().strip().strip('>').split('>')
    inp.close()
    for seq in tmp:
        lines = map(string.strip, seq.split('\n'))
        headers.append(lines[0])
        seqs.append(''.join(lines[1:]))
    return headers, seqs

def Align(params):
    headers, seqs = ReadSeq(params.fastaPath)
    if params.chip:
        aligned, idxAligned, strands = LASAGNA_ChIP(seqs, params.SCOPE)
        if params.trim:
            aligned = TrimAlignment(aligned, ICThres='mean15', covThres='mean')
    else:
        aligned, idxAligned, strands = LASAGNA(seqs, params.SCOPE)
        if params.trim:
            aligned = TrimAlignment(aligned)
    for seq, i, strand in zip(aligned, idxAligned, strands):
        print '>%s %s\n%s' % (headers[i], strand, seq)




