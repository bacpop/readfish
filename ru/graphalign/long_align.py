# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2015, Michigan State University.
# Copyright (C) 2015, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org

import sys
import string
import array
import khmer
import screed

__complementTranslation = str.maketrans('ACTGactg-=Nn', 'TGACtgac-=Nn')

REGIONSIZE = 50


def reverse_complement(s):
    """
    Build reverse complement of 's', alignment-aware.
    """
    # r = array.array('b', [ord(c) for c in s])
    # r.reverse()
    # r = str.join(r, '')
    c = s.translate(__complementTranslation)[::-1]
    # c = str.translate(r, __complementTranslation)

    return c


class GraphAlignment(object):
    def __init__(self, g, r, covs):
        assert len(g) == len(r), (g, len(g), r, len(r))
        self.g = g.upper()
        self.r = r.upper()
        self.covs = covs

    def reverse_complement(self):
        return GraphAlignment(reverse_complement(self.g),
                              reverse_complement(self.r),
                              list(reversed(self.covs)))

    rc = reverse_complement

    def __add__(self, other):
        return GraphAlignment(self.g + other.g, self.r + other.r,
                              self.covs + other.covs)

    def __len__(self):
        return len(self.g)

    def refseqlen(self):
        r = self.r
        return r.count('G') + r.count('C') + r.count('T') + r.count('A')

    def kmer_abundance(self, ct, gi):
        # note: this may differ from alignment.covs for unclear reasons.
        if self.g[gi] in '-=':
            return 0

        start = gi

        n_ch = 0
        while n_ch < ct.ksize() and gi < len(self.g):
            if self.g[gi] == '=':
                gi -= 1
                break
            elif self.g[gi] != '-':
                n_ch += 1
            gi += 1

        end = gi

        if n_ch < ct.ksize():
            gi = start
            gi -= 1
            while n_ch < ct.ksize() and gi >= 0:
                if self.g[gi] == '=':
                    break
                elif self.g[gi] != '-':
                    n_ch + 1
                gi -= 1
        kmer = self.g[start:end].replace('-', '')

        if len(kmer) != ct.ksize():
            return 0

        return ct.get(kmer)

    def __getitem__(self, i):
        if isinstance(i, slice):
            start, stop, step = i.indices(len(self.g))
            assert step == 1
            return GraphAlignment(self.g[start:stop], self.r[start:stop],
                                  self.covs[start:stop])

        return (self.g[i], self.r[i], self.covs[i])

    def compare(self):
        n = len(self.g)
        m = n

        for _, a, b in self.mismatches():
            m -= 1

        return m, n

    def mismatches(self):
        for n, (a, b) in enumerate(zip(self.g, self.r)):
            if (a in '-=' or b in '-=' or a != b):
                yield n, a, b

    def soft_clip(self):
        n = len(self.r)
        m = n

        # move through read forward, determine if start soft-clipped
        for a in self.g:
            if (a in '-='):
                m -= 1
            else:
                break

        # move through read backwards, determine if end soft-clipped
        for a in self.g[::-1]:
            if (a in '-='):
                m -= 1
            else:
                break

        return float(m/n)

    def variants(self):
        for n, (a, b) in enumerate(zip(self.g, self.r)):
            if a == '=':
                continue
            if a != b:
                yield n, a, b

    def __str__(self):
        line1 = []
        line2 = []
        line3 = []

        for (a, b) in zip(self.g, self.r):
            line1.append(a)
            line3.append(b)

            if a in '-=' or b in '-=' or a != b:
                line2.append(' ')
            else:
                line2.append('|')

        return "%s\n%s\n%s\n" % ("".join(line1), "".join(line2), "".join(line3))


def make_gap(r):
    return GraphAlignment('=' * len(r), r, [0] * len(r))


def _index_alignment(galign, freq=100):
    g_to_r = {}
    r_to_g = {}

    si = 0
    for gi, (a, b) in enumerate(zip(galign.g, galign.r)):
        if gi % freq == 0:
            g_to_r[gi] = si
        if si % freq == 0:
            r_to_g[si] = gi

        if b in 'ACGT':
            si += 1

    return g_to_r, r_to_g


class AlignmentIndex(object):
    def __init__(self, galign, freq=100):
        self.g_to_r, self.r_to_g = _index_alignment(galign, int(freq))
        self.galign = galign
        self.freq = int(freq)

    def get_gi(self, ri):
        "Return alignment coordinates cooresponding to reference seq coords."
        rpost = int(ri / self.freq) * self.freq
        gpost = self.r_to_g[rpost]

        diff = ri - rpost
        gi = gpost

        while diff > 0:
            (a, b, _) = self.galign[gi]
            if b in 'ACGT':
                diff -= 1
            gi += 1

        # make sure it's on a valid letter ;)
        while 1:
            (a, b, _) = self.galign[gi]
            if b in 'ACGT':
                break
            gi += 1

        return gi

    def get_ri(self, gi):
        "Return reference sequence coordinates from alignment coordinates."
        gpost = int(gi / self.freq) * self.freq
        ri = self.g_to_r[gpost]

        diff = gi - gpost
        while diff > 0:
            (a, b, c) = self.galign[gpost]
            diff -= 1
            gpost += 1
            if b in 'ACGT':
                ri += 1

        return ri

    def _sanityCheck(self, seq):
        alignment = self.galign

        for gi in range(len(alignment)):
            if alignment[gi][1] in 'ACGT':
                ri = self.get_ri(gi)
                assert alignment[gi][1] == seq[ri]

        for ri in range(len(seq)):
            gi = self.get_gi(ri)
            assert seq[ri] == alignment[gi][1]


def stitch(galign, K):
    """Stitch together multiple alignments, assuming all overlaps are
    off-by-one on the right end."""
    ga = []
    ra = []
    covlist = []

    len_so_far = 0

    n = 0
    for gg in galign[:-1]:
        g = gg.g
        r = gg.r
        covs = gg.covs
        assert len(g) == len(r)

        ga.append(g[:-K + 1])
        ra.append(r[:-K + 1])
        covlist.extend(covs[:-K + 1])

        n += 1

    ga.append(galign[-1].g)
    ra.append(galign[-1].r)
    covlist.extend(galign[-1].covs)

    return GraphAlignment("".join(ga), "".join(ra), covlist)


def align_segment_right(aligner, seq, next_ch=None):
    if len(seq) < 21:  # @CTB
        return 0, make_gap(seq)

    score, g, r, truncated, covs = aligner.align_forward(seq)
    if truncated:
        return 0, make_gap(seq)
    galign = GraphAlignment(g, r, covs)

    # did it fail to align across entire segment?
    if truncated:
        aligned_length = galign.refseqlen()
        unaligned_len = len(seq) - aligned_length

        # if we are given next ch, try aligning backwards from next seed
        if next_ch:
            unaligned_seq = seq[-unaligned_len:]
            sr, ggr = align_segment_left(aligner, unaligned_seq + next_ch)
            galign += ggr[:-1]
            score += sr  # need to adjust score down... @@
        else:
            # if not, just build a gap...
            galign += make_gap(seq[-unaligned_len:])

    assert galign.refseqlen() == len(seq)

    return score, galign


def align_segment_left(aligner, seq):
    seq_rc = reverse_complement(seq)
    score, galign = align_segment_right(aligner, seq_rc)
    return score, galign.rc()


def align_long(ct, aligner, sta):
    K = ct.ksize()

    # first, pick seeds for each chunk
    seeds = []
    for start in range(0, len(sta), REGIONSIZE):
        region = sta[start:start + REGIONSIZE + K - 1]
        if len(region) < K:
            assert len(sta) - start < REGIONSIZE
            break
        seed_pos = find_highest_abund_kmer(ct, region)
        seeds.append(start + seed_pos)

    # if no seeds found, return
    if not seeds:
        return 0

    assert len(seeds) == len(set(seeds))

    # then, break between-seed intervals down into regions, starting at
    # first seed.
    region_coords = []
    for i in range(len(seeds) - 1):
        seed_pos = seeds[i]
        end_seed = seeds[i + 1] - 1
        region_coords.append((seed_pos, end_seed + K))
    # account for situation where last region is too small to align.
    if len(sta) - seeds[-1] > K:
        region_coords.append((seeds[-1], len(sta)))
    else:
        (last_seed, _) = region_coords.pop()
        region_coords.append((last_seed, len(sta)))

    # start building piecewise alignments, anchored by seeds.
    alignments = []
    scores = []

    n = 0
    for (start, end) in region_coords[:-1]:
        score, galign = align_segment_right(aligner, sta[start:end],
                                            next_ch=sta[end])
        scores.append(score)
        alignments.append(galign)

        n += 1

    # deal with end (no possibility of alignment from right)
    (start, end) = region_coords[-1]
    score, galign = align_segment_right(aligner, sta[start:end])

    alignments.append(galign)
    scores.append(score)

    # deal with beginning, too: reverse align from first seed.
    leftend = sta[0:seeds[0] + K]
    score, galign = align_segment_left(aligner, leftend)

    galign = galign[:-1]  # trim off seed k-mer
    alignments.insert(0, galign)
    scores.insert(0, score)

    # stitch all the alignments together
    final = stitch(alignments, K)

    # get % matches in alignment
    clipped = final.soft_clip()

    return clipped


def find_highest_abund_kmer(ct, r):
    pos = 0
    abund = ct.get(r[:ct.ksize()])
    for kstart in range(1, len(r) - ct.ksize() + 1):
        kmer = r[kstart:kstart + ct.ksize()]
        count = ct.get(kmer)

        if count > abund:
            pos = kstart

    return pos