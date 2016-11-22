import subprocess as sp
import os


class Genome:
    def __init__(self, fasta_path=''):
        self.fasta_path = fasta_path
        self.chroms_dict = {}

    def get_sequence(self, chrom, start, stop, strand):
        """finds the sequence of a 1-indexed genomic interval '1', 387941, 388099, '+' """
        assert self.fasta_path
        pos_string = chrom + ':' + str(start) + '-' + str(stop)
        samtools_cmm = ['samtools', 'faidx', self.fasta_path, pos_string]
        samtools_out = sp.check_output(samtools_cmm).decode()
        seq = fasta2seq(samtools_out)
        if strand == '-':
            seq = comp_rev(seq)
        return seq.upper()

    @property
    def chromosomes(self):
        for chrom in self.chroms_dict.values():
            yield chrom

    @property
    def genes(self):
        for chrom in self.chromosomes:
            for gene in chrom.genes:
                yield gene

    @property
    def transcripts(self):
        for chrom in self.chromosomes:
            for trans in chrom.transcripts:
                yield trans


class Chromosome:
    def __init__(self, name, genome=None):
        self.name = name
        self.genes_dict, self.ordered_transs = {}, []
        if genome:
            self.genome = genome
            genome.chroms_dict[name] = self
        self.coverage = {'+': {}, '-': {}}

    def get_sequence(self, start, stop, strand):
        return self.genome.get_sequence(self.name, start, stop, strand)

    def ordered_transcripts_tss(self):
        for trans in self.transcripts:
            insert_pos = len(self.ordered_transs)
            for upstream_trans in self.ordered_transs[::-1]:
                if upstream_trans.tss < trans.tss:
                    self.ordered_transs[insert_pos:insert_pos] = [trans]
                    break
                insert_pos -= 1
            else:
                self.ordered_transs[0:0] = [trans]
        return self.ordered_transs

    def get_antisense_transcripts(self):
        antisenses, trans_i = [], 0
        transcripts = list(self.transcripts)
        for trans in transcripts:
            antisenses.append([trans, self._get_closest_antisense(trans, trans_i, transcripts)])
            trans_i += 1
        return antisenses

    def _get_closest_antisense(self, trans, trans_i, trans_chrom):
        strand = trans.strand
        range_downstream = range(trans_i + 1, len(trans_chrom))
        range_upstream = range(trans_i - 1, -1, -1)
        as_downstream = self._closest_in_range(range_downstream, trans_chrom, not_gene=trans.gene.id, not_strand=strand)
        as_upstream = self._closest_in_range(range_upstream, trans_chrom, not_gene=trans.gene.id, not_strand=strand)
        if abs(trans.tss - as_downstream.tss) < abs(trans.tss - as_upstream.tss):
            return as_downstream
        return as_upstream

    @staticmethod
    def _closest_in_range(trans_range, trans_chrom, not_gene, strand=None, not_strand=None):
        """returns the closest transcript from trans_chrom in the specified range and strand"""
        for trans_j in trans_range:
            if trans_chrom[trans_j] != strand and trans_chrom[trans_j].gene.id != not_gene:
                return trans_chrom[trans_j]

    @property
    def genes(self):
        for gene in self.genes_dict.values():
            yield gene

    @property
    def transcripts(self):
        for gene in self.genes:
            for trans in gene.transcripts:
                yield trans


class Gene:
    def __init__(self, id, chromosome, name='', strand=''):
        self.id = id
        self.name = name
        self.chromosome = chromosome
        chromosome.genes_dict[id] = self
        self.strand = strand
        self.trans_dict, self.exons_dict = {}, {}

    def __repr__(self):
        return self.id

    def __len__(self):
        return max([len(trans) for trans in self.transcripts] + [0])

    @property
    def exons(self):
        for exon in self.exons_dict.values():
            yield exon

    @property
    def transcripts(self):
        for trans in self.trans_dict.values():
            yield trans

    def add_exon(self, exon, start, stop):
        self.exons_dict[(start, stop)] = exon

    def get_sequence(self, start, stop):
        return self.chromosome.get_sequence(start, stop, self.strand)

    def includes(self, pos):
        for trans in self.transcripts:
            if trans.includes(pos):
                return True
        return False

    def most_upstream_trans(self):
        transs = list(self.transcripts)
        if not transs[0]:
            return
        str_mul = int(self.strand + '1')
        ups_trans = transs[0]
        ups_start = str_mul * transs[0].tss
        for trans in transs[1:]:
            trans_tss = str_mul * trans.tss
            if trans_tss < ups_start:
                ups_start = trans_tss
                ups_trans = trans
        return ups_trans


class Transcript:
    def __init__(self, id, gene=None, cds_start=None, cds_stop=None, tss=None, biotype=None):
        self.id = id
        self.gene = gene
        self.biotype = biotype
        self.gene.trans_dict[id] = self
        self.exon_dict, self.splice_sites = {}, []
        self.cds_start, self.cds_stop = fix_order(cds_start, cds_stop, self.strand)

    def __repr__(self):
        return self.id

    def __len__(self):
        return sum([len(exon) for exon in self.exon_dict.values()])

    @property
    def exons(self):
        for n in sorted(self.exon_dict):
            yield self.exon_dict[n]

    @property
    def n_splice_sites(self):
        return len(self.splice_sites)

    @property
    def strand(self):
        return self.gene.strand

    @property
    def chromosome(self):
        return self.gene.chromosome

    @property
    def three_utr_len(self):
        if not self.cds_stop:
            return 0
        return self.distance_to_end(self.cds_stop)

    @property
    def five_utr_len(self):
        if not self.cds_start:
            return 0
        return self.relative_position(self.cds_start)

    @property
    def cds_len(self):
        if not (self.cds_start and self.cds_stop):
            return 0
        return self.distance_to_end(self.cds_start) - self.three_utr_len

    @property
    def tss(self):
        if 1 in self.exon_dict:
            return self.exon_dict[1].start

    @property
    def tstop(self):
        if 1 in self.exon_dict:
            return self.exon_dict[len(self.exon_dict)].stop

    def in_cds(self, pos):
        if self.cds_start <= pos <= self.cds_stop or self.cds_start >= pos >= self.cds_stop:
            return True

    def get_sequence(self, start=None, stop=None):
        """returns the transcript sequence between start and stop, or the whole"""
        whole_seq = ''.join([exon.get_sequence() for exon in self.exons])
        rel_start = self.relative_position(start) if start else 0
        rel_stop = self.relative_position(stop) if stop else len(whole_seq)
        return whole_seq[rel_start:rel_stop]

    def get_exons_from_splice_site(self, start, stop):
        query = '_'.join([str(x) for x in fix_order(start, stop, self.strand)])
        for site_n in range(self.n_splice_sites):
            if self.splice_sites[site_n] == query:
                return self.exon_dict[site_n], self.exon_dict[site_n + 1]

    def add_exon(self, exon, number):
        self.exon_dict[number] = exon

    def add_cds(self, start, stop=None):
        if stop:
            start, stop = fix_order(start, stop, self.strand)
            self.cds_stop = fix_order(self.cds_stop, stop, self.strand)[1]
        self.cds_start = fix_order(self.cds_start, start, self.strand)[0]

    def add_splice_sites(self):
        # deprecated!
        # not executed by Gtf reader
        exons = list(self.exons)
        for exon_n in range(1, len(exons)):
            self.splice_sites.append(str(exons[exon_n - 1].stop) + '_' + str(exons[exon_n].start))

    def ups_splice_intervals(self, up_distance=200, dw_distance=200, skip_short=True):
        for exon in list(self.exons)[:-1]:
            start = self.abs_pos_upstream(exon.stop, up_distance)
            if not start:
                continue
            intervals = list(self.intervals(start, exon.stop))
            stop = move_pos(exon.stop, dw_distance, self.strand)
            if skip_short and self.includes(stop):
                continue
            intervals[-1] = (intervals[-1][0], stop)
            yield intervals

    def dws_splice_intervals(self, up_distance=200, dw_distance=200, skip_short=True):
        exons = list(self.exons)
        for exon in exons[1:]:
            stop = self.abs_pos_downstream(exon.start, dw_distance)
            if not stop:
                continue
            intervals = list(self.intervals(exon.start, stop))
            start = move_pos(exon.start, - up_distance, self.strand)
            if skip_short and self.includes(start):
                continue
            intervals[0] = (start, intervals[0][1])
            yield intervals

    def relative_position(self, pos):
        distance = 0
        for exon in self.exons:
            if exon.includes(pos):
                return distance + exon.relative_position(pos)
            distance += len(exon)
        return distance

    def distance_to_end(self, pos):
        return len(self) - self.relative_position(pos)

    def intervals(self, start, stop):
        """given chrom positions start and stop in transcript, returns the exonic intervals"""
        start, stop = fix_order(start, stop, self.strand)
        after_start = False
        for exon in self.exons:
            if after_start:
                if exon.includes(stop):
                    yield exon.start, stop
                    break
                yield exon.start, exon.stop
            elif exon.includes(start):
                after_start = True
                if exon.includes(stop):
                    yield start, stop
                    break
                else:
                    yield start, exon.stop

    @staticmethod
    def exon_n_with(pos, exons):
        """returns the number of the exon that has pos"""
        for exon_n in range(len(exons)):
            if exons[exon_n].includes(pos):
                return exon_n

    def exon_with(self, pos):
        """returns the exon that has pos"""
        for exon in self.exons:
            if exon.includes(pos):
                return exon

    def includes(self, pos):
        if self.exon_with(pos):
            return True
        return False

    def abs_pos_downstream(self, pos, length):
        """
        walks the transcript downstream, from pos for length
        returns the end position. None if span of transcript is exceeded
        """
        exons = list(self.exons)
        start_n = self.exon_n_with(pos, exons)
        path = abs(pos - exons[start_n].stop) + 1
        if path > length:
            return move_pos(pos, length, self.strand)
        length -= path
        for exon in exons[start_n + 1:]:
            if len(exon) > length:
                return move_pos(exon.start, length, self.strand)
            length -= len(exon)

    def abs_pos_upstream(self, pos, length):
        """
        walks the transcript upstream, from pos for length
        returns the end position. None if span of transcript is exceeded
        """
        exons = list(self.exons)
        start_n = self.exon_n_with(pos, exons)
        path = abs(pos - exons[start_n].start) + 1
        if path > length:
            return move_pos(pos, -length, self.strand)
        length -= path
        if start_n == 0:
            return
        for exon in exons[start_n - 1::-1]:
            if len(exon) > length:
                return move_pos(exon.stop, -length, self.strand)
            length -= len(exon)

    def get_cds_relative_starts(self, bam, min_qual=40, read_len=[]):
        """
        saves the relative position of read starts sites overlapping its cds
        :param bam: a bam object (from reader module)
        :param min_qual: default TopHat: only uniquely mapped reads
        :param read_len: only consider reads whose length is in this list
        :return: a dict with 1-based starts as key and number of reads as value
        """
        start_stop = fix_order(self.cds_start, self.cds_stop, '+')
        starts_dict = bam.get_read_starts(self.chromosome.name, start_stop[0], start_stop[1], self.strand, min_qual, read_len)
        rel_starts = {}
        for start, n_reads in starts_dict.items():
            start += 1
            if not self.in_cds(start):
                continue
            rel_pos = self.distance_from_cds_start(start)
            rel_starts[rel_pos] = n_reads
        return rel_starts

    def distance_from_cds_start(self, pos):
        if not self.cds_start:
            return
        distance = 0
        for exon in self.exons:
            if exon.includes(pos):
                if exon.includes(self.cds_start):
                    return abs(self.cds_start - pos)
                return distance + exon.relative_position(pos)
            elif exon.includes(self.cds_start):
                distance += abs(self.cds_start - exon.stop) + 1
            else:
                distance += len(exon)

    def set_cds_from_start(self, cds_start):
        genomic_cds_start, seq = 0, ''
        for exon in self.exons:
            if genomic_cds_start:
                seq += exon.get_sequence()
            elif exon.cds_start:
                genomic_cds_start = exon.cds_start
                seq = exon.get_seq_to_stop(genomic_cds_start)
        if seq:
            seq_obj = Sequence(seq)
            cds_rel_stop = seq_obj.get_orf_from(0)[1]
            if cds_rel_stop:
                self.cds_start = genomic_cds_start
                self.cds_stop = self.abs_pos_downstream(genomic_cds_start, cds_rel_stop)


class Exon:
    def __init__(self, id, num, transcript, start, stop, cds_start=None):
        self.id = id
        self.transcripts = [transcript]
        self.genomic_start = min(start, stop)
        self.genomic_stop = max(start, stop)
        self.start, self.stop = fix_order(start, stop, self.strand)
        self.cds_start = cds_start
        transcript.add_exon(self, num)

    def __repr__(self):
        return self.id

    def __len__(self):
        return self.genomic_stop - self.genomic_start + 1

    @property
    def gene(self):
        return self.transcripts[0].gene

    @property
    def strand(self):
        return self.transcripts[0].strand

    @property
    def chromosome(self):
        return self.transcripts[0].chromosome

    def coding_transcripts(self):
        for trans in self.transcripts:
            if trans.cds_stop and \
                    before_strand(self.start, trans.cds_stop, self.strand) and \
                    before_strand(trans.cds_start, self.stop, self.strand):
                yield trans

    def add_transcript(self, trans, num):
        self.transcripts.append(trans)
        trans.add_exon(self, num)

    def get_sequence(self):
        """returns the whole exon sequence"""
        return self.gene.get_sequence(self.genomic_start, self.genomic_stop)

    def get_seq_from_start(self, stop):
        """
        returns the exon sequence from start to stop,
        the whole sequence if stop is after the exon
        or '' if stop is before
        """
        stop = fix_order(self.stop, stop, self.strand)[0]
        if not self.includes(stop):
            return ''
        return self.gene.get_sequence(sorted((self.start, stop)))

    def get_seq_to_stop(self, start):
        """
        returns the exon sequence from start to stop,
        the whole sequence if start is before exon
        or '' if start is after
        """
        start = fix_order(self.start, start, self.strand)[1]
        if not self.includes(start):
            return ''
        return self.gene.get_sequence(sorted((start, self.stop)))

    def relative_position(self, pos):
        return abs(self.start - pos)

    def includes(self, pos):
        if self.genomic_start <= pos <= self.genomic_stop:
            return True
        return False

    def gtf(self):
        transs = '_'.join([x.id + '.' + self.id for x in self.transcripts])
        return '\t'.join([self.chromosome.name, 't', 'exon', str(self.genomic_start), str(self.genomic_stop), '.',
                          self.strand, 'gene_id "' + transs + '";'])

    def gtf_coding(self):
        coding_transcripts = self.coding_transcripts()
        if coding_transcripts:
            transs = '_'.join([x.id + '.' + self.id for x in coding_transcripts])
            return '\t'.join([self.chromosome.name, 't', 'exon', str(self.genomic_start), str(self.genomic_stop), '.',
                              self.strand, 'gene_id "' + transs + '";'])


class GenomicOverlay:
    def __init__(self, chromosome, n_quantiles=50):
        """for the moment it deals only with continuous intervals"""
        self.chromosome = chromosome
        self.n_quantiles = n_quantiles
        self.quantiles = [0] * n_quantiles
        self.points = {'+': {}, '-': {}}
        self.layers = []

    def add_bedgraph(self, bg):
        """adds a score (eg: coverage) at pos, strand"""
        points = self.points[bg.strand]
        for base, score in bg.iter():
            if base in points:
                points[base].score = score

    def get_point(self, pos, strand):
        """returns Point from self.points_dict, or a new Point if not present"""
        if pos in self.points[strand]:
            return self.points[strand][pos]
        point = GenomicPoint(pos, strand, 0)
        self.points[strand][pos] = point
        return point


class GenomicLayer:
    def __init__(self, start, stop, strand, overlay=None):
        self.start, self.stop = start, stop
        self.strand = strand
        if overlay:
            self.overlay = overlay
            self.overlay.layers.append(self)
            self.quantile_width = len(self) / float(overlay.n_quantiles)
            self.points = self._set_points()

    def __len__(self):
        return abs(self.stop - self.start)

    def _set_points(self):
        """sets the points, based on the overlay"""
        for pos in range(self.start, self.stop):
            n_quantile = self._position_to_quantile(pos)
            yield self.overlay.get_point(pos, self.strand, n_quantile)

    def normalize_point_score(self, divisor=None):
        if not divisor:
            divisor = max([x.score for x in self.points] + [1])
        divisor = float(divisor)
        for point in self.points:
            point.score /= divisor

    def _position_to_quantile(self, pos):
        """determines the bin number, given the position (based on absolute dis from start)"""
        return int(abs(pos - self.start) / self.quantile_width)


class GenomicPoint:
    def __init__(self, pos, score, n_quantile):
        self.pos = pos
        self.score = score
        self.n_quantile = n_quantile


class Sequence:
    def __init__(self, seq):
        self.seq = seq
        self.len = len(seq)

    def __len__(self):
        return self.len

    def get_orf_from(self, start):
        """finds the closest stop from an ATG, or the end of the seq
        :param start: ATG in the sequence
        :return: start, stop
        """
        assert self.is_start(start)
        stop = self.next_stop(start + 3)
        if stop:
            return start, stop
        return start, len(self) - 1

    def get_orfs(self):
        orfs = []
        for pos in range(len(self) - 2):
            if self.is_start(pos):
                stop = self.next_stop(pos + 3)
                if stop:
                    orfs.append([pos, stop])
        return orfs

    def previous_start(self, pos):
        while pos >= 0:
            if self.is_start(pos):
                return pos
            pos -= 3

    def next_stop(self, pos):
        while pos < len(self) - 2:
            if self.is_stop(pos):
                return pos
            pos += 3

    def is_start(self, pos):
        assert pos < len(self.seq) - 2
        return self.seq[pos:pos + 3] == 'ATG'

    def is_stop(self, pos):
        assert pos >= 0
        return self.seq[pos:pos + 3] in ['TAA', 'TGA', 'TAG']

    def run_RNAplfold(self):
        p1 = sp.Popen(['echo', self.seq], stdout=sp.PIPE, close_fds=True)
        p2 = sp.Popen(['RNAplfold', '-o'], stdin=p1.stdout, stdout=sp.PIPE, close_fds=True)
        p1.stdout.close()
        p2.communicate()
        return self._parse_RNAplfold()

    def run_RNAfold(self, skip_matches=0):
        p1 = sp.Popen(['echo', self.seq], stdout=sp.PIPE, close_fds=True)
        p2 = sp.Popen(['RNAfold', '--noPS'], stdin=p1.stdout, stdout=sp.PIPE, close_fds=True)
        p1.stdout.close()
        output = p2.communicate()
        return self._parse_RNAfold(output[0].decode(), skip_matches)

    def run_RNAduplex(self, query_seq):
        p1 = sp.Popen(['echo', self.seq, query_seq], stdout=sp.PIPE, close_fds=True)
        p2 = sp.Popen(['xargs', '-n1'], stdin=p1.stdout, stdout=sp.PIPE, close_fds=True)
        p3 = sp.Popen(['RNAduplex', '-s'], stdin=p2.stdout, stdout=sp.PIPE, close_fds=True)
        p1.stdout.close()
        p2.stdout.close()
        output = p3.communicate()
        return self._parse_RNAduplex(output[0].decode())

    def _parse_RNAplfold(self):
        tot_score, fout = 0, 'plfold_basepairs'
        if os.path.isfile(fout):
            for linea in open(fout):
                splat = [x for x in linea.split(' ') if x]
                tot_score += float(splat[2])
            os.remove(fout)
        return tot_score

    def _parse_RNAfold(self, output, skip_matches=0):
        linea = output.split('\n')[1]
        fold = linea.split(' ')[0]
        n_par, dis = 0, 0
        for char in fold:
            if char == '(':
                n_par += 1
            elif char == ')':
                n_par -= 1
            if n_par <= skip_matches:
                dis += 1
        return dis

    def _parse_RNAduplex(self, linea):
        linea = linea.rstrip('\n')
        splat = [x for x in linea.split(' ') if x]
        self_start = splat[1].split(',')[0]
        energy = float(splat[-1].lstrip('(').rstrip(')'))
        return self_start, energy


def fix_order(pos1, pos2, strand):
    """
    returns pos1, pos2 given the order of the strand
    if one of the two pos is None, returns twice the other
    """
    if not pos1:
        if pos2:
            return pos2, pos2
        else:
            return None, None
    elif not pos2:
        return pos1, pos1
    pos_sorted = sorted([pos1, pos2])
    if strand == '-':
        return pos_sorted[1], pos_sorted[0]
    return pos_sorted


def comp_rev(seq):
    """complementary reverse of a sequence (only ACTGN)"""
    seq = seq.lower()
    seq = seq[::-1]
    seq = seq.replace('a', 'T')
    seq = seq.replace('t', 'A')
    seq = seq.replace('c', 'G')
    seq = seq.replace('g', 'C')
    seq = seq.replace('n', 'N')
    return seq


def fasta2seq(fasta):
    """extracts the sequence from a multi-line fasta format"""
    assert fasta.startswith('>')
    samtools_splat = fasta.split('\n')
    return ''.join(samtools_splat[1:])


def move_pos(pos, move, strand):
    """moves the position according to the strand"""
    if strand == '+':
        return pos + move
    return pos - move


def before_strand(pos1, pos2, strand):
    if strand == '+' and pos1 <= pos2:
        return True
    if strand == '-' and pos1 >= pos2:
        return True
    return False


def intervals_to_bed12(intervals, strand):
    """
    returns 0-based bed12 fields:
    chromStart, chromStop, blockCounts, blockSizes, blockStarts
    from a 1-based set of intervals, according to strand
    the intervals are ordered depending on the strand
    """
    intervals = list(intervals)
    if strand == '-':
        intervals = [sorted(interval) for interval in intervals[::-1]]
    chromStart = intervals[0][0] - 1
    chromStop = intervals[-1][1]
    blockCounts = len(intervals)
    blockSizes, blockStarts = '', ''
    for interval in intervals:
        blockSizes += str(interval[1] - interval[0] + 1) + ','
        blockStarts += str(interval[0] - 1 - chromStart) + ','
    return chromStart, chromStop, blockCounts, blockSizes, blockStarts
