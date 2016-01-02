import subprocess as sp
import os


class Genome:

    def __init__(self, fasta_path=''):
        self.fasta_path = fasta_path
        self.chroms_dict = {}

    def get_sequence(self, chrom, start_stop, strand):
        """finds the sequence of a genomic interval '1:387941-388099', '+' """
        assert self.fasta_path
        pos_string = chrom + ':' + str(start_stop[0]) + '-' + str(start_stop[1])
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

    def get_sequence(self, start_stop, strand):
        return self.genome.get_sequence(self.name, start_stop, strand)

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
        self.trans_dict = {}

    @property
    def transcripts(self):
        for trans in self.trans_dict.values():
            yield trans

    def get_sequence(self, start_stop):
        return self.chromosome.get_sequence(start_stop, self.strand)

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

    def __init__(self, id, gene=None, cds_start=None, cds_stop=None, tss=None):
        self.id = id
        self.gene = gene
        self.gene.trans_dict[id] = self
        self.exons, self.splice_sites = [], []
        self.cds_start, self.cds_stop = fix_order(cds_start, cds_stop, self.strand)

    def __len__(self):
        return sum([len(exon) for exon in self.exons])

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
        return self.distance_from_start(self.cds_start)

    @property
    def cds_len(self):
        if not (self.cds_start and self.cds_stop):
            return 0
        return self.distance_to_end(self.cds_start) - self.three_utr_len

    @property
    def tss(self):
        if len(self.exons) > 0:
            return self.exons[0].start

    def get_sequence(self):
        """returns the whole transcript sequence"""
        return ''.join([exon.get_sequence() for exon in self.exons])

    def get_5utr_seq(self):
        if self.cds_start:
            return ''.join([exon.get_seq_from_start(self.cds_start) for exon in self.exons])

    def get_3utr_seq(self):
        if self.cds_stop:
            first_3utr = move_pos(self.cds_stop, +3, self.strand)
            return ''.join([exon.get_seq_to_stop(first_3utr) for exon in self.exons])

    def get_exons_from_site(self, start, stop):
        query = '_'.join([str(x) for x in fix_order(start, stop, self.strand)])
        for site_n in range(self.n_splice_sites):
            if self.splice_sites[site_n] == query:
                return self.exons[site_n], self.exons[site_n + 1]

    def add_exon(self, exon):
        self.exons.append(exon)
        prev_exon_i = len(self.exons) - 2
        if prev_exon_i >= 0:
            self._add_splice_site(self.exons[prev_exon_i].stop, exon.start)

    def add_cds(self, start, stop):
        start, stop = fix_order(start, stop, self.strand)
        self.cds_start = fix_order(self.cds_start, start, self.strand)[0]
        self.cds_stop = fix_order(self.cds_stop, stop, self.strand)[1]

    def get_stop_counts(self, bam):
        if self.tss:
            return bam.get_coverage(self.chromosome.name, move_pos(self.tss, -10, self.strand) - 1, strand=self.strand)

    def _add_splice_site(self, start, stop):
        self.splice_sites.append('_'.join([str(x) for x in (start, stop)]))

    def distance_from_start(self, pos):
        distance = 0
        for exon in self.exons:
            if exon.includes(pos):
                return distance + exon.distance_from_start(pos)
            distance += len(exon)
        return distance

    def distance_to_end(self, pos):
        distance = 0
        for exon in self.exons[::-1]:
            if exon.includes(pos):
                return distance + exon.distance_to_end(pos)
            distance += len(exon)

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

    def exon_n_with(self, pos):
        """returns the number of the exon that has pos"""
        for exon_n in range(len(self.exons)):
            if self.exons[exon_n].includes(pos):
                return exon_n

    def abs_pos_downstream(self, pos, length):
        """
        walks the transcript downstream, from pos for length
        returns the end position. None if span of transcript is exceeded
        """
        start_n = self.exon_n_with(pos)
        path = abs(pos - self.exons[start_n].stop) + 1
        if path > length:
            return move_pos(pos, length, self.strand)
        length -= path
        for exon in self.exons[start_n+1:]:
            if len(exon) > length:
                return move_pos(exon.start, length, self.strand)
            length -= len(exon)

    def abs_pos_upstream(self, pos, length):
        """
        walks the transcript upstream, from pos for length
        returns the end position. None if span of transcript is exceeded
        """
        start_n = self.exon_n_with(pos)
        path = abs(pos - self.exons[start_n].start) + 1
        if path > length:
            return move_pos(pos, -length, self.strand)
        length -= path
        if start_n == 0:
            return
        for exon in self.exons[start_n-1::-1]:
            if len(exon) > length:
                return move_pos(exon.stop, -length, self.strand)
            length -= len(exon)


class Exon:

    def __init__(self, number, transcript, start, stop):
        self.transcript = transcript
        self.genomic_start = min(start, stop)
        self.genomic_stop = max(start, stop)
        if len(self.transcript.exons) == number - 1:
            self.start, self.stop = fix_order(start, stop, self.strand)
            self.transcript.add_exon(self)

    @property
    def gene(self):
        return self.transcript.gene

    @property
    def strand(self):
        return self.transcript.strand

    @property
    def chromosome(self):
        return self.transcript.chromosome

    def __len__(self):
        return self.genomic_stop - self.genomic_start + 1

    def get_sequence(self):
        """returns the whole exon sequence"""
        return self.gene.get_sequence((self.genomic_start, self.genomic_stop))

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
        the whole sequence if stop is before exon
        or '' if start is after
        """
        start = fix_order(self.start, start, self.strand)[1]
        if not self.includes(start):
            return ''
        return self.gene.get_sequence(sorted((start, self.stop)))

    def distance_from_start(self, pos):
        return abs(self.start - pos)

    def distance_to_end(self, pos):
        return abs(self.stop - pos)

    def includes(self, pos):
        if self.genomic_start <= pos <= self.genomic_stop:
            return True
        return False


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

    def get_orfs(self):
        orfs = []
        for pos in range(len(self) - 2):
            if self.is_start(pos):
                stop = self.next_stop(pos + 3)
                if stop:
                    orfs.append([pos, stop])
        return orfs

    def next_stop(self, pos):
        while pos < len(self) - 2:
            if self.is_stop(pos):
                return pos
            pos += 3

    def is_start(self, pos):
        assert pos < len(self.seq) - 2
        return self.seq[pos:pos+3] == 'ATG'

    def is_stop(self, pos):
        assert pos >= 0
        return self.seq[pos:pos+3] in ['TAA', 'TGA', 'TAG']

    def run_RNAplfold(self):
        p1 = sp.Popen(['echo', self.seq], stdout=sp.PIPE, close_fds=True)
        p2 = sp.Popen(['RNAplfold', '-o'], stdin=p1.stdout, stdout=sp.PIPE, close_fds=True)
        p1.stdout.close()
        p2.communicate()
        return self._parse_RNAplfold()

    def _parse_RNAplfold(self):
        tot_score, fout = 0, 'plfold_basepairs'
        if os.path.isfile(fout):
            for linea in open(fout):
                splat = [x for x in linea.split(' ') if x]
                tot_score += float(splat[2])
            os.remove(fout)
        return tot_score


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
