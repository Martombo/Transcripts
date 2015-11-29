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
        self.genes_dict = {}
        if genome:
            self.genome = genome
            genome.chroms_dict[name] = self

    def get_sequence(self, start_stop, strand):
        return self.genome.get_sequence(self.name, start_stop, strand)

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


class Transcript:

    def __init__(self, id, gene, cds_start=None, cds_stop=None):
        self.id = id
        self.gene = gene
        self.gene.trans_dict[id] = self
        self.exons = []
        self.splice_sites = []
        self.cds_start, self.cds_stop = fix_order(cds_start, cds_stop, self.strand)

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
        after_stop = move_pos(self.cds_stop, +2, self.strand)
        return self.distance_to_end(self.cds_start) - self.distance_to_end(after_stop)

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
            n_stop_14 = bam.get_coverage(self.chromosome, move_pos(self.tss, -14, self.strand) - 1, strand=self.strand)
            n_stop_7 = bam.get_coverage(self.chromosome, move_pos(self.tss, -7, self.strand) - 1, strand=self.strand)
            n_stop = bam.get_coverage(self.chromosome, self.tss - 1, strand=self.strand)
            return max(n_stop_14, n_stop_7, n_stop)

    def _add_splice_site(self, start, stop):
        self.splice_sites.append('_'.join([str(x) for x in (start, stop)]))

    def distance_to_end(self, pos):
        distance = 0
        for exon in self.exons[::-1]:
            if exon.includes(pos):
                return distance + exon.distance_to_end(pos)
            distance += exon.len

    def distance_from_start(self, pos):
        distance = 0
        for exon in self.exons:
            if exon.includes(pos):
                return distance + exon.distance_from_start(pos)
            distance += exon.len
        return distance


class Exon:

    def __init__(self, number, transcript, start, stop):
        self.transcript = transcript
        self.genomic_start = min(start, stop)
        self.genomic_stop = max(start, stop)
        assert len(self.transcript.exons) == number - 1
        if number > 1:
            pass
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

    @property
    def len(self):
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

    def includes(self, pos1, pos2=None):
        if self.genomic_start <= pos1 <= self.genomic_stop:
            return True
        if pos2 and self.genomic_start <= pos2 <= self.genomic_stop:
            return True
        return False


class GenomicInterval:

    def __init__(self, chromosome, start, stop, name='', score='', strand='', n_quantiles=50):
        """
        positions are 1-based
        """
        self.chromosome = chromosome
        self.start, self.stop = fix_order(start, stop, strand)
        self.genomic_start = min(start, stop)
        self.genomic_stop = max(start, stop)
        self.name, self.score, self.strand = name, score, strand
        self.quantile_width, self.quantiles = self._setup_quantiles(n_quantiles)
        if strand and not name:
            self.name = 'interval'
        if strand and not score:
            self.score = 0

    @property
    def len(self):
        return self.genomic_stop - self.genomic_start + 1

    def bed(self, n_cols=6):
        if n_cols < 4:
            return '\t'.join([str(x) for x in [self.chromosome.name, self.genomic_start - 1, self.genomic_stop]])
        if self.strand and n_cols == 6:
            return '\t'.join([str(x) for x in [self.chromosome.name, self.genomic_start - 1, self.genomic_stop,
                                               self.name, self.score, self.strand]])

    def includes(self, pos1, pos2=None):
        if self.genomic_start < pos1 < self.genomic_stop:
            return True
        if pos2 and self.genomic_start < pos2 < self.genomic_stop:
            return True
        return False

    def normalize_quantiles(self, divisor=None):
        if not divisor:
            divisor = max(self.quantiles + [1])
        divisor = float(divisor)
        for i in range(len(self.quantiles)):
            self.quantiles[i] /= divisor

    def add(self, pos, n=1):
        n_bin = int(self._distance_from_start(pos) / self.quantile_width)
        self.quantiles[n_bin] += n

    def _distance_from_start(self, pos):
        return abs(pos - self.start)

    def _setup_quantiles(self, n_quantiles):
        if not n_quantiles:
            return 0, []
        quantile_width = self.len / float(n_quantiles)
        quantiles = [0] * n_quantiles
        return quantile_width, quantiles


class Sequence:

    def __init__(self, seq):
        self.seq = seq
        self.len = len(seq)

    def get_orfs(self):
        orfs = []
        for pos in range(self.len - 2):
            if self.is_start(pos):
                stop = self.next_stop(pos + 3)
                if stop:
                    orfs.append([pos, stop])
        return orfs

    def next_stop(self, pos):
        while pos < self.len - 2:
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
        p1 = sp.Popen(['echo', self.seq], stdout=sp.PIPE)
        p2 = sp.Popen(['RNAplfold', '-o'], stdin=p1.stdout, stdout=sp.PIPE)
        p1.stdout.close()
        p2.communicate()
        return self._parse_RNAplfold()

    def _parse_RNAplfold(self):
        fold_array, fout = [0] * self.len, 'plfold_basepairs'
        for linea in open(fout):
            splat = [x for x in linea.split(' ') if x]
            fold_array = self._add_fold_score(fold_array, splat[0], splat[2])
            fold_array = self._add_fold_score(fold_array, splat[1], splat[2])
        os.remove(fout)
        return [round(x, 3) for x in fold_array]

    def _add_fold_score(self, fold_array, pos_str, score):
        pos = int(pos_str) - 1
        if 0 <= pos < len(fold_array):
            fold_array[pos] += float(score)
        return fold_array


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
        if strand == '+':
            return pos + move
        return pos - move
