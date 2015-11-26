class Genome:

    def __init__(self):
        self.chroms_dict = {}

    @property
    def chromosomes(self):
        for chrom in self.chroms_dict.values():
            yield chrom


class Chromosome:

    def __init__(self, name, genome = None):
        self.name = name
        self.genes_dict = {}
        if genome:
            genome.chroms_dict[name] = self

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
        self.chromosome = chromosome.name
        chromosome.genes_dict[id] = self
        self.strand = strand
        self.trans_dict = {}

    @property
    def transcripts(self):
        for trans in self.trans_dict.values():
            yield trans


class Transcript:

    def __init__(self, id, gene, cds_start=0, cds_stop=0):
        self.id = id
        self.gene = gene
        self.gene.trans_dict[id] = self
        self.exons = []
        self.splice_sites = []
        self.n_splice_sites = 0
        self.cds_start, self.cds_stop = fix_order(cds_start, cds_stop, self.strand)

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
    def tss(self):
        return self.exons[0].start

    @property
    def five_utr_len(self):
        if not self.cds_start:
            return 0
        return self.distance_from_start(self.cds_start)

    def set_cds(self, cds_start, cds_stop):
        self.cds_start, self.cds_stop = fix_order(cds_start, cds_stop, self.strand)

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
            self.n_splice_sites += 1

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
        self.number = number
        self.transcript = transcript
        self.genomic_start = min(start, stop)
        self.genomic_stop = max(start, stop)
        self.start, self.stop = fix_order(start, stop, self.strand)
        self.transcript.add_exon(self)

    @property
    def strand(self):
        return self.transcript.strand

    @property
    def chromosome(self):
        return self.transcript.chromosome

    @property
    def len(self):
        return self.genomic_stop - self.genomic_start + 1

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
            return '\t'.join([str(x) for x in [self.chromosome, self.genomic_start - 1, self.genomic_stop]])
        if self.strand and n_cols == 6:
            return '\t'.join([str(x) for x in [self.chromosome, self.genomic_start - 1, self.genomic_stop,
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


def fix_order(pos1, pos2, strand):
    pos_sorted = sorted([pos1, pos2])
    if strand == '-':
        return pos_sorted[1], pos_sorted[0]
    return pos_sorted
