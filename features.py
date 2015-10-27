class Gene:

    def __init__(self, id, name='', chromosome='', strand='', transcripts=None):
        self.id = id
        self.name = name
        self.chromosome = chromosome
        self.strand = strand
        self.transcripts = transcripts if transcripts else {}


class Transcript:

    def __init__(self, id, gene, exons=None):
        self.id = id
        self.gene = gene
        self.gene.transcripts[id] = self
        self.exons = exons if exons else []
        self.splice_sites = []
        self.n_splice_sites = 0

    @property
    def strand(self):
        return self.gene.strand

    @property
    def chromosome(self):
        return self.gene.chromosome

    def get_splice_sites(self):
        start, splice_sites = -1, []
        for exon in self.exons:
            if start >= 0:
                stop = exon.stop
                splice_sites.append((start, stop))
            start = exon.start
        return splice_sites

    def get_exons_site(self, start, stop):
        query = '_'.join([str(x) for x in self._fix_order(start, stop)])
        for site_n in range(self.n_splice_sites):
            if self.splice_sites[site_n] == query:
                return self.exons[site_n], self.exons[site_n + 1]

    def add_exon(self, exon):
        self.exons.append(exon)
        prev_exon_i = len(self.exons) - 2
        if prev_exon_i >= 0:
            self._add_splice_site(exon.start, self.exons[prev_exon_i].stop)
            self.n_splice_sites += 1

    def _add_splice_site(self, start, stop):
        self.splice_sites.append('_'.join([str(x) for x in (start, stop)]))

    def _fix_order(self, pos1, pos2):
        reverse_it = True if self.strand == '-' else False
        reverse_it = not reverse_it if pos1 > pos2 else reverse_it
        poss = (pos2, pos1) if reverse_it else (pos1, pos2)
        return poss


class Exon:

    def __init__(self, number, transcript, start, stop):
        self.number = number
        self.start = start
        self.stop = stop
        self.transcript = transcript
        assert len(self.transcript.exons) == number - 1
        self.transcript.add_exon(self)

    @property
    def strand(self):
        return self.transcript.strand

    @property
    def chromosome(self):
        return self.transcript.chromosome
