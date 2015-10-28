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

    def get_exons_from_site(self, start, stop):
        query = '_'.join([str(x) for x in self.fix_order(start, stop)])
        for site_n in range(self.n_splice_sites):
            print(self.splice_sites[site_n], query)
            if self.splice_sites[site_n] == query:
                return self.exons[site_n], self.exons[site_n + 1]

    def add_exon(self, exon):
        self.exons.append(exon)
        prev_exon_i = len(self.exons) - 2
        if prev_exon_i >= 0:
            self._add_splice_site(self.exons[prev_exon_i].stop, exon.start)
            self.n_splice_sites += 1

    def fix_order(self, pos1, pos2):
        reverse_it = True if self.strand == '-' else False
        reverse_it = not reverse_it if pos1 > pos2 else reverse_it
        poss = (pos2, pos1) if reverse_it else (pos1, pos2)
        return poss

    def _add_splice_site(self, start, stop):
        self.splice_sites.append('_'.join([str(x) for x in (start, stop)]))


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
