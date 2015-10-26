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


class Exon:

    def __init__(self, number, transcript, start, stop):
        self.number = number
        self.start = start
        self.stop = stop
        self.transcript = transcript
        assert len(self.transcript.exons) == number - 1
        self.transcript.exons.append(self)

    @property
    def strand(self):
        return self.transcript.strand

    @property
    def chromosome(self):
        return self.transcript.chromosome
