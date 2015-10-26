import features as ft
import unittest as ut
import parser as ps


class TestGeneTransExon(ut.TestCase):

    def test_simple(self):
        gene1 = ft.Gene('g1', 'gene1', 'chr1', '+')
        trans1 = ft.Transcript('t1', gene1)
        exon1 = ft.Exon(1, trans1, 100, 200)
        splice_sites = trans1.get_splice_sites()
        self.assertEquals(0, len(splice_sites))

    def test_2exons(self):
        gene1 = ft.Gene('g1', 'gene1', 'chr1', '+')
        trans1 = ft.Transcript('t1', gene1)
        exon1 = ft.Exon(1, trans1, 100, 200)
        exon2 = ft.Exon(2, trans1, 300, 400)
        self.assertEquals(2, len(trans1.exons))

    def test_2trans(self):
        gene1 = ft.Gene('g1', 'gene1', 'chr1', '+')
        trans1 = ft.Transcript('t1', gene1)
        trans2 = ft.Transcript('t2', gene1)
        self.assertEquals(2, len(gene1.transcripts.values()))

    def test_get_genes(self):
        parser = ps.Gtf('/Users/martin/notDropbox/utils/genes/prova.gtf')
        parser.get_genes()
        for gene in parser.genes.values():
            for trans in gene.transcripts.values():
                for exon_n in range(len(trans.exons)):
                    self.assertEquals(exon_n, trans.exons[exon_n].number - 1)
