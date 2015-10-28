import unittest as ut
import reader as rd
import features as ft


class TestFeature(ut.TestCase):

    def test_simple(self):
        gene1 = ft.Gene('g1', 'gene1', 'chr1', '+')
        trans1 = ft.Transcript('t1', gene1)
        exon1 = ft.Exon(1, trans1, 100, 200)
        self.assertEquals(0, len(trans1.splice_sites))

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

    def tst_get_genes(self):
        parser = rd.Gtf('/Users/martin/notDropbox/utils/genes/prova.gtf')
        parser.get_genes()
        for gene in parser.genes.values():
            for trans in gene.transcripts.values():
                for exon_n in range(len(trans.exons)):
                    self.assertEquals(exon_n, trans.exons[exon_n].number - 1)


class TestTranscripts(ut.TestCase):

    def test_fix_order(self):
        gene = ft.Gene('', strand='+')
        trans = ft.Transcript('', gene)
        sites = trans.fix_order(10, 20)
        self.assertLess(sites[0], sites[1])

    def test_fix_order_opp(self):
        gene = ft.Gene('', strand='+')
        trans = ft.Transcript('', gene)
        sites = trans.fix_order(20, 10)
        self.assertLess(sites[0], sites[1])

    def test_fix_order_rev(self):
        gene = ft.Gene('', strand='-')
        trans = ft.Transcript('', gene)
        sites = trans.fix_order(10, 20)
        self.assertGreater(sites[0], sites[1])

    def test_fix_order_rev_opp(self):
        gene = ft.Gene('', strand='-')
        trans = ft.Transcript('', gene)
        sites = trans.fix_order(10, 20)
        self.assertGreater(sites[0], sites[1])
