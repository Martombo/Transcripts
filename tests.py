import unittest as ut
import reader as rd
import features as ft


class TestFeatureGeneral(ut.TestCase):

    def test_fix_order(self):
        sites = ft.fix_order(10, 20, '+')
        self.assertLess(sites[0], sites[1])

    def test_fix_order_opp(self):
        sites = ft.fix_order(20, 10, '+')
        self.assertLess(sites[0], sites[1])

    def test_fix_order_rev(self):
        sites = ft.fix_order(10, 20, '-')
        self.assertGreater(sites[0], sites[1])

    def test_fix_order_rev_opp(self):
        sites = ft.fix_order(20, 10, '-')
        self.assertGreater(sites[0], sites[1])


class TestTranscripts(ut.TestCase):

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


class TestGenomicInterval(ut.TestCase):

    gi = ft.GenomicInterval(1, 1, 100)

    def test_len(self):
        self.assertEquals(100, self.gi._len())

    def test_setup_quantiles_None(self):
        gi = ft.GenomicInterval(1, 1, 100)
        quant_width, quants = gi._setup_quantiles(None)
        self.assertEquals(0, quant_width)
        self.assertEquals([], quants)

    def test_setup_quantiles(self):
        gi = ft.GenomicInterval(1, 1, 100)
        quant_width, quants = gi._setup_quantiles(5)
        self.assertEquals(20, quant_width)
        self.assertEquals(5, len(quants))
        for k in quants:
            self.assertEquals(0,k)

    def test_distance_from_start(self):
        self.assertEquals(0, self.gi._distance_from_start(1))
        self.assertEquals(100, self.gi._distance_from_start(101))

    def test_distance_from_start_rev(self):
        gi = ft.GenomicInterval(1, 1, 100, strand = '-')
        self.assertEquals(1, gi._distance_from_start(99))
        self.assertEquals(100, gi._distance_from_start(0))

    def test_quantiles_add(self):
        gi = ft.GenomicInterval(1, 1, 100)
        for k in range(1,101):
            gi.add(k)
        for j in gi.quantiles:
            self.assertEquals(2, j)

    def test_quantiles_normalize_empty(self):
        gi = ft.GenomicInterval(1, 1, 100)
        gi.normalize_quantiles()
        self.assertEquals(0, sum(self.gi.quantiles))

    def test_quantiles_normalize_divisor(self):
        gi = ft.GenomicInterval(1, 1, 100, n_quantiles=10)
        gi.add(1)
        gi.normalize_quantiles(divisor = 10)
        self.assertEquals(0.1, gi.quantiles[0])
        self.assertEquals(0.1, sum(gi.quantiles))

    def test_quantiles_normalize(self):
        gi = ft.GenomicInterval(1, 1, 100, n_quantiles=50)
        gi.add(1)
        gi.add(1)
        gi.add(3)
        gi.normalize_quantiles()
        self.assertEquals(1, gi.quantiles[0])
        self.assertEquals(0.5, gi.quantiles[1])

    def test_includes_not(self):
        self.assertFalse(self.gi.includes(0))

    def test_includes(self):
        self.assertTrue(self.gi.includes(24))

    def test_includes_pos2(self):
        self.assertTrue(self.gi.includes(0,2))

    def test_includes_pos2_not(self):
        self.assertTrue(self.gi.includes(23,66))
