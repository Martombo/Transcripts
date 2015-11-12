import unittest as ut
import features as ft
import groups as gr
import reader as rd


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


class TestExons(ut.TestCase):

    gene1 = ft.Gene('g1', 'gene1', 'chr1', '+')
    gene1_rev = ft.Gene('g1', 'gene1', 'chr1', '-')

    def test_distance_from_start(self):
        trans1 = ft.Transcript('t1', self.gene1)
        exon1 = ft.Exon(1, trans1, 100, 200)
        self.assertEqual(10, exon1.distance_from_start(110))

    def test_distance_from_start_rev(self):
        trans1 = ft.Transcript('t1', self.gene1_rev)
        exon1 = ft.Exon(1, trans1, 100, 200)
        self.assertEqual(70, exon1.distance_from_start(130))


class TranscriptsGroup(ut.TestCase):

    gene = ft.Gene('g1', 'gene1', 'chr1', '+')
    trans_100 = ft.Transcript('t1', gene)
    exon_100 = ft.Exon(1, trans_100, 100, 200)
    trans_130 = ft.Transcript('t2', gene)
    exon_130 = ft.Exon(1, trans_130, 130, 230)
    gene_rev = ft.Gene('g1', 'gene2', 'chr1', '-')
    trans_rev_250 = ft.Transcript('t3', gene_rev)
    exon_rev_250 = ft.Exon(1, trans_rev_250, 100, 250)
    trans_rev_30 = ft.Transcript('t4', gene_rev)
    exon_30 = ft.Exon(1, trans_rev_30, 10, 30)

    def test_add_1tss(self):
        trans_group = gr.TranscriptsGroup()
        trans_group._add_tss(self.trans_100)
        self.assertEquals(1, len(trans_group.transcripts['chr1']))

    def test_add_2tss(self):
        trans_group = gr.TranscriptsGroup()
        trans_group._add_tss(self.trans_100)
        trans_group._add_tss(self.trans_130)
        self.assertEquals(2, len(trans_group.transcripts['chr1']))
        self.assertEquals(100, trans_group.transcripts['chr1'][0].tss)
        self.assertEquals(130, trans_group.transcripts['chr1'][1].tss)

    def test_add_2tss_order(self):
        trans_group = gr.TranscriptsGroup()
        trans_group._add_tss(self.trans_130)
        trans_group._add_tss(self.trans_100)
        self.assertEquals(100, trans_group.transcripts['chr1'][0].tss)
        self.assertEquals(130, trans_group.transcripts['chr1'][1].tss)

    def test_add_3tss_rev(self):
        trans_group = gr.TranscriptsGroup()
        trans_group._add_tss(self.trans_100)
        trans_group._add_tss(self.trans_rev_250)
        trans_group._add_tss(self.trans_rev_30)
        self.assertEquals(3, len(trans_group.transcripts['chr1']))
        self.assertEquals(30, trans_group.transcripts['chr1'][0].tss)
        self.assertEquals(100, trans_group.transcripts['chr1'][1].tss)
        self.assertEquals(250, trans_group.transcripts['chr1'][2].tss)


    def test_get_antisense_in_range(self):
        trans_group = gr.TranscriptsGroup()
        trans_group._add_tss(self.trans_100)
        trans_group._add_tss(self.trans_rev_250)
        antisense = trans_group._get_antisense_in_range(range(1, 2), trans_group.transcripts['chr1'], '+')
        self.assertEquals(self.trans_rev_250, antisense)

    def test_get_closest_antisense(self):
        trans_group = gr.TranscriptsGroup()
        trans_group._add_tss(self.trans_100)
        trans_group._add_tss(self.trans_rev_250)
        antisense = trans_group._get_closest_antisense(self.trans_100, 0, trans_group.transcripts['chr1'])
        self.assertEquals(self.trans_rev_250, antisense)

    def test_get_closest_antisense2(self):
        trans_group = gr.TranscriptsGroup()
        trans_group._add_tss(self.trans_100)
        trans_group._add_tss(self.trans_rev_250)
        trans_group._add_tss(self.trans_rev_30)
        antisense = trans_group._get_closest_antisense(self.trans_100, 1, trans_group.transcripts['chr1'])
        self.assertEquals(self.trans_rev_30, antisense)


class TestGenomicInterval(ut.TestCase):

    gi = ft.GenomicInterval(1, 1, 100)

    def test_len(self):
        self.assertEquals(100, self.gi.len)

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
