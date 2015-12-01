import unittest as ut
import features as ft
import reader as rd
import pysam


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

    def test_move_pos(self):
        pos = ft.move_pos(100,+10, '+')
        self.assertEqual(110, pos)

    def test_move_pos2(self):
        pos = ft.move_pos(100,+10, '-')
        self.assertEqual(90, pos)


class TestExons(ut.TestCase):

    chr1 = ft.Chromosome('chr1')
    gene1 = ft.Gene('g1', chr1, 'gene1', '+')
    gene1_rev = ft.Gene('g1', chr1, 'gene1', '-')

    def test_distance_from_start(self):
        trans1 = ft.Transcript('t1', self.gene1)
        exon1 = ft.Exon(1, trans1, 100, 200)
        self.assertEqual(10, exon1.distance_from_start(110))

    def test_distance_from_start_rev(self):
        trans1 = ft.Transcript('t1', self.gene1_rev)
        exon1 = ft.Exon(1, trans1, 100, 200)
        self.assertEqual(70, exon1.distance_from_start(130))

    def test_includes(self):
        trans1 = ft.Transcript('t1', self.gene1_rev)
        exon1 = ft.Exon(1, trans1, 100, 200)
        self.assertTrue(exon1.includes(exon1.start))
        self.assertTrue(exon1.includes(exon1.genomic_stop))


class TestTranscripts(ut.TestCase):

    chr1 = ft.Chromosome('chr1')
    gene1 = ft.Gene('g1', chr1, 'gene1', '+')
    gene2 = ft.Gene('g1', chr1, 'gene1', '-')

    def test_simple(self):
        trans1 = ft.Transcript('t1', self.gene1)
        exon1 = ft.Exon(1, trans1, 100, 200)
        self.assertEquals(0, len(trans1.splice_sites))

    def test_2exons(self):
        trans1 = ft.Transcript('t1', self.gene1)
        exon1 = ft.Exon(1, trans1, 100, 200)
        exon2 = ft.Exon(2, trans1, 300, 400)
        self.assertEquals(2, len(trans1.exons))

    def test_2trans(self):
        trans1 = ft.Transcript('t1', self.gene1)
        trans2 = ft.Transcript('t2', self.gene1)
        self.assertEquals(2, len(self.gene1.trans_dict.values()))

    def test_add_cds(self):
        trans = ft.Transcript('t', self.gene1)
        trans.add_cds(100, 200)
        self.assertEqual(100, trans.cds_start)
        self.assertEqual(200, trans.cds_stop)

    def test_add_cds_rev(self):
        trans = ft.Transcript('t', self.gene2)
        trans.add_cds(100, 200)
        self.assertEqual(200, trans.cds_start)
        self.assertEqual(100, trans.cds_stop)

    def test_add_2cds(self):
        trans = ft.Transcript('t', self.gene1)
        trans.add_cds(100, 200)
        trans.add_cds(500, 600)
        self.assertEqual(100, trans.cds_start)
        self.assertEqual(600, trans.cds_stop)

    def test_add_2cds_rev(self):
        trans = ft.Transcript('t', self.gene2)
        trans.add_cds(500, 600)
        trans.add_cds(100, 200)
        self.assertEqual(600, trans.cds_start)
        self.assertEqual(100, trans.cds_stop)


class TestGene(ut.TestCase):

    chr1 = ft.Chromosome('chr1')
    gene1 = ft.Gene('g1', chr1, 'gene1', '+')
    trans1 = ft.Transcript('t1', gene1)

    def test_transcripts(self):
        n = 0
        for trans in self.gene1.transcripts:
            n += 1
        self.assertEqual(1, n)


class TestChromosome(ut.TestCase):

    chr1 = ft.Chromosome('chr1')
    gene1 = ft.Gene('g1', chr1)
    trans1 = ft.Transcript('t1', gene1)
    trans2 = ft.Transcript('t2', gene1)
    gene2 = ft.Gene('g2', chr1)
    trans3 = ft.Transcript('t3', gene2)

    def test_transcripts(self):
        n = 0
        for k in self.chr1.transcripts:
            n += 1
        self.assertEqual(3, n)


# class TestTranscriptsGroup(ut.TestCase):
#
#     chr1 = ft.Chromosome('chr1')
#     gene = ft.Gene('g1', chr1, 'gene1', '+')
#     trans_100 = ft.Transcript('t1', gene)
#     exon_100 = ft.Exon(1, trans_100, 100, 200)
#     trans_130 = ft.Transcript('t2', gene)
#     exon_130 = ft.Exon(1, trans_130, 130, 230)
#     gene_rev = ft.Gene('g2', chr1, 'gene2', '-')
#     trans_rev_250 = ft.Transcript('t3', gene_rev)
#     exon_rev_250 = ft.Exon(1, trans_rev_250, 100, 250)
#     trans_rev_30 = ft.Transcript('t4', gene_rev)
#     exon_30 = ft.Exon(1, trans_rev_30, 10, 30)
#
#     def test_add_1tss(self):
#         trans_group = gr.TranscriptsGroup()
#         trans_group._add_tss(self.trans_100)
#         self.assertEquals(1, len(trans_group.transcripts['chr1']))
#
#     def test_add_2tss(self):
#         trans_group = gr.TranscriptsGroup()
#         trans_group._add_tss(self.trans_100)
#         trans_group._add_tss(self.trans_130)
#         self.assertEquals(2, len(trans_group.transcripts['chr1']))
#         self.assertEquals(100, trans_group.transcripts['chr1'][0].tss)
#         self.assertEquals(130, trans_group.transcripts['chr1'][1].tss)
#
#     def test_add_2tss_order(self):
#         trans_group = gr.TranscriptsGroup()
#         trans_group._add_tss(self.trans_130)
#         trans_group._add_tss(self.trans_100)
#         self.assertEquals(100, trans_group.transcripts['chr1'][0].tss)
#         self.assertEquals(130, trans_group.transcripts['chr1'][1].tss)
#
#     def test_add_3tss_rev(self):
#         trans_group = gr.TranscriptsGroup()
#         trans_group._add_tss(self.trans_100)
#         trans_group._add_tss(self.trans_rev_250)
#         trans_group._add_tss(self.trans_rev_30)
#         self.assertEquals(3, len(trans_group.transcripts['chr1']))
#         self.assertEquals(30, trans_group.transcripts['chr1'][0].tss)
#         self.assertEquals(100, trans_group.transcripts['chr1'][1].tss)
#         self.assertEquals(250, trans_group.transcripts['chr1'][2].tss)
#
#     def test_get_antisense_in_range(self):
#         trans_group = gr.TranscriptsGroup()
#         trans_group._add_tss(self.trans_100)
#         trans_group._add_tss(self.trans_rev_250)
#         antisense = trans_group._as_in_range(range(1, 2), trans_group.transcripts['chr1'], '+', 'lulli')
#         self.assertEquals(self.trans_rev_250, antisense)
#
#     def test_get_closest_antisense(self):
#         trans_group = gr.TranscriptsGroup()
#         trans_group._add_tss(self.trans_100)
#         trans_group._add_tss(self.trans_rev_250)
#         antisense = trans_group._get_closest_antisense(self.trans_100, 0, trans_group.transcripts['chr1'])
#         self.assertEquals(self.trans_rev_250, antisense)
#
#     def test_get_closest_antisense2(self):
#         trans_group = gr.TranscriptsGroup()
#         trans_group._add_tss(self.trans_100)
#         trans_group._add_tss(self.trans_rev_250)
#         trans_group._add_tss(self.trans_rev_30)
#         antisense = trans_group._get_closest_antisense(self.trans_100, 1, trans_group.transcripts['chr1'])
#         self.assertEquals(self.trans_rev_30, antisense)


class TestGenomicInterval(ut.TestCase):

    chr1 = ft.Chromosome('1')
    gi = ft.GenomicInterval(chr1, 1, 100)

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


class TestSequence(ut.TestCase):

    def test_is_start(self):
        seq = ft.Sequence('AAAAAATGCCCCCC')
        self.assertFalse(seq.is_start(0))
        self.assertFalse(seq.is_start(11))
        self.assertTrue(seq.is_start(5))

    def test_is_stop(self):
        seq = ft.Sequence('AAAATAGCCGCCCC')
        self.assertFalse(seq.is_stop(0))
        self.assertFalse(seq.is_stop(10))
        self.assertTrue(seq.is_stop(4))

    def test_next_stop(self):
        seq = ft.Sequence('ACGACGTGACCGCCCC')
        self.assertEqual(6, seq.next_stop(0))
        self.assertIsNone(seq.next_stop(1))
        self.assertIsNone(seq.next_stop(7))

    def test_get_orfs(self):
        seq = ft.Sequence('AATGCCCTGACCC')
        orfs = seq.get_orfs()
        self.assertEqual(1, len(orfs))
        self.assertEqual(1, orfs[0][0])
        self.assertEqual(7, orfs[0][1])

    def test_get_2orfs(self):
        seq = ft.Sequence('AATGCCCTGACCCGATGTAA')
        orfs = seq.get_orfs()
        self.assertEqual(2, len(orfs))
        self.assertEqual(1, orfs[0][0])
        self.assertEqual(7, orfs[0][1])
        self.assertEqual(14, orfs[1][0])
        self.assertEqual(17, orfs[1][1])

    def test_rnaplfold(self):
        seq = ft.Sequence('aaaaaaaaagggggggggttttttttt')
        folds = seq.run_RNAplfold()
        self.assertGreater(sum(folds), 0)


class TestReaderGtf(ut.TestCase):

    attr = 'gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_id "ENST00000249857"; exon_number "1";\n'
    mocking = '\t'.join(['1', 'havana', 'exon', '11869', '14409', '.', '+', '.', attr])

    def test_parse_line(self):
        gtf = rd.Gtf('', test=True)
        line_dic = gtf._parse_line(self.mocking)
        self.assertIsNotNone(line_dic)
        gene = gtf._get_gene(line_dic)
        self.assertIsNotNone(gene)

    def test_it(self):
        gtf = rd.Gtf('/Users/martin/notDropbox/utils/genes/Homo_sapiens.GRCh38.81_head.gtf')
        genome = gtf.get_genome()


class TestReaderBam(ut.TestCase):

    read = pysam.AlignedSegment()
    parser_bam = rd.Bam(test=True)

    def test_det_strand(self):
        self.read.is_reverse = False
        self.read.is_read2 = False
        strand = self.parser_bam.determine_strand(self.read)
        self.assertEquals('+', strand)

    def test_det_strand_rev(self):
        self.read.is_reverse = True
        self.read.is_read2 = False
        self.assertNotEqual('reverse',self.parser_bam.reads_orientation)
        strand = self.parser_bam.determine_strand(self.read)
        self.assertEquals('-', strand)

    def test_det_strand_rev_mate(self):
        self.read.is_reverse = True
        self.read.is_read2 = True
        strand = self.parser_bam.determine_strand(self.read)
        self.assertEquals('+', strand)

    def test_det_strand_opp(self):
        self.read.is_reverse = False
        self.read.is_read2 = False
        self.parser_bam.reads_orientation = 'reverse'
        strand = self.parser_bam.determine_strand(self.read)
        self.parser_bam.reads_orientation = 'forward'
        self.assertEquals('-', strand)
