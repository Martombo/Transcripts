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

    def test_before_strand(self):
        self.assertTrue(ft.before_strand(100, 200, '+'))

    def test_before_strand_rev(self):
        self.assertTrue(ft.before_strand(200, 100, '-'))

    def test_before_strand_false(self):
        self.assertFalse(ft.before_strand(200, 100, '+'))


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
        self.assertEquals(2, len(trans1.exon_dict))

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

    def test_abs_pos_dw_1exon(self):
        trans = ft.Transcript('t', self.gene1)
        exon1 = ft.Exon(1, trans, 100, 200)
        pos = trans.abs_pos_downstream(150, 10)
        self.assertEqual(160, pos)

    def test_abs_pos_dw_2exons(self):
        trans = ft.Transcript('t', self.gene1)
        exon1 = ft.Exon(1, trans, 100, 200)
        exon2 = ft.Exon(2, trans, 300, 400)
        pos = trans.abs_pos_downstream(150, 100)
        self.assertEqual(349, pos)

    def test_abs_pos_dw_3exons(self):
        trans = ft.Transcript('t', self.gene1)
        exon1 = ft.Exon(1, trans, 100, 200)
        exon2 = ft.Exon(2, trans, 300, 400)
        exon2 = ft.Exon(3, trans, 800, 900)
        pos = trans.abs_pos_downstream(150, 200)
        self.assertEqual(848, pos)

    def test_abs_pos_dw_1exon_rev(self):
        trans = ft.Transcript('t', self.gene2)
        exon1 = ft.Exon(1, trans, 100, 200)
        pos = trans.abs_pos_downstream(150, 10)
        self.assertEqual(140, pos)

    def test_abs_pos_dw_2exons_rev(self):
        trans = ft.Transcript('t', self.gene2)
        exon1 = ft.Exon(1, trans, 300, 400)
        exon2 = ft.Exon(2, trans, 100, 200)
        pos = trans.abs_pos_downstream(350, 100)
        self.assertEqual(151, pos)

    def test_abs_pos_dw_3exons_rev(self):
        trans = ft.Transcript('t', self.gene2)
        exon1 = ft.Exon(1, trans, 600, 700)
        exon2 = ft.Exon(2, trans, 300, 400)
        exon3 = ft.Exon(3, trans, 100, 200)
        pos = trans.abs_pos_downstream(650, 200)
        self.assertEqual(152, pos)

    def test_abs_pos_up_1exon(self):
        trans = ft.Transcript('t', self.gene1)
        exon1 = ft.Exon(1, trans, 100, 200)
        pos = trans.abs_pos_upstream(150, 10)
        self.assertEqual(140, pos)

    def test_abs_pos_up_2exons(self):
        trans = ft.Transcript('t', self.gene1)
        exon1 = ft.Exon(1, trans, 100, 200)
        exon2 = ft.Exon(2, trans, 300, 400)
        pos = trans.abs_pos_upstream(340, 100)
        self.assertEqual(141, pos)

    def test_abs_pos_up_3exons(self):
        trans = ft.Transcript('t', self.gene1)
        exon1 = ft.Exon(1, trans, 100, 200)
        exon2 = ft.Exon(2, trans, 300, 400)
        exon2 = ft.Exon(3, trans, 800, 923)
        pos = trans.abs_pos_upstream(820, 210)
        self.assertEqual(112, pos)

    def test_abs_pos_up_1exon_rev(self):
        trans = ft.Transcript('t', self.gene2)
        exon1 = ft.Exon(1, trans, 100, 200)
        pos = trans.abs_pos_upstream(150, 10)
        self.assertEqual(160, pos)

    def test_abs_pos_up_2exons_rev(self):
        trans = ft.Transcript('t', self.gene2)
        exon1 = ft.Exon(1, trans, 400, 500)
        exon2 = ft.Exon(2, trans, 100, 200)
        pos = trans.abs_pos_upstream(120, 110)
        self.assertEqual(429, pos)

    def test_intervals_1exon(self):
        trans = ft.Transcript('t', self.gene1)
        exon1 = ft.Exon(1, trans, 100, 200)
        inter = list(trans.intervals(140, 160))[0]
        self.assertEqual(140, inter[0])
        self.assertEqual(160, inter[1])

    def test_intervals_2exon(self):
        trans = ft.Transcript('t', self.gene1)
        exon1 = ft.Exon(1, trans, 100, 200)
        exon2 = ft.Exon(2, trans, 300, 400)
        inter = list(trans.intervals(140, 360))
        self.assertEqual(140, inter[0][0])
        self.assertEqual(200, inter[0][1])
        self.assertEqual(300, inter[1][0])
        self.assertEqual(360, inter[1][1])

    def test_intervals_2exon_rev(self):
        trans = ft.Transcript('t', self.gene2)
        exon1 = ft.Exon(1, trans, 300, 400)
        exon2 = ft.Exon(2, trans, 100, 200)
        inter = list(trans.intervals(360, 140))
        self.assertEqual(360, inter[0][0])
        self.assertEqual(300, inter[0][1])
        self.assertEqual(200, inter[1][0])
        self.assertEqual(140, inter[1][1])

    def test_distance_from_tss_1exon(self):
        trans = ft.Transcript('t', self.gene1)
        trans.cds_start = 110
        exon1 = ft.Exon(1, trans, 100, 200)
        self.assertEqual(10, trans.distance_from_cds_start(120))

    def test_distance_from_tss_2exons(self):
        trans = ft.Transcript('t', self.gene1)
        trans.cds_start = 110
        exon1 = ft.Exon(1, trans, 100, 200)
        exon2 = ft.Exon(2, trans, 300, 400)
        self.assertEqual(91, trans.distance_from_cds_start(300))


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
        self.assertGreater(folds, 0)


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
