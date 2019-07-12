from unittest import TestCase
from numpy import array
import Validate_NMR_Restraints as V


class TestValidate_NMR_Restraints(TestCase):
    # def test_generate_json(self):
    #     self.fail()
    #
    # def test_restraints_type_statistics(self):
    #     self.fail()

    def test_get_coordinates(self):
        VR = V.Validate_NMR_Restraints()
        models, maps = VR.get_coordinates('nef_examples/2m2e.cif')
        self.assertEqual(len(models), 20)
        self.assertEqual(len(models[1]), 1196)
        self.assertEqual(len(maps), 20)
        self.assertEqual(len(maps[1]), 1196)

    def test_get_restraints(self):
        VR = V.Validate_NMR_Restraints()
        distance, angle, chain_dict = VR.get_restraints('nef_examples/2m2e.str')
        self.assertEqual(len(distance), 1539)
        self.assertEqual(len(angle), 125)
        self.assertEqual(len(chain_dict), 0)

    def test_get_distance(self):
        VR = V.Validate_NMR_Restraints()
        self.assertAlmostEqual(VR.get_distance(array([2.092, - 21.669, 10.607]), array([-2.563, - 20.933, 7.711])),
                               5.531504045013436)
        self.assertAlmostEqual(VR.get_distance(array([2.284, - 21.812, 10.579]), array([-3.656, - 20.662, 6.527])),
                               7.281813235726387)
        self.assertAlmostEqual(VR.get_distance(array([2.284, - 21.812, 10.579]), array([-3.483, - 19.822, 8.07])),
                               6.596474058161679)
        self.assertAlmostEqual(VR.get_distance(array([2.284, - 21.812, 10.579]), array([-2.487, - 21.234, 7.717])),
                               5.593529207933038)
        self.assertAlmostEqual(VR.get_distance(array([2.122, - 21.675, 10.452]), array([-2.512, - 21.129, 7.424])),
                               5.562450539105944)

    def test_get_dihedral_angle(self):
        VR = V.Validate_NMR_Restraints()
        self.assertAlmostEqual(
            VR.get_dihedral_angle(array([-39.343, - 18.738, 25.524]), array([-38.452, - 19.489, 26.2]),
                                  array([-38.82, - 20.666, 27.017]),
                                  array([-39.206, - 20.263, 28.455])), -86.3831)
        self.assertAlmostEqual(
            VR.get_dihedral_angle(array([-38.301, - 15.653, 28.525]), array([-37.546, - 16.49, 29.271]),
                                  array([-38.124, - 17.512, 30.168]), array([-38.817, - 16.826, 31.368])), -69.0679)
        self.assertAlmostEqual(
            VR.get_dihedral_angle(array([-38.401, - 18.131, 26.899]), array([-38.627, - 19.387, 27.656]),
                                  array([-39.485, - 19.101, 28.901]), array([-39.113, - 18.045, 29.638])), -46.9763)
        self.assertAlmostEqual(
            VR.get_dihedral_angle(array([-37.639, - 20.237, 26.025]), array([-38.122, - 21.486, 26.662]),
                                  array([-38.755, - 21.184, 28.03]), array([-38.373, - 20.041, 28.618])), -22.4749)
        self.assertAlmostEqual(
            VR.get_dihedral_angle(array([-36.636, - 22.659, 26.068]), array([-36.611, - 23.994, 26.692]),
                                  array([-37.297, - 23.953, 28.069]), array([-37.269, - 22.769, 28.711])), -27.4804)
        self.assertAlmostEqual(
            VR.get_dihedral_angle(array([-37.452, - 19.287, 28.582]), array([-37.772, - 20.674, 28.987]),
                                  array([-37.757, - 20.825, 30.526]), array([-36.981, - 19.957, 31.217])), -25.954)

    def test_r6sum(self):
        VR = V.Validate_NMR_Restraints()
        self.assertAlmostEqual(
            VR.r6sum([3.7383778567715753, 2.5503327233912043, 3.8017359192874038, 3.0633733040555153]),
            2.377152042171228)
        self.assertAlmostEqual(
            VR.r6sum([3.1577354227357315, 2.6224320772900893, 3.7876963975482538, 3.7942507824338696]),
            2.4384378568386813)
        self.assertAlmostEqual(
            VR.r6sum([2.9301226936768354, 3.803568193157576, 2.5916317253807493, 3.7542907452673391]),
            2.3753017323279115)
        self.assertAlmostEqual(VR.r6sum(
            [6.19694150367744, 5.2189753783669079, 5.7625093492331976, 5.2472198353032624, 3.7908882864046562,
             3.9508109294168956]), 3.323332779586526)
        self.assertAlmostEqual(VR.r6sum(
            [5.3662990039691216, 4.7131832130737283, 6.0063677043617627, 5.2097326226976355, 5.6391735210046523,
             6.2487019452042993]), 3.978124439298335)
        self.assertAlmostEqual(VR.r6sum(
            [5.729679659457414, 5.5690237923715156, 4.7813129995849462, 5.4157897854329606, 5.5754325392744182,
             4.4649329222285079]), 3.7766787060218725)
        self.assertAlmostEqual(VR.r6sum(
            [5.7141960064386979, 6.6101051428853959, 5.3858857210304771, 6.2696460027660228, 5.0519800078780968,
             4.9115820261907466]), 4.041359228884005)
        self.assertAlmostEqual(VR.r6sum(
            [5.7753563526417988, 5.0417998770280441, 6.1301254473297657, 4.975516857573691, 5.7589886264864223,
             6.0961632195996804]), 4.071618867159647)
        self.assertAlmostEqual(VR.r6sum(
            [4.1302409130703266, 5.476898483631043, 4.1014237771778719, 6.2213707492802586, 5.6205676759558729,
             5.1053908763188751]), 3.4893696605967657)
        self.assertAlmostEqual(VR.r6sum(
            [6.0595155746973681, 5.7729182395041745, 4.967711545571059, 5.7677892645276136, 6.0508190354694964,
             4.9451896829140924]), 4.04095191120266)
        self.assertAlmostEqual(VR.r6sum(
            [7.6136000026268791, 8.4906253008833232, 6.8902264839408582, 8.042190870652103, 7.6417561463318115,
             8.5295367986778778, 5.8465535146785426, 6.5694326239029213, 4.2585158212691905, 4.8715276864655124,
             4.7346136061985025, 5.8463897406861269]), 3.650949040820588)
        self.assertAlmostEqual(VR.r6sum(
            [7.6549069883310814, 8.7528249725445804, 6.7277271793674878, 8.060614182554577, 7.5010949200766683,
             8.6124757764535964, 4.8805365483725272, 6.0826333935229089, 6.1925426926263514, 7.0841752519259433,
             4.6521592836015415, 5.3865690378941578]), 3.88988867366983)
        self.assertAlmostEqual(VR.r6sum(
            [5.2693563174262525, 4.2433343021732304, 6.7166828122221149, 5.4186812048689479, 6.9056374072202775,
             5.9027740935936226, 8.3091696937780721, 7.0431365881970498, 8.0815244230281191, 6.5802338104356135,
             7.9136819496363371, 6.4942862579347373]), 3.774507721666569)


    # def test_calculate_distance_violations(self):
    #     self.fail()
    #
    # def test_calculate_angle_violations(self):
    #     self.fail()
    #
    # def test_write_xml(self):
    #     self.fail()
    #
    # def test_write_xml_simple(self):
    #     self.fail()
    #
    # def test_calculate_violation_statistics(self):
    #     self.fail()
    #
    # def test_bin_distance_violations(self):
    #     self.fail()
    #
    # def test_bin_angle_violations(self):
    #     self.fail()
