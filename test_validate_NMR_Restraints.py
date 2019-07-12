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


    #
    # def test_r6sum(self):
    #     self.fail()
    #
    # def test_r6average(self):
    #     self.fail()
    #
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
