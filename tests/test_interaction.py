#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `Interaction` module."""


import unittest

from pexmd import interaction
import numpy as np

class TestInteraction(unittest.TestCase):
  """Tests for `Interaction` module."""

  def setUp(self):
    """Set up test fixtures, if any."""
    self.four_by3 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                              [-1.0, 0.0, 0.0], [0.0, 1.0, 8.0]], dtype=np.float32)
    self.three_by3 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                               [-1.0, 0.0, 0.0]], dtype=np.float32)
    self.four_type = np.array([1, 1, 1, 1], dtype=np.int32)
    self.two_two_type = np.array([1, 1, 2, 2], dtype=np.int32)

  def tearDown(self):
    """Tear down test fixtures, if any."""
    pass

  def test_create_interaction(self):
    i = interaction.Interaction()
    f, e = i.forces(self.four_by3, self.four_by3, self.four_type)
    np.testing.assert_array_almost_equal(f, np.zeros_like(self.four_by3))


  def test_create_shortrange(self):
    interaction.ShortRange(5.4, "None")

  def test_create_lennardjones(self):
    interaction.LennardJones(5.4, 1.0, 1.0, "None")

  def test_lj_two_noshift(self):
    lj = interaction.LennardJones(5.4, 1.0, 1.0, "None")
    f = lj.pair_force(np.array([2.0**(1.0/6), 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    e = lj.pair_energ(np.array([2.0**(1.0/6), 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    np.testing.assert_array_almost_equal(f, np.zeros(3))
    np.testing.assert_almost_equal(e, -1.0)
    f = lj.pair_force(np.array([0.0, 1.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    e = lj.pair_energ(np.array([0.0, 1.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    np.testing.assert_array_almost_equal(f, np.array([0, 24, 0]))
    np.testing.assert_almost_equal(e, 0.0)
    f = lj.pair_force(np.array([0.0, 7.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    e = lj.pair_energ(np.array([0.0, 7.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    np.testing.assert_array_almost_equal(f, np.zeros(3))
    np.testing.assert_almost_equal(e, 0.0)

  def test_lj_two_displace(self):
    lj = interaction.LennardJones(5.4, 1.0, 1.0, "Displace")
    vcut = -0.0001613169181702531

    f = lj.pair_force(np.array([2.0**(1.0/6), 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    e = lj.pair_energ(np.array([2.0**(1.0/6), 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    np.testing.assert_almost_equal(e, -1.0 - vcut)
    np.testing.assert_array_almost_equal(f, np.zeros(3))
    f = lj.pair_force(np.array([0.0, 1.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    e = lj.pair_energ(np.array([0.0, 1.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    np.testing.assert_array_almost_equal(f, np.array([0, 24, 0]))
    np.testing.assert_almost_equal(e, 0.0 - vcut)
    f = lj.pair_force(np.array([0.0, 7.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    e = lj.pair_energ(np.array([0.0, 7.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    np.testing.assert_array_almost_equal(f, np.zeros(3))
    np.testing.assert_almost_equal(e, 0.0)

  def test_lj_forces_equal(self):
    lj = interaction.LennardJones(5.4, 1.0, 1.0, "None")
    f, e = lj.forces(self.four_by3, self.four_by3)
    force_by_hand = np.array([[0.0, 0.0, 0.0], [23.818359, 0.0, 0.0],
                              [-23.818359, 0.0, 0.0], [0.0, 0.0, 0.0]])
    np.testing.assert_array_almost_equal(f, force_by_hand)

  def test_lj_forces_diff(self):
    lj = interaction.LennardJones(5.4, 1.0, 1.0, "None")
    f, e = lj.forces(self.four_by3, self.four_by3, pairs=np.array([[0, 2], [0, 3], [1, 2], [1, 3]], dtype=np.int64))
    force_by_hand = np.array([[24.0, 0.0, 0.0], [-0.181641, 0.0, 0.0],
                              [-23.818359, 0.0, 0.0], [0.0, 0.0, 0.0]])
    np.testing.assert_array_almost_equal(f, force_by_hand)

#-------------------------- MORSE ----------------------------------------------

  def test_create_morse(self):
    interaction.Morse(5.4, 1.0, 1.0, 2.0, "None")

  def test_morse_two_noshift(self):
    mor = interaction.Morse(5.4, 1.0, 1.0, 2.0, "None")
    f = mor.pair_force(np.array([2.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    e = mor.pair_energ(np.array([2.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    np.testing.assert_array_almost_equal(f, np.zeros(3))
    np.testing.assert_almost_equal(e, 0.0)
    f = mor.pair_force(np.array([0.0, np.log(2)+2, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    e = mor.pair_energ(np.array([0.0, np.log(2)+2, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    np.testing.assert_array_almost_equal(f, np.array([0, -0.5, 0]))
    np.testing.assert_almost_equal(e, 0.25)
    f = mor.pair_force(np.array([0.0, 7.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    e = mor.pair_energ(np.array([0.0, 7.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    np.testing.assert_array_almost_equal(f, np.zeros(3))
    np.testing.assert_almost_equal(e, 0.0)

  def test_morse_two_displace(self):
    mor = interaction.Morse(5.4, 1.0, 1.0, 2.0, "Displace")
    vcut = 0.93436723522719264411398895669

    f = mor.pair_force(np.array([2.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    e = mor.pair_energ(np.array([2.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    np.testing.assert_array_almost_equal(f, np.zeros(3))
    np.testing.assert_almost_equal(e, 0.0 - vcut)
    f = mor.pair_force(np.array([0.0, np.log(2)+2, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    e = mor.pair_energ(np.array([0.0, np.log(2)+2, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    np.testing.assert_array_almost_equal(f, np.array([0, -0.5, 0]))
    np.testing.assert_almost_equal(e, 0.25 - vcut)
    f = mor.pair_force(np.array([0.0, 7.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    e = mor.pair_energ(np.array([0.0, 7.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    np.testing.assert_array_almost_equal(f, np.zeros(3))
    np.testing.assert_almost_equal(e, 0.0)

  def test_morse_forces_equal(self):
    mor = interaction.Morse(5.4, 1.0, 1.0, 2.5, "None")
    f, e = mor.forces(self.four_by3, self.four_by3)
    force_by_hand = np.array([[0.0, 0.0, 0.0], [31.2076957+2.13912211, 0.0, 0.0],
                              [-31.2076957-2.13912211, 0.0, 0.0], [0.0, 0.0, 0.0]])
    np.testing.assert_array_almost_equal(f, force_by_hand)

  def test_morse_forces_diff(self):
    mor = interaction.Morse(5.4, 1.0, 1.0, 2.5, "None")
    f, e = mor.forces(self.four_by3, self.four_by3, pairs=np.array([[0, 2], [0, 3], [1, 2], [1, 3]], dtype=np.int64))
    force_by_hand = np.array([[31.2076957, 0.0, 0.0], [2.13912211, 0.0, 0.0],
                              [-2.13912211-31.2076957, 0.0, 0.0], [0.0, 0.0, 0.0]])
    np.testing.assert_array_almost_equal(f, force_by_hand)

#-------------------------- PAULI ----------------------------------------------

  def test_create_pauli(self):
    interaction.Pauli(5.4, 1.0, 1.0, 2.0, "None")

  def test_pauli_two_noshift(self):
    pauli = interaction.Pauli(5.4, 1.0, 1.0, 2.0, "None")
    f = pauli.pair_force(np.array([2.0, 0.0, 0.0]),
                      np.array([4.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    g = pauli.pair_gorce(np.array([2.0, 0.0, 0.0]),
                      np.array([4.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    e = pauli.pair_energ(np.array([2.0, 0.0, 0.0]),
                      np.array([4.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    np.testing.assert_array_almost_equal(f, np.array([0.036631277777468357, 0, 0]))
    np.testing.assert_array_almost_equal(g, np.array([0.018315638888734179, 0, 0]))
    np.testing.assert_almost_equal(e, 0.018315638888734179)
    f = pauli.pair_force(np.array([0.0, np.sqrt(np.log(2)), 0.0]),
                      np.array([0.0, 2*np.sqrt(np.log(2)), 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    g = pauli.pair_gorce(np.array([0.0, np.sqrt(np.log(2)), 0.0]),
                      np.array([0.0, 2*np.sqrt(np.log(2)), 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    e = pauli.pair_energ(np.array([0.0, np.sqrt(np.log(2)), 0.0]),
                      np.array([0.0, 2*np.sqrt(np.log(2)), 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    np.testing.assert_array_almost_equal(f, np.array([0, 0.41627730557884884, 0]))
    np.testing.assert_array_almost_equal(g, np.array([0, 0.20813865278942442, 0]))
    np.testing.assert_almost_equal(e, 0.5)
    f = pauli.pair_force(np.array([0.0, 7.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    g = pauli.pair_gorce(np.array([0.0, 7.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    e = pauli.pair_energ(np.array([0.0, 7.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    np.testing.assert_array_almost_equal(f, np.zeros(3))
    np.testing.assert_array_almost_equal(g, np.zeros(3))
    np.testing.assert_almost_equal(e, 0.0)

  def test_pauli_two_displace(self):
    pauli = interaction.Pauli(5.4, 1.0, 1.0, 2.0, "Displace")
    vcut = 4.6557157157830781e-07

    f = pauli.pair_force(np.array([2.0, 0.0, 0.0]),
                      np.array([4.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    g = pauli.pair_gorce(np.array([2.0, 0.0, 0.0]),
                      np.array([4.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    e = pauli.pair_energ(np.array([2.0, 0.0, 0.0]),
                      np.array([4.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    np.testing.assert_array_almost_equal(f, np.array([0.036631277777468357, 0, 0]))
    np.testing.assert_array_almost_equal(g, np.array([0.018315638888734179, 0, 0]))
    np.testing.assert_almost_equal(e, 0.018315638888734179-vcut)
    f = pauli.pair_force(np.array([0.0, np.sqrt(np.log(2)), 0.0]),
                      np.array([0.0, 2*np.sqrt(np.log(2)), 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    g = pauli.pair_gorce(np.array([0.0, np.sqrt(np.log(2)), 0.0]),
                      np.array([0.0, 2*np.sqrt(np.log(2)), 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    e = pauli.pair_energ(np.array([0.0, np.sqrt(np.log(2)), 0.0]),
                      np.array([0.0, 2*np.sqrt(np.log(2)), 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    np.testing.assert_array_almost_equal(f, np.array([0, 0.41627730557884884, 0]))
    np.testing.assert_array_almost_equal(g, np.array([0, 0.20813865278942442, 0]))
    np.testing.assert_almost_equal(e, 0.5-vcut)
    f = pauli.pair_force(np.array([0.0, 7.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    g = pauli.pair_gorce(np.array([0.0, 7.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    e = pauli.pair_energ(np.array([0.0, 7.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]),
                      np.array([0.0, 0.0, 0.0]))
    np.testing.assert_array_almost_equal(f, np.zeros(3))
    np.testing.assert_array_almost_equal(g, np.zeros(3))
    np.testing.assert_almost_equal(e, 0.0)


  def test_pauli_forces_equal(self):
    pauli = interaction.Pauli(5.4, 1.0, 1.0, 2.0, "None")
    f, e = pauli.forces(self.four_by3, self.four_by3)
    g, e = pauli.gorces(self.four_by3, self.four_by3)
    force_by_hand = 4*np.array([[0.0, 0.0, 0.0], [0.17485785644169696, 0.0, 0.0],
                              [-0.17485785644169696, 0.0, 0.0], [0.0, 0.0, 0.0]])
    gorce_by_hand = np.array([[0.0, 0.0, 0.0], [0.17485785644169696, 0.0, 0.0],
                              [-0.17485785644169696, 0.0, 0.0], [0.0, 0.0, 0.0]])
    np.testing.assert_array_almost_equal(f, force_by_hand)
    np.testing.assert_array_almost_equal(g, gorce_by_hand)


  def test_pauli_forces_diff(self):
    pauli = interaction.Pauli(5.4, 1.0, 1.0, 2.0, "None")
    f, e = pauli.forces(self.four_by3, self.four_by3, pairs=np.array([[0, 2], [0, 3], [1, 2], [1, 3]], dtype=np.int64))
    g, e = pauli.gorces(self.four_by3, self.four_by3, pairs=np.array([[0, 2], [0, 3], [1, 2], [1, 3]], dtype=np.int64))
    force_by_hand = 4*np.array([[0.13381535712974757, 0.0, 0.0], [0.0410424993119494, 0.0, 0.0],
                              [-0.17485785644169696, 0.0, 0.0], [0.0, 0.0, 0.0]])
    gorce_by_hand = np.array([[0.13381535712974757, 0.0, 0.0], [0.0410424993119494, 0.0, 0.0],
                              [-0.17485785644169696, 0.0, 0.0], [0.0, 0.0, 0.0]])
    np.testing.assert_array_almost_equal(f, force_by_hand)
    np.testing.assert_array_almost_equal(g, gorce_by_hand)