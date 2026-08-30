#!/usr/bin/env python3
"""Unit tests for the NDT continuity-gap diagnostic helpers."""

import unittest

from run_interval_aki_primary import Spell
from run_ndt_continuity_gap_sensitivity import maximum_adjacent_gap, persistence_support_chain


class ContinuityGapTest(unittest.TestCase):
    def setUp(self):
        self.spell = Spell("synthetic", end=7 * 24 * 60)

    def test_maximum_gap_empty_or_singleton(self):
        self.assertEqual(maximum_adjacent_gap([]), 0)
        self.assertEqual(maximum_adjacent_gap([(0, 2.0)]), 0)

    def test_maximum_gap_chain(self):
        self.assertEqual(maximum_adjacent_gap([(0, 2.0), (24 * 60, 2.1), (49 * 60, 2.2)]), 25 * 60)

    def test_chain_stops_at_observed_support_end(self):
        values = [(0, 1.6), (24 * 60, 1.7), (49 * 60, 1.8), (60 * 60, 0.9), (72 * 60, 2.0)]
        chain = persistence_support_chain(self.spell, values, onset_upper=0, baseline=1.0, recovery_lower=49 * 60)
        self.assertEqual([time for time, _ in chain], [0, 24 * 60, 49 * 60])

    def test_subthreshold_value_is_not_positive_state_support(self):
        values = [(0, 1.6), (12 * 60, 1.2), (24 * 60, 1.7)]
        chain = persistence_support_chain(self.spell, values, onset_upper=0, baseline=1.0, recovery_lower=24 * 60)
        self.assertEqual([time for time, _ in chain], [0, 24 * 60])


if __name__ == "__main__":
    unittest.main()
