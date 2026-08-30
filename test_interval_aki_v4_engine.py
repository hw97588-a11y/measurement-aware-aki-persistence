"""Synthetic boundary tests for the v4 interval-censored phenotype engine."""

from __future__ import annotations

import unittest

from interval_aki_v4_engine import M48, classify_first_episode
from run_interval_aki_primary import Spell


def spell(end_minutes: float = 7 * 24 * 60) -> Spell:
    return Spell("synthetic", end_minutes, extra={"cluster_id": "synthetic"})


class IntervalAkiV4Tests(unittest.TestCase):
    def classify(self, values, end=7 * 24 * 60, policy="first_recovery"):
        episode = classify_first_episode("synthetic", spell(end), values, policy)
        self.assertIsNotNone(episode)
        return episode

    def test_recovery_just_before_or_at_48_hours_is_transient(self):
        before = self.classify([(0, 1.0), (60, 1.5), (M48 - 1, 1.0)])
        exact = self.classify([(0, 1.0), (60, 1.5), (M48, 1.0)])
        self.assertEqual(before.category, "definite_transient")
        self.assertEqual(exact.category, "definite_transient")

    def test_recovery_just_after_48_hours_is_indeterminate(self):
        episode = self.classify([(0, 1.0), (60, 1.5), (M48 + 1, 1.0)])
        self.assertEqual(episode.category, "interval_indeterminate")

    def test_lower_duration_strictly_over_48_is_persistent(self):
        # Persistence requires an observed AKI-positive value more than 48 h
        # after the latest possible onset.  A late first recovery alone leaves
        # the recovery interval open at its lower end and is indeterminate.
        episode = self.classify([(0, 1.0), (60, 1.5), (60 + M48 + 1, 1.5), (60 + M48 + 2, 1.0)])
        self.assertEqual(episode.category, "definite_persistent")

    def test_late_first_recovery_without_late_positive_is_indeterminate(self):
        episode = self.classify([(0, 1.0), (60, 1.5), (60 + M48 + 1, 1.0)])
        self.assertEqual(episode.category, "interval_indeterminate")
        self.assertEqual(episode.duration_lower, 0.0)

    def test_onset_and_recovery_intervals_can_jointly_cross_boundary(self):
        episode = self.classify([(0, 1.0), (6 * 60, 1.5), (50 * 60, 1.0)])
        self.assertEqual(episode.category, "interval_indeterminate")
        self.assertLessEqual(episode.duration_lower, M48)
        self.assertGreater(episode.duration_upper, M48)

    def test_absent_recovery_is_not_automatically_persistent(self):
        episode = self.classify([(0, 1.0), (60, 1.5)], end=60 + 47 * 60)
        self.assertEqual(episode.category, "right_censored_unresolved")
        self.assertFalse(episode.coverage_48h)
        self.assertTrue(episode.structural_coverage_censored)

    def test_no_recheck_after_49_hours_is_monitoring_not_structural_censoring(self):
        episode = self.classify([(0, 1.0), (60, 1.5)], end=60 + 49 * 60)
        self.assertEqual(episode.category, "right_censored_unresolved")
        self.assertTrue(episode.coverage_48h)
        self.assertFalse(episode.structural_coverage_censored)

    def test_duplicate_values_are_median_deduplicated(self):
        episode = self.classify([(0, 1.0), (60, 1.5), (M48, 0.9), (M48, 1.1)])
        self.assertEqual(episode.category, "definite_transient")
        self.assertEqual(episode.recovery_value, 1.0)

    def test_confirmed_recovery_requires_six_hour_gap(self):
        no_gap = self.classify([(0, 1.0), (60, 1.5), (180, 1.0), (180 + 5 * 60, 1.0)], policy="two_recoveries_6h")
        six_hour_gap = self.classify([(0, 1.0), (60, 1.5), (180, 1.0), (180 + 6 * 60, 1.0)], policy="two_recoveries_6h")
        self.assertEqual(no_gap.category, "right_censored_unresolved")
        self.assertEqual(six_hour_gap.category, "definite_transient")

    def test_recurrent_aki_blocks_recovery_confirmation_until_two_recovered_values(self):
        episode = self.classify([
            (0, 1.0), (60, 1.5), (180, 1.0), (240, 1.5),
            (360, 1.0), (360 + 6 * 60, 1.0),
        ], policy="two_recoveries_6h")
        self.assertEqual(episode.recovery_upper, 360)
        self.assertEqual(episode.category, "definite_transient")

    def test_first_48_hours_aki_history_is_detected(self):
        episode = self.classify([(0, 1.0), (30, 1.5), (M48 + 60, 1.5)])
        self.assertFalse(episode.no_identifiable_aki_first48h)


if __name__ == "__main__":
    unittest.main(verbosity=2)
