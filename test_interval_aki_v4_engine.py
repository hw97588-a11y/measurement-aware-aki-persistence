"""Synthetic boundary tests for the v4 interval-censored phenotype engine."""

from __future__ import annotations

import unittest
from datetime import datetime

import numpy as np

from interval_aki_v4_engine import M48, classify_first_episode
from run_interval_aki_primary import (
    Spell,
    eicu_age_years,
    mimic_admission_age,
    sicdb_icu_covered_end_minutes,
    sicdb_sex,
)
from run_v4_controlled_thinning import draw_phases_by_identifier, index_episode_match_state


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

    def test_late_day6_onset_uses_followup_beyond_icu_day7(self):
        episode = self.classify([
            (5 * 24 * 60, 1.0),
            (6 * 24 * 60, 1.5),
            (8 * 24 * 60 + 1, 1.5),
            (9 * 24 * 60, 1.0),
        ], end=10 * 24 * 60)
        self.assertEqual(episode.category, "definite_persistent")
        self.assertEqual(episode.recovery_upper, 9 * 24 * 60)

    def test_index_search_end_is_distinct_from_followup_end(self):
        indexed_spell = Spell(
            "synthetic",
            10 * 24 * 60,
            extra={"cluster_id": "synthetic", "index_search_end_minutes": 2 * 24 * 60},
        )
        late_only = classify_first_episode(
            "synthetic",
            indexed_spell,
            [(0, 1.0), (3 * 24 * 60, 1.5), (4 * 24 * 60, 1.0)],
        )
        self.assertIsNone(late_only)
        early_with_late_recovery = classify_first_episode(
            "synthetic",
            indexed_spell,
            [(0, 1.0), (24 * 60, 1.5), (4 * 24 * 60, 1.0)],
        )
        self.assertIsNotNone(early_with_late_recovery)
        self.assertEqual(early_with_late_recovery.recovery_upper, 4 * 24 * 60)

    def test_sicdb_coverage_end_uses_same_icuoffset_origin(self):
        self.assertEqual(sicdb_icu_covered_end_minutes(10 * 3600, 2 * 3600), 8 * 60)
        self.assertIsNone(sicdb_icu_covered_end_minutes(2 * 3600, 2 * 3600))

    def test_mimic_admission_age_uses_anchor_year_offset(self):
        self.assertEqual(mimic_admission_age("17", "2200", datetime(2202, 1, 1)), 19.0)
        self.assertEqual(mimic_admission_age("52", "2180", datetime(2180, 5, 6)), 52.0)
        self.assertIsNone(mimic_admission_age("invalid", "2180", datetime(2180, 5, 6)))

    def test_eicu_protected_age_is_retained_as_a_top_coded_lower_bound(self):
        self.assertEqual(eicu_age_years("> 89"), 90.0)
        self.assertEqual(eicu_age_years(">89"), 90.0)
        self.assertEqual(eicu_age_years("67"), 67.0)
        self.assertIsNone(eicu_age_years("unknown"))

    def test_sicdb_documented_sex_reference_mapping(self):
        self.assertEqual(sicdb_sex("735"), "Male")
        self.assertEqual(sicdb_sex("736"), "Female")
        self.assertEqual(sicdb_sex("737"), "Unknown")
        self.assertEqual(sicdb_sex(""), "Unknown")

    def test_later_recurrence_cannot_replace_missed_index_episode(self):
        original = self.classify([
            (0, 1.0), (60, 1.5), (12 * 60, 1.0),
            (72 * 60, 1.0), (73 * 60, 1.5),
        ])
        scheduled = self.classify([
            (0, 1.0), (12 * 60, 1.0),
            (72 * 60, 1.0), (73 * 60, 1.5),
        ])
        self.assertEqual(
            index_episode_match_state(original, scheduled),
            "index_not_retained_later_recurrence_detected",
        )

    def test_patient_specific_phase_is_shared_within_patient_only(self):
        spells = {
            "stay-a": Spell("stay-a", 7 * 24 * 60, extra={"uniquepid": "patient-1"}),
            "stay-b": Spell("stay-b", 7 * 24 * 60, extra={"uniquepid": "patient-1"}),
            "stay-c": Spell("stay-c", 7 * 24 * 60, extra={"uniquepid": "patient-2"}),
        }
        phases, units = draw_phases_by_identifier(
            ["stay-a", "stay-b", "stay-c"], spells, np.random.default_rng(123), 24 * 60, "patient-specific",
        )
        self.assertEqual(units, 2)
        self.assertEqual(phases["stay-a"], phases["stay-b"])
        self.assertNotEqual(phases["stay-a"], phases["stay-c"])

    def test_global_phase_is_shared_across_reference_episodes(self):
        spells = {
            "stay-a": Spell("stay-a", 7 * 24 * 60, extra={"uniquepid": "patient-1"}),
            "stay-b": Spell("stay-b", 7 * 24 * 60, extra={"uniquepid": "patient-2"}),
        }
        phases, units = draw_phases_by_identifier(
            ["stay-a", "stay-b"], spells, np.random.default_rng(123), 24 * 60, "global",
        )
        self.assertEqual(units, 1)
        self.assertEqual(phases["stay-a"], phases["stay-b"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
