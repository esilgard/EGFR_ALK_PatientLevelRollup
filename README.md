# EGFR_ALK_PatientLevelRollup
basic logic for assigning patient level EGFR and ALK status labels based on report level labels

- default to "unknown" by "none" (the test was not done/not mentioned)

- if the report level label is in the trumping list (ordered in patient_level_rollup.json) - give it the appropriate rank

- if it is not a feasible label and not in the trumping list (e.g. "Positive by None") and there is a test result (i.e. positive or negative), then use that result and default the method to the most common method (and give it the appropriate rank)
