import json
from pathlib import Path


class GCContentTestData:
    @classmethod
    def get_cases(cls):
        cases = []
        with open(
            "tests/test_cases/find_most_likely_common_ancestor/find_most_likely_common_ancestor_test_cases.json",
            "r",
        ) as file:
            data = json.load(file)

            for key, value in data.items():
                input_fasta = data.get(key, {}).get("input_fasta")
                output_consensus = data.get(key, {}).get("output_consensus")
                output_profile = data.get(key, {}).get("output_profile")
                cases.append((input_fasta, (output_consensus, output_profile)))
        return cases


n = GCContentTestData()
n.get_cases()
