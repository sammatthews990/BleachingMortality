import unittest
from pathlib import Path

from src.pipeline.run_pipeline import (
    load_pipeline,
    project_root,
    select_stages,
    stage_index,
)


class PipelineRunnerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = project_root(Path(__file__).parent)
        cls.config = load_pipeline(cls.root)

    def test_stage_names_are_unique(self):
        indexed = stage_index(self.config)
        self.assertEqual(len(indexed), len(self.config['stages']))

    def test_current_profile_has_selected_model_before_diagnostics(self):
        stages = select_stages(self.config, 'current', [])
        names = [stage['name'] for stage in stages]
        self.assertLess(names.index('selected-model'), names.index('diagnostics'))
        self.assertLess(names.index('diagnostics'), names.index('current-reports'))

    def test_all_commands_point_to_existing_sources(self):
        for stage in self.config['stages']:
            for command in stage.get('commands', []):
                parts = command.split()
                if parts[0] in {'python', 'Rscript'}:
                    self.assertTrue((self.root / parts[1]).is_file(), command)


if __name__ == '__main__':
    unittest.main()
