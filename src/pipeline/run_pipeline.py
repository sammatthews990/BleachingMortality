#!/usr/bin/env python3
'''Readable runner for the ordered commands in config/pipeline.toml.'''

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError as exc:  # pragma: no cover
    raise SystemExit('Python 3.11 or newer is required.') from exc


def project_root(start: Path | None = None) -> Path:
    candidate = (start or Path.cwd()).resolve()
    for directory in (candidate, *candidate.parents):
        if (directory / 'config' / 'pipeline.toml').is_file():
            return directory
    raise FileNotFoundError('Could not find config/pipeline.toml.')


def load_pipeline(root: Path) -> dict:
    with (root / 'config' / 'pipeline.toml').open('rb') as handle:
        return tomllib.load(handle)


def stage_index(config: dict) -> dict[str, dict]:
    stages = config.get('stages', [])
    indexed = {stage['name']: stage for stage in stages}
    if len(indexed) != len(stages):
        raise ValueError('Pipeline stage names must be unique.')
    return indexed


def select_stages(config: dict, profile: str | None, names: list[str]) -> list[dict]:
    indexed = stage_index(config)
    selected_names = list(names)
    if profile:
        profiles = config.get('profiles', {})
        if profile not in profiles:
            raise KeyError(f'Unknown profile: {profile}')
        selected_names = list(profiles[profile]) + selected_names
    if not selected_names:
        raise ValueError('Choose --profile, --stage or --list.')
    ordered = []
    for name in selected_names:
        if name not in indexed:
            raise KeyError(f'Unknown stage: {name}')
        if name not in {stage['name'] for stage in ordered}:
            ordered.append(indexed[name])
    return ordered


def print_pipeline(config: dict) -> None:
    print('Profiles:')
    for name, stages in config.get('profiles', {}).items():
        route = ' -> '.join(stages)
        print(f'  {name}: {route}')
    print('\nStages:')
    for stage in config.get('stages', []):
        name = stage['name']
        description = stage['description']
        print(f'  {name}: {description}')
        for command in stage.get('commands', []):
            print(f'    {command}')


def missing_inputs(root: Path, stages: list[dict]) -> list[tuple[str, str]]:
    missing = []
    for stage in stages:
        for relative in stage.get('inputs', []):
            if not (root / relative).exists():
                missing.append((stage['name'], relative))
    return missing


def run(stages: list[dict], root: Path, dry_run: bool) -> None:
    for stage in stages:
        name = stage['name']
        description = stage['description']
        print(f'\n[{name}] {description}', flush=True)
        for command in stage.get('commands', []):
            print(f'> {command}', flush=True)
            if not dry_run:
                subprocess.run(command, cwd=root, shell=True, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--list', action='store_true', help='Print profiles and commands.')
    parser.add_argument('--profile', help='Run a named ordered profile.')
    parser.add_argument('--stage', action='append', default=[], help='Run one stage; repeatable.')
    parser.add_argument('--check-inputs', action='store_true', help='Report missing declared inputs and exit.')
    parser.add_argument('--dry-run', action='store_true', help='Print commands without executing them.')
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = project_root()
    config = load_pipeline(root)
    if args.list:
        print_pipeline(config)
        return 0
    try:
        stages = select_stages(config, args.profile, args.stage)
    except (KeyError, ValueError) as exc:
        print(exc, file=sys.stderr)
        return 2
    if args.check_inputs:
        missing = missing_inputs(root, stages)
        if missing:
            print('Missing declared inputs:')
            for stage, path in missing:
                print(f'  [{stage}] {path}')
            return 1
        print('All declared inputs are present.')
        return 0
    run(stages, root, args.dry_run)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
