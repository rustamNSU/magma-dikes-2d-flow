import os
import sys
import json
import runpy
from pathlib import Path

# Set repository and simulation paths
repository_dir = Path.cwd()
simulations_dir = repository_dir / "simulations"
sys.path.append(str(repository_dir))

# Find all sim<N>.py files in the current directory
current_dir = Path(__file__).parent
sim_files = sorted(current_dir.glob("sim*.py"), key=lambda p: int(p.stem[3:]))

input_paths = []

# Run each simulation script and collect input.json paths
for sim_file in sim_files:
    print(f"Running {sim_file.name}")
    runpy.run_path(str(sim_file))  # Execute the simulation script

    sim_id = sim_file.stem[3:]  # extract <N> from sim<N>.py
    sim_dir = simulations_dir / f"simID{sim_id}"
    input_json_path = sim_dir / "input.json"

    if input_json_path.exists():
        input_paths.append(str(input_json_path.resolve()))
    else:
        print(f"Warning: {input_json_path} does not exist")

# Write runner.json file with all input paths
runner_json_path = repository_dir / "runner.json"
with open(runner_json_path, "w") as f:
    json.dump({"input_paths": input_paths}, f, indent=4)
print(f"\nrunner.json written with {len(input_paths)} entries")
