from exporttomongo import load_trajectory_data
from pathlib import Path

trajectories = [{"structurefile": Path(__file__).parent / "TRAJshort.xyz",
                 "references": [
                     {"Author": "Johan Bergsma", "year": "2010", "note": "Any further remarks", "id": "Bergsma2010"}],
                 "custom_fields": {"program": "CPMD",
                                   "dimension_types": [1, 1, 1],
                                   "nperiodic_dimensions": 3},
                 "first_frame":2,
                 "last_frame":10},
                {"structurefile": Path(__file__).parent / "output.pdb",
                 "trajectoryfiles": [Path(__file__).parent / "output.dcd"]}
                ]

for trajectory in trajectories:
    load_trajectory_data(trajectory["structurefile"],
                         trajectory.get("trajectoryfiles", None),
                         trajectory.get("references", None),
                         last_frame=trajectory.get("last_frame", None),
                         frame_step=trajectory.get("frame_step", 1),
                         first_frame=trajectory.get("first_frame", 0),
                         custom_fields=trajectory.get("custom_fields", {}))
