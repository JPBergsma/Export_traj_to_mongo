from exporttomongo import load_trajectory_data
from pathlib import Path

trajectories = [{"structurefile": Path(__file__).parent/"TRAJshort.xyz",
                 "references":[{"Author": "Johan Bergsma", "year": "2010", "note": "Any further remarks", "id": "Bergsma2010"}]},
               ]

for trajectory in trajectories:
    load_trajectory_data(trajectory["structurefile"],
                         trajectory.get("trajectoryfiles", None),
                         trajectory.get("references", None),
                         last_frame=trajectory.get("last_frame", None),
                         frame_step=1,
                         first_frame=trajectory.get("first_frame", 0),
                         storage_dir=Path(__file__).parent/"hdf5_files/")
