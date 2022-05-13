from exporttomongo import load_trajectory_data

trajectories = [#{"structurefile": "/home/kwibus/Documents/Cecam/testfiles/OBr2.xyz"},
                #{"structurefile": "/home/kwibus/Documents/Cecam/testfiles/OF2.xyz"},
                {"structurefile": "/home/kwibus/Documents/Cecam/testfiles/water_rotating.xyz"}
                #{"structurefile": "/home/kwibus/Documents/Cecam/testfiles/ala_small_traj.pdb"},
               #{"structurefile": "/home/kwibus/Documents/Cecam/testfiles/TRAJ.xyz", "last_frame": 99, "first_frame": 0},
                #{"structurefile": "/home/kwibus/Documents/Cecam/testfiles/1KX5.pdb"}
                 # {"structurefile": "/home/kwibus/Documents/Cecam/testfiles/systemstripped.pdb",
                 #  "trajectoryfiles": [
                 #       "/home/kwibus/Documents/Cecam/testfiles/sarscov2-10875754-no-water-zinc-glueCA-0000.dcd"],
                 #  "references": [{"title": "Creative Commons Attribution 4.0 International Public License",
                 #                 "url": "https://creativecommons.org/licenses/by/4.0/legalcode",
                 #                 "note": "Entries with a reference to this licence were and are shared under this licence.",
                 #                 "organization": "Creative Commons",
                 #                 "id": "CCBy4"},
                 #                {"title": "Molecular Dynamics Simulations Related to SARS-CoV-2",
                 #                 "organization": "D. E. Shaw Research",
                 #                 "year": "2020",
                 #                 "note": "Only a small fraction of trajectory DESRES-ANTON-10875754 is shared here as an example.",
                 #                 "url": "http://www.deshawresearch.com/resources_sarscov2.html",
                 #                 "id": "Shaw2020a"}], "last_frame": 100, "first_frame": 0}
               ]

for trajectory in trajectories:
    load_trajectory_data(trajectory["structurefile"],
                         trajectory.get("trajectoryfiles", None),
                         trajectory.get("references", None),
                         last_frame=trajectory.get("last_frame", None),
                         frame_step=1,
                         first_frame=trajectory.get("first_frame", 0),
                         storage_dir="/home/kwibus/Documents/Cecam/testfiles/")
