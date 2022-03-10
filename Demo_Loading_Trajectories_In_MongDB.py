from exporttomongo import load_trajectory_data

trajectories = [{"structurefile": "/home/kwibus/Documents/Cecam/testfiles/OBr2.xyz"},
                {"structurefile": "/home/kwibus/Documents/Cecam/testfiles/OF2.xyz"},
                {"structurefile": "/home/kwibus/Documents/Cecam/testfiles/water_rotating.xyz"}]
#{"structurefile": "/home/kwibus/Documents/Cecam/testfiles/ala_small_traj.pdb"},
                #{"structurefile": "/home/kwibus/Documents/Cecam/testfiles/TRAJshort.xyz"},
                # {"structurefile": "/home/kwibus/Documents/Cecam/testfiles/systemstripped.pdb",
                #  "trajectoryfiles": [
                #      "/home/kwibus/Documents/Cecam/testfiles/sarscov2-10875754-no-water-zinc-glueCA-0000.dcd"]}
               #]

for trajectory in trajectories:
    load_trajectory_data(trajectory["structurefile"],
                         trajectory.get("trajectoryfiles", None),
                         storage_dir="/home/kwibus/Documents/Cecam/testfiles/")
