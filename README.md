# Export_traj_to_mongo

This python script is used for reading Molecular dynamics trajectories and extracting the OPTIMADE fields to store it in an OPTIMADE database that is complient with the optimade python tools .
The cartesian_site_positions, i.e. the particle positions are stored in a hdf5 file because that field would be to large to store in a mongo DB document.  
The OPTIMADE standard can be found at https://www.optimade.org/
The proposed modification of the standard to support trajectory data can be found here: https://github.com/Materials-Consortia/OPTIMADE/pull/377
The optimade python tools that are in devellopment in parrallel can be found here: https://github.com/Materials-Consortia/optimade-python-tools/pull/1065
This package depends on the Optimade python tools and will store the data in the mongoDB trajectories_collection which is defined in the config file of the optimade python tools. 
At the moment this script is still under devellopment and there is no guarantee that it works.

An example on how it can be used is given in Demo_Loading_Trajectories_In_MongDB.py


load_trajectory_data(structure_file,
                         trajectory.get("trajectoryfiles", None),
                         trajectory.get("references", None),
                         last_frame=trajectory.get("last_frame", None),
                         frame_step=1,
                         first_frame=trajectory.get("first_frame", 0),
                         storage_dir="/home/kwibus/Documents/Cecam/testfiles/")
    arguments: structure_file: string: A file path to a file containing the structural information about the compounds in the trajectory.
                               If there is no separate trajectory file it can also contain trajectory information.
               trajectory_files: list of strings [Optional]: A list of paths to files containing the trajectory.
               references: dictionary [Optional]: A list of dictionaries containing references belonging to this trajectory.
                           Valid fields are all the fields defined by the bibtech standard.
                           Each referene should have an id field(string) that is unique within the database.
               first_frame: integer [Optional]: In case only a part of the trajectory should be used this indicates the first frame that should be stored in the data base. Default = 1.
               frame_step: integer [Optional]: Only 1 out of every frame step will be stored on the server. Default = 1
               last_frame: integer [Optional]: The last frame which should be stored in the database. The default value is equal to the number of frames in the trajectory.
               !Warning: The slicing parameter first_frame, frame_step and last_frame have not tested this slicing properly yet.               
               traj_id: string [Optional]: An id for this trajectory that is unique within the database. If no id is provided the id from mongo DB will be used.
               reference_frame: integer [Optional]: This indicates which frame will be used to generate the reference structure. Default = 1.

