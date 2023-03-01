# Export_traj_to_mongo

#### Introduction

This python script is used for extracting the OPTIMADE fields from Molecular dynamics or Monte Carlo trajectories and to store them in an MongoDB database that is compatible with the optimade python tools.
At the moment this script is still under development and there is no guarantee that it works.
Currently, only the cartesion_site_positions are stored as an indexable field. 
The other properties are stored as a constant field.
The cartesian_site_positions, i.e. the particle positions are stored in an hdf5 file because that field would be too large to store in a mongo DB document.  
I do hope to add support for gridFS in the future.

The OPTIMADE standard can be found at [https://www.optimade.org](https://www.optimade.org)
The proposed changes to the OPTIMADE standard matching this version of the script can be found on: [https://github.com/JPBergsma/OPTIMADE/tree/Trajectory_proposal_v0.1](https://github.com/JPBergsma/OPTIMADE/tree/Trajectory_proposal_v0.1)
The version of the optimade python tools that matches this version of the proposal can be found here: [https://github.com/JPBergsma/optimade-python-tools/tree/optimade_python_tools_trajectory_0.1](https://github.com/JPBergsma/optimade-python-tools/tree/optimade_python_tools_trajectory_0.1)
The discussion on the modification of the standard to support trajectory data can be found here: [https://github.com/Materials-Consortia/OPTIMADE/pull/377](https://github.com/Materials-Consortia/OPTIMADE/pull/377)
The latest version of the proposal can be found here: [https://github.com/JPBergsma/OPTIMADE/tree/JPBergsma_add_Trajectories](https://github.com/JPBergsma/OPTIMADE/tree/JPBergsma_add_Trajectories)
The latest development version of the optimade python tools that are in development in parallel can be found here: https://github.com/JPBergsma/optimade-python-tools/tree/JPBergsma_add_trajectory


#### installation

This script can be cloned with:

`git clone https://github.com/JPBergsma/https://github.com/JPBergsma/Export_traj_to_mongo.git`

or installed as a library with:

`pip install git+https://github.com/JPBergsma/Export_traj_to_mongo@master`

#### Usage

This package will store the metadata of the trajectory in a mongoDB collection.
The name of this collection is read from the `"trajectories_collection"` field in the configuration file of the optimade-python-tools, i.e. `.optimade.json`
The default location for this configuration file is: `~/.optimade.json`
For more information about the configuration file and how to specify a different location, see the [configuration.md](https://github.com/JPBergsma/optimade-python-tools/blob/optimade_python_tools_trajectory_0.1/docs/configuration.md) file. 

To load a trajectory into the mongoDB collection you use the load_trajectory_data function.
The indexing for the first_frame, last_frame and reference_frame arguments is zero based. i.e. the first frame is frame 0. (Indexing in OPTIMADE is in contrast 1 based.)
It takes as arguments:
* structure_file: 
  * Description: A string or Path containing the location of the file with structural information about the compounds in the trajectory.
  If there is no separate trajectory file it can also contain trajectory information.
  It is also possible to use a stringIO object. In that case the file type should be specified with `file_format`.
  * Type: String, Path, StringIO
  * Optional: False
  
* trajectory_files: 
  * Description: a list of strings or paths containing the location of the files containing the trajectory.
  * Type: List of Strings or Paths
  * Optional: True

* references:
  * Description: A list with a separate dictionary for each reference. 
    These dictionaries can contain all the fields defined by [BibTeX](https://www.bibtex.com/format/).
  * Type: List of Dictionaries
  * Optional: True
  
* first_frame: 
  * Description: The first frame for which data should be stored. For example to exclude the equilibration part of the trajectory. 
  * Type: Integer
  * Optional: True
  * Default: 0

* frame_step: 
  * Description: Only 1 out of every frame_step frames will be shared.  
  * Type: Integer
  * Optional: True
  * Default: 1

* last_frame
  * Description: Indicates the last frame to be stored.
    It is exclusive, so if `last_frame = 100` the last frame that will be included in the database is the frame with index 99. 
  * Type: Integer
  * Optional: True
  * Default: The number of frames in the trajectory.

* storage_dir
  * Description: The location at which the HDF5 files containing the particle positions should be stored.
  * Type: Path or string
  * Optional: False

* traj_id:
  * Description: An id for this trajectory that is unique within the database. 
    If no id is provided the id from mongo DB will be used. 
  * Type: string
  * Optional: True
  * Default: The mongo DB ID

* custom_fields: A dictionary with extra database specific fields. 
  If a standard optimade field is within custom_fields the value from the trajectory file will be overwritten.
  * Type: Dictionary
  * Optional: True
  * Default: An empty dictionary

* reference_frame:
  * Description: The frame that should be used to generate the reference structure.
  * Type: integer
  * Optional: True
  * Default: the last_frame

An example on how it can be used is given in Demo_Loading_Trajectories_In_MongDB.py in the Demo folder.

The program uses the MDAnalysis package to read the trajectory files. 
The file types that it supports can be found here: [https://userguide.mdanalysis.org/stable/formats/index.html](https://userguide.mdanalysis.org/stable/formats/index.html)
