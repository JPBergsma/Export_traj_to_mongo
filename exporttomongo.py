import re
import MDAnalysis as mda
import pymatgen.core
import datetime
import h5py
import warnings
import sys

sys.path.append("/home/kwibus/PycharmProjects/Optimade/optimade-python-tools/")
from optimade.server.config import CONFIG
from optimade.server.entry_collections import create_collection
from optimade.models.trajectories import TrajectoryResource
from optimade.models.references import ReferenceResource
from optimade.models.utils import ANONYMOUS_ELEMENTS
from optimade.server.mappers import TrajectoryMapper, ReferenceMapper


def load_trajectory_data(
        structure_file,
        trajectory_files=None,
        references=None,
        first_frame=0,
        last_frame=None,
        storage_dir=None,
        frame_step=1,
        id=None,
        reference_frame=0,
):
    # step 1: load file(s) with trajectory data
    # TODO It is probably better to use paths from the Path library instead of strings for the filenames
    # TODO Determine whether the file is small enough to store in memory in that case set "in_memory=True"
    if trajectory_files:
        traj = mda.Universe(structure_file, trajectory_files)
    else:
        traj = mda.Universe(structure_file)

    # Step2 generate pymatgen structure from MDAnalysis structure so we can also use the methods of pymatgen
    struct = generate_pymatgen_from_mdtraj(traj, reference_frame)

    # step 3: generate all the neccesary OPTIMADE fields from the data
    # TODO add automatically reading the units from the data if present

    if not last_frame:
        last_frame = len(traj.trajectory)
    n_frames = 1 + (((last_frame-1) - first_frame) // frame_step)

    # TODO it would be nice to also allow adding structures in the same way we add trajectories
    # type
    #if n_frames > 1:
    type = "trajectories"
    #else:
    #    type = "structures"

    # immutable_id (OPTIONAL)
    # last_modified handled by the function current_time()
    # elements
    elements = sorted(struct.symbol_set)
    # nelements
    nelements = len(elements)
    # elements_ratios
    elements_ratios = [
        float(i)
        for i in re.sub(
            "[A-Z][a-z]*",
            "",
            struct.composition.fractional_composition.alphabetical_formula,
        ).split()
    ]
    # chemical_formula_descriptive
    chemical_formula_descriptive = struct.composition.alphabetical_formula.replace(
        " ", ""
    )
    # chemical_formula_reduced and chemical_formula_anonymous
    chemical_formula_reduced, chemical_formula_anonymous = get_formula_reduced_and_anonymous(struct)

    # chemical_formula_hill(OPTIONAL) Not yet implemented

    # dimension_types
    dimension_types = None  # TODO: The files I use for testing do not explicitly store this information. Most likely it is periodic but to be sure I should probably check whether: 1. particles move through the periodic boundaries 2. There are chemical bonds that go through the periodic boundary 3. particles that are outside the unitcell

    # nperiodic_dimensions
    if dimension_types:
        nperiodic_dimensions = sum(dimension_types)
    else:
        nperiodic_dimensions = None

    # lattice_vectors
    if hasattr(traj.trajectory[reference_frame], "dimensions"):
        lattice_vectors = struct.lattice.matrix.tolist()
    else:
        lattice_vectors = None

    # cartesian_site_positions
    cartesian_site_positions = struct.cart_coords.tolist()

    # nsites
    nsites = traj.atoms.n_atoms

    # species_at_sites
    species_at_sites = traj.atoms.names.tolist()

    # species  # TODO the atom names/labels may not be unique enough in some cases. In that case extra descriptors such as the number of attached hydrogens or the charge have to be added.
    species = []
    for specie in set(species_at_sites):
        index = species_at_sites.index(specie)
        species.append(
            {
                "name": specie,
                "chemical_symbols": [traj.atoms.elements[index]],
                "concentration": [getoccu(traj, index)],
                "mass": [traj.atoms.masses[index]],
                "_biomol_atom_name": specie
            }
        )

    # assemblies(OPTIONAL)

    # structure_features
    # TODO make more checks to see which properties should be set here.
    structure_features = []
    for specie in species:
        if len(specie["chemical_symbols"]) > 1:
            structure_features.append("disorder")
            break

    available_properties = {
        "cartesian_site_positions": {
            "frame_serialization_format": "explicit",
            "nvalues": n_frames,
            # TODO find examples of trajectories where the number of coordinates sets and the number of frames does not match.
        },
        "lattice_vectors": {"frame_serialization_format": "constant"},
        "species": {"frame_serialization_format": "constant"},
        "dimension_types": {"frame_serialization_format": "constant"},
        "species_at_sites": {"frame_serialization_format": "constant"},
    }

    # MDAnalysis throws a warning when the timestep is not specified but does return a reasonable value of 1.0 ps. This can be confusing, so I therefore choose to catch this warning.
    # We do not want this warning to be displayed to the user so we temporarely allow only errors to be reported.
    warnings.filterwarnings("error")
    try:
        dt = traj.trajectory[0].dt * frame_step
        time_present = True
    except UserWarning:
        time_present = False
    warnings.filterwarnings("default")

    # TODO find examples of trajectories where the number of coordinates sets and the number of frames does not match.

    if time_present:  # if the time step is not zero or none
        available_properties["time"] = {"frame_serialization_format": "linear"}

    reference_structure = {
        "elements": elements,
        "nelements": nelements,
        "elements_ratios": elements_ratios,
        "chemical_formula_descriptive": chemical_formula_descriptive,
        "chemical_formula_reduced": chemical_formula_reduced,
        "chemical_formula_anonymous": chemical_formula_anonymous,
        "dimension_types": dimension_types,
        "nperiodic_dimensions": nperiodic_dimensions,
        "lattice_vectors": lattice_vectors,
        "cartesian_site_positions": cartesian_site_positions,
        "nsites": nsites,
        "species_at_sites": species_at_sites,
        "species": species,
        "structure_features": structure_features,
    }

    entry = {
        "reference_structure": reference_structure,
        "nframes": n_frames,
        "reference_frame": reference_frame,
        "available_properties": available_properties,
        "type": type,
        "last_modified": last_modified(),
        "cartesian_site_positions": {
            "frame_serialization_format": "explicit",
            "nvalues": n_frames,
            "_storage_location": "file",
        },
        "lattice_vectors": {
            "frame_serialization_format": "constant",
            "_storage_location": "mongo",
            "values": [lattice_vectors],
        },
        "species": {
            "frame_serialization_format": "constant",
            "_storage_location": "mongo",
            "values": [species],
        },
        "dimension_types": {
            "frame_serialization_format": "constant",
            "_storage_location": "mongo",
            "values": [dimension_types],
        },
        "species_at_sites": {
            "frame_serialization_format": "constant",
            "_storage_location": "mongo",
            "values": [species_at_sites],
        },
    }
    if time_present:
        entry["_exmpl_time"] = {
            "_storage_location": "mongo",
            "frame_serialization_format": "linear",
            "offset_linear": traj.trajectory[0].time,
            "step_size_linear": dt,
        }

    # step 3.5 add references
    if references:
        entry["relationships"] = generate_relationships(references)

    # Generate biomolecular fields:

    if hasattr(traj.residues, "icodes"): # TODO need a more thorough way to determine if it is a biomolecular simulation.
        sequences, residues, chains = get_biomol_fields(traj)
        entry["sequences"] = sequences
        entry["residues"] = residues
        entry["chains"] = chains


    # Step 4: store trajectory data

    trajectories_coll = create_collection(
        name=CONFIG.trajectories_collection,
        resource_cls=TrajectoryResource,
        resource_mapper=TrajectoryMapper,
    )
    mongoid = trajectories_coll.insert([entry]).inserted_ids[0]
    if not id:
        id = str(mongoid)
    fields_to_add = {"id": id}

    if (type == "trajectories"):
        # Write trajectory data in HDF5 format
        if n_frames * nsites * 3 * 4 > 16 * 1024:  # If the trajectory is larger than about 16 kb store it in hdf5 file
            hdf5path = storage_dir + id + ".hdf5"
            fields_to_add["_hdf5file_path"] = hdf5path

            # TODO It would be better to use a try and except around storing the data. If writing the data to the hdf5 file failes the corresponding entry should be removed from the mongo DB.
            with h5py.File(hdf5path, "w") as hdf5file:
                #TODO It would be nice as we could store all the trajectory data in the HDF5 file So we should still add the storing of the other relevant trajectory info here as well.

                if traj.trajectory[
                    reference_frame
                ].has_positions:  # TODO check whether there are more properties that can be stored such as force and velocities
                    arr = hdf5file.create_dataset(
                        "cartesian_site_positions/values",
                        (n_frames, nsites, 3),
                        chunks=True,
                        dtype=traj.trajectory[0].positions[0][0].dtype,
                    )  # TODO allow for varying number of particles
                    for i in range(first_frame, last_frame,
                                   frame_step):  # TODO offer the option to place only a part of the frames in the database
                        arr[(i - first_frame) // frame_step] = traj.trajectory[i].positions

        else:  # If the trajectory is small it can be stored locally
            positions = []

            if "F" in elements:  # To cover all the different cases for testing I encode the information about the trajectory in a different ways here.
                import random
                frames = []
                for i in range(n_frames):
                    if random.randrange(0, 10) >= 5:
                        frames.append(i)
                        positions.append(traj.trajectory[i].positions.tolist())
                fields_to_add.update({"available_properties.cartesian_site_positions.frame_serialization_format": "explicit_custom_sparse",
                    "available_properties.cartesian_site_positions.nvalues": len(frames),
                    "cartesian_site_positions.frames": frames,
                    "cartesian_site_positions.frame_serialization_format": "explicit_custom_sparse",
                    "cartesian_site_positions.nvalues": len(frames)})

            elif "Br" in elements:
                for i in range(0, n_frames, 2):
                    positions.append(traj.trajectory[i].positions.tolist())
                fields_to_add.update({"cartesian_site_positions.step_size_sparse": 2,
                    "cartesian_site_positions.offset_sparse": 0,
                    "cartesian_site_positions.frame_serialization_format": "explicit_regular_sparse",
                    "cartesian_site_positions.nvalues": len(positions),
                    "available_properties.cartesian_site_positions.frame_serialization_format": "explicit_regular_sparse",
                    "available_properties.cartesian_site_positions.nvalues": len(positions)})

            else:
                for i in range(first_frame, last_frame, frame_step):
                    positions.append(traj.trajectory[i].positions.tolist())
            fields_to_add.update({"cartesian_site_positions._storage_location": "mongo",
                                            "cartesian_site_positions.values": positions})
    trajectories_coll.collection.update_one(
        {"_id": mongoid}, {"$set": fields_to_add}
    )


def flip_chem_form_anon(chemical_formula_anonymous: str) -> str:
    """Converts an anonymous chemical formula with the most numerous element in the last position to an
    anonymous chemical formula with the most numerous element in the first position and vice versa."""
    numbers = [n for n in re.split(r"[A-Z][a-z]*", chemical_formula_anonymous)]
    anon_elem = [
        n
        for n in re.findall(
            "[A-Z][^A-Z]*", re.sub("[0-9]+", "", chemical_formula_anonymous)
        )
    ]
    return "".join(
        [
            y
            for x in range(len(anon_elem))
            for y in [anon_elem[x], numbers[len(anon_elem) - x]]
        ]
    )


def getoccu(traj, index):
    if hasattr(traj.atoms, "occupancies"):
        return traj.atoms.occupancies[index]
    else:
        return 1.0


def generate_relationships(references):
    references_coll = create_collection(
        name=CONFIG.references_collection,
        resource_cls=ReferenceResource,
        resource_mapper=ReferenceMapper,
    )
    list_references = []
    for reference in references:
        list_references.append({"type": "references", "id": reference["id"]})
        reference["last_modified"] = last_modified()
        references_coll.insert([
            reference])  # TODO check whether id is unique if it is already in the data base it would need a postfix if not identical.

    return {"references": {"data": list_references}}


def last_modified():
    return datetime.datetime.utcnow().replace(
        microsecond=0
    )  # MongeDB does not accept microseconds


def get_res_type(resname):
    if resname in AMINOACID_DICT:
        return "amino_acid"
    if resname in DNA_DICT:
        return "DNA"
    if resname in RNA_LIST:
        return "RNA"
    if resname in SOLVENTS:
        return "solvent"
    if resname in IONS:
        return "ion"
    else:
        return "other"


def get_chain_indixes(seg_id_set, chain_id_to_index):
    chains = []
    for chain_id in seg_id_set:
        if chain_id not in chain_id_to_index:
            chain_id_to_index[chain_id] = len(chain_id_to_index)
        chains.append(chain_id_to_index[chain_id])
    return chains, chain_id_to_index


def generate_pymatgen_from_mdtraj(traj, frame):
    """Generates a pymatgen structure object from a frame of a MDanalysis trajectory object.

    arguments: traj MDanalysis trajectory object.
               frame  the frame of the trajectory for which the pymatgen structure object should be generated.
    """

    nsites = traj.atoms.n_atoms
    if traj.trajectory[frame].dimensions is not None:
        boxdim = traj.trajectory[frame].dimensions[:3]
        dictstruct = {
            "lattice": {
                "a": boxdim[0],
                "b": boxdim[1],
                "c": boxdim[2],
                "alpha": traj.trajectory[frame].dimensions[3],
                "beta": traj.trajectory[frame].dimensions[4],
                "gamma": traj.trajectory[frame].dimensions[5],
            }
        }
    else:
        boxdim = [1.0, 1.0, 1.0]
        dictstruct = {
            "lattice": {
                "a": boxdim[0],
                "b": boxdim[1],
                "c": boxdim[2],
                "alpha": 90.0,
                "beta": 90.0,
                "gamma": 90.0,
            }
        }

    # dictstruct["charge"] = TODO MDTRAJ does not read charges from pdb files yet so we cannot implement this yet.
    sites = []

    single_frame = (len(traj.trajectory) == 1)
    if single_frame:  # somehow reading a structure as a trajectory is very slow so instead we read the positions from the atoms.
        for i in range(nsites):
            site = {
                "species": [{"element": traj.atoms.elements[i], "occu": getoccu(traj, i)}],
                "abc": traj.atoms.positions[i] / boxdim[:3],
            }
            sites.append(site)
    else:
        for i in range(nsites):
            site = {
                "species": [{"element": traj.atoms.elements[i], "occu": getoccu(traj, i)}],
                "abc": traj.trajectory[frame].positions[i] / boxdim[:3],
            }
            sites.append(site)

    dictstruct["sites"] = sites
    return pymatgen.core.structure.IStructure.from_dict(dictstruct)


def get_formula_reduced_and_anonymous(struct):

    def sort_second(val):
        return val[1]

    # TODO find a better method for rounding down the elemental composition. if one element is
    el = struct.composition.alphabetical_formula
    formula_components = re.split("([A-Z][a-z]?)", el)
    formula_pairs = []
    for i in range(1, len(formula_components), 2):
        if formula_components[i + 1] == '':
            formula_components[i + 1] = '1'
        formula_pairs.append((formula_components[i], round(float(formula_components[i + 1]))))
    formula_reduced = "".join(["".join([pair[0], re.sub("^[1]$", "", str(pair[1]))]) for pair in formula_pairs])
    formula_pairs.sort(key=sort_second, reverse=True)
    formula_anonymous = "".join(
        ["".join([ANONYMOUS_ELEMENTS[i], re.sub("^[1]$", "", str(pair[1]))]) for i, pair in enumerate(formula_pairs)])

    return formula_reduced, formula_anonymous


def get_biomol_fields(traj):
    biomol_residues = []
    biomol_chains = []
    chain_id_to_index = {}
    for res in traj.residues:
        if res.icode == '':
            insertion_code = None
        else:
            insertion_code = res.icode
        residue_dict = {"name": res.resname, "number": int(res.resnum), "insertion_code": insertion_code, "sites": res.atoms.indices.tolist()}
        biomol_residues.append(residue_dict)

        # MDanalysis does not have a datatype equivalent to the chain in a pdb file so we have to construct it based on the information in the residues and the segments
        chain_index, chain_id_to_index = get_chain_indixes([res.segid], chain_id_to_index)
        type = get_res_type(res.resname)
        if chain_index[0] >= len(biomol_chains):  # a new chain has been found
            chain_dict = {"name": res.segid, "residues": [int(res.resindex)], "types": [type], "sequences": [],
                          "sequence_types": []}
            biomol_chains.append(chain_dict)
        else:  # case existing chain
            biomol_chains[chain_index[0]]["residues"].append(int(res.resindex))
            if type not in biomol_chains[chain_index[0]]["types"]:
                biomol_chains[chain_index[0]]["types"].append(type)

    # _biomol_sequences

    biomol_sequences = []

    for seg in traj.segments:
        residue_name_list = []
        residue_index_list = []
        seg_id_set = set()
        current_type = "other"
        for res in seg.residues:
            res_type = get_res_type(res.resname)
            if (current_type in MONOMER_TYPES) and res_type != current_type:
                # reached the end of a sequence so make a new sequence
                sequence = ''.join(residue_name_list)
                chains, chain_id_to_index = get_chain_indixes(seg_id_set, chain_id_to_index)
                biomol_sequences.append({"sequence": sequence, "type": current_type, "chains": chains,
                                         "residues": residue_index_list})  # TODO In the specifictaion written by Dani https://github.com/Materials-Consortia/OPTIMADE/pull/400/files a se
                current_type = res_type
                residue_index_list = []
                residue_name_list = []
                seg_id_set = set()
            if res_type in MONOMER_TYPES:  # No sequence is generated for segment types other than RNA, DNA and aminoacids
                # sequence continues
                if res_type == "DNA":
                    resname = DNA_DICT[res.resname]
                elif res_type == "amino_acid":
                    resname = AMINOACID_DICT[res.resname]
                else:  # case RNA
                    resname = res.resname
                residue_index_list.append(int(res.resindex))
                residue_name_list.append(resname)
                seg_id_set.add(res.segid)
                current_type = res_type
        # finally at the end of the sequence store the final sequence
        if current_type in MONOMER_TYPES:
            sequence = ''.join(residue_name_list)
            chains, chain_id_to_index = get_chain_indixes(seg_id_set, chain_id_to_index)
            biomol_sequences.append(
                {"sequence": sequence, "type": current_type, "chains": chains, "residues": residue_index_list})

    # Add sequence info to chain
    for sequence in biomol_sequences:
        for chain_index in sequence["chains"]:
            biomol_chains[chain_index]["sequences"].append(sequence["sequence"])
            biomol_chains[chain_index]["sequence_types"].append(sequence["type"])

    return biomol_sequences, biomol_residues, biomol_chains


SOLVENTS = ["HOH"]
IONS = ["CL", "MN"]
MONOMER_TYPES = ("RNA", "DNA", "amino_acid")
DNA_DICT = {"DA": "A", "DC": "C", "DG": "G", "DT": "T", "DI": "I"}
RNA_LIST = ["C", "G", "A", "U", "I"]
AMINOACID_DICT = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "ASX": "B", "CYS": "C",
                  "GLU": "E", "GLN": "Q", "GLX": "Z", "GLY": "G", "HIS": "H", "ILE": "I",
                  "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
                  "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O",
                  "XLE": "J", "XAA": "X"}
