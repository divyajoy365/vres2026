from simulariumio import TrajectoryConverter, TrajectoryData, AgentData
import numpy as np
from simulariumio.physicell import PhysicellConverter, PhysicellData
from simulariumio import UnitData, MetaData, DisplayData, DISPLAY_TYPE, ModelMetaData

box_size = 1500.0

example_data = PhysicellData(
    timestep=6.0,
    path_to_output_dir="/mnt/c/VRES_2025/PhysiCell-1.14.2/output",
    meta_data=MetaData(
        box_size=np.array([box_size, box_size, box_size]),
        scale_factor=0.01,
        trajectory_title="Some parameter set",
        model_meta_data=ModelMetaData(
            title="Some agent-based model",
            version="8.1",
            authors="A Modeler",
            description=(
                "An agent-based model run with some parameter set"
            ),
            doi="10.1016/j.bpj.2016.02.002",
            source_code_url="https://github.com/simularium/simulariumio",
            source_code_license_url="https://github.com/simularium/simulariumio/blob/main/LICENSE",
            input_data_url="https://allencell.org/path/to/native/engine/input/files",
            raw_output_data_url="https://allencell.org/path/to/native/engine/output/files",
        ),
    ),
    nth_timestep_to_read=1,
    display_data={
        0: DisplayData(
            name="T cell",
            display_type=DISPLAY_TYPE.OBJ,
            url="tumor_cell.obj",
            color="#FFD700",
        ),
        1: DisplayData(   # Rod cells
        name="Rod cell",
        display_type=DISPLAY_TYPE.OBJ,
        url="tumor_cell.obj",
        color="#1E8F3A",  # green #1E8F3A

    ),
    },
    phase_names={
        0: {4: "G0G1"},
        1: {4: "G0G1"},
    },
    time_units=UnitData("s"),  # seconds
)

PhysicellConverter(example_data).save("Check_runs/base4.gif")