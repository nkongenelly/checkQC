from pathlib import Path
import json
import pytest

from checkQC.app import App
from checkQC.app import run_new_checkqc


@pytest.fixture
def bcl2fastq_runfolder_path():
    return str(Path(__file__).parent / "resources/bcl2fastq/170726_D00118_0303_BCB1TVANXX/")


@pytest.fixture
def bclconvert_runfolder_path():
    return str(Path(__file__).parent / "resources/bclconvert/200624_A00834_0183_BHMTFYTINY/")


def test_run(bcl2fastq_runfolder_path):
    app = App(runfolder=bcl2fastq_runfolder_path)
    # The test data contains fatal qc errors
    assert app.run() == 1


def test_run_json_mode(bcl2fastq_runfolder_path):
    app = App(runfolder=bcl2fastq_runfolder_path, json_mode=True)
    # The test data contains fatal qc errors
    assert app.run() == 1


def test_run_use_closest_read_length(bcl2fastq_runfolder_path):
    config_file = Path("tests") / "resources/read_length_not_in_config.yaml"
    app = App(
        runfolder=bcl2fastq_runfolder_path,
        config_file=config_file,
        use_closest_read_length=True,
    )
    # The test data contains fatal qc errors
    assert app.run() == 1


def test_run_downgrade_error(bcl2fastq_runfolder_path):
    app = App(
        runfolder=bcl2fastq_runfolder_path,
        downgrade_errors_for="ReadsPerSampleHandler",
    )
    # Test data should not produce fatal qc errors anymore
    assert app.run() == 0


def test_run_new_checkqc(bclconvert_runfolder_path):
    exit_status, reports = run_new_checkqc(
        None,
        bclconvert_runfolder_path,
        downgrade_errors_for="ReadsPerSampleHandler",
        use_closest_read_length=True,
        demultiplexer="bclconvert",
    )

    expected_data = (
        Path(bclconvert_runfolder_path) / "expected_qc.json").read_text()

    assert exit_status == 1
    assert reports + '\n' == expected_data


def test_bclconvert_and_bcl2fastq_qc(bclconvert_runfolder_path):
    bclconvert_path = Path(bclconvert_runfolder_path) / "expected_qc.json"
    bcl2fastq_tiny_runfolder = (
        Path(__file__).parent
        / "resources/bcl2fastq/200624_A00834_0183_BHMTFYTINY"
    )
    bclconvert_data = json.loads(bclconvert_path.read_text())
    bcl2fastq_data = json.loads(
        (bcl2fastq_tiny_runfolder / "expected_qc.json").read_text()
    )

    qc_title_conversion = {
        "reads_per_sample": "ReadsPerSampleHandler",
        "undetermined_percentage": "UndeterminedPercentageHandler",
    }

    bcl2fastq_index = {} # format: {handler_name: {lane_sample_key: data}}
    for handler_name, results in bcl2fastq_data.items():
        if handler_name in {"exit_status", "run_summary"}:
            continue

        lane_map = {}
        for result in results:
            data = result["data"]
            lane = int(data["lane"])
            if "sample_name" in data:
                lane_map[(lane, data["sample_name"])] = data
            else:
                lane_map[lane] = data
        bcl2fastq_index[handler_name] = lane_map

    for lane_no, lane_data in bclconvert_data["lane_reports"].items():
        lane = int(lane_no)
        for handler_name, handler_results in lane_data.items():
            if handler_name == "error_rate": # not in our test bcl2fastq results
                continue

            bcl2fastq_handler = qc_title_conversion[handler_name]
            lookup = bcl2fastq_index[bcl2fastq_handler]

            for result in handler_results:
                result_data = result["data"]
                if handler_name == "reads_per_sample":
                    expected_data = lookup[(lane, result_data["sample_id"])]
                    assert expected_data["sample_reads"] == result_data["sample_reads"]
                if handler_name == "undetermined_percentage":
                    expected_data = lookup[lane]
                    assert expected_data["percentage_undetermined"] == result_data["percentage_undetermined"]

