from pathlib import Path

import pytest

from checkQC.parsers.illumina import (
    _read_interop_summary,
    _read_demultiplexing_metrics,
    _read_run_metadata,
    _read_samplesheet,
)


@pytest.fixture
def runfolder_path():
    return (
        Path(__file__).parent.parent
        / "resources/bclconvert/200624_A00834_0183_BHMTFYTINY"
    )


def test_read_interop_summary(runfolder_path):
    run_summary, index_summary, _ = _read_interop_summary(runfolder_path)

    total_reads_pf = run_summary.at(0).at(0).reads_pf()
    assert total_reads_pf == 532464327

    sample_id = index_summary.at(0).at(0).sample_id()
    assert sample_id == "14574-qiagen-indexset1-sp-lane1"


def test_read_quality_metrics(runfolder_path):
    quality_metrics = _read_demultiplexing_metrics(
            runfolder_path / "Reports/Quality_Metrics.csv")

    assert len(quality_metrics) == 6
    assert quality_metrics[0] == {
        '% Q30': '0.96',
        'Lane': '1',
        'Mean Quality Score (PF)': '36.43',
        'QualityScoreSum': '2066796745',
        'ReadNumber': '1',
        'SampleID': '14574-qiagen-indexset1-sp-lane1',
        # 'Sample_Project': 'AB-1234',
        'Yield': '56739672',
        'YieldQ30': '54699221',
        'index': 'GAACTGAGCG',
        'index2': 'TCGTGGAGCG',
    }


def test_read_to_unknown_barcodes(runfolder_path):
    top_unknown_barcodes = _read_demultiplexing_metrics(
            runfolder_path / "Reports/Top_Unknown_Barcodes.csv")

    assert len(top_unknown_barcodes) == 2000
    assert top_unknown_barcodes[:3] == [
        {
            '# Reads': '2008558',
            '% of All Reads': '0.003772',
            '% of Unknown Barcodes': '0.003787',
            'Lane': '1',
            'index': 'ATATCTGCTT',
            'index2': 'TAGACAATCT',
        },
        {
            '# Reads': '1934726',
            '% of All Reads': '0.003634',
            '% of Unknown Barcodes': '0.003648',
            'Lane': '1',
            'index': 'CACCTCTCTT',
            'index2': 'CTCGACTCCT',
        },
        {
            '# Reads': '1905527',
            '% of All Reads': '0.003579',
            '% of Unknown Barcodes': '0.003593',
            'Lane': '1',
            'index': 'ATGTAACGTT',
            'index2': 'ACGATTGCTG',
        },
    ]


# TODO add tests with paired end reads
def test_read_run_metadata(runfolder_path):
    instrument, read_length = _read_run_metadata(runfolder_path)
    assert instrument == "novaseq_SP"
    assert read_length == 36


def test_read_samplesheet(runfolder_path):
    samplesheet = _read_samplesheet(runfolder_path)

    assert len(samplesheet) == 4
    assert all(
        key.lower() == key
        for row in samplesheet
        for key in row
    )
    assert all(
        " " not in row["index"] and " " not in row["index2"]
        for row in samplesheet
    )
    assert samplesheet[0] == {
        'index': 'GAACTGAGCG',
        'index2': 'TCGTGGAGCG',
        'lane': 1,
        'sample_id': '14574-qiagen-indexset1-sp-lane1',
        'sample_project': 'AB-1234',
        "overridecycles": "Y36;I10;I10",
        'custom_description': 'LIBRARY_NAME:test',
        'custom_sample_name': '14574-qiagen-indexset1-sp-lane1',
    }
