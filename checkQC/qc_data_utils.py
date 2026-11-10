import numpy as np
from checkQC.parsers.illumina import _read_interop_summary


def bclconvert_test_runfolder(qc_data, runfolder_path):
    _, _, run_info = _read_interop_summary(runfolder_path)
    flowcell_id = run_info.flowcell_id()
    if "HMTFYDRXX" in flowcell_id:
        return {
            "qc_data": qc_data,
            "expected_instrument": "novaseq_SP",
            "expected_read_length": 36,
            "expected_samplesheet": {
                "len": 4,
                "head": [
                    {
                        "lane": 1,
                        "sample_id": "14574-qiagen-indexset1-sp-lane1",
                        "index": "GAACTGAGCG",
                        "index2": "TCGTGGAGCG",
                        "sample_project": "AB-1234",
                        "overridecycles": "Y36;I10;I10",
                        "custom_description": "LIBRARY_NAME:test",
                        "custom_sample_name": "14574-qiagen-indexset1-sp-lane1",
                    },
                    {
                        "lane": 1,
                        "sample_id": "CD-5678-Zymo-IndexSet1-SP-Lane1",
                        "index": "GTCCGTAAGC",
                        "index2": "CGTGTATCTT",
                        "sample_project": "CD-5678",
                        "overridecycles": "Y36;I10;I10",
                        "custom_description": "LIBRARY_NAME:test",
                        "custom_sample_name": "CD-5678-Zymo-IndexSet1-SP-Lane1",
                    },
                    {
                        "lane": 2,
                        "sample_id": "14574-Qiagen-IndexSet1-SP-Lane2",
                        "index": "GAACTGAGCG",
                        "index2": "TCGTGGAGCG",
                        "sample_project": "AB-1234",
                        "overridecycles": "Y36;I10;I10",
                        "custom_description": "LIBRARY_NAME:test",
                        "custom_sample_name": "14574-Qiagen-IndexSet1-SP-Lane2",
                    },
                    {
                        "lane": 2,
                        "sample_id": "CD-5678-Zymo-IndexSet1-SP-Lane2",
                        "index": "GTCCGTAAGC",
                        "index2": "CGTGTATCTT",
                        "sample_project": "CD-5678",
                        "overridecycles": "Y36;I10;I10",
                        "custom_description": "LIBRARY_NAME:test",
                        "custom_sample_name": "CD-5678-Zymo-IndexSet1-SP-Lane2",
                    },
                ],
            },
            "expected_sequencing_metrics": {
                1: {
                    "total_reads_pf": 532_464_327,
                    "total_reads": 638_337_024,
                    "raw_density": 2_961_270.5,
                    "pf_density": 2_470_118.25,
                    "yield": 19_168_715_772,
                    "yield_undetermined": 19_091_771_604,
                    "top_unknown_barcodes": {
                        "len": 1000,
                        "head": [
                            {
                                'index': 'ATATCTGCTT', 'index2': 'TAGACAATCT',
                                'count': 2008558,
                            },
                            {
                                'index': 'CACCTCTCTT', 'index2': 'CTCGACTCCT',
                                'count': 1934726,
                            },
                            {
                                'index': 'ATGTAACGTT', 'index2': 'ACGATTGCTG',
                                'count': 1905527,
                            },
                            {
                                'index': 'TCACCGCGCT', 'index2': 'TCGATACTAG',
                                'count': 1819959,
                            },
                            {
                                'index': 'TTCGGTGTGA', 'index2': 'GAACAAGTAT',
                                'count': 1818376,
                            },
                        ],
                    },
                    "reads": {
                        1: {
                            "mean_error_rate": np.nan,
                            "percent_q30": 95.70932006835938,
                            "is_index": False,
                            "mean_percent_phix_aligned": 0.,
                        },
                        2: {
                            "mean_error_rate": np.nan,
                            "percent_q30": 92.57965850830078,
                            "is_index": True,
                            "mean_percent_phix_aligned": np.nan,
                        },
                        3: {
                            "mean_error_rate": np.nan,
                            "percent_q30": 90.3790283203125,
                            "is_index": True,
                            "mean_percent_phix_aligned": np.nan,
                        },
                    },
                    "reads_per_sample": [
                        {
                            "sample_id": "14574-qiagen-indexset1-sp-lane1",
                            "reads_pf":1576102,
                            "percent_of_lane": 0.3,
                            "percent_perfect_index_reads": 97.93,
                            "mean_quality_score": 36.43,
                            "percent_q30": 96.0,
                        },
                        {
                            "sample_id": "CD-5678-Zymo-IndexSet1-SP-Lane1",
                            "reads_pf":561236,
                            "percent_of_lane": 0.11,
                            "percent_perfect_index_reads": 97.97,
                            "mean_quality_score": 36.43,
                            "percent_q30": 96.0,
                        },
                    ],
                },
                2: {
                    "total_reads_pf": 530_917_565,
                    "total_reads": 638_337_024,
                    "raw_density": 2_961_270.5,
                    "pf_density": 2_462_942.5,
                    "yield": 19_113_032_340,
                    "yield_undetermined": 19_036_322_532,
                    "top_unknown_barcodes": {
                        "len": 1000,
                        "head": [
                            {
                                'index': 'ATATCTGCTT', 'index2': 'TAGACAATCT',
                                'count': 2006759,
                            },
                            {
                                'index': 'CACCTCTCTT', 'index2': 'CTCGACTCCT',
                                'count': 1933013,
                            },
                            {
                                'index': 'ATGTAACGTT', 'index2': 'ACGATTGCTG',
                                'count': 1898240,
                            },
                            {
                                'index': 'TTCGGTGTGA', 'index2': 'GAACAAGTAT',
                                'count': 1818840,
                            },
                            {
                                'index': 'TCACCGCGCT', 'index2': 'TCGATACTAG',
                                'count': 1800407,
                            },
                        ],
                    },
                    "reads": {
                        1: {
                            "mean_error_rate": np.nan,
                            "percent_q30": 95.75276184082031,
                            "is_index": False,
                            "mean_percent_phix_aligned": 0.,
                        },
                        2: {
                            "mean_error_rate": np.nan,
                            "percent_q30": 92.60448455810547,
                            "is_index": True,
                            "mean_percent_phix_aligned": np.nan,
                        },
                        3: {
                            "mean_error_rate": np.nan,
                            "percent_q30": 90.2811050415039,
                            "is_index": True,
                            "mean_percent_phix_aligned": np.nan,
                        },
                    },
                    "reads_per_sample": [
                        {
                            "sample_id": "14574-Qiagen-IndexSet1-SP-Lane2",
                            "reads_pf":1568445,
                            "percent_of_lane": 0.3,
                            "percent_perfect_index_reads": 97.82,
                            "mean_quality_score": 36.44,
                            "percent_q30": 96.0,
                        },
                        {
                            "sample_id": "CD-5678-Zymo-IndexSet1-SP-Lane2",
                            "reads_pf":562383,
                            "percent_of_lane": 0.11,
                            "percent_perfect_index_reads": 97.98,
                            "mean_quality_score": 36.45,
                            "percent_q30": 97.0,
                        },
                    ],
                },
            },
        }
    else:
        raise Exception(
            "This function is only compatible with the run with flowcell_id: 'HMTFYDRXX', "
                f"the supplied runfolder has flowcell_id: {flowcell_id}"
        )

