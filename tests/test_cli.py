
from lsvtool.cli.init import init
from lsvtool.cli.run import run

from pathlib import Path


TESTDATA_PATH = Path("tests/testdata")
DATA_10X_NA28385 = TESTDATA_PATH / "10x_NA24385.GRCh38.dels.reformatted.vcf.gz"
DATA_PACBIO_NA28385 = TESTDATA_PATH / "PacBio_HG002_GRCh38.pbsv.vcf.gz"


def test_init(tmp_path):
    init(
        vcfs=[DATA_10X_NA28385, DATA_PACBIO_NA28385],
        sample_names=["10x", "PacBio"],
        directory=tmp_path / "results"
        )


def test_run(tmp_path):
    outdir = tmp_path / "results"
    init(
        vcfs=[DATA_10X_NA28385, DATA_PACBIO_NA28385],
        sample_names=["10x", "PacBio"],
        directory=outdir,
        )
    run(workdir=outdir)
