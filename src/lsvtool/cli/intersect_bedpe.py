"""
Filter VCF with SV call based on overlap with BEDPE breakpoints.
"""
import logging
from collections import defaultdict
from itertools import chain

import vcf
from lsvtool.utils import Summary

logger = logging.getLogger(__name__)


def add_arguments(parser):
    parser.add_argument("vcf", help="Input sorted VCF with SV calls.")
    parser.add_argument("bedpe", help="Input BEDPE with breakpoints")
    parser.add_argument("-o", "--output-vcf", default="/dev/stdout",
                        help="Output VCF. Default to stdout.")
    parser.add_argument("-d", "--dist", type=int, default=0,
                        help="Min distance from breakpoint. Default: %(default)s")


def main(args):
    run_intersect(
        input_vcf=args.vcf,
        bedpe=args.bedpe,
        output_vcf=args.output_vcf,
        dist=args.dist,
    )


def run_intersect(
    input_vcf: str,
    bedpe: str,
    output_vcf: str,
    dist: int = 0,
):
    summary = Summary()

    # Parse BEDPE
    chrom_breakpoints = defaultdict(list)
    with open(bedpe) as f:
        # Check first line to confirm format is BEDPE
        first_line = next(f)

        # Skip header with comments if exists
        while first_line.startswith("#"):
            first_line = next(f)

        els = first_line.strip().split("\t")

        # Need at least 6 columns
        assert len(els) > 5

        # Columns 2,3,5,6 should be coordinates
        assert all(els[i].isdigit() for i in [1, 2, 4, 5])

        # Join first with rest
        f = chain([first_line], f)

        for line in f:
            summary["BEDPE records"] += 1
            els = line.strip().split("\t")
            chrom = els[0]

            # TODO handle breakpoints on different chromosomes
            if els[3] != chrom:
                summary["BEDPE records on different chroms"] += 1
                continue

            breakpoint1 = (max(0, int(els[1]) - dist), int(els[2]) + dist)
            breakpoint2 = (max(0, int(els[4]) - dist), int(els[5]) + dist)
            if breakpoint1[0] > breakpoint2[0]:
                breakpoint1, breakpoint2 = breakpoint2, breakpoint1

            chrom_breakpoints[chrom].append((breakpoint1, breakpoint2))

    # Parse VCF and write output
    vcf_reader = vcf.Reader(filename=input_vcf)
    if "END" not in vcf_reader.infos:
        logger.warning("Missing END information in VCF header")

    vcf_writer = vcf.Writer(stream=open(output_vcf, "w"), template=vcf_reader)
    prev_chrom = None
    breakpoints = None
    for record in vcf_reader:
        summary["VCF records read"] += 1
        chrom = record.CHROM
        start = record.start
        try:
            end = record.INFO["END"]
            assert start < end
        except KeyError:
            logger.info(f"Record missing END: {record}")
            end = -1

        if chrom != prev_chrom or prev_chrom is None:
            prev_chrom = chrom
            breakpoints = chrom_breakpoints.get(chrom, [])
            breakpoints.sort()

        intersect = False
        for bp1, bp2 in breakpoints:
            # Skip if starts before first breakpoint
            if start < bp1[0]:
                continue

            # Skip if starts after first breakpoint
            if start > bp1[1]:
                continue

            # Skip if ends before first breakpoint
            if end < bp2[0]:
                continue

            # Break if ends after first breakpoint, no match possible
            if end > bp2[1]:
                break

            intersect = True
            summary["VCF record matching breakpoints"] += 1
            logger.debug(f"Matched breakpoint on {chrom}.\n"
                         f" {bp1[0]:,} < start = {start:,} < {bp1[1]:,}\n"
                         f" {bp2[0]:,} <   end = {end:,} < {bp2[1]:,}")
            break

        if not intersect:
            summary["VCF record written"] += 1
            vcf_writer.write_record(record)

    summary.print_stats(name=__name__)
