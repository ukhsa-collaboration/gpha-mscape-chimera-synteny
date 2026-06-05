#!/usr/bin/env python3
from importlib import resources
import argparse
import datetime
import os
import shlex
import shutil
import subprocess
import tempfile
import uuid
from pathlib import Path

import pandas as pd
import plotly.express as px
from Bio import Entrez
from onyx import OnyxClient, OnyxConfig, OnyxEnv, OnyxField
from taxaplease import TaxaPlease
from tqdm import tqdm
from functional import pseq

__version__ = "0.0.2"

###########
# Globals #
###########

TESTING = True
SAMTOOLS_CMD = "conda run -n samtools samtools"

## harmless pandas options
pd.options.display.max_columns = 0
pd.options.mode.copy_on_write = True
pd.set_option("display.max_colwidth", 0)

## set up onyx config
config = OnyxConfig(
    domain=os.environ[OnyxEnv.DOMAIN],
    token=os.environ[OnyxEnv.TOKEN],
)

## init taxaplease as none - to overwrite in main
tp = None


#############
# Functions #
#############


def init_argparser():
    parser = argparse.ArgumentParser(
        prog="chimera-synteny",
        description="Generates a HTML report with synteny plots, given CLIMB IDs as input",
    )

    parser.add_argument("input", help="CLIMB ID", nargs="+")
    parser.add_argument(
        "--email", help="Email address to supply for Entrez queries", required=True
    )
    parser.add_argument(
        "--outdir", help="Output directory. 'stats' by default", default="stats"
    )
    parser.add_argument(
        "--allow_ttv",
        help="Allow Anelloviridae in output (Torque Teno Virus)",
        action="store_true",
    )
    parser.add_argument(
        "--allow_singleton",
        help="Allow figures with only a single taxon",
        action="store_true",
    )
    parser.add_argument(
        "--allow_phage", help="Allow phages in output", action="store_true"
    )
    parser.add_argument(
        "--allow_tmv",
        help="Allow Virgaviridae in output (Tobacco Mosaic Virus)",
        action="store_true",
    )
    parser.add_argument(
        "--onyx-project",
        help="Speicfy Onyx project to use (mSCAPE/synthscape)",
        default="mscape"
    )
    parser.add_argument(
        "--database",
        help="Path to a TaxaPlease database",
        default=None,
    )

    return parser


def resolve_samtools_path():
    samtools_path = shutil.which("samtools")
    if samtools_path:
        return shlex.quote(samtools_path)

    return SAMTOOLS_CMD


def get_file_from_s3(*, input_s3_uri=None, output_folder=None):
    if not output_folder:
        output_folder = os.getcwd()

    outfile_path = Path(output_folder, os.path.basename(input_s3_uri))

    cmd = ["s3cmd", "get", "--quiet", "--skip-existing", input_s3_uri, outfile_path]

    proc = subprocess.run(cmd)

    return str(outfile_path)


def get_chimera_bam_uri_by_climb_id(input_climb_id_list, onyx_project):
    """
    Takes in a list of climb ids
    Returns a dataframe with the columns climb_id and chimera_bam
    The latter containing a S3 URI to download the chimera bam
    file corresponding to the climb id
    """
    ## secretly turn single string inputs
    ## into a list despite the fact the argument
    ## clearly wants a list
    if isinstance(input_climb_id_list, str):
        input_climb_id_list = [input_climb_id_list]

    with OnyxClient(config) as client:
        return_data = pd.DataFrame(
            client.query(
                project=onyx_project,
                query=(OnyxField(climb_id__in=input_climb_id_list)),
                include=("climb_id", "chimera_bam"),
            )
        )

    return return_data


def process_bam_uri_to_pileup_gz(input_bam_uri, *, outdir=None):
    if not outdir:
        outdir = "stats"

    outfile_name = f"{os.path.basename(input_bam_uri)}.pileup.gz"
    outfilepath = str(os.path.join(outdir, outfile_name))

    if os.path.isfile(outfilepath):
        return outfilepath

    with tempfile.TemporaryDirectory() as tempdir:
        if Path(input_bam_uri).is_file():
            bam_file = input_bam_uri
        else:
            bam_file = get_file_from_s3(input_s3_uri=input_bam_uri, output_folder=tempdir)

        if not bam_file:
            print(f"Failed to retrieve {input_bam_uri}")

        ## this will prefer samtools on PATH, otherwise uses a conda fallback
        samtools_cmd = resolve_samtools_path()
        cmd = f"{samtools_cmd} mpileup -a {shlex.quote(bam_file)} | cut -f 1,2,4 | gzip > {shlex.quote(outfilepath)}"

        os.system(cmd)

    return outfilepath


def get_taxid_from_ncbi_accession(ncbi_accession_list):
    Entrez.email = ENTREZ_EMAIL
    ## get handle
    handle = Entrez.efetch(
        db="nuccore", id=",".join(ncbi_accession_list), retmode="xml"
    )
    ## parse XML into dict
    result = Entrez.parse(handle)
    ## it's a generator, so we force it into a list
    result = list(result)
    return_result = []

    ## get the bits we need
    for r in result:
        return_result.append(
            (
                r.get("GBSeq_accession-version"),
                r.get("GBSeq_primary-accession"),
                [
                    x.get("GBQualifier_value", [[]])
                    for x in r.get("GBSeq_feature-table", [[]])[0].get(
                        "GBFeature_quals", []
                    )
                    if (x.get("GBQualifier_name") == "db_xref")
                ][0],
                (
                    [
                        x.get("GBQualifier_value", [[]])
                        for x in r.get("GBSeq_feature-table", [[]])[0].get(
                            "GBFeature_quals", []
                        )
                        if (x.get("GBQualifier_name") == "segment")
                    ]
                    or [""]
                )[0],
            )
        )

    return return_result


def update_lookup_in_place(pileup_df):
    try:
        ncbi_accession_set = set(pileup_df["accession"])

        accession_lookup_dict = pd.read_csv(
            get_lookup_path(), sep="\t", names=["accession", "_", "taxid", "segment"]
        )
        accession_lookup_dict["taxid"] = accession_lookup_dict["taxid"].apply(
            lambda x: x.split(":")[-1] if ":" in x else x
        )
        accession_lookup_dict = dict(
            accession_lookup_dict[["accession", "taxid"]].to_records(index=False)
        )

        stuff_to_lookup = [
            x for x in ncbi_accession_set if x not in accession_lookup_dict.keys()
        ]
        if not len(stuff_to_lookup):
            return

        result = get_taxid_from_ncbi_accession(stuff_to_lookup)

        with open(get_lookup_path(), "a") as outfile:
            for record in result:
                outfile.write("\t".join(record))
                outfile.write("\n")
    except Exception as e:
        print(e)
        raise e


def process_df(input_df):
    max_pos = input_df["position"].max()
    n_bins = 100
    input_df["pos_ratio"] = input_df["position"] / max_pos
    input_df["pos_ratio_bin"] = input_df["pos_ratio"].apply(
        lambda x: round(x * n_bins)
    )  ## dodgy

    return_df = (
        input_df[["accession", "pos_ratio_bin", "depth"]]
        .groupby(["accession", "pos_ratio_bin"])
        .mean()
        .reset_index()
    )

    return return_df


def process_pileup_into_figure_array(
    pileup,
    *,
    suppress_ttv=False,
    suppress_singleton_families=False,
    suppress_phage=False,
    suppress_tmv=False,
):
    """
    Relies on some side-effects and code that self-modifies the package.
    Be careful!

    Takes in a path to a samtools mpileup dump
    Figures out what references are in the pileup and makes per ref
    dotplots showing coverage and depth.

    Args
    ----
    pileup
        Path to a samtools mpileup dump
    suppress_ttv
        Suppress torque teno virus and other Anelloviridae
        from the output
    suppress_singleton_families
        Suppress plots where only one accession exists in the dataset
        for that family
    suppress_phage
        Suppress plots for accessions that correspond to a phage
    suppress_tmv
        Suppress tobacco mosaic virus and other Virgaviridae

    Returns
    -------
    list[plotly.express.scatter]
        list of plotly figs
    """
    ## read in the pileup
    pileup_df = pd.read_csv(pileup, sep="\t", names=["accession", "position", "depth"])

    ## update installed lookup.txt
    update_lookup_in_place(pileup_df)

    ## group by accession, add in metadata
    ## then reconcatenate the data
    try:
        concat_df = pd.concat(
            [process_df(df) for accession, df in pileup_df.groupby("accession")]
        )
    except Exception as e:
        ## most likely excpetion is that there is nothing to concatenate
        ## so we return an empty list to allow things to progress
        ## printing the error as we go so it isn't completely silent
        print(f"Non-fatal error processing {pileup}")
        print(e)
        return []

    ## get accession to taxid mappings
    accession_lookup_dict = pd.read_csv(
        get_lookup_path(), sep="\t", names=["accession", "_", "taxid", "segment"]
    )
    accession_lookup_dict["taxid"] = accession_lookup_dict["taxid"].apply(
        lambda x: x.split(":")[-1] if ":" in x else x
    )
    accession_lookup_dict = dict(
        accession_lookup_dict[["accession", "taxid"]].to_records(index=False)
    )
    ## get accession to segment mappings
    segment_lookup_dict = pd.read_csv(
        get_lookup_path(), sep="\t", names=["accession", "_", "taxid", "segment"]
    )
    segment_lookup_dict = dict(
        segment_lookup_dict[["accession", "segment"]].to_records(index=False)
    )

    ## add in taxon metadata
    concat_df["taxid"] = concat_df["accession"].apply(
        lambda x: accession_lookup_dict.get(x)
    )
    concat_df["taxon"] = concat_df["taxid"].apply(
        lambda x: (tp.get_record(x) or {}).get("name")
    )
    concat_df["taxon_accession"] = concat_df["taxon"] + "_" + concat_df["accession"]
    concat_df["genus"] = (
        concat_df["taxid"]
        .apply(lambda x: (tp.get_specified_rank_record(x, "genus") or {}).get("name"))
        .fillna("missing genus lookup")
    )
    concat_df["family"] = (
        concat_df["taxid"]
        .apply(lambda x: (tp.get_specified_rank_record(x, "family") or {}).get("name"))
        .fillna("missing family lookup")
    )
    concat_df["is_phage"] = concat_df["taxid"].apply(lambda x: tp.isPhage(x))
    concat_df["segment"] = (
        concat_df["accession"].apply(lambda x: segment_lookup_dict.get(x, 0)).fillna(0)
    )

    ## variable for plotly figures
    figlist = []

    ## create a plot for each family in the dataset
    ## effectively a dotplot with colour indicating depth
    for rank, temp_df in concat_df[concat_df["depth"] > 0].groupby("family"):
        if suppress_ttv and (rank == "Anelloviridae"):
            continue
        if suppress_singleton_families and (len(set(temp_df["accession"])) == 1):
            continue
        if suppress_phage and list(temp_df["is_phage"])[0]:
            continue
        if suppress_tmv and (rank == "Virgaviridae"):
            continue

        try:
            temp_df = temp_df.sort_values("segment")
        except:
            print(f"Non-fatal error: cannot sort segments for {pileup}")

        # print(rank)
        fig = px.scatter(
            temp_df,
            x="pos_ratio_bin",
            y="taxon",
            color="depth",
            height=350,  ##+(len(temp_df) * 1),
            title=f"<b>Position bin vs reference accession (family: {rank})</b><br>Coloured by average read depth. Generated: {datetime.datetime.now()}. Input: {os.path.basename(pileup)}",
            facet_col="segment",
        )

        fig.update_xaxes(range=(-0.9, 100.9))

        fig.update_yaxes(tickmode="linear")

        figlist.append(fig)

    return figlist


def figure_array_to_html(
    list_of_figures,
    *,
    outhtml_filename=None,
    outhtml_folder=None,
):
    """
    Takes in a list of plotly.express figures,
    a filename and a path to write a HTML report to.

    Chucks all the figures into the HTML and
    calls it a good job.
    """
    ## default to current working dir
    outhtml_folder = os.getcwd() if not outhtml_folder else outhtml_folder
    ## make up a filename if unspecified
    outhtml_filename = (
        f"{uuid.uuid4()}.html" if not outhtml_filename else outhtml_filename
    )
    ## make sure the filename ends in html
    outhtml_filename = (
        f"{outhtml_filename}.html"
        if not outhtml_filename.endswith(".html")
        else outhtml_filename
    )

    outfile_full_path = os.path.join(outhtml_folder, outhtml_filename)

    with open(outfile_full_path, "wb") as outhtml:
        for idx, fig in enumerate(list_of_figures):
            outhtml.write(
                fig.to_html(full_html=False, include_plotlyjs=bool(not idx)).encode(
                    "utf-8"
                )
            )

    return outfile_full_path


def get_js_text():
    return resources.files("chimera_synteny.data").joinpath("main.js").read_text()


def get_css_text():
    return resources.files("chimera_synteny.data").joinpath("main.css").read_text()


def get_lookup_path():
    return str(
        resources.files("chimera_synteny.data").joinpath("lookup.txt")
    )


def generate_report(
    input,
    *,
    outdir=None,
    suppress_ttv=True,
    suppress_singleton_families=True,
    suppress_phage=True,
    suppress_tmv=True,
    onyx_project
):
    ## make an output directory
    Path(outdir).mkdir(exist_ok=True)

    local_bams = [x for x in input if Path(x).is_file()]
    climb_ids  = list(set(input) - set(local_bams))

    s3_uris = []
    if climb_ids:
        print(f"{datetime.datetime.now()} Querying Onyx...")
        climb_id_to_chimera_uri_df = get_chimera_bam_uri_by_climb_id(climb_ids, onyx_project)
        s3_uris = list(climb_id_to_chimera_uri_df["chimera_bam"])

    ## download (if needed) and process bams in parallel
    print(f"{datetime.datetime.now()} Retrieving bams and generating pileups")
    pileup_filepath_list = pseq(local_bams + s3_uris).map(
        lambda x: process_bam_uri_to_pileup_gz(x, outdir=outdir)
    )

    ## make figures and write them out to a single HTML file
    outreportpath = os.path.join(
        outdir, f"synteny_run.{'test' if TESTING else uuid.uuid4()}.html"
    )

    with open(outreportpath, "w") as outhtml:
        plotly_js_iswritten = False

        outhtml.write(
            "\n".join(
                [
                    "<html>",
                    "<head>",
                    "<style>",
                    get_css_text(),
                    "</style>",
                    "</head>",
                    "<body>",
                ]
            )
        )

        print(f"{datetime.datetime.now()} Processing pileups into figures")

        for pileup in tqdm(
            pileup_filepath_list, total=len(input)
        ):
            outhtml.write("<div class='figureAndHeaderContainer'>\n")
            outhtml.write(f"<h2>{os.path.basename(pileup)}</h2>\n")
            outhtml.write("<div class='figureContainer'>\n")

            figs = process_pileup_into_figure_array(
                pileup,
                suppress_ttv=suppress_ttv,
                suppress_singleton_families=suppress_singleton_families,
                suppress_phage=suppress_phage,
                suppress_tmv=suppress_tmv,
            )

            for fig in figs:
                outhtml.write(
                    fig.to_html(
                        full_html=False, include_plotlyjs=bool(not plotly_js_iswritten)
                    )
                )

                plotly_js_iswritten = True

            if not len(figs):
                outhtml.write("No figures to report.")

            outhtml.write(f"</div>\n</div>\n")
        outhtml.write(
            "\n".join(
                [
                    "</body>",
                    "<script>",
                    get_js_text(),
                    "</script>",
                    "</html>",
                ]
            )
        )

    print(f"{datetime.datetime.now()} Report written to {outreportpath}")


def main():
    global ENTREZ_EMAIL, tp

    args = init_argparser().parse_args()

    ENTREZ_EMAIL = args.email
    tp = TaxaPlease(database=args.database) if args.database else TaxaPlease()

    generate_report(
        args.input,
        outdir=args.outdir,
        suppress_ttv=not args.allow_tmv,
        suppress_singleton_families=not args.allow_singleton,
        suppress_phage=not args.allow_phage,
        suppress_tmv=not args.allow_tmv,
        onyx_project=args.onyx_project
    )


if __name__ == "__main__":
    main()
