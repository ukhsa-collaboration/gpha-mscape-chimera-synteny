# chimera-synteny

## What is this?

An argparse tool that generates a HTML report with synteny plots, given CLIMB IDs as input.

## Installation?

```bash
## pip install straight from the github repo
pip install git+https://github.com/ukhsa-collaboration/gpha-mscape-chimera-synteny.git
```

## How do I use this?

After pip installing, you can run it directly from the commandline:

```bash
chimera-synteny C-00000000 C-111111111 C-22222222 C-3333333 --email "example@example.com"
```

Replacing the CLIMB IDs with real IDs, and the email address as appropriate (used for querying Entrez)

You can additionally specify the path to a [taxaPlease](https://github.com/ukhsa-collaboration/gpha-mscape-taxaplease)
database using the `--database` argument, or provide local bam files as an argument instead of CLIMB IDs.

## What is Chimera?

That would be this: <https://github.com/CLIMB-TRE/chimera>

## What is synteny?

The notion that gene order remains mostly conserved, both intraspecies and possibly interspecies to some degree.

With that as an assumption, we can scale all genomes within a taxonomic family to 100% length, plot the
depth across the genome based on alignments from Chimera, and this should show us if the coverage switches
across multiple references.

This might be particularly useful for something like Influenza A, where good coverage and depth might be seen
for segments from one accession, but other segments might align to a different accession. This would potentially
indicate a reassortment event.

Such plots would also potentially be useful for identifying off-target alignments, mixed infection and recombination.
