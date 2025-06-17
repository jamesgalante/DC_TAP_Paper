import csv
import os
import click
from jinja2 import Environment, FileSystemLoader

TEMPLATE_MAP = {
    "crispr": "./templates/dc_tap_seq_crispr_seqspec.yaml.j2",
    "rna": "./templates/dc_tap_seq_rna_seqspec.yaml.j2"
}

@click.command()
@click.option('--tsv', required=True, type=click.Path(exists=True), help='Input TSV file with variables for the template.')
@click.option('--modality', required=True, type=click.Choice(['crispr', 'rna']), help='Modality to use: crispr or rna.')
@click.option('--output-dir', required=False, type=click.Path(), help='Directory to write output YAML files. If not set, output to stdout.')
def fill_template(tsv, modality, output_dir):
    """
    Fill a Jinja YAML template using variables from a TSV file.

    The template is selected based on the modality (crispr or rna).
    Each row in the TSV is rendered as a separate YAML file.
    Output files are named using the value of 'read1_name' or row number if not present.
    """
    # Select template path
    template_relpath = TEMPLATE_MAP.get(modality)
    if not template_relpath:
        raise click.ClickException(f"Unknown modality: {modality}")
    template_path = os.path.join(os.path.dirname(__file__), template_relpath)
    if not os.path.exists(template_path):
        raise click.ClickException(f"Template file not found: {template_path}")

    # Prepare Jinja environment
    template_dir, template_file = os.path.split(os.path.abspath(template_path))
    env = Environment(loader=FileSystemLoader(template_dir))
    tmpl = env.get_template(template_file)

    # Read TSV
    with open(tsv, newline='') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for idx, row in enumerate(reader):
            rendered = tmpl.render(**row)
            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
                outname = row.get('read1_name', f'output_{idx+1}')
                outfile = os.path.join(output_dir, f'{outname}.yaml')
                with open(outfile, 'w') as f:
                    f.write(rendered)
                click.echo(f'Wrote {outfile}')
            else:
                click.echo(rendered)
                click.echo('---')


if __name__ == '__main__':
    fill_template()
