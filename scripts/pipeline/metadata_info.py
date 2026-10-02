import json, os


def load_assembly_reports(report_path):
    """
    Loads an NCBI assembly_data_report.jsonl file into a dictionary.

    Args
    ------
    report_path (str): The path to the assembly_data_report.jsonl file
                       (one JSON record per line, one record per genome).

    Return
    ------
    reports (dict): Dictionary {genome accession (str): record (dict)}.
    """
    reports = {}
    with open(report_path) as f:
        for line in f:
            line = line.strip()
            if line:
                record = json.loads(line)
                reports[record["accession"]] = record
    return reports

def get_nested(data, path):
    """
    Gets a value from a nested JSON record using a dotted path.
    If a list is found along the path, every element is explored and the
    results are flattened (e.g. several bioprojects for one genome).

    Args
    ------
    data (dict): The JSON record (one genome).
    path (str): Dotted path to the value,
                e.g. "assemblyInfo.bioprojectLineage.bioprojects.accession".

    Return
    ------
    value (str): The value found. If several unique values are found, they
                 are joined with ";". "NA" if the path does not exist.
    """
    values = [data]
    for key in path.split("."):
        next_values = []
        for v in values:
            if isinstance(v, dict) and key in v:
                child = v[key]
                if isinstance(child, list):
                    next_values.extend(child)
                else:
                    next_values.append(child)
        values = next_values

    if not values:
        return "NA"
    # dict.fromkeys removes duplicates keeping the order
    return ";".join(dict.fromkeys(str(v) for v in values))

def add_genome_metadata(hits_df, extra_columns = None):
    """
    Adds the species and, optionally, extra columns taken from the
    assembly_data_report.jsonl of the genome each hit belongs to.

    Expected directory structure:
    .../{species}/{species}/ncbi_dataset/data/{accession}/ (genome_dir)
    .../{species}/{species}/ncbi_dataset/data/assembly_data_report.jsonl

    Args
    ------
    hits_df (pd.DataFrame): Hits table with the columns genome and genome_dir.
    extra_columns (dict, optional): Extra columns to add, as
                                    {column name: dotted path in the JSON}.
                                    Defaults to None (only species is added).

    Return
    ------
    hits_df (pd.DataFrame): A copy of the input with the "species" column
                            (name of the species folder) and one column per
                            entry of extra_columns ("NA" if not found).
    """
    hits_df = hits_df.copy()
    extra_columns = extra_columns or {}
    reports_cache = {} # Each jsonl is read only once
    species = []
    extras = {col: [] for col in extra_columns}

    for genome_dir, genome in zip(hits_df["genome_dir"], hits_df["genome"]):
        data_dir = os.path.dirname(genome_dir) # .../ncbi_dataset/data
        species.append(os.path.basename(os.path.dirname(os.path.dirname(data_dir))))

        if extra_columns:
            if data_dir not in reports_cache:
                reports_cache[data_dir] = load_assembly_reports(
                    os.path.join(data_dir, "assembly_data_report.jsonl"))
            record = reports_cache[data_dir].get(genome, {})
            for col, path in extra_columns.items():
                extras[col].append(get_nested(record, path))

    hits_df["species"] = species
    for col, values in extras.items():
        hits_df[col] = values
    return hits_df