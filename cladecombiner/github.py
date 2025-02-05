import csv
import datetime
import re
import string
from typing import Collection, Optional

import github.ContentFile
import github.Repository
from github import Github


def _get_gh_sha_as_of(
    repo: github.Repository.Repository,
    file_path: str,
    as_of: datetime.datetime,
):
    """
    Get the SHA of the last commit in a GitHub repository prior to as_of.

    Parameters
    ---------
    repo : github.Repository.Repository
        PyGithub representation of the repo.
    file_path : str
        Relative path to file from repo root.
    as_of: datetime.datetime
        We will get the sha of the last commit before this datetime.

    Returns
    -------
    str
        The sha of the commit.

    Raises
    -------
    RuntimeError
        If the as-of date is prior to the earliest commit.
    """
    commits = repo.get_commits(path=file_path, until=as_of)
    try:
        commit = commits[0]
    # Intercept and provide a more interpretable error
    except IndexError as e:
        raise RuntimeError(
            f"There are no commits prior to as-of date {as_of}."
        ) from e

    time = datetime.datetime.strptime(
        commit.commit.raw_data["committer"]["date"], "%Y-%m-%dT%H:%M:%S%z"
    )

    assert time <= as_of

    return commit.sha


def get_gh_file_contents_as_of(
    repo_name: str, file_path: str, as_of: Optional[datetime.date]
) -> str:
    """
    Get the contents of a file in a GitHub repository as of the last commit prior to `as_of`.

    All times are taken to be in UTC.

    Parameters
    ---------
    repo_name : str
        The username/repo combination, e.g. CDCGov/cladecombiner for this
        package's repo.
    file_path : str
        Relative path to file from repo root, e.g. cladecombiner/utils.py for
        the file where this function is defined.
    as_of: Optional[datetime.date]
        The as-of date for getting the file. If None, current contents are
        retrieved. Otherwise, contents are obtained for the state at the
        end of the as_of date. That is, all changes on that day are included.

    Returns
    -------
    str
        The sha of the commit.

    Raises
    -------
    RuntimeError
        If the as-of date is prior to the earliest commit.
    """

    g = Github()

    repo = g.get_repo(repo_name)

    if as_of is None:
        content_file = repo.get_contents(file_path)
    else:
        as_of_dt = datetime.datetime(
            as_of.year,
            as_of.month,
            as_of.day + 1,
            tzinfo=datetime.timezone.utc,
        )
        sha = _get_gh_sha_as_of(repo, file_path, as_of_dt)
        content_file = repo.get_contents(file_path, ref=sha)

    assert isinstance(content_file, github.ContentFile.ContentFile)
    assert isinstance(content_file.decoded_content, bytes)

    return content_file.decoded_content.decode("utf-8")


def _nextstrain_sc2_extractor(
    file_content: str, paranoid: bool = True
) -> Collection[str]:
    """
    For parsing Nextstrain clades listed in nextstrain/ncov/defaults/clades.tsv
    """
    clades_reader = csv.DictReader(file_content.split("\n"), delimiter="\t")
    taxa = list(set(row["clade"] for row in clades_reader))

    # Handle clade names like 21L (Omicron)
    clade_regex = re.compile(r"\d\d[A-Z]")
    for i in range(len(taxa)):
        matches = clade_regex.match(taxa[i])
        assert matches is not None
        taxa[i] = matches.group(0)

    if paranoid:
        taxa.sort()
        assert taxa[0] == "19A"
        by_year = {}
        for taxon in taxa:
            year = int(taxon[:2])
            letter = taxon[2:]
            if year in by_year:
                by_year[year].append(letter)
            else:
                by_year[year] = [letter]

        assert set(range(min(by_year), max(by_year) + 1)) == set(by_year)
        for letters in by_year.values():
            highest = max(letters)
            for letter in string.ascii_uppercase:
                if letter <= highest:
                    assert letter in letters
                else:
                    break
    assert all(isinstance(taxon, str) for taxon in taxa)
    return set(taxa)


def _pango_sc2_extractor(file_content: str) -> Collection[str]:
    """
    For parsing Pango lineages listed in cov-lineages/pango-designation/lineage_notes.txt
    """
    lineages_reader = csv.DictReader(file_content.split("\n"), delimiter="\t")
    # TODO: probably ought to "coax" these names, requiring move of code elsewhere to avoid circular imports
    return set(
        row["Lineage"] for row in lineages_reader if row["Lineage"][0] != r"*"
    ) | set([""])
