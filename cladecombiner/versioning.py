import csv
import datetime
import re
import string
from typing import Iterable, TypeAlias, Union

import github.ContentFile
import github.Repository
from github import Github

Datelike: TypeAlias = Union[datetime.datetime, str]


def _twentythree_fiftynine(date: str) -> datetime.datetime:
    """
    Take YYYY-MM-DD date string and return a datetime object for 23:59:59:999999 UTC on that date
    """
    y, m, d = date.split("-")
    return datetime.datetime(
        int(y), int(m), int(d), 23, 59, 59, 999999, datetime.timezone.utc
    )


def _parse_date(as_of: Datelike) -> datetime.datetime:
    """
    as_of: Datelike
        The as-of date for getting the file. If providing a str, must be
        YYYY-MM-DD, and the time be taken to be 23:59:59:999999 such that a
        commit on the specified date will be used.
    """
    if isinstance(as_of, str):
        return _twentythree_fiftynine(as_of)
    elif isinstance(as_of, datetime.date):
        return as_of
    else:
        raise TypeError(f"Cannot parse date {as_of} of type {type(as_of)}")


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
        The as-of date for getting the file.

    Returns
    -------
    str
        The sha of the commit.

    Raises
    -------
    RuntimeError
        If the as-of date is prior to the earliest commit.
    """
    commits = repo.get_commits(path=file_path)
    sha = None
    # TODO: even for small histories this takes a noticeable amount of time, a smarter search would be good
    for commit in commits:
        time = datetime.datetime.strptime(
            commit.commit.raw_data["committer"]["date"], "%Y-%m-%dT%H:%M:%S%z"
        )
        if time <= as_of:
            sha = commit.sha
            break
    if sha is None:
        raise RuntimeError(
            f"Earliest commit in repo, at {time}, is after as-of date {as_of}."
        )
    return sha


def get_gh_file_contents_as_of(
    repo_name: str, file_path: str, as_of: Datelike
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
    as_of: Datelike
        The as-of date for getting the file. If providing a str, must be
        YYYY-MM-DD, and the time be taken to be 23:59:59:999999 such that a
        commit on the specified date will be used.

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

    sha = _get_gh_sha_as_of(repo, file_path, _parse_date(as_of))

    content_file = repo.get_contents(file_path, ref=sha)

    assert isinstance(content_file, github.ContentFile.ContentFile)
    assert isinstance(content_file.decoded_content, bytes)

    return content_file.decoded_content.decode("utf-8")


def _nextstrain_sc2_extractor(
    file_content: str, paranoid: bool = True
) -> Iterable[str]:
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


def _pango_sc2_extractor(file_content: str) -> Iterable[str]:
    """
    For parsing Pango lineages listed in cov-lineages/pango-designation/lineage_notes.txt
    """
    lineages_reader = csv.DictReader(file_content.split("\n"), delimiter="\t")
    return set(
        row["Lineage"] for row in lineages_reader if row["Lineage"][0] != r"*"
    )
