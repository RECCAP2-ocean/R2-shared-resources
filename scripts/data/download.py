"""
This is an extract of a package I've written to download and process data. 
It's probably a bit overkill, but since I already have it, I thought why not :)
That being said, there might be some bugs... 
"""

import logging
import warnings

import fsspec
import pooch

warnings.filterwarnings("ignore", category=RuntimeWarning)
logging.basicConfig(
    level=15, 
    format="%(asctime)s [RECCAP]  %(message)s", 
    datefmt="%Y-%m-%d %H:%M:%S"
)

prefix = "info_"
logfile = f"{prefix}logging.log"
readme = f"{prefix}readme.txt"
cache = f"{prefix}filelist.cache"


def read_catalog(catalog_name):
    """
    Used to read YAML files that contain download information. 
    YAML file entries require:
        url: remote path to file/s. Can contain *
        dest: where the file/s will be stored
        meta:
            doi: url to the data source
            description: info about the data
            citation: how to cite this dataset
    """
    from dotenv import find_dotenv
    from envyaml import EnvYAML

    catalog_raw = EnvYAML(
        yaml_file=catalog_name, env_file=find_dotenv(), strict=True
    ).export()

    catalog = {}
    for k in catalog_raw:
        if ("." not in k) and (isinstance(catalog_raw[k], dict)):
            catalog[k] = catalog_raw[k]

    return catalog


def download(
    url="",
    dest="./",
    name="",
    n_jobs=8,
    login={},
    use_cache=True,
    verbose=True,
    **kwargs,
):
    """
    A high level function to retrieve data from a url with a wildcard.
    
    Parameters
    ----------
    url: str
        a url with wildcards (*) formatter. Python string formatters that 
        match kwarg entries will be replaced. 
    dest: str
        where the files will be saved to. String formatting supported (as with url)
    use_cache: bool
        if set to True, will use cached url list instead of fetching a
        new list. This is useful for updating data
    n_jobs: int
        the number of parallel downloads. Will not show progress bar 
        when n_jobs > 1
    login: dict
        required if username and passwords are required for protocol
    name: str
        used to keep track in logging. can set this to the data source
    verbose: bool / int
        if verbose is False, logging level set to ERROR (40)
        if verbose is True, logging level set to 15
        if verbose is intiger, then sets logging level directly
    **kwargs:
        not used.
    """

    import ftplib
    from collections import defaultdict
    
    if isinstance(verbose, bool):
        logging_level = 15 if verbose else 40
    elif isinstance(verbose, int):
        logging_level = verbose
    else:
        logging_level = 15
    logging.getLogger().setLevel(logging_level)

    kwargs = {**get_kwargs(), **kwargs}
    dest = dest.format_map(kwargs)
    
    log_fname = f"{dest}/{logfile}"
    log_to_file(log_fname)

    kwargs.update({'download_logging': log_fname})
    create_download_readme(**kwargs)

    urls = get_url_list(
        url=url.format_map(kwargs),
        raise_on_empty=False,
        cache_path=f"{dest}/{cache}",
        use_cache=use_cache,
        **login,
    )

    logging.info(
        f"{name.upper(): <15s} {len(urls): >3} " f"files {url.format_map(kwargs)}"
    )

    logging.log(15, f"{len(urls)} files exist for {url}".format_map(kwargs))
    logging.log(20, f"Files will be saved to {dest}")

    if len(urls) == 0:
        return []

    downloader = choose_downloader(url)

    flist = None
    flist = download_urls(
        urls,
        n_jobs=n_jobs,
        dest_path=dest,
        downloader=downloader(progressbar=True, **login)
    )
    if flist is None:
        raise ValueError("Files could not be downloaded")

    return flatten_list(flist)


def get_url_list(
    url,
    username=None,
    password=None,
    cache_path=None,
    use_cache=True,
    raise_on_empty=True,
):
    """
    If a url has a wildcard (*) value, remote files will be searched for.
    Leverages off the `fsspec` package. This doesnt work for all HTTP urls.
    Parameters
    ----------
    username: str
        if required for given url and protocol (e.g. FTP)
    password: str
        if required for given url and protocol (e.g. FTP)
    cache_path: str
        the path where the cached files will be stored
    use_cache: bool
        if there is a file with cached remote urls, then those
        values will be returned as a list
    raise_on_empty: bool
        if there are no files, raise an error or silently pass
    Returns
    -------
    a sorted list of urls
    """
    from pathlib import Path as posixpath
    from urllib.parse import urlparse

    from aiohttp import ClientResponseError
    from pandas import Series, read_csv

    if cache_path is None:
        cache_path = get_cache_path(url)
    else:
        cache_path = posixpath(cache_path)

    if cache_path.is_file() and use_cache:
        flist = read_csv(str(cache_path), index_col=False).iloc[:, 0].to_list()
        logging.log(15, f"Fetched {len(flist)} files from flist cache: {cache_path}")
        logging.debug(flist)

        return sorted(flist)

    purl = urlparse(url)
    protocol = purl.scheme
    host = purl.netloc
    path = purl.path

    logging.log(15, f"Fetching filenames from {url}")

    props = {"protocol": protocol}
    if not protocol.startswith("http"):
        props.update({"host": host})
    if username is not None:
        props["username"] = username
    if password is not None:
        props["password"] = password

    fs = fsspec.filesystem(**props)
    if protocol.startswith("http"):
        path = f"{protocol}://{host}/{path}"
        try:
            flist = fs.glob(path)
        except ClientResponseError:
            if raise_on_empty:
                raise ValueError(f"No files could be found for the url: {url}")
            else:
                return []
    else:
        flist = [f"{protocol}://{host}{f}" for f in fs.glob(path)]

    no_files = len(flist) == 0
    if no_files and raise_on_empty:
        raise ValueError(f"No files could be found for the url: {url}")

    if no_files and not use_cache:
        return flist

    cache_path.parent.mkdir(exist_ok=True, parents=True)
    # writing url list to cache file
    Series(flist, dtype="str").to_csv(str(cache_path), index=False)
    logging.log(15, f"Cached {len(flist)} urls to: {cache_path}")
    logging.debug(flist)

    return sorted(flist)


def download_urls(
    urls, 
    n_jobs=8, 
    dest_path="./{t:%Y}/{t:%m}", 
    date_format="%Y%m%d", 
    **kwargs
):
    """
    Downloads the given list of urls to a specified destination path using
    the `pooch` package in Python.
    NOTE: `fsspec` is not used as it fails for some FTP and SFTP protocols.
    Parameters
    ----------
    urls: list
        the list of URLS to download - may not contain wildcards
    dest_path: str
        the location where the files will be downloaded to. May contain
        date formatters that are labelled with "{t:%fmt} to create subfolders
    date_format: str
        the format of the date in the urls that will be used to fill in the
        date formatters in `dest_path` kwarg. Matches limited to 1970s to 2020s
    kwargs: key=value
        will be passed to pooch.retrieve. Can be used to set the downloader
        with username and password and the processor for unzipping. See
        `choose_downloader` for more info.
    Returns
    -------
    file names of downloaded urls
    """
    def pooch_retrieve_handling(kwargs):

        pooch.get_logger().setLevel(1000)
        url = kwargs.get('url')

        try:
            logging.log(15, f"retrieving {url}")
            return 0, pooch.retrieve(**kwargs)
        except:
            pass

        try:
            # this is for when the server does not allow the file size to be fetched
            kwargs['downloader'].progressbar = False
            return 0, pooch.retrieve(**kwargs)
        except:
            pass

        # this will raise the error
        try:
            pooch.retrieve(**kwargs)
        except Exception as e:
            if '550' in str(e):
                message = f'ERROR: Check file permissions: {url}. '
                logging.log(20, message)
            return 1, url

    import pandas as pd
    from datetime_matcher import DatetimeMatcher
    from pooch import Unzip

    re_date = DatetimeMatcher()
    # format will limit between 1970s and 2020s (with exceptions)
    re_date.format_code_to_regex_map["Y"] = "[12][90][789012][0-9]"

    download_args = []
    for url in urls:
        if "{t:" in dest_path:
            date = re_date.extract_datetime(date_format, url)
            if date is not None:
                date = pd.to_datetime(date, format=date_format)
                fpath = dest_path.format(t=date)
            else:
                raise ValueError("No date found in the dest_path")
        else:
            fpath = dest_path
        
        download_args += dict(
            url=url, 
            known_hash=None, 
            fname=url.split("/")[-1], 
            path=fpath, 
            processor=choose_processor(url), 
            **kwargs),
    
    n_jobs = min([n_jobs, len(download_args)])
    if n_jobs == 1:
        flist = [pooch_retrieve_handling(d) for d in download_args]
    elif 1 < n_jobs <= 8:
        from joblib import Parallel, delayed
        [setattr(d['downloader'], 'progressbar', False) for d in download_args]
        flist = Parallel(n_jobs=n_jobs, prefer="threads")(
            delayed(pooch_retrieve_handling)(d) for d in download_args)
    else:
        raise Exception('n_jobs must be between 1 and 8 to avoid too many requests')
    
    failed = [f for o, f in flist if o > 0]
    passed = [f for o, f in flist if o == 0]
    logging.info(
        f'SUMMARY: Retrieved={len(passed)}, Failed={len(failed)} listing failed below: \n' + 
        '\n'.join(failed)
    )
    
    return passed


def get_cache_path(url, cache_dir=None):
    """
    Creates the path for the cache used to store remote file names
    Saves time in updating the
    """
    import hashlib
    import tempfile
    from pathlib import Path as posixpath

    if cache_dir is None:
        cache_dir = tempfile.gettempdir()

    cache_fname = hashlib.md5(str(url).encode()).hexdigest()
    cache_path = posixpath(f"{cache_dir}/{cache_fname}")

    return cache_path

    
def choose_processor(url):
    """
    chooses the processor to uncompress if required
    """
    known_processors = {
        pooch.Decompress(): (".gz2", ".gz"),
        pooch.Untar(): (".tar", ".tgz", ".tar.gz"),
        pooch.Unzip(): (".zip",),
        None: "*"
    }
    
    chosen = None
    for processor, extensions in known_processors.items():
        for ext in extensions:
            if ext in url:
                chosen = processor
    return chosen


def choose_downloader(url):
    """
    Will automatically select the correct downloader for the given url.
    Pass result to pooch.retrieve(downloader=downloader(**kwargs))
    Parameters
    ----------
    url: str
        the path of a url
    Returns
    -------
    pooch.Downloader as a function. Resulting function Can be called with
    (username, password, progressbar) options.
    """
    from urllib.parse import urlparse as parse_url
    import pooch

    known_downloaders = {
        "ftp": pooch.FTPDownloader,
        "https": pooch.HTTPDownloader,
        "http": pooch.HTTPDownloader,
    }

    parsed_url = parse_url(url)
    if parsed_url.scheme not in known_downloaders:
        raise ValueError(
            f"Unrecognized URL protocol '{parsed_url.scheme}' in '{url}'. "
            f"Must be one of {known_downloaders.keys()}."
        )
    downloader = known_downloaders[parsed_url.scheme]

    return downloader


def flatten_list(list_of_lists):
    if len(list_of_lists) == 0:
        return list_of_lists
    if isinstance(list_of_lists[0], list):
        return flatten_list(list_of_lists[0]) + flatten_list(list_of_lists[1:])
    return list_of_lists[:1] + flatten_list(list_of_lists[1:])


def create_download_readme(**source_dict):
    import inspect
    from pathlib import Path as posixpath
    import logging

    dest = source_dict.get("dest").format_map(source_dict)
    cache_fname = f"{source_dict.get('dest')}/{cache}"
    manipulation = inspect.cleandoc(
        f"""
    Data has been downloaded directly from the server shown in URL.
    There has been no modification to the original files.
    There may be a data cache located in the annual subfolders of each
    with the format {cache_fname.replace('//', '/')}
    """
    )

    args = [
        source_dict.get("name", ''),
        source_dict.get("meta", {}).get("doi", None),
        source_dict.get("url", None),
        source_dict.get("meta", {}).get("citation", None),
        source_dict.get("meta", {}).get("description", None),
        source_dict.get("variables", []),
        manipulation,
    ]

    readme_fname = posixpath(f"{dest}/{readme}")
    readme_fname.parent.mkdir(parents=True, exist_ok=True)

    email = source_dict.get("email", None)
    logging = source_dict.get("download_logging", "None")

    readme_text = make_readme_file(*args, email=email, download_logging=logging)

    with open(readme_fname, "w") as file:
        file.write(readme_text)
        
        
def make_readme_file(
    dataset_name,
    doi_or_link,
    url,
    citation,
    description,
    variables,
    manipulation,
    download_logging="None",
    email=None,
):
    """
    Adheres to the UP group's (ETHZ) readme prerequisites.

    Parameters
    ----------
    dataset_name: str
        The name of the dataset that will be at the top of the file
    doi_or_link: str
        The link to where the dataset can be downloaded or more info
        can be fetched
    url: str
        The url used to download the data - may be useful for other
        downloaders. May contain wildcards and placeholders.
    citation: str
        Just making it easier for people to cite the dataset. If there
        is more than one doi, please cite that too. This may be the case
        where a dataset is published alongside a paper.
    description: str
        The description of the data - can be copied directly from the website
        of the data.
    manipulation: str
        Any manipulation or changes you make to the data before saving.
    variables: list
        A list of the names of the variables that are downloaded
    download_loggin: str
        The path to the download logging. Defaults to None
    email: str
        Defaults to USER@ethz.ch if not provided.

    """
    import inspect
    import os
    import pwd
    from textwrap import wrap

    import pandas as pd

    if email is None:
        username = pwd.getpwuid(os.getuid())[0]
        email = f"{username}@ethz.ch"
    today = pd.Timestamp.today().strftime("%Y-%m-%d")

    w = "\n" + " " * 4
    if variables == []:
        variables = ""
    elif isinstance(variables, list):
        variables = f"{w}Variables:{w}" + f"{w}".join(["- " + v for v in variables])
    else:
        variables = ""
    
    citation = w.join(wrap(citation.replace("\n", " "), 80))
    description = w.join(wrap(description.replace("\n", " "), 80))
    manipulation = w.join(wrap(manipulation.replace("\n", " "), 80))

    readme_text = inspect.cleandoc(
        f"""
    {'='*len(dataset_name)}
    {dataset_name}
    {'='*len(dataset_name)}

    Contact: {email}
    Date:    {today}
    Source:  {doi_or_link}
    URL:     {url}
    Logging: {download_logging}


    ------------
    Dataset info
    ------------
    {citation}

    {description}
    {variables}

    ------------------
    Dataset processing
    ------------------
    {manipulation}



    readme file was automatically created using netCDFdownloader tool
    https://github.com/lukegre/netCDF-Downloader

    """
    )

    return readme_text


def get_kwargs():
    import inspect

    frame = inspect.currentframe().f_back
    keys, _, _, values = inspect.getargvalues(frame)
    kwargs = {}
    for key in keys:
        if key != "self":
            kwargs[key] = values[key]
    return kwargs


def log_to_file(fname):
    import logging
    from pathlib import Path as posixpath

    fname = posixpath(fname)
    fname.parent.mkdir(exist_ok=True, parents=True)

    rootLogger = logging.getLogger()

    # remove existing file handlers
    for handler in rootLogger.handlers:
        if isinstance(handler, logging.FileHandler):
            rootLogger.handlers.remove(handler)

    # add the new logger with the formatting
    logFormatter = logging.Formatter(
        "%(asctime)s [DOWNLOAD]  %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    fileHandler = logging.FileHandler(fname)
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)

    logging.info("=" * 80 + "\n" * 2)
    logging.info("Start of logging session")
    logging.info("-" * 80)

