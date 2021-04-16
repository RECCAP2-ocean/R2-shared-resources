"""
High level functions to read in data. Note that these functions rely the data
catalog information that is kept in `.data.reccap2_data.yaml`

RECCAP data is loaded with a special function that allows to browse 
the RECCAP data interactively. 


"""
from munch import Munch as _munch


def obs_surface_co2(
    data_dict, 
    specs=[    
        { # variables
            'fgco2_reg': "fgco2_reg",
            'fgco2_glob': "fgco2_glob",
            'fgco2': '(fgco2)\w(?!glob|reg)',
            'spco2': "spco2",
            'pco2atm': 'p.?co2.?atm',
            'Kw': "Kw",
            'alpha': "(alpha)\w(?!.*skin)",
            'alpha_skin': "alpha_skin",
            'alpha_subskin': "alpha_subskin",
            'fice': 'siconc|fice',
            'talk': "talk",
            'ph': "(?!al)(ph)(?!a)",
            'tos': "tos",
            'sos': "sos",
            'dissic': "dissic",
            'area': "area"},
        { # data products
            'UOEX': 'UOEX_Wat20',
            'NIES': 'NIES',
            'ETHZ': 'OceanSODAETHZ',
            'CSIR': 'CSIRML6',
            'CMEMS': 'CMEMS-LSCE-FFNN'}], 
    verbose=True, 
    **kwargs
):
    """
    kwargs are passed to data.downloads.retrieve_files
    """
    from .download import download 
    
    flist = download(**data_dict, verbose=verbose, **kwargs)
    obj = _read_reccap2_products(flist, specs)
    
    with open(data_dict['dest'] + '/info_summary.txt', 'w') as f:
        f.write(obj.__repr__())
    
    return obj


def model_surface_co2(
    data_dict, 
    specs=[
        { # variables
            'fgco2_reg': "fgco2_reg",
            'fgco2_glob': "fgco2_glob",
            'fgco2': '(fgco2)\w(?!glob|reg)',
            'spco2': 'spco2',
            'dpco2': 'dpco2',
            'alpha': 'alpha',
            'Kw': 'Kw',
            'pco2atm': 'pco2atm',
            'dissicos': 'dissicos',
            'talkos': 'talkos',
            'tos': 'tos',
            'sos': 'sos',
            'no3os': 'no3os',
            'sios': 'sios',
            'po4os': 'po4os',
            'o2os': 'o2os',
            'mld': 'mld',
            'fice': 'fice',
            'siconc': 'siconc',
            'icfriver': 'icfriver',
            'ocfriver': 'ocfriver',
            'mlotst': 'mlotst'},
        { # models - code breaks with more than two levels
            'BSOSE': 'BSOSE',
            'CNRM_A': 'CNRM-ESM2-1_A',
            'CNRM_B': 'CNRM-ESM2-1_B',
            'CNRM_C': 'CNRM-ESM2-1_C',
            'CNRM_D': 'CNRM-ESM2-1_D',
            'ECCO_A': 'ECCO-Darwin_A',
            'ETHZ_A': 'CESM-ETHZ_A',
            'ETHZ_B': 'CESM-ETHZ_B',
            'ETHZ_C': 'CESM-ETHZ_C',
            'ETHZ_D': 'CESM-ETHZ_D',
            'FESOM_A': 'FESOM_REcoM_LR_A', 
            'FESOM_B': 'FESOM_REcoM_LR_B', 
            'FESOM_C': 'FESOM_REcoM_LR_C', 
            'FESOM_D': 'FESOM_REcoM_LR_D', 
            'GEOMAR_A': 'ORCA025-GEOMAR_A',
            'GEOMAR_B': 'ORCA025-GEOMAR_B',
            'GEOMAR_C': 'ORCA025-GEOMAR_C',
            'GEOMAR_D': 'ORCA025-GEOMAR_D',
            'GOA': 'GOA-COBALT',
            'MRI': 'MRI-ESM2',
            'NorESM_A': 'NorESM-OC1.2_A',
            'NorESM_B': 'NorESM-OC1.2_B',
            'NorESM_C': 'NorESM-OC1.2_C',
            'NorESM_D': 'NorESM-OC1.2_D',
        }
    ], 
    verbose=True, 
    **kwargs
):
    """
    kwargs are passed to data.downloads.retrieve_files
    """
    from .download import download 
    
    flist = download(**data_dict, verbose=verbose, **kwargs)
    obj = _read_reccap2_products(flist, specs)
    
    with open(data_dict['dest'] + '/info_summary.txt', 'w') as f:
        f.write(obj.__repr__())
    
    return obj


def socat2020(
    data_dict, 
    variables=['fco2_ave_unwtd', 'fco2_std_unwtd'], 
    processed_dest='../data/processed/socat2020_monthlygridded.nc', 
    verbose=True, 
    **kwargs
):
    """kwargs are passed to scripts.data.download.retrieve_files"""
    import xarray as xr
    from .preprocess import preprocess 
    from .download import download, flatten_list as _flatten, read_catalog
    from pathlib import Path as path
    from xarray import open_dataset
    
    if path(processed_dest).is_file():
        return open_dataset(processed_dest)
    
    socat = xr.open_mfdataset(
        download(**data_dict, verbose=verbose, **kwargs), 
        preprocess=preprocess(decode_times=False)
    )[variables]
    
    socat.to_netcdf(processed_dest, encoding={k: dict(zlib=True, complevel=4) for k in socat})
    
    return socat


def soccom_float_liar(
    data_dict, 
    processed_dest='../data/processed/soccom_float_liar.nc', 
    verbose=True, 
    **kwargs
):
    from .non_reccap_data import grid_soccom_argo_float
    from .preprocess import preprocess 
    from .download import download, flatten_list as _flatten, read_catalog
    from pathlib import Path as path
    from xarray import open_dataset
    
    if path(processed_dest).is_file():
        return open_dataset(processed_dest)
    
    flist = _flatten(download(**data_dict, verbose=verbose, **kwargs))
    flist = [f for f in flist if f.endswith('.nc')]
    
    variables = [
        'Temperature', 'Salinity', 
        'pCO2_LIAR', 'TALK_LIAR', 'DIC_LIAR', 'TALK_LIAR_QFA', 'DIC_LIAR_QFA']

    liar = grid_soccom_argo_float(flist, keep_vars=variables, agg_funcs=['mean'])

    qc = (liar.DIC_LIAR_QFA == 0) & (liar.TALK_LIAR_QFA == 0)
    liar = liar.where(qc).drop(['DIC_LIAR_QFA', 'TALK_LIAR_QFA'])

    liar = preprocess()(liar)
    
    liar = liar.rename(
        Salinity='so',
        Temperature='thetao',
        pCO2_LIAR='pco2',
        TALK_LIAR='talk',
        DIC_LIAR='dissic',
    )
    
    liar.to_netcdf(processed_dest, encoding={k: dict(zlib=True, complevel=4) for k in liar})

    return liar


class _RECCAP_dict(_munch):
        
    @property
    def data(self, dim_name='variable'):
        if hasattr(self, '_data'):
            return self._data
        
        from xarray import concat, merge, DataArray
        from pandas import Index
        is_array = [isinstance(self[k], DataArray) for k in self]
        is_munch = [isinstance(self[k], self.__class__) for k in self]
        if all(is_array):
            idx = Index(self.keys(), name=dim_name)
            try:
                xda = concat([self[k] for k in self], dim=idx)
            except:
                xda = merge([self[k] for k in self])
            self._data = xda
            return self._data
        else:
            dataarrays = []
            for k in self:
                if isinstance(self[k], self.__class__):
                    dataarrays += self[k].data,
            xds = merge(dataarrays)
            self._data = xds
            return self._data
        
    def __repr__(self):
        def get_info(obj):
            unit = obj.attrs.get('units', '-')
            dims = ".".join([f"{k}({obj[k].size})" for k in obj.dims])
            hist = obj.attrs.get('processing', '')
            
            unit = f"[UNITS] {unit} "
            dims = f"[DIMS] {dims}"
            hist = (
                f"[PREPROCESSING] {hist}"
                .replace(": ", "=")
                .replace("{", "[")
                .replace("}", "]"))
            info = f"{unit: <25}{dims: <40}{hist}"
            return info
        def walk_through_dictionary(d):
            new_dict = {}
            for k, v in d.items():
                if isinstance(v, dict):
                    new_dict[f"{k} {'-'*(100 - len(k))}"] = walk_through_dictionary(v)
                elif isinstance(v, xr.DataArray):
                    new_dict[f"{k: <10}"] = get_info(v)
            return new_dict
        import json
        import xarray as xr
        import re
        import os
        
        printable = walk_through_dictionary(self)
        pretty = f"{str(self.__class__)}"
        pretty += json.dumps(printable, indent=2, sort_keys=True)
        pretty = re.sub('[",:{}]', '', pretty)
        
        m = self.get('not_matched', [])
        if m != []:
            s = "\n    "
            not_matched = s + s.join([os.path.relpath(f) for f in m])
            pretty += f"  not_matched {'='*69}{not_matched}"

        return pretty

        
def _read_reccap2_products(flist, fname_specs=[]):
    def flatten_list(list_of_lists):
        if len(list_of_lists) == 0:
            return list_of_lists
        if isinstance(list_of_lists[0], list):
            return flatten_list(list_of_lists[0]) + flatten_list(list_of_lists[1:])
        return list_of_lists[:1] + flatten_list(list_of_lists[1:])
    
    def build_tree(tree_list):
        if tree_list:
            if len(tree_list) > 2:
                return {tree_list[0]: build_tree(tree_list[1:])}
            else:
                try:
                    xds = xr.open_mfdataset(tree_list[1], 
                                            decode_times=False, 
                                            preprocess=preprocess(decode_times=True, center_months=True))
                except OSError: 
                    return {}
                if len(xds.data_vars) == 1:
                    xds = xds[list(xds.data_vars.keys())[0]]
                elif len(xds.data_vars) > 1:
                    for key in xds:
                        if key in tree_list[1]:
                            xds = xds[key].assign_attrs(processing=xds.attrs.get('processing', ''))
                            break
                return {tree_list[0]: xds}
        return {}
    
    import re
    import xarray as xr
    import mergedeep
    from munch import Munch
    from .preprocess import preprocess
    
    flist = flatten_list(flist)

    output = []
    failed = []
    for f in flist:
        if f.endswith('.nc'):
            matches = []
            for spec in fname_specs:
                matches += [m for m in spec if re.findall(spec[m] if spec[m] else m, f)]
            if len(matches) == len(fname_specs):
                tree = [re.sub("[^0-9a-zA-Z]+", "_", m) for m in matches] + [f]
                output += build_tree(tree),
            else:
                failed += f,
    merged = mergedeep.merge({}, *output)            
    merged['not_matched'] = failed
    obj = _RECCAP_dict.fromDict(merged)
    
    return obj



