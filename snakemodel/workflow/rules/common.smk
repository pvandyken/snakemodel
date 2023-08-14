def get_entity_filters(filters):
    return {
        key: val for key, val in filters.items()
        if key not in {"scope", "extension"} and "_" not in key
    }

def log_transform(*args):
    def _log_transform(x):
        return np.log10(x)

    return preproc.FunctionTransformer(_log_transform)


def impute_transform(df, field, strategy, fill=None):
    return impute.SimpleImputer(
        strategy=strategy,
        fill_value=(
            df[field].dtype.type(fill) if fill is not None else None
        )
    )



@ft.lru_cache(None)
def _get_participants():
    filters = {
        filt[0]: filt[1:] for filt in config.get("filter_participants", []) or []
    }
    df = filter_participants(
        Path(config['bids_dir'], 'participants.tsv'),
        **filters,
    )
    df = df[
        df["participant_id"].isin(
            map(lambda s: f"sub-{s}", inputs[comp].entities["subject"])
        )
    ]

    if config["transform"]:
        transforms = {
            "log": log_transform,
            "impute": impute_transform
        }
        for col in config["transform"]:
            txfs = [shlex.split(txf) for txf in col[1:]]
            pp = pipeline.make_pipeline(
                *(transforms[txf[0]](df, col[0], *txf[1:]) for txf in txfs)
            )
            transformed = compose.make_column_transformer(
                (pp, [col[0]])
            ).fit_transform(df)
            df[col[0]] = transformed[:, 0]


    # if config["impute"]:
    #     for imp in config["impute"]:
    #         fill = itx.nth(imp, 2, None)
    #         imputed = compose.make_column_transformer(
    #             (
    #                 impute.SimpleImputer(
    #                     strategy=imp[1],
    #                     fill_value=(
    #                         df[imp[0]].dtype.type(fill) if fill is not None else None
    #                     )
    #                 ),
    #                 [imp[0]],
    #             )
    #         ).fit_transform(df)
    #         df[imp[0]] = imputed[:, 0]

    # terms = _get_terms(df)
    # print(terms)
    # if config["skip_nan"]:
    #     return df[~df.isna().any(axis=1)]

    # for term in terms:
    #     if np.any(pd.isnull(df.get(term, []))):
    #         raise ConfigError(
    #             f"Found nan in column '{term}'. NaNs must either be explicitly "
    #             "filtered using --skip-nan, or an imputation strategy must be used"
    #         )
    return df.sort_values(
        by='participant_id',
    )

@ft.lru_cache(None)
def _unpack_contrasts():
    contrasts, *_params = it.zip_longest(*config["ttest"], fillvalue=None)
    params = [
        dict(
            param.split("=", 1) for param in row if param is not None
        )
        for row in zip(*_params)
    ]
    return contrasts, params or [{} for _ in range(len(contrasts))]

def _get_terms(df):
    return dmatrix(config["design"], df).design_info.term_names

def _expand_subjects_from_group(path, allow_missing=False, **kwargs):
    def inner(wcards):
        filtered = inputs[comp].filter(
            subject=_get_design_ids().map(lambda s: s[4:])
        )
        if len(filtered.zip_lists["subject"]) != len(filtered.entities["subject"]):
            raise ConfigError(
                "Some subjects have multiple data files associated. Use filters to "
                "ensure each subject only appears once"
            )
        if _get_design_matrix().shape[0] != len(filtered.expand()):
            raise ConfigError(
                "Number of data files different than number of rows in design matrix. "
                "You may need to add additional filters"
            )
        return expand(
            expand(
                path,
                zip,
                **pd.DataFrame(filtered.zip_lists).sort_values(by="subject"),
            ),
            allow_missing=allow_missing,
            **kwargs
        )
    return inner


def _get_design_matrix():
    df = _get_participants()
    na_action = 'drop' if config['skip_nan'] else 'raise'
    return dmatrix(config["design"], df, NA_action=na_action)

def _get_design_ids():
    df = _get_participants().set_index('participant_id')
    na_action = 'drop' if config['skip_nan'] else 'raise'
    return dmatrix(
        config["design"], df, NA_action=na_action, return_type="dataframe"
    ).index

def _get_contrast_matrix(dm):
    contrasts, _ = _unpack_contrasts()
    contrast = ", ".join(contrasts)
    return dm.design_info.linear_constraint(contrast).coefs

class Counter:
    def __init__(self):
        self.data = {
            "": 0
        }
    
    def next(self, item):
        if item not in self.data:
            self.data[item] = 0
            return item
        self.data[item] += 1
        return f"{item}{self.data[item]}"

def _get_contrast_labels(dm):
    labels = []
    counter = Counter()
    for contrast, params in zip(*_unpack_contrasts()):
        num_contrasts = dm.design_info.linear_constraint(contrast).coefs.shape[0]
        name = params.get("name", "")
        for _ in range(num_contrasts):
            labels.append(counter.next(name))
    return labels


rule test_patsy_dm:
    input:
        Path(config['bids_dir'], 'participants.tsv'),
    run:
        dm = _get_design_matrix()
        print(dm.design_info.column_names)
        print(dm)
        con = _get_contrast_matrix(dm)
        print("Contrast")
        print(repr(con))

rule make_patsy_dm:
    input:
        Path(config['bids_dir'], 'participants.tsv'),
    output:
        mat=tempout("make_patsy_design_matrix", None, ".mat.txt"),
        con=tempout("make_patsy_design_matrix", None, ".con.txt"),
    run:
        dm = _get_design_matrix()
        np.savetxt(output.mat, dm, fmt="%i")
        np.savetxt(output.con, _get_contrast_matrix(dm), fmt="%i")


rule split_contrast_matrix:
    input:
        rules.make_patsy_dm.output.con,
    output:
        tempout("split_contrast_matrix", None, ".con.txt", "contrast_label"),
    params:
        contrast_idx=lambda wcards: (
            _get_contrast_labels(_get_design_matrix()).index(wcards["contrast_label"])
            + 1
        )
    shell:
        "sed '{params.contrast_idx}q;d' {input} > {output}"


rule convert_design_matrices_to_fsl:
    input:
        mat=rules.make_patsy_dm.output["mat"],
        con=rules.split_contrast_matrix.output[0]
    output:
        mat=tempout("convert_design_matrices_to_fsl", None, ".mat", "contrast_label"),
        con=tempout("convert_design_matrices_to_fsl", None, ".con", "contrast_label"),
    log: log(f"convert_design_matrices_to_fsl/{outdir}", None, "contrast_label")
    benchmark: benchmark(f"convert_design_matrices_to_fsl/{outdir}", None, "contrast_label")
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "fsl/6.0.4"
        
    threads: 1
    shell:
        """
        Text2Vest {input.mat} {output.mat}
        Text2Vest {input.con} {output.con}
        """
