rule create_average_image:
    input:
        _expand_subjects_from_group(
            inputs[fa].path
        ),
    output:
        mean=source / outdir / "mean.nii.gz",
        mask=source / outdir / "mask.nii.gz",
    log: log("create_average_image", None)
    benchmark: benchmark("create_average_image", None)
    group: 'tbss'
    threads: 1
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "fsl/6.0.4",
    resources:
        mem_mb=6000,
        runtime=5,
    shadow: 'minimal'
    shell:
        """
        fslmerge -t stack.nii.gz {input}
        fslmaths stack.nii.gz -max 0 -Tmin -bin mask.nii.gz -odt char
        fslmaths stack.nii.gz -mas mask.nii.gz stack.nii.gz
        fslmaths stack.nii.gz -Tmean {output.mean}
        fslmaths mask.nii.gz -Tmean {output.mask}
        """

rule skeletonize_average_image:
    input:
        rules.create_average_image.output.mean
    output:
        tempout("skeletonize_average_image", None, ".nii.gz")
    log: log("skeletonize_average_image", None)
    benchmark: benchmark("skeletonize_average_image", None)
    group: 'tbss'
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "fsl/6.0.4",
    threads: 1
    resources:
        mem_mb=1000,
        runtime=1,
    shell:
        "tbss_skeleton -i {input} -o {output}"


rule threshold_skeletonized_image:
    input:
        rules.skeletonize_average_image.output
    output:
        tempout("threshold_skeletonized_image", None, ".nii.gz")
    log: log("threshold_skeletonized_image", None)
    benchmark: benchmark("threshold_skeletonized_image", None)
    group: "tbss"
    threads: 1
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "fsl/6.0.4",
    resources:
        mem_mb=1000,
        runtime=1,
    params:
        threshold=0.2
    shell:
        "fslmaths {input} -thr {params.threshold} -bin {output}"

rule get_distance_map:
    input:
        skeleton=rules.threshold_skeletonized_image.output,
        mask=rules.create_average_image.output.mask,
    output:
        tempout("get_distance_map", None, ".nii.gz")
    log: log("get_distance_map", None)
    benchmark: benchmark("get_distance_map", None)
    group: 'tbss'
    threads: 1
    resources:
        mem_mb=1000,
        runtime=1,
    shadow: 'minimal'
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "fsl/6.0.4",
    shell:
        """
        fslmaths {input.mask} -mul -1 -add 1 -add {input.skeleton} \
            dist_mask.nii.gz
        distancemap -i dist_mask.nii.gz -o {output}
        """

_skeleton_label = config.get(
    "desc",
    config["pybids_inputs"][comp]
    .get('filters', {})
    .get('suffix', "data")
)

rule project_onto_skeleton:
    input:
        fa=inputs[fa].path,
        distance=rules.get_distance_map.output,
        mean=rules.create_average_image.output.mean,
        mask=rules.create_average_image.output.mask,
        **(
            {"data": inputs[comp].path} if comp is not fa else {}
        )
    output:
        bids(
            source,
            **{
                **get_entity_filters(config['pybids_inputs'][comp]["filters"]),
                **inputs[comp].wildcards,
                "desc": "skeletonized",
                "suffix": _skeleton_label  + ''.join(Path(inputs[comp].path).suffixes)
            }
        )
    log: log(f"project_onto_skeleton/{_skeleton_label}", inputs[comp])
    benchmark: benchmark(f"project_onto_skeleton{_skeleton_label}", inputs[comp])
    group: 'tbss'
    threads: 1
    resources:
        mem_mb=500,
        runtime=1,
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "fsl/6.0.4",
    params:
        threshold=0.2,
        param_map=lambda wcards, input: str(input.data) if comp is not fa else ""
    shadow: 'minimal'
    shell:
        """
        fslmaths "{input.fa}" -mas "{input.mask}" fa_masked.nii.gz
        if [[ -n "{params.param_map}" ]]; then
            fslmaths "{params.param_map}" -mas "{input.mask}" map_masked.nii.gz
            param_map="-a map_masked.nii.gz"
        else
            param_map=
        fi

        tbss_skeleton -i {input.mean} -p {params.threshold} {input.distance} \\
            ${{FSLDIR}}/data/standard/LowerCingulum_1mm fa_masked.nii.gz \\
            {output} $param_map
        """


rule skeleton_glm:
    input:
        data=_expand_subjects_from_group(
            rules.project_onto_skeleton.output,
        ),
        design=rules.convert_design_matrices_to_fsl.output.mat,
        contrast=rules.convert_design_matrices_to_fsl.output.con,
        mask=rules.create_average_image.output.mask,
    output:
        tstat=bids(
            output / outdir,
            label="{contrast_label}",
            model="tbss",
            desc="{parameter}",
            suffix="tstat.nii.gz",
        ),
        corrp=bids(
            output / outdir,
            label="{contrast_label}",
            model="tbss",
            desc="{parameter}",
            suffix="corrp.nii.gz",
        ),
    log: log(f"skeleton_glm/{outdir}", None, "contrast_label", "parameter")
    benchmark: benchmark(f"skeleton_glm/{outdir}", None, "contrast_label", "parameter")
    params:
        permutations=10000,
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "fsl/6.0.4",
    # Need shallow so that log file copies over
    shadow: 'shallow'
    shell:
        """
        fslmerge -t merged.nii.gz {input.data}

        randomise -i merged.nii.gz -o out \\
            -d {input.design} -t {input.contrast} -m {input.mask} \\
            -n {params.permutations} --T2 |
            tee {log} |
            awk -v n=100 'NR%n==5 {{ print "[{wildcards.contrast_label}] " $0 }}'

        mv "out_tstat1.nii.gz" {output.tstat}
        mv "out_tfce_corrp_tstat1.nii.gz" {output.corrp}
        """

rule tbss_fill:
    input:
        corrp=rules.skeleton_glm.output["corrp"],
        mean=rules.create_average_image.output["mean"],
    output:
        filled=bids(
            output / outdir,
            label="{contrast_label}",
            model="tbss",
            desc="{parameter}",
            suffix="filled.nii.gz",
        ),
    log: log(f"tbss_fill/{outdir}", None, "contrast_label", "parameter")
    benchmark: benchmark(f"tbss_fill/{outdir}", None, "contrast_label", "parameter")
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "fsl/6.0.4",
    shell:
        "tbss_fill {input.corrp} 0.95 {input.mean} {output}"
    