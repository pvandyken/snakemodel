rule prepare_nbs_connectomes:
    input:
        connectomes=_expand_subjects_from_group(
            inputs["connectome"].path,
        ),
        metadata=resource("atlases/atlas-brainnetome246/labels.tsv"),
    output:
        tempout(f"prepare_nbs_connectomes/{outdir}", None, ".mat")
    log: log(f"prepare_nbs_connectomes/{outdir}", None)
    benchmark: benchmark(f"prepare_nbs_connectomes/{outdir}", None)
    envmodules:
        'python/3.10'
    group: 'nbs'
    threads: 1
    resources:
        mem_mb=1000,
        runtime=2,
    params:
        hemi=lambda wcards: wcards.get('hemi', '-')
    shell:
        boost(
            design_matrix_env.script,
            pyscript(
                "scripts/nbsprep.py",
                input=["connectomes", "metadata"],
                params=["hemi"]
            ),
        )


rule run_nbs:
    input:
        connectomes=rules.prepare_nbs_connectomes.output,
        mat=rules.make_patsy_dm.output["mat"],
        con=rules.split_contrast_matrix.output[0],
    output:
        matrix=temp(Path(config['output_dir'], outdir, "nbs_{contrast_label}.mat"))
    log: log(f"run_nbs/{outdir}/{config.get('desc', '')}", None, "contrast_label")
    benchmark: benchmark(f"run_nbs/{outdir}/{config.get('desc', '')}", None, "contrast_label")
    container: "docker://pvandyken/nbs:1.2"
    group: 'nbs'
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
    params:
        test="t-test",
        thresh=3.0,
        perms=10_000,
        alpha = 0.05,
    shell:
        """
        /app/.venv/bin/nbs --test {params.test} --thresh {params.thresh} \\
            --perms {params.perms} --alpha {params.alpha} \\
            --contrast {input.con} --design {input.mat} \\
            --matrices {input.connectomes} {output} \\
            > {log}
        """

rule post_nbs:
    input: 
        matrix=rules.run_nbs.output[0]
    output:
        matrix=bids(
            output / outdir,
            model="nbs",
            label="{contrast_label}",
            **(
                {"desc": config["desc"]}
                if "desc" in config
                else {}
            ),
            suffix="nbs.hdf5",
        )
    shell:
        boost(
            design_matrix_env.script,
            pyscript(
                "scripts/post_nbs.py",
                input=["matrix"],
                output=["matrix"],
            ),
        )