# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 14:00:20 2021

@author: cow082
"""

BASEDIR = workflow.basedir # this is full path to the snakefile itself


rule annotation_all:
    input:
        expand("02_irma_assembly/{sample}/annotation/", sample=config["samples"])
        

rule graph_workflow:
    input:
        snakefile="%s/workflow.snk" % BASEDIR
    output:
        imgfile="graphics/rulegraph.png"
    run:
        cmd = utils.graph_workflow(input.snakefile, output.imgfile, config)
        shell(cmd)


rule run_annotation:
    message:
        """
        ** annotation **
        Blasting and finding cleavage site
        """
    input:
       run_IRMA = "02_irma_assembly/{sample}/IRMA_COMPLETE",
    params:
       IRMA_dir = "02_irma_assembly/{sample}/irma_output/",
       email = config["blast_email"],
       database = config["blast_flu"],
    output:
       output_dir = directory("02_irma_assembly/{sample}/annotation/")
    shell:
        """
		mkdir 02_irma_assembly/{wildcards.sample}/annotation
		
        python {config[program_dir]}annotation/rules/analyse.py \
            -input_dir {params.IRMA_dir} \
            -db {params.database} \
            -output_dir {output.output_dir} \
            -email {params.email}
        """