# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 14:00:20 2021

@author: cow082
"""

BASEDIR = workflow.basedir # this is full path to the snakefile itself


#rule annotation_all:
#    input:
#        expand("03_annotation/{sample}/ANNOTATION_COMPLETE", sample=config["samples"])
        

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
       known_cleavage_sites = config["cleavage_sites"],
       output_dir = "03_annotation/{sample}/"
    output:
       "03_annotation/{sample}/ANNOTATION_COMPLETE"
    shell:
        """
        # need to make a directory for the python script to put the output
        # -p makes the command create multiple directory levels at once
        mkdir -p 03_annotation/{wildcards.sample}/
        
        python {config[program_dir]}annotation/rules/analyse.py \
            -input_dir {params.IRMA_dir} \
            -db {params.database} \
            -known_sites {params.known_cleavage_sites} \
            -output_dir {params.output_dir} \
            -email {params.email}
            
        if [ "$?" == "0" ]; then
            touch 03_annotation/{wildcards.sample}/ANNOTATION_COMPLETE
        fi
        """
        








        
