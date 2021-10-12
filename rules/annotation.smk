# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 14:00:20 2021

@author: cow082
"""

BASEDIR = workflow.basedir # this is full path to the snakefile itself


rule all:
    input:
        "graphics/rulegraph.png"
        

rule graph_workflow:
    input:
        snakefile="%s/workflow.snk" % BASEDIR
    output:
        imgfile="graphics/rulegraph.png"
    run:
        cmd = utils.graph_workflow(input.snakefile, output.imgfile, config)
        shell(cmd)

