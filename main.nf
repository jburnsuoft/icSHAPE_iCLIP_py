#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters
params.clipfile        = 'ZRANB2.bam.bed'
params.shape_range     = 25
params.cluster_distance= 0
params.plus_shape      = 'HEK293_NP_VIVO_PLUS_HG19.bed'
params.minus_shape     = 'HEK293_NP_VIVO_MINUS_HG19.bed'
params.out_fasta       = 'ZRANB2_NEW25.bed_VIVO_sequences_fasta.fa'
params.seq_folder      = 'ZRANB2_NEW25.bed_SEQ_VIVO_20_0/'
params.shape_folder    = 'ZRANB2_NEW25.bed_FULL_SHAPE_VIVO_20_0/'
params.ct_folder       = 'ZRANB2_NEW25.bed_CT_VIVO_20_0/'
params.db_folder       = 'ZRANB2_NEW25.bed_FULL_DOT_VIVO_20_0/'
params.pname           = 'example_pname'
params.vitroShapeFolder= 'example_vitroShapeFolder'
params.vivoShapeFolder = 'example_vivoShapeFolder'
params.plFolder        = 'example_plFolder'
params.plshuffle       = 'example_plshuffle'

// Process to run CLIPSHAPE
process runCLIPSHAPE {
    output:
        file("${params.out_fasta}") into out_fasta

    script:
    """
    python /workspaces/icSHAPE_iCLIP/CLIPSHAPE_PY/CLIPSHAPE \\
      ${params.clipfile} \\
      ${params.shape_range} \\
      ${params.cluster_distance} \\
      ${params.plus_shape} \\
      ${params.minus_shape} \\
      ${params.out_fasta} \\
      ${params.seq_folder} \\
      ${params.shape_folder} \\
      ${params.ct_folder} \\
      ${params.db_folder}
    """
}

// Process to run MAPSHAPE
process runMAPSHAPE {
    input:
    val pname
    val vitroShapeFolder
    val vivoShapeFolder
    val plFolder
    val plshuffle

    output:
    file("${pname}_median_plot.pdf")

    script:
    """
    python MAPSHAPE.py --pname $pname --vitroShapeFolder $vitroShapeFolder --vivoShapeFolder $vivoShapeFolder --plFolder $plFolder --plshuffle $plshuffle
    """
}

// Workflow definition
workflow {
    runCLIPSHAPE()
    runMAPSHAPE(params.pname, params.vitroShapeFolder, params.vivoShapeFolder, params.plFolder, params.plshuffle)
}
