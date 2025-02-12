nextflow.enable.dsl=2

params.pname             = "MyProtein"           // modify as needed
params.vitroShapeFolder  = "/path/to/vitro"      // modify as needed
params.vivoShapeFolder   = "/path/to/vivo"       // modify as needed
params.plFolder          = "/path/to/pl"         // modify as needed
params.plshuffle         = "/path/to/plshuffle"  // modify as needed

process runCLIPSHAPE {
    output:
        file "${params.pname}_median_plot.pdf"

    script:
        """
        python CLIPSHAPE.py --pname ${params.pname} \
            --vitroShapeFolder ${params.vitroShapeFolder} \
            --vivoShapeFolder ${params.vivoShapeFolder} \
            --plFolder ${params.plFolder} \
            --plshuffle ${params.plshuffle}
        """
}

workflow {
    runCLIPSHAPE()
}
