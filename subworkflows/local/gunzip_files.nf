
include { GUNZIP_LOCAL } from '../../modules/local/gunzip_local'



workflow GUNZIP_FILES {
    take:
    input_file

    

    main:

    GUNZIP_LOCAL( input_file )
    GUNZIP_LOCAL.out.gunzip.view()

    // emit:
    // temp
    
}




