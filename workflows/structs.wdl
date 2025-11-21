version 1.0

struct PGENFileset {
    String chr
    File pgen
    File psam
    File pvar
}

struct PLINKFileset {
    String chr
    File bed
    File bim
    File fam
}

struct RegenieStep1Files {
    Array[File] phenotypes_loco
    File list
}