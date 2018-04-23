setClass("genealogy",representation(history="matrix", N="integer", Ngen="integer", model="character"))
setClass("samplegenealogy", representation(ancestry="matrix", samplesize="integer", sample="integer", MRCA="integer"),contains="genealogy")
