from ..Geoline.seqTemplate import *

def sampleMapAttr(assay: str) -> Dict:
    samples = samplesSection(title=f'{assay}:tube_name',
                             sourceName=f'{assay}:biospecimen',
                             organism='subject:name',
                             characteristics=characteristics(addAnother, {}),
                             molecule='',
                             description=f'{assay}:notes',
                             processedDataFile=f'{assay}:gene_expression',
                             rawFile={
                                 'raw file': f'{assay}:raw_fastq',
                                 'raw file ': f'{assay}:raw_fastqs',
                             })
    return samples

def processedFilesMapAttr(assay: str) -> Dict:
    pass
    #processedFiles = processedFilesSection(fileName=f'{assay}:',
