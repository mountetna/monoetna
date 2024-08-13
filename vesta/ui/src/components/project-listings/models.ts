import { ThemeData } from "../themes/models";


export enum ProjectStatus {
    team = 'Team',
    community = 'Community',
}

export enum ProjectType {
    coproject = 'CoProject',
    consortium = 'Consortium',
    copilot = 'CoPilot',
    collab = 'Collaboration',
    external = 'External',
}

export enum DataType {
    cytof = 'CyTOF',
    bulkRnaSeq = 'Bulk RNA-Seq',
    citeSeq = 'CITE-Seq',
    tcrBcrRepertoire = 'TCR-BCR Repertoire',
    zeNith = 'ZeNITH',
    hlaSequencing = 'HLA Sequencing',
    phIpSeq = 'PhIP-Seq',
    wholeExomeSeq = 'Whole-exome seq',
    atacSeq = 'ATAC-seq',
    visiumSpatialTranscriptomics = 'Visium Spatial Transcriptomics',
    organoidExpansion = 'Organoid expansion',
    vdjSeq = 'V(D)Jseq',
    microbiome = 'Microbiome',
    tcrSeq = 'TCR-seq',
}

export interface PrincipalInvestigator {
    name: string
    title?: string
    imageUrl?: string
    color: string
    altColor: string
}

export interface Project {
    name: string
    fullName: string
    heading?: string
    description: string
    fundingSource: string
    principalInvestigators: PrincipalInvestigator[]
    status: ProjectStatus
    type: ProjectType
    dataTypes: DataType[]
    hasSamples: boolean
    hasAssays: boolean
    hasClinicalData: string,
    species: string
    startDate: Date
    dataCollectionComplete: boolean
    userCount: number

    theme: ThemeData
}

export enum ExternalProjectStatus {
    init = 'Initial Launch',
    sampling = 'Sampling',
    analysis = 'Data Analysis',
    public = 'Public',
}

export function getExternalProjectStatus(proj: Project): ExternalProjectStatus {
    if (proj.status === ProjectStatus.community) {
        return ExternalProjectStatus.public
    } else if (proj.hasAssays) {
        return ExternalProjectStatus.analysis
    } else if (proj.hasSamples) {
        return ExternalProjectStatus.sampling
    } else {
        return ExternalProjectStatus.init
    }
}

export enum ProjectHeadingInfoSet {
    default = 'Default',
    dataTypes = 'Data Types',
    pis = 'Principal Investigators',
}

export const PROJECTS_SEARCH_PARAMS_KEY = 'projects'

export interface ProjectsSearchParamsControls {
    viewSet?: string
    page?: number
}

export interface ProjectsSearchParamsState {
    filters?: Record<string, string[]>
    controls?: ProjectsSearchParamsControls
}