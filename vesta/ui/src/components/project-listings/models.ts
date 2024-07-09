import { ThemeData } from "../themes/models";


export enum ProjectStatus {
    team = 'Team',
    community = 'Community',
}

export enum ExternalProjectStatus {
    init = 'Initial Launch',
    sampling = 'Sampling',
    analysis = 'Analysis',
    public = 'Public',
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
}

export interface PrincipalInvestigator {
    name: string
    title?: string
    imageUrl?: string
    color: string
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
    species: string
    startDate: Date
    dataCollectionComplete: boolean
    userCount: number

    theme: ThemeData
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