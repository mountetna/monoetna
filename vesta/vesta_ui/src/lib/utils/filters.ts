import { DATA_TYPES } from '@/lib/fixtures';
import { Project, ProjectDataType } from '@/components/project-explorer/models';

export const projectDataTypes = (project:Project):ProjectDataType[] => [ ...new Set(
  project.dataTypes.map(
    (modelName:string) => Object.keys(DATA_TYPES).find(
      (dt:string) => DATA_TYPES[dt as ProjectDataType].includes(modelName)
    )
  ).filter((_:any)=>_) as ProjectDataType[]
) ]
