export interface UserFlagsInterface {
  email: string;
  name: string;
  projects: string[];
  flags: string[] | null;
}

export interface Project {
  project_name: string;
  project_name_full: string;
}
