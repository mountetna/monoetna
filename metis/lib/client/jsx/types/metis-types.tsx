export type Folder = {
  folder_name: string;
  folder_path: string;
  bucket_name: string;
  project_name: string;
  read_only: boolean;
};

export type File = {
  file_name: string;
  file_path: string;
  bucket_name: string;
  project_name: string;
  read_only: boolean;
  download_url: string;
};

export type UiControlItem = {
  label: string;
  callback: () => void;
  role: string;
};
