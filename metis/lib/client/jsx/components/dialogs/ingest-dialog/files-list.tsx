import React from 'react';

import FormControlLabel from '@material-ui/core/FormControlLabel';
import Checkbox from '@material-ui/core/Checkbox';
import Typography from '@material-ui/core/Typography';
import {CheckBoxOutlineBlankSharp} from '@material-ui/icons';

export type IngestFile = {
  name: string;
  host: string;
};

const FilesList = ({
  allFiles,
  selectedFileNames,
  onChangeFileNames
}: {
  allFiles: IngestFile[];
  selectedFileNames: string[];
  onChangeFileNames: (newFileNames: string[]) => void;
}) => {
  function onSelectAllClick() {
    onChangeFileNames(
      selectedFileNames.length === allFiles.length
        ? []
        : allFiles.map((f) => f.name)
    );
  }

  function onClickSingleFile(fileName: string) {
    let updatedFiles = [...selectedFileNames];
    if (updatedFiles.includes(fileName)) {
      updatedFiles = updatedFiles.filter((f) => f !== fileName);
    } else {
      updatedFiles.push(fileName);
    }

    onChangeFileNames(updatedFiles);
  }

  if (0 === allFiles.length) return null;

  return (
    <React.Fragment>
      <FormControlLabel
        control={
          <Checkbox
            indeterminate={
              selectedFileNames.length > 0 &&
              selectedFileNames.length < allFiles.length
            }
            checked={
              selectedFileNames.length > 0 &&
              selectedFileNames.length === allFiles.length
            }
            onChange={onSelectAllClick}
            inputProps={{'aria-label': 'select all files'}}
          />
        }
        label='All'
      />
      {allFiles.map((file: IngestFile, index: number) => (
        <FormControlLabel
          key={index}
          control={
            <Checkbox
              checked={selectedFileNames.includes(file.name)}
              onChange={() => {
                onClickSingleFile(file.name);
              }}
            />
          }
          label={file.name}
        />
      ))}
    </React.Fragment>
  );
};

export default FilesList;
