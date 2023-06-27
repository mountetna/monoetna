import {json_post} from 'etna-js/utils/fetch';
import { metisPath } from 'etna-js/api/metis_api';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogContentText from '@material-ui/core/DialogContentText';
import DialogTitle from '@material-ui/core/DialogTitle';

import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import { dateFormat, authorFormat } from 'etna-js/utils/format';
import {makeStyles, Theme} from '@material-ui/core/styles';
import React, {
  useState,
  useEffect,
  useCallback,
  useContext,
  useMemo
} from 'react';

const useStyles = makeStyles((theme: Theme) => ({
  container: {
    maxHeight: '50vh'
  }
}));

type File = {
  file_name: string;
  file_path: string;
  bucket_name: string;
  project_name: string;
  read_only: boolean;
  download_url: string;
};

const TestFileMatch = ({projectName,bucketName,script}:{
  projectName: string;
  bucketName: string;
  script: Script;
}) => {
  const [ files, setFiles ] = useState<File[]>([]);
  const [ showFiles, setShowFiles ] = useState<boolean>(false);

  const classes = useStyles();

  const testMatch = () => {
    // with the find api we can list files matching the given glob
    json_post(metisPath(`${projectName}/find/${bucketName}`), {
      params: [{
        attribute: 'name',
        predicate: 'glob',
        type: 'file',
        value: script.folder_path + '/**/' + script.file_match }]
    }).then(({files}) => setFiles(files)); 
    setShowFiles(true);
  }
  const touchFiles = () => {
    json_post(metisPath(`/${projectName}/tail/${bucketName}`)); 
  }
  return <Grid item xs={2}>
    <Button onClick={ testMatch }>Test</Button>
    <Dialog fullWidth maxWidth='lg' open={ showFiles } onClose={ () => setShowFiles(false) }>
      <DialogTitle>Matching Files</DialogTitle>
      <DialogContent>
        <TableContainer className={classes.container}>
	  <Table stickyHeader>
	    <TableHead>
	      <TableRow>
		<TableCell>File Path</TableCell>
		<TableCell align="right">Updated At</TableCell>
	      </TableRow>
	    </TableHead>
	    <TableBody>
	      {
		files.map(
		  file => <TableRow key={file.id}>
                    <TableCell>{file.file_path}</TableCell>
                    <TableCell align="right">{dateFormat(file.updated_at)}</TableCell>
                  </TableRow>
                )
              }
	    </TableBody>
	  </Table>
	</TableContainer>
      </DialogContent>
      <DialogActions>
        <Button onClick={() => setShowFiles(false)}>
          Cancel
        </Button>
        <Button onClick={touchFiles}>
          Touch
        </Button>
      </DialogActions>
    </Dialog>
  </Grid>
};

export default TestFileMatch;
