import {json_post,json_error} from 'etna-js/utils/fetch';
import { metisPath } from 'etna-js/api/metis_api';
import { MetisFile, Script } from '../polyphemus';
import Typography from '@material-ui/core/Typography';
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

const TestFileMatch = ({projectName,bucketName,script,className}:{
  projectName: string;
  bucketName: string;
  script: Script;
  className: string;
}) => {
  const [ files, setFiles ] = useState<MetisFile[]>([]);
  const [ showFiles, setShowFiles ] = useState<boolean>(false);
  const [ confirmTouch, setConfirmTouch ] = useState<boolean>(false);
  const [ error, setError ] = useState<string|null>(null);

  const hideFiles = () => {
    setShowFiles(false);
    setError(null);
  }

  const classes = useStyles();

  const touchPath = script.folder_path + '/**/' + script.file_match;

  const testMatch = () => {
    // with the find api we can list files matching the given glob
    json_post(metisPath(`${projectName}/find/${bucketName}`), {
      params: [{
        attribute: 'name',
        predicate: 'glob',
        type: 'file',
        value: touchPath
      }]
    }).then(
      ({files}) => setFiles(files)
    ).catch(
      json_error( error => setError(error))
    );
    setShowFiles(true);
  }

  const touchFiles = useCallback(() => {
    json_post(
      metisPath(`/${projectName}/file/touch/${bucketName}`),
      { file_paths: files.map(f => f.file_path) }
    ).then(
      ({files}) => setFiles(files)
    ).catch(
      json_error( error => setError(error))
    );
    setConfirmTouch(false);
  }, [ files ]);

  return <Grid className={className} item xs={2}>
    <Button onClick={ testMatch }>Match</Button>
    <Dialog fullWidth maxWidth='lg' open={ showFiles } onClose={ hideFiles }>
      <DialogTitle>Matching Files</DialogTitle>
      <DialogContent>
        {
          error && <DialogContentText><Typography color="error">{error}</Typography></DialogContentText>
        }
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
        <Button onClick={() => setConfirmTouch(true)} color="primary">
          Touch
        </Button>
        <Button onClick={hideFiles} color="secondary">
          Cancel
        </Button>
      </DialogActions>
    </Dialog>
    <Dialog
      open={confirmTouch}
      onClose={() => setConfirmTouch(false)}
    >
      <DialogTitle>Touch files</DialogTitle>
      <DialogContent>
        <DialogContentText>
	  Set the updated_at timestamp for files matching <em>{touchPath}</em>?
	</DialogContentText>
      </DialogContent>
      <DialogActions>
	<Button onClick={touchFiles} color="primary">
	  Update
	</Button>
	<Button disabled={files.length==0} onClick={() => setConfirmTouch(false)} color="secondary" autoFocus>
	  Cancel
	</Button>
      </DialogActions>
    </Dialog>
  </Grid>
};

export default TestFileMatch;
