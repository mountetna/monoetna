import React, { useState, useEffect, useCallback, useRef } from 'react';
import { json_get } from 'etna-js/utils/fetch';
import { getType } from 'etna-js/utils/types';

import Typography from '@material-ui/core/Typography';
import Card from '@material-ui/core/Card';
import CardContent from '@material-ui/core/CardContent';
import {makeStyles} from '@material-ui/core/styles';
import CardActions from '@material-ui/core/CardActions';
import Button from '@material-ui/core/Button';
import Box from '@material-ui/core/Box';
import Collapse from '@material-ui/core/Collapse';


import { JSONEditor, defaults } from '@json-editor/json-editor';

JSONEditor.defaults.options.collapsed = true;

const Config = ({config, schema}) =>  {
  const editorRef = useRef();
  useEffect(() => {
    const editor = new JSONEditor( editorRef.current, {
      schema, startval: config,
      display_required_only: true,
      theme: 'barebones',
      iconlib: 'fontawesome5',
      remove_button_labels: true,
      compact: true,
      disable_array_add: false,
      disable_array_delete: false,
      disable_array_delete_all_rows: true,
      disable_array_delete_last_row: true,
      disable_array_reorder: true,
      disable_edit_json: true,
      disable_properties: false
    });
  });
  return <div ref={editorRef}/>
}

const PolyphemusMain = ({project_name}) => {
  const [ etls, setEtls ] = useState(null);
  const [ jobs, setJobs ] = useState(null);
  const [ expanded, setExpanded ] = useState(false);

  const toggleExpanded = () => setExpanded(!expanded);

  useEffect( () => {
    json_get(`/api/${project_name}/etl/configs`).then(setEtls)
    json_get(`/api/${project_name}/etl/jobs`).then(setJobs)
  }, [] );
  return <div id='polyphemus-main'>
    {project_name}
    {
      !(jobs && etls) ? null :
        <div className='etls'>
          <span className='title'>ETLs</span>
          {
            etls.map( ({etl,name,config}) => <Card className='etl' key={etl}>
              <CardContent>
                <Typography component="h5" variant="h5">
                  { name }
                </Typography>

                <CardActions>
                  <Button variant="text" onClick={toggleExpanded} size="small">Configure</Button>
                </CardActions>
                <Collapse in={expanded} timeout="auto" unmountOnExit>
                  <CardContent>
                    <Config config={ config } schema={ jobs.find( j => j.name == etl).schema } />
                  </CardContent>
                </Collapse>
              </CardContent>
            </Card> )
          }
        </div>
    }
  </div>;
}

export default PolyphemusMain;
