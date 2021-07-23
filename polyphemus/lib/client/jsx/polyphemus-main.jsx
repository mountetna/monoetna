import React, { useState, useEffect, useCallback } from 'react';
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

const Config = ({config, grammar}) =>  {
  // the base term is config
  //
  const term = grammar.config;

  return <ConfigBlock term={term} name={'config'} config={config}/>;
}

const ConfigList =({config}) => {
  return <React.Fragment> [
    {
      config.map( (k,i) => <Config key={i} config={k}/> ).reduce((p,c) => [ p, ', ', c ])
    } ]
  </React.Fragment>
}

const getBlockTerm = (term, key, block) => {
  if (key in term) return term[key];

  Object.keys(term).forEach(term_key => {
    if (matchTemplate(term_key, key)) {
      return term[term_key];
    }
  })
}

const ConfigBlockItem = ({config, term, name}) => {
  return <Box className='config'>
    <Typography component="h1" variant="h5">
      { name }
      {
      }
    </Typography>
  </Box>
}

const matchProduction = (production, entry) => {
  let match = production.match(/^[\.\*\+] ([^#]*) (# .*)?$/)
}

const matchTermName = (item_name, term_name) => (
  term_name == item_name
    || (term_name[0] == '?' && term_name.slice(1) == item_name)
    || matchProduction(term_name, item_name)
);

const validateBlock = (config, term) => {
  // find a term for each item in config. If there are no matching terms for an item, generate an error.
  Object.keys(config).forEach( item_name => {
    Object.keys(term).find( term_name => matchTermName(item_name, term_name) )
  })
  // If a term has too many matching items, generate an error
  // Return possible terms that may be added to this block
}

const ConfigBlock=({config, name, term}) => {
  // the block shows a hash
  //
  return <Box>
    <Typography>
      { name }
    </Typography>
    {
      Object.keys(config).map( k => <ConfigBlockItem key={k} config={config} term={term} name={k}/>)
    }
  </Box>
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
                    <Config config={ config } grammar={ jobs.find( j => j.name == etl).grammar } />
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
