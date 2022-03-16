import React from 'react';

import Icon from 'etna-js/components/icon';

import {workflowName} from '../../selectors/workflow_selectors';
import ImageMemo from './image_memo';
import Tag from '../tag';
import Typography from '@material-ui/core/Typography';

// To get webpack to pick up the files.
require('../../../img/umap.png');
require('../../../img/barplot.png');
require('../../../img/scatter.png');
require('../../../img/yplot.png');
require('../../../img/default.png');
require('../../../img/add_integers.png');

function MultiLineOutput(iterator) {
  return iterator.map((item, ind) => (
    <div key={ind} className='sub_item'>
      {item}
    </div>
  ));
}
import {makeStyles} from '@material-ui/core/styles';

const useStyles = makeStyles((theme) => ({
  card: {
    margin: '2rem',
    border: '1px solid lightgray',
    boxShadow: '0 0 0 15px #eee, 0 0 4px 15px #aaa',
    width: '260px',
    height: '300px',
    cursor: 'pointer',
    display: 'flex',
    flexDirection: 'column'
  },
  selectedCard: {
    margin: '2rem',
    border: '1px solid lightgray',
    boxShadow: '0 0 0 15px #ffc060, 0 0 15px 13px #fcc',
    width: '260px',
    height: '300px',
    cursor: 'pointer',
    display: 'flex',
    flexDirection: 'column'
  },
  image: {
    height: '60%',
    overflow: 'hidden',
    margin: '0',
    borderTopLeftRadius: '2px',
    borderTopRightRadius: '2px',
    '& img': {
      display: 'block',
      width: '150%',
      margin: '-25%'
    }
  },
  description: {
    display: 'flex',
    flexDirection: 'column',
    fontSize: '12px',
    flex: 1,
    borderTop: '1px solid lightgray'
  },
  row: {
    flex: 1,
    display: 'flex'
  },
  label: {
    cursor: 'default',
    borderRight: '1px solid forestgreen',
    display: 'inline-block',
    flex: '0 0 140px',
    paddingRight: '8px',
    paddingTop: '4px',
    textAlign: 'right',
    color: 'darkblue'
  },
  cardtags: {
    fontSize: '1em'
  },
  value: {
    display: 'inline-block',
    padding: '4px 8px',
    boxSizing: 'border-box',
    flex: 1
  },
  sub_item: {
    borderBottom: '1px solid #eee'
  }
}));

export default function Card({workflow, onClick, selected}) {
  const classes = useStyles();
  return (
    <div
      className={selected ? classes.selectedCard : classes.card}
      onClick={onClick}
    >
      <figure className={classes.image}>
        <ImageMemo
          src={`/images/${workflow.image || 'default.png'}`}
          alt='Workflow image'
        />
      </figure>
      <div className={classes.description}>
        <div className={classes.row}>
          <div className={classes.value}>
            <Typography variant='subtitle1'>
            {workflow.displayName || workflowName(workflow)}
            </Typography>
          </div>
        </div>
        <div className={classes.row}>
          <div className={classes.value}>
            { workflow.authors.join(', ') }
          </div>
        </div>
        <div className={classes.row}>
          <div className={classes.value}>{workflow.lastModified}</div>
        </div>
        <div className={classes.row}>
          <div className={classes.value}>{
            workflow.tags.map( t => <Tag key={t} label={t}/>)
          }</div>
        </div>
      </div>
    </div>
  );
}
