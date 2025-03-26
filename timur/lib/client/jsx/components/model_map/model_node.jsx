import React from 'react';
import { makeStyles } from '@material-ui/core/styles';

const useStyles = makeStyles((theme) => ({
  model_node: {
    cursor: 'pointer',
    position: 'absolute',

    fontSize: '0.75em',

    maxWidth: 'var(--model-width)',
    overflow: 'hidden',
    textOverflow: 'ellipsis',

    borderRadius: '3px',
    border: '2px solid #aaa',
    background: '#f2f2f2',
    boxShadow: '0 0 2px 0 #ccc',
    padding: '5px',

    transform: 'translate(-50%, -50%)',
    '&:hover': {
      background: '#cec',
      maxWidth: 'unset',
      overflow: 'unset',
      zIndex: 1000
    }
  },
  selected: {
    border: '2px solid #0c0',
    background: '#cfc',
  },
  disabled: {
    background: '#eee',
    border: 'none',
    color: '#ccc'
  }
}));

const ModelNode = ({model_name, center, size, selected, handler, disabled, numSiblings}) => {

  const classes = useStyles();
  if (!center) return null;
  return <div
      className={`${classes.model_node} ${selected ? classes.selected : ''} ${
        disabled ? classes.disabled : ''
      }`}

      style={{
        top: center.y,
        left: center.x,
        maxWidth: `${600 / (numSiblings + 1) - 15}px`
      }}

      onClick={() => {
        if (!disabled) handler(model_name);
      }}
    >
      {model_name}
  </div>
}

export default ModelNode;
