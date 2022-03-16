import React from 'react';
import Chip from '@material-ui/core/Chip';
import hash from 'object-hash';

//const colors = ['#f77189', '#bb9832', '#50b131', '#36ada4', '#3ba3ec', '#e866f4'];
const colors = ['#f77189', '#ce9032', '#97a431', '#32b166', '#36ada4', '#39a7d0', '#a48cf4', '#f561dd']

//const colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3']


const Tag = (props:any) => {
  const backgroundColor = colors[ parseInt( hash(props.label).slice(0,8), 16) % colors.length ]
  return <Chip style={{backgroundColor}} {...props}/>;
}

export default Tag;
