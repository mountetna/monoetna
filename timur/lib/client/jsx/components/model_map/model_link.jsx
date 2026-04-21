import React from 'react';

const dist = (p1, p2) =>
  Math.sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));

const ModelLink = ({center, parent, size, type}) => {
  if (!parent || !center) return null;

  let scale = size / dist(center, parent);

  let p1 = parent;
  let p2 = { x: center.x, y: center.y - size / 2 };
  let s = dist(p1, p2);
  let tan = { x: (p1.y - p2.y)/s, y: (p2.x - p1.x)/s };
  if (p1.x < p2.x) {
    tan.x = - tan.x;
    tan.y = - tan.y;
  }
  let mid = { x: (p1.x + p2.x) / 2 + tan.x * s / 3, y: (p1.y + p2.y) / 2 + tan.y * s  / 3}

  let element = type == 'link'
    ? <path
        d={`M ${p1.x} ${p1.y}
          Q ${mid.x} ${mid.y}
            ${p2.x} ${p2.y}`}
        stroke="goldenrod"
        strokeWidth="3"
        strokeOpacity="0.1"
        fill='none'
        markerEnd='url(#path_arrow)'
      />
    : <line
        x2={center.x + (parent.x - center.x) * scale}
        y2={center.y - size / 2}
        x1={parent.x}
        y1={parent.y}
        stroke="#aca"
        strokeWidth="3"
        markerEnd='url(#line_arrow)'
      />

  return (
    <g className='model_link'>
      {element}
    </g>
  );
};

export const Arrowhead = (props) => (
  <marker
    {...props}
    markerWidth='3'
    markerHeight='3'
    refX='0'
    refY='1'
    orient='auto'
    markerUnits='strokeWidth'
  >
    <path d='M0,0 L0,2 L3,1 z' />
  </marker>
);

export default ModelLink;
