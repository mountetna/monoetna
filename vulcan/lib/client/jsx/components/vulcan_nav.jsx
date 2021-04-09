// Framework libraries.
import React, {useContext, useEffect, useState, useCallback} from 'react';
import {connect} from 'react-redux';

import {VulcanContext} from '../contexts/vulcan_context';
import {workflowName} from '../selectors/workflow_selectors';
import Nav from 'etna-js/components/Nav';
import Link from 'etna-js/components/link';
import {selectUser} from 'etna-js/selectors/user-selector';

const {sin, cos, PI, random, max, min, pow, abs, sqrt} = Math;

const RAD = PI / 180;

const JITTERS = Array(300).fill().map(random);

const tween = (a, b, pct) => a + (b - a) * pct;
const roll = (a, b) => random() * (b - a) + a;

const x = (r, deg) => cos(RAD * deg) * r;
const y = (r, deg) => sin(RAD * deg) * r;

const Segment = function ({sx, sy, radius}) {
  const numsegs = 20;
  const tx = (r, deg) => x(r, deg) + radius;
  const ty = (r, deg) => y(r, deg) + radius;

  let counter = 0;

  const start = roll(10, 270);
  const stop = roll(10, 270);
  const lifetime = roll(10, 30);
  const offset = roll(0, 360);
  const velocity = (stop - start) / lifetime;
  const amp = max(roll(velocity / 3, velocity), 6);
  const freq = roll(0, 2);
  const pulse = roll(0, 6);

  this.color = '#ffe';
  this.width = min(3, 2 / pow(abs(velocity), 0.5));

  this.update = () => {
    counter = counter + 1;

    return counter < lifetime;
  };

  this.path = () => {
    let deg = (start + counter * velocity) % 360;
    let ox = tx(radius, deg);
    let oy = ty(radius, deg);
    let m = (sx - ox) / (oy - sy);
    let dx = sqrt(1 / (m ** 2 + 1));
    let dy = m * dx;

    return (
      `M${sx} ${sy} ` +
      Array(numsegs + 1)
        .fill()
        .map((_, i) => {
          let ix = tween(sx, ox, i / numsegs);
          let iy = tween(sy, oy, i / numsegs);

          let jitter =
            1 *
              // ease
              (-(cos(PI * (1 - 2 * abs(i / numsegs - 0.5))) - 1) / 2) *
              // basic sin
              sin(counter * pulse + (freq * PI * i) / numsegs) *
              amp +
            // random jitter
            JITTERS[(counter + i) % 300] * 3 -
            1.5;

          return `L${ix + jitter * dx} ${iy + jitter * dy}`;
        })
        .join(' ')
    );
  };
};

const Halo = ({radius}) => {
  let [segments, updateSegments] = useState([]);
  let sr = 0.85 * radius;
  let sth = -15;

  let sx = x(sr, sth) + radius;
  let sy = y(sr, sth) + radius;

  let updatePath = useCallback(() => {
    let runningSegments = segments.filter((s) => s.update());
    let needsSegment = random() < 1 / (runningSegments.length + 1) ** 2;
    updateSegments([
      ...runningSegments,
      ...(needsSegment ? [new Segment({sx, sy, radius})] : [])
    ]);
  }, [segments]);

  useEffect(() => {
    const timer = setInterval(updatePath, 50);
    return () => clearInterval(timer);
  }, [segments]);

  return (
    <div className='halo'>
      <svg>
        <defs>
          <radialGradient id='glow'>
            <stop offset='0%' stopColor='#fdd' stopOpacity='100%' />
            <stop offset='20%' stopColor='#fdc' stopOpacity='90%' />
            <stop offset='100%' stopColor='orange' stopOpacity='0%' />
          </radialGradient>
        </defs>
        {segments.map((seg, i) => (
          <path
            key={i}
            d={seg.path()}
            style={{
              stroke: seg.color,
              strokeWidth: seg.width,
              filter: `drop-shadow(0px 0px 2px ${seg.color})`
            }}
          />
        ))}
      </svg>
    </div>
  );
};

function Logo() {
  let {state} = useContext(VulcanContext);
  // TODO: Change this to look for running steps.
  const {calculating} = state;

  return (
    <div id='vulcan-logo'>
      <div className='image' />
      {calculating && <Halo radius={25} />}
    </div>
  );
}

const getTabs = (workflow) => ({
  workflow: ROUTES.workflow(workflowName(workflow)),
  help: 'https://mountetna.github.io/vulcan.html'
});

const ModeBar = ({mode, workflow}) => (
  <div id='nav'>
    {Object.entries(getTabs(workflow)).map(([tab_name, route]) => (
      <div
        key={tab_name}
        className={`nav_tab ${mode == tab_name ? 'selected' : ''}`}
      >
        <Link link={route}>{tab_name}</Link>
      </div>
    ))}
  </div>
);

const VulcanNav = ({mode, user}) => {
  let {state} = useContext(VulcanContext);
  let {workflow} = state;

  return (
    <Nav user={user} logo={Logo} app='vulcan'>
      {mode !== 'home' && <ModeBar mode={mode} workflow={workflow} />}
    </Nav>
  );
};

export default connect((state) => ({
  user: selectUser(state)
}))(VulcanNav);
