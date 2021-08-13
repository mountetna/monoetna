import React, { useEffect } from 'react';
import Nav from 'etna-js/components/Nav';
import {selectUser} from 'etna-js/selectors/user-selector';
import {useReduxState} from 'etna-js/hooks/useReduxState';

import { ReactSVG } from 'react-svg';

import img from 'etna-js/images/polyphemus.svg';

const dist = (p1,p2) => Math.sqrt( (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) );

const norm = (p) => {
  let n = dist(p, {x:0,y:0});
  return { x: p.x / n, y: p.y/n };
}

const moveEye = (e) => {
  const iris = document.querySelector('#iris');

  if (!iris) return;

  const center = { x: 55, y: 55 };
  const maxdist = dist(center, { x: window.innerWidth, y: window.innerHeight});
  const maxradius = 20;

  const point = { x: e.clientX, y: e.clientY };

  const radius = maxradius * Math.sqrt(dist(center, point) / maxdist);

  const vector = norm({ x: point.x - center.x, y: point.y - center.y });

  iris.setAttribute('transform',
      `translate(${vector.x * radius},${vector.y*radius})`)
}

const Logo = () => {

  useEffect( () => {
    document.addEventListener('mousemove', moveEye );
    return () => document.removeEventListener('mousemove', moveEye);
  }, []);

  return <div id='logo'>
      <ReactSVG src="/images/polyphemus.svg" 
        renumerateIRIElements={false}
        beforeInjection={(svg) => {
          svg.setAttribute('style', 'width: 50px; height: 50px')
        }}
      />
    </div>
}

const NavBar = ({user}) => <div id='nav'></div>

const PolyphemusNav = () => {
  let user = useReduxState( state => selectUser(state) );
  return <Nav logo={Logo} user={user} app='polyphemus'>
    <NavBar user={user}/>
  </Nav>;
}

export default PolyphemusNav;
