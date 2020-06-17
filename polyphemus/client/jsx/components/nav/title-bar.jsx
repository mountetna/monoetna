import * as React from 'react'

export default class TitleBar extends React.Component{

  constructor(){

    super();
  }

  render(){
    
    return (

      <div id='title-menu'>

        <button className='title-menu-btn'>
            
          { 'Polyphemus' }
          <br />
          <span>

            { 'NETWORK MANAGER' }
          </span>
        </button>
        <img id='ucsf-logo' src='/img/ucsf_logo_black.png' alt='' />
      </div>
    );
  }
}